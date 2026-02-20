use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::io::Cursor;
use std::ops::Deref;
use std::str::FromStr;

use pyo3::basic::CompareOp;
use pyo3::create_exception;
use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyList, PyType};
use serde_json::Value;

use crate::base::InternalIndex;
use crate::json::JsonFormat;
use crate::rank::TaxRank;
use crate::Taxonomy as TaxonomyTrait;
use crate::{gtdb, json, ncbi, newick, phyloxml, prune_away, prune_to, GeneralTaxonomy};

create_exception!(taxonomy, TaxonomyError, pyo3::exceptions::PyException);

// Avoid some boilerplate with the error handling
macro_rules! py_try {
    ($call:expr) => {
        $call.map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?
    };
    ($call:expr, $msg:expr) => {
        $call.map_err(|_| PyErr::new::<TaxonomyError, _>($msg.to_owned()))?
    };
}

fn json_value_to_pyobject(val: &Value) -> PyObject {
    Python::with_gil(|py| match val {
        Value::Null => py.None(),
        Value::Bool(b) => b.to_object(py),
        Value::Number(n) => {
            if let Some(n1) = n.as_i64() {
                return n1.to_object(py);
            }
            n.as_f64().unwrap().to_object(py)
        }
        Value::String(s) => s.to_object(py),
        Value::Array(arr) => {
            let pylist = PyList::empty(py);
            for v in arr {
                pylist
                    .append(json_value_to_pyobject(v))
                    .expect("can add items to list");
            }
            pylist.to_object(py)
        }
        Value::Object(obj) => {
            let pydict = PyDict::new(py);
            for (key, val) in obj.iter() {
                pydict
                    .set_item(key, json_value_to_pyobject(val))
                    .expect("can add items to dict");
            }
            pydict.to_object(py)
        }
    })
}

/// The data returned when looking up a taxonomy by id or by name
#[pyclass]
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct TaxonomyNode {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    parent: Option<String>,
    #[pyo3(get)]
    rank: String,
    // Ideally this would be private
    extra: HashMap<String, Value>,
}

#[pymethods]
impl TaxonomyNode {
    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.id.hash(&mut hasher);
        self.name.hash(&mut hasher);
        self.parent.hash(&mut hasher);
        self.rank.hash(&mut hasher);
        for key in self.extra.keys() {
            key.hash(&mut hasher);
        }
        hasher.finish()
    }

    fn __richcmp__(&self, other: PyRef<TaxonomyNode>, op: CompareOp) -> Py<PyAny> {
        let py = other.py();
        match op {
            CompareOp::Eq => (self == other.deref()).into_py(py),
            CompareOp::Ne => (self != other.deref()).into_py(py),
            _ => py.NotImplemented(),
        }
    }

    fn __getitem__(&self, obj: &PyAny, py: Python<'_>) -> PyResult<PyObject> {
        let key: &str = obj.extract()?;
        match key {
            "id" => Ok(self.id.to_object(py)),
            "name" => Ok(self.name.to_object(py)),
            "parent" => Ok(self.parent.to_object(py)),
            "rank" => Ok(self.rank.to_object(py)),
            _ => {
                if self.extra.contains_key(key) {
                    Ok(json_value_to_pyobject(self.extra.get(key).unwrap()))
                } else {
                    return Err(PyKeyError::new_err(format!("Key {} not found", key)));
                }
            }
        }
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "<TaxonomyNode (id=\"{}\" rank=\"{}\" name=\"{}\")>",
            self.id, self.rank, self.name
        ))
    }

    /// get_data(self)
    /// --
    ///
    /// Get all extra data fields as a dictionary.
    ///
    /// Returns:
    ///     dict: Dictionary containing all additional fields from the taxonomy
    fn get_data(&self, py: Python<'_>) -> PyResult<PyObject> {
        let pydict = PyDict::new(py);
        for (key, val) in self.extra.iter() {
            pydict.set_item(key, json_value_to_pyobject(val))?;
        }
        Ok(pydict.to_object(py))
    }

    /// get_data_keys(self)
    /// --
    ///
    /// Get list of all available data field keys.
    ///
    /// Returns:
    ///     list: List of field names available in the extra data
    fn get_data_keys(&self, py: Python<'_>) -> PyResult<PyObject> {
        let pylist = PyList::empty(py);
        for key in self.extra.keys() {
            pylist.append(key)?;
        }
        Ok(pylist.to_object(py))
    }

    /// get(self, key: str, default=None)
    /// --
    ///
    /// Get a data field value with optional default.
    ///
    /// Args:
    ///     key: Field name to retrieve
    ///     default: Default value if field doesn't exist (default: None)
    ///
    /// Returns:
    ///     Value of the field or default if not found
    fn get(&self, key: &str, default: Option<&PyAny>, py: Python<'_>) -> PyResult<PyObject> {
        if self.extra.contains_key(key) {
            Ok(json_value_to_pyobject(self.extra.get(key).unwrap()))
        } else if let Some(def) = default {
            Ok(def.to_object(py))
        } else {
            Ok(py.None())
        }
    }

    /// Convenience properties for common NCBI fields
    #[getter]
    fn genetic_code_id(&self, py: Python<'_>) -> PyResult<PyObject> {
        self.get("genetic_code_id", None, py)
    }

    #[getter]
    fn embl_code(&self, py: Python<'_>) -> PyResult<PyObject> {
        self.get("embl_code", None, py)
    }

    #[getter]
    fn division_id(&self, py: Python<'_>) -> PyResult<PyObject> {
        self.get("division_id", None, py)
    }

    #[getter]
    fn mitochondrial_genetic_code_id(&self, py: Python<'_>) -> PyResult<PyObject> {
        self.get("mitochondrial_genetic_code_id", None, py)
    }
}

/// The Taxonomy object provides the primary interface for exploring a
/// biological taxonomy.
#[pyclass]
#[derive(Debug, Clone)]
pub struct Taxonomy {
    tax: GeneralTaxonomy,
}

/// Some private fns we don't want to share in Python but that make the Python code easier to
/// write
impl Taxonomy {
    pub(crate) fn get_name<'t>(&'t self, tax_id: &'t str) -> PyResult<&'t str> {
        let name = py_try!(self.tax.name(tax_id));
        Ok(name)
    }

    pub(crate) fn get_rank(&self, tax_id: &str) -> PyResult<String> {
        let rank = py_try!(self.tax.rank(tax_id)).to_string();
        Ok(rank)
    }

    pub(crate) fn as_node(&self, tax_id: &str) -> PyResult<TaxonomyNode> {
        let name = self.get_name(tax_id)?;
        let rank = self.get_rank(tax_id)?;
        let parent = py_try!(self.tax.parent(tax_id)).map(|(p, _)| p.to_string());
        let extra = py_try!(self.tax.data(tax_id));

        Ok(TaxonomyNode {
            id: tax_id.to_string(),
            name: name.to_string(),
            rank,
            extra: (*extra).to_owned(),
            parent,
        })
    }
}

#[pymethods]
impl Taxonomy {
    /// from_gtdb(cls, value: str)
    /// --
    ///
    /// Load a Taxonomy from a GTDB-encoded string.
    #[classmethod]
    fn from_gtdb(_cls: &PyType, value: &str) -> PyResult<Taxonomy> {
        let mut c = Cursor::new(value);
        let tax = py_try!(gtdb::load(&mut c));
        Ok(Taxonomy { tax })
    }

    /// from_json(cls, value: str, /, json_pointer: str)
    /// --
    ///
    /// Load a Taxonomy from a JSON-encoded string. The format can either be
    /// of the tree or node_link_data types and will be automatically detected.
    /// If `path` is specified, the JSON will be traversed to that sub-object
    /// before being parsed as a taxonomy. `path` has to be a valid JSON path string.
    #[classmethod]
    fn from_json(_cls: &PyType, value: &str, json_pointer: Option<&str>) -> PyResult<Taxonomy> {
        let mut c = Cursor::new(value);
        let tax = py_try!(json::load(&mut c, json_pointer));
        Ok(Taxonomy { tax })
    }

    /// from_newick(cls, value: str)
    /// --
    ///
    /// Load a Taxonomy from a Newick-encoded string.
    #[classmethod]
    fn from_newick(_cls: &PyType, value: &str) -> PyResult<Taxonomy> {
        let mut c = Cursor::new(value);
        let tax = py_try!(newick::load(&mut c));
        Ok(Taxonomy { tax })
    }

    /// from_ncbi(cls, dump_dir: str)
    /// --
    ///
    /// Load a Taxonomy from a directory.
    /// The directory must contain the `nodes.dmp` and `names.dmp` files.
    ///
    /// All fields from both nodes.dmp and names.dmp are loaded and accessible
    /// via the node's data methods or dict-like interface.
    ///
    /// Args:
    ///     dump_dir: Path to directory containing NCBI taxonomy dump files
    ///
    /// Returns:
    ///     Taxonomy: Loaded taxonomy with all NCBI fields accessible
    #[classmethod]
    fn from_ncbi(_cls: &PyType, dump_dir: &str) -> PyResult<Taxonomy> {
        let tax = py_try!(ncbi::load(dump_dir));
        Ok(Taxonomy { tax })
    }

    /// from_phyloxml(cls, value: str)
    /// --
    ///
    /// Load a Taxonomy from a PhyloXML-encoded string.
    ///
    /// Experimental.
    #[classmethod]
    fn from_phyloxml(_cls: &PyType, value: &str) -> PyResult<Taxonomy> {
        let mut c = Cursor::new(value);
        let tax = py_try!(phyloxml::load(&mut c));
        Ok(Taxonomy { tax })
    }

    /// clone(self)
    /// --
    ///
    /// Clone the current taxonomy
    pub fn clone(&self) -> Taxonomy {
        Clone::clone(self)
    }

    /// to_json_tree(self)
    /// --
    ///
    /// Export a Taxonomy as a JSON-encoded byte string in a tree format
    fn to_json_tree(&self, py: Python<'_>) -> PyResult<PyObject> {
        let mut bytes = Vec::new();
        py_try!(json::save::<_, &str, _>(
            &mut bytes,
            &self.tax,
            JsonFormat::Tree,
            None
        ));
        Ok(PyBytes::new(py, &bytes).into())
    }

    /// to_json_node_links(self)
    /// --
    ///
    /// Export a Taxonomy as a JSON-encoded byte string in a node link format
    fn to_json_node_links(&self, py: Python<'_>) -> PyResult<PyObject> {
        let mut bytes = Vec::new();
        py_try!(json::save::<_, &str, _>(
            &mut bytes,
            &self.tax,
            JsonFormat::NodeLink,
            None
        ));
        Ok(PyBytes::new(py, &bytes).into())
    }

    /// to_ncbi(self, output_dir: str)
    /// --
    ///
    /// Export a Taxonomy to NCBI format files (nodes.dmp and names.dmp).
    /// The output directory will be created if it doesn't exist.
    ///
    /// Args:
    ///     output_dir: Path to the directory where nodes.dmp and names.dmp will be written
    fn to_ncbi(&self, output_dir: &str) -> PyResult<()> {
        py_try!(ncbi::save::<&str, _, _>(&self.tax, output_dir));
        Ok(())
    }

    /// to_newick(self)
    /// --
    ///
    /// Export a Taxonomy as a Newick-encoded byte string.
    fn to_newick(&self, py: Python<'_>) -> PyResult<PyObject> {
        let mut bytes = Vec::new();
        py_try!(newick::save(
            &mut bytes,
            &self.tax,
            Some(TaxonomyTrait::<&str>::root(&self.tax))
        ));
        Ok(PyBytes::new(py, &bytes).into())
    }

    /// node(self, tax_id: str) -> Optional[TaxonomyNode]
    /// --
    ///
    /// Find a node by its id. Returns `None` if not found
    fn node(&self, tax_id: &str) -> Option<TaxonomyNode> {
        self.as_node(tax_id).ok()
    }

    /// find_all_by_name(self, name: str) -> List[TaxonomyNode]
    /// --
    ///
    /// Find a node by its name, Raises an exception if not found.
    fn find_all_by_name(&self, name: &str) -> PyResult<Vec<TaxonomyNode>> {
        let res = self
            .tax
            .find_all_by_name(name)
            .into_iter()
            .map(|tax_id| self.as_node(tax_id))
            .collect::<PyResult<Vec<TaxonomyNode>>>()?;
        Ok(res)
    }

    /// parent_with_distance(self, tax_id: str, /, at_rank: str)
    /// --
    ///
    /// Return the immediate parent taxonomy node of the node id provided and the distance to it
    ///
    /// If `at_rank` is provided, scan all the nodes in the node's lineage and return
    /// the parent id at that rank.
    fn parent_with_distance(
        &self,
        tax_id: &str,
        at_rank: Option<&str>,
    ) -> PyResult<(Option<TaxonomyNode>, Option<f32>)> {
        let parent_res = if let Some(rank) = at_rank {
            if let Ok(rank) = TaxRank::from_str(rank) {
                self.tax.parent_at_rank(tax_id, rank)
            } else {
                return Err(PyErr::new::<TaxonomyError, _>(format!(
                    "Rank {} could not be understood",
                    rank
                )));
            }
        } else {
            self.tax.parent(tax_id)
        };

        if let Ok(Some((id, distance))) = parent_res {
            Ok((self.as_node(id).ok(), Some(distance)))
        } else {
            Ok((None, None))
        }
    }

    /// parent(self, tax_id: str, /, at_rank: str)
    /// --
    ///
    /// Return the immediate parent taxonomy node of the node id provided.
    ///
    /// If `at_rank` is provided, scan all the nodes in the node's lineage and return
    /// the parent id at that rank.
    fn parent(&self, tax_id: &str, at_rank: Option<&str>) -> PyResult<Option<TaxonomyNode>> {
        let (node, _) = self.parent_with_distance(tax_id, at_rank)?;
        Ok(node)
    }

    /// children(self, tax_id: str)
    /// --
    ///
    /// Return a list of direct child taxonomy nodes from the node id provided.
    fn children(&self, tax_id: &str) -> PyResult<Vec<TaxonomyNode>> {
        let res = py_try!(self.tax.children(tax_id))
            .into_iter()
            .map(|tax_id| self.as_node(tax_id))
            .collect::<PyResult<Vec<TaxonomyNode>>>()?;
        Ok(res)
    }

    /// descendants(self, tax_id: str)
    /// --
    ///
    /// Return a list of all child taxonomy nodes from the node id provided.
    fn descendants(&self, tax_id: &str) -> PyResult<Vec<TaxonomyNode>> {
        let res = py_try!(self.tax.descendants(tax_id))
            .into_iter()
            .map(|tax_id| self.as_node(tax_id))
            .collect::<PyResult<Vec<TaxonomyNode>>>()?;
        Ok(res)
    }

    /// lineage(self, tax_id: str)
    /// --
    ///
    /// Return a list of all the parent taxonomy nodes of the node id provided
    /// (including that node itself).
    fn lineage(&self, tax_id: &str) -> PyResult<Vec<TaxonomyNode>> {
        let res = py_try!(self.tax.lineage(tax_id))
            .into_iter()
            .map(|tax_id| self.as_node(tax_id))
            .collect::<PyResult<Vec<TaxonomyNode>>>()?;
        Ok(res)
    }

    /// data(self, tax_id: str)
    /// --
    ///
    /// Get all extra data fields for a taxonomy node as a dictionary.
    ///
    /// Args:
    ///     tax_id: The taxonomy ID to look up
    ///
    /// Returns:
    ///     dict: Dictionary containing all additional fields from the taxonomy
    fn data(&self, tax_id: &str, py: Python<'_>) -> PyResult<PyObject> {
        let data = py_try!(self.tax.data(tax_id));
        let pydict = PyDict::new(py);
        for (key, val) in data.iter() {
            pydict.set_item(key, json_value_to_pyobject(val))?;
        }
        Ok(pydict.to_object(py))
    }

    /// get_field(self, tax_id: str, field: str, default=None)
    /// --
    ///
    /// Get a specific data field for a taxonomy node.
    ///
    /// Args:
    ///     tax_id: The taxonomy ID to look up
    ///     field: Field name to retrieve
    ///     default: Default value if field doesn't exist (default: None)
    ///
    /// Returns:
    ///     Value of the field or default if not found
    fn get_field(&self, tax_id: &str, field: &str, default: Option<&PyAny>, py: Python<'_>) -> PyResult<PyObject> {
        let data = py_try!(self.tax.data(tax_id));
        if let Some(value) = data.get(field) {
            Ok(json_value_to_pyobject(value))
        } else if let Some(def) = default {
            Ok(def.to_object(py))
        } else {
            Ok(py.None())
        }
    }

    /// internal_index(self, tax_id: str)
    /// --
    ///
    /// Return the internal integer ID generated by the taxonomy library
    fn internal_index(&self, tax_id: &str) -> PyResult<usize> {
        let internal_index = py_try!(self.tax.to_internal_index(tax_id));
        Ok(internal_index)
    }

    /// parents(self, tax_id: str)
    /// --
    ///
    /// Return a list of all the parent taxonomy nodes of the node id provided.
    /// It is equivalent to `lineage` except it doesn't include itself
    fn parents(&self, tax_id: &str) -> PyResult<Vec<TaxonomyNode>> {
        let mut lineage = self.lineage(tax_id)?;
        lineage.drain(..1);
        Ok(lineage)
    }

    /// lca(self, id1: str, id2: str)
    /// --
    ///
    /// Return the lowest common ancestor of two taxonomy nodes.
    fn lca(&self, id1: &str, id2: &str) -> PyResult<Option<TaxonomyNode>> {
        let lca_id = py_try!(self.tax.lca(id1, id2));
        Ok(self.node(lca_id))
    }

    /// prune(self, keep: List[str], remove: List[str])
    /// --
    ///
    /// Return a copy of the taxonomy containing:
    ///  - only the nodes in `keep` and their parents if provided
    ///  - all of the nodes except those in remove and their children if provided
    fn prune(&self, keep: Option<Vec<&str>>, remove: Option<Vec<&str>>) -> PyResult<Taxonomy> {
        let mut tax = self.tax.clone();
        if let Some(k) = keep {
            tax = py_try!(prune_to(&tax, &k, false));
        }
        if let Some(r) = remove {
            tax = py_try!(prune_away(&tax, &r));
        }
        Ok(Taxonomy { tax })
    }

    /// remove_node(self, tax_id: str)
    /// --
    ///
    /// Remove the node from the tree.
    fn remove_node(&mut self, tax_id: &str) -> PyResult<()> {
        py_try!(self.tax.remove(tax_id));
        Ok(())
    }

    /// add_node(self, parent_id: str, tax_id: str, name: str, rank: str)
    /// --
    ///
    /// Add a new node to the tree at the parent provided.
    fn add_node(&mut self, parent_id: &str, tax_id: &str, name: &str, rank: &str) -> PyResult<()> {
        if self.node(tax_id).is_some() {
            return Err(PyErr::new::<TaxonomyError, _>(format!(
                "A node with tax id {} already exists",
                tax_id
            )));
        }
        py_try!(self.tax.add(parent_id, tax_id));
        py_try!(self.edit_node(tax_id, Some(name), Some(rank), None, None));
        Ok(())
    }

    /// edit_node(self, tax_id: str, /, name: str, rank: str, parent_id: str, parent_dist: float)
    /// --
    ///
    /// Edit properties on a taxonomy node.
    fn edit_node(
        &mut self,
        tax_id: &str,
        name: Option<&str>,
        rank: Option<&str>,
        parent_id: Option<&str>,
        parent_distance: Option<f32>,
    ) -> PyResult<()> {
        let idx = py_try!(self.tax.to_internal_index(tax_id));

        if let Some(r) = rank {
            self.tax.ranks[idx] = py_try!(TaxRank::from_str(r), "Rank could not be understood");
        }
        if let Some(n) = name {
            self.tax.names[idx] = n.to_string();
        }
        if let Some(p) = parent_id {
            if tax_id == TaxonomyTrait::<&str>::root(&self.tax) {
                return Err(PyErr::new::<TaxonomyError, _>("Root cannot have a parent"));
            }
            let lineage = py_try!(self.tax.lineage(p), "New parent has bad lineage?");
            if lineage.contains(&tax_id) {
                return Err(PyErr::new::<TaxonomyError, _>(
                    "Node can not be moved to its child",
                ));
            }

            let old_parent_idx = self.tax.parent_ids[idx];
            let new_parent_idx = py_try!(self.tax.to_internal_index(p));

            // update id to parent_id
            self.tax.parent_ids[idx] = py_try!(self.tax.to_internal_index(p));

            // remove current node from children lookup for old parent
            let removal_index = self.tax.children_lookup[old_parent_idx]
                .binary_search(&idx)
                .unwrap();
            self.tax.children_lookup[old_parent_idx].remove(removal_index);

            self.tax.children_lookup[old_parent_idx].sort_unstable();

            // add idx as a child of new parent idx
            self.tax.children_lookup[new_parent_idx].push(idx);
            self.tax.children_lookup[new_parent_idx].sort_unstable();
        }

        if let Some(p) = parent_distance {
            if tax_id == TaxonomyTrait::<&str>::root(&self.tax) {
                return Err(PyErr::new::<TaxonomyError, _>("Root cannot have a parent"));
            }
            self.tax.parent_distances[idx] = p;
        }

        Ok(())
    }

    #[getter]
    fn root(&self) -> TaxonomyNode {
        let key: &str = self.tax.root();
        self.as_node(key).unwrap()
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "<Taxonomy ({} nodes)>",
            TaxonomyTrait::<InternalIndex>::len(&self.tax)
        ))
    }

    fn __len__(&self) -> PyResult<usize> {
        Ok(TaxonomyTrait::<InternalIndex>::len(&self.tax))
    }

    // TODO: way to set name and rank
    // TODO: also way to set parent and parent distance?

    fn __getitem__(&self, tax_id: &str) -> PyResult<TaxonomyNode> {
        self.as_node(tax_id)
    }

    fn __delitem__(&mut self, tax_id: &str) -> PyResult<()> {
        Ok(py_try!(self.tax.remove(tax_id)))
    }

    fn __contains__(&self, tax_id: &str) -> PyResult<bool> {
        Ok(self.tax.to_internal_index(tax_id).is_ok())
    }

    fn __iter__(slf: PyRefMut<Self>, py: Python<'_>) -> PyResult<TaxonomyIterator> {
        let root = slf.tax.root();
        let root_idx = slf.tax.to_internal_index(root).unwrap();
        Ok(TaxonomyIterator {
            t: slf.into_py(py),
            nodes_left: vec![root_idx],
            visited_nodes: Vec::new(),
        })
    }
}

#[pyclass]
pub struct TaxonomyIterator {
    t: PyObject,
    visited_nodes: Vec<usize>,
    nodes_left: Vec<usize>,
}

#[pymethods]
impl TaxonomyIterator {
    fn __next__(mut slf: PyRefMut<Self>, py: Python<'_>) -> PyResult<Option<String>> {
        let traverse_preorder = true;
        loop {
            if slf.nodes_left.is_empty() {
                return Ok(None);
            }

            let cur_node = *slf.nodes_left.last().unwrap();
            let node_visited = {
                let last_visited = slf.visited_nodes.last();
                Some(&cur_node) == last_visited
            };
            let node = if node_visited {
                slf.visited_nodes.pop();
                slf.nodes_left.pop().unwrap() // postorder
            } else {
                slf.visited_nodes.push(cur_node);
                let tax: PyRef<Taxonomy> = slf.t.extract(py)?;
                let cur_node_str = tax.tax.from_internal_index(cur_node).unwrap();
                let children = tax
                    .tax
                    .children(cur_node_str)
                    .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
                let children: Vec<_> = children
                    .into_iter()
                    .map(|c| tax.tax.to_internal_index(c).unwrap())
                    .collect();
                // drop the reference to tax, we don't need it anymore
                drop(tax);
                if !children.is_empty() {
                    slf.nodes_left.extend(children);
                }
                cur_node // preorder
            };
            if node_visited != traverse_preorder {
                let tax: PyRef<Taxonomy> = slf.t.extract(py)?;
                return Ok(Some(
                    tax.tax
                        .from_internal_index(node)
                        .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?
                        .to_string(),
                ));
            }
        }
    }
}

/// The taxonomy module
#[pymodule]
fn taxonomy(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Taxonomy>()?;
    m.add("TaxonomyError", py.get_type::<TaxonomyError>())?;

    Ok(())
}
