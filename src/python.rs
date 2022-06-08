use std::collections::HashMap;
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
use crate::{json, ncbi, newick, phyloxml, prune_away, prune_to, GeneralTaxonomy};

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
    let gil_guard = Python::acquire_gil();
    let py = gil_guard.python();
    match val {
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
    }
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
    fn __richcmp__(&self, other: PyRef<TaxonomyNode>, op: CompareOp) -> Py<PyAny> {
        let py = other.py();
        match op {
            CompareOp::Eq => (self == other.deref()).into_py(py),
            CompareOp::Ne => (self != other.deref()).into_py(py),
            _ => py.NotImplemented(),
        }
    }

    fn __getitem__(&self, obj: &PyAny) -> PyResult<PyObject> {
        let key: &str = obj.extract()?;
        let gil_guard = Python::acquire_gil();
        let py = gil_guard.python();
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
}

/// The Taxonomy object provides the primary interface for exploring a
/// biological taxonomy.
#[pyclass]
#[derive(Debug)]
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

    pub(crate) fn get_ncbi_rank(&self, tax_id: &str) -> PyResult<&str> {
        let rank = py_try!(self.tax.rank(tax_id)).to_ncbi_rank();
        Ok(rank)
    }

    pub(crate) fn as_node(&self, tax_id: &str) -> PyResult<TaxonomyNode> {
        let name = self.get_name(tax_id)?;
        let rank = self.get_ncbi_rank(tax_id)?;
        let parent = py_try!(self.tax.parent(tax_id)).map(|(p, _)| p.to_string());
        let extra = py_try!(self.tax.data(tax_id));

        Ok(TaxonomyNode {
            id: tax_id.to_string(),
            name: name.to_string(),
            rank: rank.to_string(),
            extra: (*extra).to_owned(),
            parent,
        })
    }
}

#[pymethods]
impl Taxonomy {
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

    /// to_json_tree(self)
    /// --
    ///
    /// Export a Taxonomy as a JSON-encoded byte string in a tree format
    fn to_json_tree(&self) -> PyResult<PyObject> {
        let mut bytes = Vec::new();
        py_try!(json::save::<_, &str, _>(
            &mut bytes,
            &self.tax,
            JsonFormat::Tree,
            None
        ));
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyBytes::new(py, &bytes).into())
    }

    /// to_json_node_links(self)
    /// --
    ///
    /// Export a Taxonomy as a JSON-encoded byte string in a node link format
    fn to_json_node_links(&self) -> PyResult<PyObject> {
        let mut bytes = Vec::new();
        py_try!(json::save::<_, &str, _>(
            &mut bytes,
            &self.tax,
            JsonFormat::NodeLink,
            None
        ));
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyBytes::new(py, &bytes).into())
    }

    /// to_newick(self)
    /// --
    ///
    /// Export a Taxonomy as a Newick-encoded byte string.
    fn to_newick(&self) -> PyResult<PyObject> {
        let mut bytes = Vec::new();
        py_try!(newick::save(
            &mut bytes,
            &self.tax,
            Some(TaxonomyTrait::<InternalIndex>::root(&self.tax))
        ));
        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyBytes::new(py, &bytes).into())
    }

    /// node(self, tax_id: str) -> Optional[TaxonomyNode]
    /// --
    ///
    /// Find a node by its id. Returns `None` if not found
    fn node(&self, tax_id: &str) -> Option<TaxonomyNode> {
        self.as_node(tax_id).ok()
    }

    /// find_by_name(self, name: str) -> TaxonomyNode
    /// --
    ///
    /// Find a node by its name, Raises an exception if not found.
    fn find_by_name(&self, name: &str) -> Option<TaxonomyNode> {
        if let Some(tax_id) = self.tax.find_by_name(name) {
            self.as_node(tax_id).ok()
        } else {
            None
        }
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
    /// Return a list of child taxonomy nodes from the node id provided.
    fn children(&self, tax_id: &str) -> PyResult<Vec<TaxonomyNode>> {
        let children: Vec<&str> = py_try!(self.tax.children(tax_id));
        let mut res = Vec::with_capacity(children.len());
        for key in children {
            let child = self.node(key);
            if let Some(c) = child {
                res.push(c);
            } else {
                return Err(PyErr::new::<TaxonomyError, _>(format!(
                    "Node {} is missing in children",
                    key
                )));
            }
        }

        Ok(res)
    }

    /// lineage(self, tax_id: str)
    /// --
    ///
    /// Return a list of all the parent taxonomy nodes of the node id provided
    /// (including that node itself).
    fn lineage(&self, tax_id: &str) -> PyResult<Vec<TaxonomyNode>> {
        let lineage: Vec<&str> = py_try!(self.tax.lineage(tax_id));
        let mut res = Vec::with_capacity(lineage.len());
        for key in lineage {
            let ancestor = self.node(key);
            if let Some(a) = ancestor {
                res.push(a);
            } else {
                return Err(PyErr::new::<TaxonomyError, _>(format!(
                    "Node {} is missing in lineage",
                    key
                )));
            }
        }

        Ok(res)
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

    /// add_node(self, parent_id: str, tax_id: str)
    /// --
    ///
    /// Add a new node to the tree at the parent provided.
    fn add_node(&mut self, parent_id: &str, tax_id: &str) -> PyResult<()> {
        py_try!(self.tax.add(parent_id, tax_id));
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
            let parent = py_try!(self.tax.parent(p))
                .ok_or_else(|| PyErr::new::<PyKeyError, _>("New parent ID is not in taxonomy"))?
                .0;
            let lineage = py_try!(self.tax.lineage(parent), "New parent has bad lineage?");
            if lineage.contains(&tax_id) {
                return Err(PyErr::new::<TaxonomyError, _>(
                    "Node can not be moved to its child",
                ));
            }
            self.tax.parent_ids[idx] = py_try!(self.tax.to_internal_index(parent));
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

    fn __iter__(slf: PyRefMut<Self>) -> PyResult<TaxonomyIterator> {
        let root = slf.tax.root();
        let root_idx = slf.tax.to_internal_index(root).unwrap();
        let gil_guard = Python::acquire_gil();
        let py = gil_guard.python();
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
    fn __next__(mut slf: PyRefMut<Self>) -> PyResult<Option<String>> {
        let gil_guard = Python::acquire_gil();
        let py = gil_guard.python();
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
