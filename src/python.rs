//! Implementation of the Python API.
//!
//! Only enabled when `cargo build --features python` is run or when the
//! Python API is build with `python setup.py develop`.
use std::collections::HashMap;
use std::fs::File;
use std::io::Cursor;
use std::str::FromStr;

use pyo3::class::{PyIterProtocol, PyMappingProtocol, PyObjectProtocol, PySequenceProtocol};
use pyo3::exceptions::KeyError;
use pyo3::prelude::*;
use pyo3::types::{IntoPyDict, PyBytes, PyDict, PyType};
use pyo3::{create_exception, wrap_pyfunction, wrap_pymodule};

use crate::base::{GeneralTaxonomy, IntTaxId};
use crate::edit::{prune_away, prune_to};
use crate::formats::json::{load_json, save_json};
use crate::formats::ncbi::load_ncbi;
use crate::formats::newick::{load_newick, save_newick};
use crate::formats::phyloxml::load_phyloxml;
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy as TaxTrait;
use crate::weights::{maximum_weighted_path as tax_mwp, rollup_weights as tax_rw};

create_exception!(taxonomy, TaxonomyError, pyo3::exceptions::Exception);

/// The data returned when looking up a taxonomy by id or by name
#[pyclass]
#[derive(Debug, Clone)]
pub struct TaxonomyNode {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    parent: Option<String>,
    #[pyo3(get)]
    rank: String,
}

#[pyproto]
impl PyObjectProtocol for TaxonomyNode {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "<TaxonomyNode (id={} rank=\"{}\" name=\"{}\"))>",
            self.id, self.rank, self.name
        ))
    }
}

/// The Taxonomy object provides the primary interface for exploring a
/// biological taxonomy. Iterating over the Taxonomy returns all the taxonomy
/// ids in pre-order fashion and indexing into the Taxonomy object with a
/// taxonomy id will return a tuple of the name and rank of that id.
#[pyclass]
#[derive(Debug)]
pub struct Taxonomy {
    t: GeneralTaxonomy,
}

// Avoid some boilerplate with the error handling
macro_rules! py_try {
    ($call:expr) => {
        $call.map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?
    };
}

/// Some private fns we don't want to share in Python but that make the Python code easier to
/// write
impl Taxonomy {
    pub(crate) fn get_int_id(&self, key: &str) -> PyResult<usize> {
        self.t
            .to_internal_id(key)
            .map_err(|_| PyErr::new::<KeyError, _>("Tax ID is not in taxonomy"))
    }

    pub(crate) fn get_name(&self, key: &str) -> PyResult<&str> {
        let name = py_try!(self.t.name(key));
        Ok(name)
    }

    pub(crate) fn get_rank(&self, key: &str) -> PyResult<&str> {
        let rank = py_try!(self.t.rank(key)).to_ncbi_rank();
        Ok(rank)
    }

    pub(crate) fn as_node(&self, key: &str) -> PyResult<TaxonomyNode> {
        let name = self.get_name(key)?;
        let rank = self.get_rank(key)?;
        let parent = py_try!(self.t.parent(key)).map(|(p, _)| p.to_string());

        Ok(TaxonomyNode {
            id: key.to_string(),
            name: name.to_string(),
            rank: rank.to_string(),
            parent,
        })
    }
}

#[pyproto]
impl PyObjectProtocol for Taxonomy {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "<Taxonomy ({} nodes)>",
            TaxTrait::<IntTaxId, _>::len(&self.t)
        ))
    }
}

#[pymethods]
impl Taxonomy {
    /// from_json(cls, value: str, /, path: List[str])
    /// --
    ///
    /// Load a Taxonomy from a JSON-encoded string. The format can either be
    /// of the tree or node_link_data types and will be automatically detected.
    /// If `path` is specified, the JSON will be traversed to that sub-object
    /// before being parsed as a taxonomy.
    #[classmethod]
    fn from_json(_cls: &PyType, value: &str, path: Option<Vec<&str>>) -> PyResult<Taxonomy> {
        let mut c = Cursor::new(value);
        let t = py_try!(load_json(&mut c, path.as_deref()));
        Ok(Taxonomy { t })
    }

    /// from_newick(cls, value: str)
    /// --
    ///
    /// Load a Taxonomy from a Newick-encoded string.
    #[classmethod]
    fn from_newick(_cls: &PyType, value: &str) -> PyResult<Taxonomy> {
        let mut c = Cursor::new(value);
        let t = py_try!(load_newick(&mut c));
        Ok(Taxonomy { t })
    }

    /// from_ncbi(cls, nodes_path: str, names_path: str)
    /// --
    ///
    /// Load a Taxonomy from a pair of NCBI dump files. The paths specified are
    /// to the individual files in the NCBI taxonomy directory (e.g. nodes.dmp
    /// and names.dmp).
    #[classmethod]
    fn from_ncbi(_cls: &PyType, nodes_path: &str, names_path: &str) -> PyResult<Taxonomy> {
        let nodes_file = File::open(nodes_path)?;
        let names_file = File::open(names_path)?;
        let t = py_try!(load_ncbi(nodes_file, names_file));
        Ok(Taxonomy { t })
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
        let t = py_try!(load_phyloxml(&mut c));
        Ok(Taxonomy { t })
    }

    /// to_json(self, /, as_node_link_data: bool)
    /// --
    ///
    /// Export a Taxonomy as a JSON-encoded byte string. By default, the JSON format
    /// is a tree format unless the `as_node_link_data` parameter is set to True.
    #[args(as_node_link_data = false)]
    fn to_json(&self, as_node_link_data: bool) -> PyResult<PyObject> {
        let mut s = Vec::new();
        py_try!(save_json::<&str, _, _, _>(
            &self.t,
            &mut s,
            None,
            as_node_link_data
        ));

        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyBytes::new(py, &s).into())
    }

    /// to_newick(self)
    /// --
    ///
    /// Export a Taxonomy as a Newick-encoded byte string.
    fn to_newick(&self) -> PyResult<PyObject> {
        let mut s = Vec::new();
        py_try!(save_newick::<&str, _, _, _>(&self.t, &mut s, None));

        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyBytes::new(py, &s).into())
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
        if let Ok(tax_id) = self.t.find_by_name(name) {
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
                self.t.parent_at_rank(tax_id, rank)
            } else {
                return Err(PyErr::new::<TaxonomyError, _>(format!(
                    "Rank {} could not be understood",
                    rank
                )));
            }
        } else {
            self.t.parent(tax_id)
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
        let children: Vec<&str> = py_try!(self.t.children(tax_id));
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
        let lineage: Vec<&str> = py_try!(self.t.lineage(tax_id));
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
        let lca_id = py_try!(self.t.lca(id1, id2));
        Ok(self.node(lca_id))
    }

    /// prune(self, keep: List[str], remove: List[str])
    /// --
    ///
    /// Return a copy of the taxonomy containing:
    ///  - only the nodes in `keep` and their parents if provided
    ///  - all of the nodes except those in remove and their children if provided
    fn prune(&self, keep: Option<Vec<&str>>, remove: Option<Vec<&str>>) -> PyResult<Taxonomy> {
        let mut tax = self.t.clone();
        if let Some(k) = keep {
            tax = py_try!(prune_to(&tax, &k, false));
        }
        if let Some(r) = remove {
            tax = py_try!(prune_away(&tax, &r));
        }
        Ok(Taxonomy { t: tax })
    }

    /// remove_node(self, tax_id: str)
    /// --
    ///
    /// Remove the node from the tree.
    fn remove_node(&mut self, tax_id: &str) -> PyResult<()> {
        let int_id = self.get_int_id(tax_id)?;
        py_try!(self.t.remove(int_id));
        Ok(())
    }

    /// add_node(self, parent_id: str, tax_id: str)
    /// --
    ///
    /// Add a new node to the tree at the parent provided.
    fn add_node(&mut self, parent_id: &str, tax_id: &str) -> PyResult<()> {
        let int_id = self.get_int_id(parent_id)?;
        py_try!(self.t.add(int_id, tax_id));
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
        let int_id = self.get_int_id(tax_id)?;

        if let Some(r) = rank {
            self.t.ranks[int_id] = TaxRank::from_str(r)
                .map_err(|_| PyErr::new::<TaxonomyError, _>("Rank could not be understood"))?;
        }
        if let Some(n) = name {
            self.t.names[int_id] = n.to_string();
        }
        if let Some(p) = parent_id {
            if int_id == TaxTrait::<IntTaxId, _>::root(&self.t) {
                return Err(PyErr::new::<TaxonomyError, _>("Root has no parent"));
            }
            let parent = self
                .t
                .to_internal_id(p)
                .map_err(|_| PyErr::new::<KeyError, _>("New parent ID is not in taxonomy"))?;
            if TaxTrait::<IntTaxId, _>::lineage(&self.t, parent)
                .map_err(|_| PyErr::new::<TaxonomyError, _>("New parent has bad lineage?"))?
                .contains(&int_id)
            {
                return Err(PyErr::new::<TaxonomyError, _>(
                    "Node can not be moved to its child",
                ));
            }
            self.t.parent_ids[int_id] = parent;
        }
        if let Some(p) = parent_distance {
            if int_id == TaxTrait::<IntTaxId, _>::root(&self.t) {
                return Err(PyErr::new::<TaxonomyError, _>("Root has no parent"));
            }
            self.t.parent_dists[int_id] = p;
        }
        Ok(())
    }

    #[getter]
    fn root(&self) -> TaxonomyNode {
        let key: &str = self.t.root();
        self.as_node(key).unwrap()
    }
}

#[pyproto]
impl PyMappingProtocol for Taxonomy {
    fn __len__(&self) -> PyResult<usize> {
        Ok(TaxTrait::<IntTaxId, _>::len(&self.t))
    }

    // TODO: way to set name and rank
    // TODO: also way to set parent and parent distance?

    fn __getitem__(&self, key: &str) -> PyResult<TaxonomyNode> {
        self.as_node(key)
    }

    fn __delitem__(&mut self, key: &str) -> PyResult<()> {
        let int_id = self.get_int_id(key)?;
        py_try!(self.t.remove(int_id));
        Ok(())
    }
}

#[pyproto]
impl PySequenceProtocol for Taxonomy {
    fn __contains__(&self, key: &str) -> PyResult<bool> {
        Ok(self.get_int_id(key).is_ok())
    }
}

#[pyproto]
impl PyIterProtocol for Taxonomy {
    fn __iter__(slf: PyRefMut<Self>) -> PyResult<TaxonomyIterator> {
        let root = TaxTrait::<IntTaxId, _>::root(&slf.t);
        let gil_guard = Python::acquire_gil();
        let py = gil_guard.python();
        Ok(TaxonomyIterator {
            t: slf.into_py(py),
            nodes_left: vec![root],
            visited_nodes: Vec::new(),
        })
    }
}

#[pyclass]
pub struct TaxonomyIterator {
    t: PyObject,
    visited_nodes: Vec<IntTaxId>,
    nodes_left: Vec<IntTaxId>,
}

#[pyproto]
impl PyIterProtocol for TaxonomyIterator {
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
                slf.visited_nodes.push(cur_node.clone());
                let tax: PyRef<Taxonomy> = slf.t.extract(py)?;
                let children = tax
                    .t
                    .children(cur_node)
                    .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
                // drop the reference to tax, we don't need it anymore
                drop(tax);
                if !children.is_empty() {
                    slf.nodes_left.extend(children);
                }
                cur_node // preorder
            };
            if node_visited == !traverse_preorder {
                let tax: PyRef<Taxonomy> = slf.t.extract(py)?;
                return Ok(Some(
                    tax.t
                        .from_internal_id(node)
                        .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?
                        .to_string(),
                ));
            }
        }
    }
}

/// maximum_weighted_path(tax: Taxonomy, weights: Dict[str, float])
/// --
///
#[pyfunction]
fn maximum_weighted_path(tax: &Taxonomy, weights: &PyDict) -> PyResult<Option<(String, f32)>> {
    // TODO: remove this once https://github.com/PyO3/pyo3/pull/702 lands in a release
    let mut hash_weights: HashMap<&str, f32> = HashMap::new();
    for (k, v) in weights.iter() {
        hash_weights.insert(k.extract()?, v.extract()?);
    }
    let weights = hash_weights;

    let mwp = tax_mwp(&tax.t, &weights, false)
        .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;

    Ok(mwp.map(|x| (x.0.to_string(), x.1)))
}

/// rollup_weights(tax: Taxonomy, weights: Dict[str, float])
/// --
///
#[pyfunction]
fn rollup_weights(tax: &Taxonomy, weights: &PyDict) -> PyResult<PyObject> {
    // TODO: remove this once https://github.com/PyO3/pyo3/pull/702 lands in a release
    let mut hash_weights: HashMap<&str, f32> = HashMap::new();
    for (k, v) in weights.iter() {
        hash_weights.insert(k.extract()?, v.extract()?);
    }
    let weights = hash_weights;

    let gil = Python::acquire_gil();
    let py = gil.python();

    let rolled: Vec<(&str, f32)> =
        tax_rw(&tax.t, &weights).map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
    Ok(rolled.into_py_dict(py).into())
}

/// Functions related to calculations using node weights
#[pymodule]
fn weights(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(maximum_weighted_path))?;
    m.add_wrapped(wrap_pyfunction!(rollup_weights))?;

    Ok(())
}

/// The taxonomy module
#[pymodule]
fn taxonomy(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Taxonomy>()?;
    m.add("TaxonomyError", py.get_type::<TaxonomyError>())?;

    let m_weights = wrap_pymodule!(weights)(py);
    m.add("weights", &m_weights)?;
    // this allows `from taxonomy.weights import ...`
    // see https://github.com/PyO3/pyo3/issues/471
    py.import("sys")?
        .dict()
        .get_item("modules")
        .unwrap()
        .set_item("taxonomy.weights", m_weights)?;

    Ok(())
}

#[cfg(test)]
mod test {
    use pyo3::prelude::*;
    use pyo3::types::PyDict;

    use crate::base::test::create_example;
    use crate::taxonomy::Taxonomy as TaxTrait;

    use super::{weights, Taxonomy, TaxonomyNode};
    use crate::TaxRank;

    fn setup_ctx(py: Python) -> PyResult<Option<&PyDict>> {
        let tax = Py::new(
            py,
            Taxonomy {
                t: create_example(),
            },
        )?;
        let ctx = PyDict::new(py);
        ctx.set_item("tax", tax)?;
        ctx.set_item("Taxonomy", py.get_type::<Taxonomy>())?;
        Ok(Some(ctx))
    }

    #[test]
    fn test_get() {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py).unwrap();

        let node: TaxonomyNode = py
            .eval(r#"tax.node("53452")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(node.id, "53452");
        assert_eq!(node.parent, Some("1046".to_string()));
        assert_eq!(node.name, "Lamprocystis");

        let node: Option<TaxonomyNode> = py
            .eval(r#"tax.node("unknown")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert!(node.is_none());
    }

    #[test]
    fn test_find_by_name() {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py).unwrap();

        let node: TaxonomyNode = py
            .eval(r#"tax.find_by_name("Lamprocystis")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(node.id, "53452");
        assert_eq!(node.parent, Some("1046".to_string()));
        assert_eq!(node.name, "Lamprocystis");

        let node: Option<TaxonomyNode> = py
            .eval(r#"tax.find_by_name("unknown")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert!(node.is_none());
    }

    #[test]
    fn test_parent() {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py).unwrap();

        let parent: TaxonomyNode = py
            .eval(r#"tax.parent("53452")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(parent.id, "1046");

        let parent: TaxonomyNode = py
            .eval(r#"tax.parent("53452", at_rank="phylum")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(parent.id, "1224");

        let parent: Option<TaxonomyNode> = py
            .eval(r#"tax.parent("bad id")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert!(parent.is_none());

        let parent: Option<TaxonomyNode> = py
            .eval(r#"tax.parent("bad id", at_rank="species")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert!(parent.is_none());

        // A bad rank does raise
        assert!(py
            .eval(r#"tax.parent("53452", at_rank="bad rank")"#, None, ctx)
            .is_err());
    }

    #[test]
    fn test_children() {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py).unwrap();

        let children: Vec<TaxonomyNode> = py
            .eval(r#"tax.children("53452")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(
            children.iter().map(|c| &c.id).collect::<Vec<_>>(),
            vec!["61598"]
        );
    }

    #[test]
    fn test_lineage() {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py).unwrap();

        let lineage: Vec<TaxonomyNode> = py
            .eval(r#"tax.lineage("1224")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(
            lineage.iter().map(|c| &c.id).collect::<Vec<_>>(),
            vec!["1224", "2", "131567", "1"]
        );
    }

    #[test]
    fn test_parents() {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py).unwrap();

        let lineage: Vec<TaxonomyNode> = py
            .eval(r#"tax.parents("1224")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(
            lineage.iter().map(|c| &c.id).collect::<Vec<_>>(),
            vec!["2", "131567", "1"]
        );
    }

    #[test]
    fn test_lca() {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py).unwrap();

        let lca: TaxonomyNode = py
            .eval(r#"tax.lca("56812", "765909")"#, None, ctx)
            .unwrap()
            .extract()
            .unwrap();
        assert_eq!(lca.id, "1236");
    }

    #[test]
    fn test_info_methods() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;

        let root: TaxonomyNode = py.eval(r#"tax.root"#, None, ctx)?.extract()?;
        assert_eq!(root.id, "1");

        let node: TaxonomyNode = py.eval(r#"tax["1224"]"#, None, ctx)?.extract()?;
        assert_eq!(node.rank, "phylum");
        assert_eq!(node.name, "Proteobacteria");
        Ok(())
    }

    #[test]
    fn test_repr() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;

        let repr: String = py.eval(r#"repr(tax)"#, None, ctx)?.extract()?;
        println!("{}", repr);
        assert_eq!(repr, "<Taxonomy (14 nodes)>");
        Ok(())
    }

    #[test]
    fn test_editing() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;

        let tax: PyRef<Taxonomy> = ctx.unwrap().get_item("tax").unwrap().extract()?;
        let n_tax_items = TaxTrait::<&str, _>::len(&tax.t);
        drop(tax);

        assert!(py
            .eval(r#"tax.add_node("bad id", "new id")"#, None, ctx)
            .is_err());
        py.eval(r#"tax.add_node("1236", "91347")"#, None, ctx)?;
        let tax: PyRef<Taxonomy> = ctx.unwrap().get_item("tax").unwrap().extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), n_tax_items + 1);
        drop(tax);

        assert!(py.eval(r#"tax.remove_node("bad id")"#, None, ctx).is_err());
        py.eval(r#"tax.remove_node("135622")"#, None, ctx)?;
        let tax: PyRef<Taxonomy> = ctx.unwrap().get_item("tax").unwrap().extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), n_tax_items);
        drop(tax);

        assert!(py.eval(r#"del tax["bad id"])"#, None, ctx).is_err());
        py.run(r#"del tax["1046"]"#, None, ctx)?;
        let tax: PyRef<Taxonomy> = ctx.unwrap().get_item("tax").unwrap().extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), n_tax_items - 1);
        drop(tax);

        // FIXME: make these fail?
        // assert!(py.eval(r#"tax.prune(keep=["bad id"])"#, None, ctx).is_err());
        // assert!(py.eval(r#"tax.prune(remove=["bad id"])"#, None, ctx).is_err());
        py.eval(r#"tax.prune(keep=["22"])"#, None, ctx)?;
        py.eval(r#"tax.prune(remove=["22"])"#, None, ctx)?;

        Ok(())
    }

    #[test]
    fn test_edit_node() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py).expect("setup context");

        assert!(py
            .eval(r#"tax.edit_node("bad id", name="New Name")"#, None, ctx)
            .is_err());
        py.eval(r#"tax.edit_node("135622", name="New Name")"#, None, ctx)?;
        let tax: PyRef<Taxonomy> = ctx.unwrap().get_item("tax").unwrap().extract()?;
        assert_eq!(tax.t.name("135622").unwrap(), "New Name");
        drop(tax);

        py.eval(r#"tax.edit_node("135622", rank="genus")"#, None, ctx)?;
        let tax: PyRef<Taxonomy> = ctx.unwrap().get_item("tax").unwrap().extract()?;
        assert_eq!(tax.t.rank("135622").unwrap(), TaxRank::Genus);
        drop(tax);

        assert!(py
            .eval(r#"tax.edit_node("2", parent_id="135622")"#, None, ctx,)
            .is_err());
        assert!(py
            .eval(r#"tax.edit_node("1", parent_id="131567")"#, None, ctx,)
            .is_err());
        assert!(py
            .eval(r#"tax.edit_node("1", parent_distance=3.5)"#, None, ctx,)
            .is_err());
        py.eval(
            r#"tax.edit_node("135622", parent_id="2", parent_distance=3.5)"#,
            None,
            ctx,
        )?;
        let tax: PyRef<Taxonomy> = ctx.unwrap().get_item("tax").unwrap().extract()?;
        assert_eq!(tax.t.parent("135622").unwrap(), Some(("2", 3.5)));
        drop(tax);

        Ok(())
    }

    #[test]
    fn test_iteration() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;

        let n_nodes: u64 = py.eval(r#"len(tax)"#, None, ctx)?.extract()?;
        assert_eq!(n_nodes, 14);

        let n_nodes: u64 = py.eval(r#"len([t for t in tax])"#, None, ctx)?.extract()?;
        assert_eq!(n_nodes, 14);

        let contained: bool = py.eval(r#""1224" in tax"#, None, ctx)?.extract()?;
        assert_eq!(contained, true);
        let contained: bool = py.eval(r#""nonexistant" in tax"#, None, ctx)?.extract()?;
        assert_eq!(contained, false);
        Ok(())
    }

    #[test]
    fn test_import() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;

        let tax: PyRef<Taxonomy> = py
            .eval(r#"Taxonomy.from_newick("(A,B)C;")"#, None, ctx)?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 3);

        let tax: PyRef<Taxonomy> = py
            .eval(
                r#"Taxonomy.from_json("{\"nodes\": [], \"links\": []}")"#,
                None,
                ctx,
            )?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 0);

        let tax: PyRef<Taxonomy> = py
            .eval(r#"Taxonomy.from_json("{\"id\": \"1\"}")"#, None, ctx)?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 1);

        let tax: PyRef<Taxonomy> = py
            .eval(
                r#"Taxonomy.from_json("{\"test\": {\"id\": \"1\"}}", path=["test"])"#,
                None,
                ctx,
            )?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 1);

        let tax: PyRef<Taxonomy> = py
            .eval(r#"Taxonomy.from_phyloxml("<phylogeny rooted=\"true\"><clade><id>root</id></clade></phylogeny>")"#, None, ctx)?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 1);

        let tax: PyRef<Taxonomy> = py
            .eval(r#"Taxonomy.from_ncbi("tests/data/ncbi_subset_tax.nodes.dmp", "tests/data/ncbi_subset_tax.names.dmp")"#, None, ctx)?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 9);

        Ok(())
    }

    #[test]
    fn test_export() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;

        let json: &[u8] = py.eval(r#"tax.to_json()"#, None, ctx)?.extract()?;
        assert_eq!(json[0], b'{');

        let nwk: Vec<u8> = py.eval(r#"tax.to_newick()"#, None, ctx)?.extract()?;
        assert_eq!(nwk[0], b'(');
        Ok(())
    }

    #[test]
    fn test_weights() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;
        ctx.map(|ctx| {
            let w = PyModule::new(py, "weights").unwrap();
            weights(py, w).unwrap();
            ctx.set_item("weights", w).unwrap();
        });

        let rolled: &PyDict = py
            .eval(
                r#"weights.rollup_weights(tax, {"1": 1, "2": 4})"#,
                None,
                ctx,
            )?
            .extract()?;
        assert_eq!(rolled.get_item("1").unwrap().extract::<f32>()?, 5.);
        assert_eq!(rolled.get_item("131567").unwrap().extract::<f32>()?, 4.);

        let mwp: (&str, f32) = py
            .eval(
                r#"weights.maximum_weighted_path(tax, {"1": 1, "2": 4})"#,
                None,
                ctx,
            )?
            .extract()?;
        assert_eq!(mwp.0, "2");
        assert_eq!(mwp.1, 5.);
        Ok(())
    }
}
