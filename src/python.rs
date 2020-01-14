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

use crate::base::{GeneralTaxonomy, IntTaxID};
use crate::edit::{prune_away, prune_to};
use crate::formats::json::{load_json, save_json};
use crate::formats::ncbi::load_ncbi;
use crate::formats::newick::{load_newick, save_newick};
use crate::formats::phyloxml::load_phyloxml;
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy as TaxTrait;
use crate::weights::{maximum_weighted_path as tax_mwp, rollup_weights as tax_rw};

create_exception!(taxonomy, TaxonomyError, pyo3::exceptions::Exception);

/// The Taxonomy object provides the primary interface for exploring a
/// biological taxonomy. Iterating over the Taxonomy returns all the taxonomy
/// ids in pre-order fashion and indexing into the Taxonomy object with a
/// taxonomy id will return a tuple of the name and rank of that id.
#[pyclass]
pub struct Taxonomy {
    t: GeneralTaxonomy,
    visited_nodes: Vec<IntTaxID>,
    nodes_left: Vec<IntTaxID>,
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
        let t = load_json(&mut c, path.as_deref())
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(Taxonomy {
            t,
            nodes_left: Vec::new(),
            visited_nodes: Vec::new(),
        })
    }

    /// from_newick(cls, value: str)
    /// --
    ///
    /// Load a Taxonomy from a Newick-encoded string.
    #[classmethod]
    fn from_newick(_cls: &PyType, value: &str) -> PyResult<Taxonomy> {
        let mut c = Cursor::new(value);
        let t =
            load_newick(&mut c).map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(Taxonomy {
            t,
            nodes_left: Vec::new(),
            visited_nodes: Vec::new(),
        })
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
        let t = load_ncbi(nodes_file, names_file)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(Taxonomy {
            t,
            nodes_left: Vec::new(),
            visited_nodes: Vec::new(),
        })
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
        let t =
            load_phyloxml(&mut c).map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(Taxonomy {
            t,
            nodes_left: Vec::new(),
            visited_nodes: Vec::new(),
        })
    }

    /// to_json(self, /, as_node_link: bool)
    /// --
    ///
    /// Export a Taxonomy as a JSON-encoded byte string. By default, the JSON format
    /// is a tree format unless the `as_node_link` parameter is set to True.
    #[args(as_node_link = false)]
    fn to_json(&self, as_node_link: bool) -> PyResult<PyObject> {
        let mut s = Vec::new();
        save_json::<&str, _, _, _>(&self.t, &mut s, None, as_node_link)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;

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
        save_newick::<&str, _, _, _>(&self.t, &mut s, None)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;

        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyBytes::new(py, &s).into())
    }

    /// parent(self, id: str, /, at_rank: str, include_dist: bool)
    /// --
    ///
    /// Return the immediate parent taxonomy node of the node id provided.
    ///
    /// If `at_rank` is provided, scan all the nodes in the node's lineage and return
    /// the parent id at that rank.
    ///
    /// If `include_dist` is provided, return the distance from the node id provided
    /// to the parent returned.
    #[args(include_dist = false)]
    fn parent(
        &self,
        id: &str,
        at_rank: Option<&str>,
        include_dist: bool,
    ) -> PyResult<Option<PyObject>> {
        let gil = Python::acquire_gil();
        let py = gil.python();

        if let Some(rank) = at_rank {
            if let Ok(rank) = TaxRank::from_str(rank) {
                self.t
                    .parent_at_rank(id, rank)
                    .map(|o| {
                        o.map(|(i, d)| {
                            if include_dist {
                                (i.to_string(), d).to_object(py)
                            } else {
                                i.to_string().to_object(py)
                            }
                        })
                    })
                    .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
            } else {
                Err(PyErr::new::<TaxonomyError, _>(format!(
                    "Rank {} could not be understood",
                    rank
                )))
            }
        } else {
            self.t
                .parent(id)
                .map(|o| {
                    o.map(|(i, d)| {
                        if include_dist {
                            (i.to_string(), d).to_object(py)
                        } else {
                            i.to_string().to_object(py)
                        }
                    })
                })
                .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
        }
    }

    /// children(self, id: str)
    /// --
    ///
    /// Return a list of child taxonomy nodes from the node id provided.
    fn children(&self, id: &str) -> PyResult<Vec<String>> {
        self.t
            .children(id)
            .map(|v| v.into_iter().map(|i| i.to_string()).collect())
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
    }

    /// lineage(self, id: str)
    /// --
    ///
    /// Return a list of all the parent taxonomy nodes of the node id provided
    /// (including that node itself).
    fn lineage(&self, id: &str) -> PyResult<Vec<String>> {
        self.t
            .lineage(id)
            .map(|v| v.into_iter().map(|i| i.to_string()).collect())
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
    }

    /// lca(self, id1: str, id2: str)
    /// --
    ///
    /// Return the lowest common ancestor of two taxonomy nodes.
    fn lca(&self, id1: &str, id2: &str) -> PyResult<String> {
        self.t
            .lca(id1, id2)
            .map(|i| i.to_string())
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
    }

    /// name(self, id: str)
    /// --
    ///
    /// Return the name of the node id provided.
    fn name(&self, id: &str) -> PyResult<String> {
        let name = self
            .t
            .name(id)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(name.to_string())
    }

    /// rank(self, id: str)
    /// --
    ///
    /// Return the taxonomic rank of the node id provided.
    fn rank(&self, id: &str) -> PyResult<String> {
        let rank = self
            .t
            .rank(id)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(rank
            .map(|x| x.to_ncbi_rank())
            .unwrap_or("no rank")
            .to_string())
    }

    /// prune(self, keep: List[str], remove: List[str])
    /// --
    ///
    /// Return a copy of the taxonomy containing:
    ///  - only the nodes in `keep` and their parents if provided
    ///  - all of the nodes except those in remove and their children if provided
    fn prune(&self, keep: Option<Vec<&str>>, remove: Option<Vec<&str>>) -> PyResult<Taxonomy> {
        let mut t = self.t.clone();
        if let Some(k) = keep {
            t = prune_to(&t, &k, false)
                .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        }
        if let Some(r) = remove {
            t = prune_away(&t, &r).map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        }
        Ok(Taxonomy {
            t,
            nodes_left: Vec::new(),
            visited_nodes: Vec::new(),
        })
    }

    /// remove(self, tax_id: str)
    /// --
    ///
    /// Remove the node from the tree.
    fn remove(&mut self, tax_id: &str) -> PyResult<()> {
        let int_id = self
            .t
            .to_internal_id(tax_id)
            .map_err(|_| PyErr::new::<KeyError, _>("Tax ID is not in taxonomy"))?;
        self.t
            .remove(int_id)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(())
    }

    /// add(self, parent_id: str, tax_id: str)
    /// --
    ///
    /// Add a new node to the tree at the parent provided.
    fn add(&mut self, parent_id: &str, tax_id: &str) -> PyResult<()> {
        let int_id = self
            .t
            .to_internal_id(parent_id)
            .map_err(|_| PyErr::new::<KeyError, _>("Parent tax ID is not in taxonomy"))?;
        self.t
            .add(int_id, tax_id)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(())
    }

    #[getter]
    fn get_root(&self) -> PyResult<String> {
        let root: &str = self.t.root();
        Ok(root.to_string())
    }
}

#[pyproto]
impl PyObjectProtocol for Taxonomy {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "<Taxonomy ({} nodes)>",
            TaxTrait::<IntTaxID, _>::len(&self.t)
        ))
    }
}

#[pyproto]
impl PyMappingProtocol for Taxonomy {
    fn __getitem__(&self, key: &str) -> PyResult<(String, Option<String>)> {
        // TODO: return namedtuple?

        // let gil = Python::acquire_gil();
        // let py = gil.python();
        // let namedtuple = py.import("collections")?.get("namedtuple")?;

        let int_id = self
            .t
            .to_internal_id(key)
            .map_err(|_| PyErr::new::<KeyError, _>("Tax ID is not in taxonomy"))?;
        let name = self
            .t
            .name(int_id)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        let rank = self
            .t
            .rank(int_id)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok((name.to_string(), rank.map(|x| x.to_ncbi_rank().to_string())))
    }

    // TODO: way to set name and rank
    // TODO: also way to set parent and parent distance?

    fn __delitem__(&mut self, key: &str) -> PyResult<()> {
        let int_id = self
            .t
            .to_internal_id(key)
            .map_err(|_| PyErr::new::<KeyError, _>("Tax ID is not in taxonomy"))?;
        self.t
            .remove(int_id)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(())
    }

    fn __len__(&self) -> PyResult<usize> {
        Ok(TaxTrait::<IntTaxID, _>::len(&self.t))
    }
}

#[pyproto]
impl PySequenceProtocol for Taxonomy {
    fn __contains__(&self, key: &str) -> PyResult<bool> {
        Ok(self.t.to_internal_id(key).is_ok())
    }
}

#[pyproto]
impl PyIterProtocol for Taxonomy {
    fn __iter__(slf: PyRefMut<Self>) -> PyResult<Self> {
        Ok(Taxonomy {
            t: slf.t.clone(),
            nodes_left: vec![TaxTrait::<IntTaxID, _>::root(&slf.t)],
            visited_nodes: Vec::new(),
        })

        // TODO: this would save a lot of memory if we could just
        // use the GC to make another copy

        // let gil_guard = Python::acquire_gil();
        // let py = gil_guard.python();
        // let inst = pyo3::AsPyRef::as_ref(slf, py);
        // Ok(inst)
    }

    fn __next__(mut slf: PyRefMut<Self>) -> PyResult<Option<String>> {
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
                let children = slf
                    .t
                    .children(cur_node)
                    .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
                if !children.is_empty() {
                    slf.nodes_left.extend(children);
                }
                cur_node // preorder
            };
            if node_visited == !traverse_preorder {
                return Ok(Some(
                    slf.t
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

    let m_weights = wrap_pymodule!(weights)(py);
    m.add("weights", &m_weights)?;
    // this allows `from taxonomy.weights import ...`
    // see https://github.com/PyO3/pyo3/issues/471
    py.import("sys")?
        .dict()
        .get_item("modules")
        .unwrap()
        .downcast_mut::<PyDict>()?
        .set_item("taxonomy.weights", m_weights)?;

    Ok(())
}

#[cfg(test)]
mod test {
    use pyo3::prelude::*;
    use pyo3::types::PyDict;

    use crate::base::test::create_example;
    use crate::taxonomy::Taxonomy as TaxTrait;

    use super::{weights, Taxonomy};

    fn setup_ctx(py: Python) -> PyResult<Option<&PyDict>> {
        let tax = Py::new(
            py,
            Taxonomy {
                t: create_example(),
                visited_nodes: vec![],
                nodes_left: vec![],
            },
        )?;
        let ctx = PyDict::new(py);
        ctx.set_item("tax", tax)?;
        ctx.set_item("Taxonomy", py.get_type::<Taxonomy>())?;
        Ok(Some(ctx))
    }

    #[test]
    fn test_parent() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;

        let parent: String = py.eval(r#"tax.parent("53452")"#, None, ctx)?.extract()?;
        assert_eq!(parent, "1046");

        let parent: String = py
            .eval(r#"tax.parent("53452", at_rank="phylum")"#, None, ctx)?
            .extract()?;
        assert_eq!(parent, "1224");

        assert!(py.eval(r#"tax.parent("bad id")"#, None, ctx).is_err());
        assert!(py
            .eval(r#"tax.parent("bad id", at_rank="species")"#, None, ctx)
            .is_err());
        assert!(py
            .eval(r#"tax.parent("53452", at_rank="bad rank")"#, None, ctx)
            .is_err());
        Ok(())
    }

    #[test]
    fn test_tree_methods() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;

        let children: Vec<String> = py.eval(r#"tax.children("53452")"#, None, ctx)?.extract()?;
        assert_eq!(children, vec!["61598"]);

        let lineage: Vec<String> = py.eval(r#"tax.lineage("1224")"#, None, ctx)?.extract()?;
        assert_eq!(lineage, vec!["1224", "2", "131567", "1"]);

        let lca: String = py
            .eval(r#"tax.lca("56812", "765909")"#, None, ctx)?
            .extract()?;
        assert_eq!(lca, "1236");
        Ok(())
    }

    #[test]
    fn test_info_methods() -> PyResult<()> {
        let gil = Python::acquire_gil();
        let py = gil.python();
        let ctx = setup_ctx(py)?;

        let name: String = py.eval(r#"tax.name("1224")"#, None, ctx)?.extract()?;
        assert_eq!(name, "Proteobacteria");

        let rank: String = py.eval(r#"tax.rank("1224")"#, None, ctx)?.extract()?;
        assert_eq!(rank, "phylum");
        let rank: String = py.eval(r#"tax.rank("1")"#, None, ctx)?.extract()?;
        assert_eq!(rank, "no rank");

        let root: String = py.eval(r#"tax.root"#, None, ctx)?.extract()?;
        assert_eq!(root, "1");

        let rank: (String, String) = py.eval(r#"tax["1224"]"#, None, ctx)?.extract()?;
        assert_eq!(rank, ("Proteobacteria".to_string(), "phylum".to_string()));
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

        let tax: &Taxonomy = ctx.unwrap().get_item("tax").unwrap().extract()?;
        let n_tax_items = TaxTrait::<&str, _>::len(&tax.t);

        assert!(py
            .eval(r#"tax.add("bad id", "new id")"#, None, ctx)
            .is_err());
        py.eval(r#"tax.add("1236", "91347")"#, None, ctx)?;
        let tax: &Taxonomy = ctx.unwrap().get_item("tax").unwrap().extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), n_tax_items + 1);

        assert!(py.eval(r#"tax.remove("bad id")"#, None, ctx).is_err());
        py.eval(r#"tax.remove("135622")"#, None, ctx)?;
        let tax: &Taxonomy = ctx.unwrap().get_item("tax").unwrap().extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), n_tax_items);

        assert!(py.eval(r#"del tax["bad id"])"#, None, ctx).is_err());
        py.run(r#"del tax["1046"]"#, None, ctx)?;
        let tax: &Taxonomy = ctx.unwrap().get_item("tax").unwrap().extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), n_tax_items - 1);

        // FIXME: make these fail?
        // assert!(py.eval(r#"tax.prune(keep=["bad id"])"#, None, ctx).is_err());
        // assert!(py.eval(r#"tax.prune(remove=["bad id"])"#, None, ctx).is_err());
        py.eval(r#"tax.prune(keep=["22"])"#, None, ctx)?;
        py.eval(r#"tax.prune(remove=["22"])"#, None, ctx)?;

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

        let tax: &Taxonomy = py
            .eval(r#"Taxonomy.from_newick("(A,B)C;")"#, None, ctx)?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 3);

        let tax: &Taxonomy = py
            .eval(
                r#"Taxonomy.from_json("{\"nodes\": [], \"links\": []}")"#,
                None,
                ctx,
            )?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 0);

        let tax: &Taxonomy = py
            .eval(r#"Taxonomy.from_json("{\"id\": \"1\"}")"#, None, ctx)?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 1);

        let tax: &Taxonomy = py
            .eval(
                r#"Taxonomy.from_json("{\"test\": {\"id\": \"1\"}}", path=["test"])"#,
                None,
                ctx,
            )?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 1);

        let tax: &Taxonomy = py
            .eval(r#"Taxonomy.from_phyloxml("<phylogeny rooted=\"true\"><clade><id>root</id></clade></phylogeny>")"#, None, ctx)?
            .extract()?;
        assert_eq!(TaxTrait::<&str, _>::len(&tax.t), 1);

        let tax: &Taxonomy = py
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
