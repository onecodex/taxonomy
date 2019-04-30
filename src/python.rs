//! Implementation of the Python API.
//!
//! Only enabled when `cargo build --features python` is run or when the
//! Python API is build with `python setup.py develop`.
use std::collections::HashMap;
use std::fs::File;
use std::io::Cursor;
use std::str::FromStr;

use pyo3::class::*;
use pyo3::prelude::*;
use pyo3::types::exceptions::KeyError;
use pyo3::types::{PyBytes, PyType};
use pyo3::AsPyRef;

use crate::base::{GeneralTaxonomy, IntTaxID};
use crate::edit::{prune_away, prune_to};
use crate::formats::json::{load_json, save_json};
use crate::formats::ncbi::load_ncbi;
use crate::formats::newick::{load_newick, save_newick};
use crate::formats::phyloxml::load_phyloxml;
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy as TaxTrait;

py_exception!(taxonomy, TaxonomyError, pyo3::exceptions::Exception);

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
    /// from_json(cls, value)
    /// --
    ///
    /// Load a Taxonomy from a JSON-encoded string. The format can either be
    /// of the tree or node_link_data types.
    #[classmethod]
    fn from_json(_cls: &PyType, value: &str) -> PyResult<Taxonomy> {
        let mut c = Cursor::new(value);
        let t = load_json(&mut c, None)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok(Taxonomy {
            t,
            nodes_left: Vec::new(),
            visited_nodes: Vec::new(),
        })
    }

    /// from_newick(cls, value)
    /// --
    ///
    /// Load a Taxonomy from a Newick-encoded string.
    ///
    /// Experimental.
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

    /// from_ncbi(cls, nodes_path, names_path)
    /// --
    ///
    /// Load a Taxonomy from a pair of NCBI dump files (nodes.dmp and names.dmp).
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

    /// from_phyloxml(cls, value)
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

    /// to_json(self, /, as_node_link)
    /// --
    ///
    /// Export a Taxonomy as a JSON-encoded byte string. By default, the JSON format
    /// is a tree format unless the `as_node_link` parameter is set to True.
    #[args(as_node_link = false)]
    fn to_json(&self, as_node_link: bool) -> PyResult<Py<PyBytes>> {
        let mut s = Vec::new();
        save_json::<&str, _, _, _>(&self.t, &mut s, None, as_node_link)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;

        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyBytes::new(py, &s))
    }

    /// to_newick(self)
    /// --
    ///
    /// Export a Taxonomy as a Newick-encoded byte string.
    ///
    /// Experimental.
    fn to_newick(&self) -> PyResult<Py<PyBytes>> {
        let mut s = Vec::new();
        save_newick::<&str, _, _, _>(&self.t, &mut s, None)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;

        let gil = Python::acquire_gil();
        let py = gil.python();
        Ok(PyBytes::new(py, &s))
    }

    /// parent(self, id, /, at_rank, include_dist)
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
            if let Some(rank) = TaxRank::from_str(rank).ok() {
                self.t
                    .parent_at_rank(id, rank)
                    .map(|o| {
                        o.map(|(i, d)| {
                            if include_dist {
                                (i.to_string(), d).into_object(py)
                            } else {
                                i.to_string().into_object(py)
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
                            (i.to_string(), d).into_object(py)
                        } else {
                            i.to_string().into_object(py)
                        }
                    })
                })
                .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
        }
    }

    /// children(self, id)
    /// --
    ///
    /// Return a list of child taxonomy nodes from the node id provided.
    fn children(&self, id: &str) -> PyResult<Vec<String>> {
        self.t
            .children(id)
            .map(|v| v.iter().map(|i| i.to_string()).collect())
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
    }

    /// lineage(self, id)
    /// --
    ///
    /// Return a list of all the parent taxonomy nodes of the node id provided
    /// (including that node itself).
    fn lineage(&self, id: &str) -> PyResult<Vec<String>> {
        self.t
            .lineage(id)
            .map(|v| v.iter().map(|i| i.to_string()).collect())
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
    }

    /// lca(self, id1, id2)
    /// --
    ///
    /// Return the lowest common ancestor of two taxonomy nodes.
    fn lca(&self, id1: &str, id2: &str) -> PyResult<String> {
        self.t
            .lca(id1, id2)
            .map(|i| i.to_string())
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
    }

    /// name(self, id)
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

    /// rank(self, id)
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

    /// prune(self, keep, remove)
    /// --
    ///
    /// Return a copy of the taxonomy containing:
    ///  - only the nodes in `keep` and their parents if provided
    ///  - all of the nodes except those in remove and their children if provided
    fn prune(&self, keep: Option<Vec<&str>>, remove: Option<Vec<&str>>) -> PyResult<Taxonomy> {
        let mut t = self.t.clone();
        if let Some(k) = keep {
            let keep_ids: Result<Vec<_>, _> = k
                .iter()
                .map(|i| {
                    self.t
                        .to_internal_id(i)
                        .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
                })
                .collect();
            t = prune_to(&t, &keep_ids?, false)
                .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        }
        if let Some(r) = remove {
            let remove_ids: Result<Vec<_>, _> = r
                .iter()
                .map(|i| {
                    self.t
                        .to_internal_id(i)
                        .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))
                })
                .collect();
            t = prune_away(&t, &remove_ids?)
                .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        }
        Ok(Taxonomy {
            t,
            nodes_left: Vec::new(),
            visited_nodes: Vec::new(),
        })
    }

    fn maximum_weighted_path(&self, weights: &PyObject) -> PyResult<(String, f64)> {
        let weights: &HashMap<&str, f32> = PyObjectRef::extract(weights)?;
        Ok(("".to_string(), 0.))
    }

    #[getter]
    fn get_root(&self) -> PyResult<String> {
        let root: IntTaxID = self.t.root();
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
        // TODO: use KeyError instead of TaxonomyError for "no key found" situations
        let name = self
            .t
            .name(key)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        let rank = self
            .t
            .rank(key)
            .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
        Ok((name.to_string(), rank.map(|x| x.to_ncbi_rank().to_string())))
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
    fn __iter__(&mut self) -> PyResult<Self> {
        Ok(Taxonomy {
            t: self.t.clone(),
            nodes_left: vec![TaxTrait::<IntTaxID, _>::root(&self.t)],
            visited_nodes: Vec::new(),
        })

        // TODO: this would save a lot of memory if we could just
        // use the GC to make another copy

        // let gil_guard = Python::acquire_gil();
        // let py = gil_guard.python();
        // let inst = pyo3::AsPyRef::as_ref(self, py);
        // Ok(inst)
    }

    fn __next__(&mut self) -> PyResult<Option<String>> {
        let traverse_preorder = true;
        loop {
            if self.nodes_left.is_empty() {
                return Ok(None);
            }

            let cur_node = *self.nodes_left.last().unwrap();
            let node_visited = {
                let last_visited = self.visited_nodes.last();
                Some(&cur_node) == last_visited
            };
            let node = if node_visited {
                self.visited_nodes.pop();
                self.nodes_left.pop().unwrap() // postorder
            } else {
                self.visited_nodes.push(cur_node.clone());
                let children = self
                    .t
                    .children(cur_node)
                    .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?;
                if !children.is_empty() {
                    self.nodes_left.extend(children);
                }
                cur_node // preorder
            };
            if node_visited == !traverse_preorder {
                return Ok(Some(
                    self.t
                        .from_internal_id(node)
                        .map_err(|e| PyErr::new::<TaxonomyError, _>(format!("{}", e)))?
                        .to_string(),
                ));
            }
        }
    }
}

#[pymodinit]
fn taxonomy(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Taxonomy>()?;

    Ok(())
}
