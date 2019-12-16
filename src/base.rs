//! The base implementation of the Taxonomy
use std::cmp::PartialOrd;
use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::iter::Sum;

use serde_derive::{Deserialize, Serialize};

use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy;
use crate::{Result, TaxonomyError};

pub type IntTaxID = usize;

/// A concrete implementation of the Taxonomy trait suitable for interconversions
/// of different formats.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct GeneralTaxonomy {
    pub(crate) tax_ids: Vec<String>,
    pub(crate) parent_ids: Vec<IntTaxID>,
    pub(crate) parent_dists: Vec<f32>,
    pub(crate) ranks: Vec<Option<TaxRank>>,
    pub(crate) names: Vec<String>,

    // these are optional lookup tables to speed up some operations
    pub(crate) tax_id_lookup: Option<HashMap<String, IntTaxID>>,
    pub(crate) children_lookup: Option<Vec<Vec<IntTaxID>>>,
}

impl GeneralTaxonomy {
    /// Generates indices to speed up searching by tax_id and child
    /// retrieval. We almost always want these (so `new` always
    /// calls this), but it's possible that a user might want to save
    /// memory/start-up time so it's theoretically possible to use this
    /// struct without running this.
    fn index(mut self, tax_id: bool, children: bool) -> Self {
        if tax_id && self.tax_id_lookup == None {
            let mut tax_id_lookup = HashMap::with_capacity(self.tax_ids.len());
            for (ix, tax_id) in self.tax_ids.iter().enumerate() {
                tax_id_lookup.insert(tax_id.clone(), ix);
            }

            self.tax_id_lookup = Some(tax_id_lookup);
        }
        if children && self.children_lookup == None {
            let mut children_lookup = vec![Vec::new(); self.tax_ids.len()];
            for (ix, parent_ix) in self.parent_ids.iter().enumerate() {
                if ix != 0 {
                    children_lookup[*parent_ix].push(ix);
                }
            }

            self.children_lookup = Some(children_lookup);
        }
        self
    }

    /// Initializer for a new Taxonomy.
    ///
    /// All `Vec`s must be the same length or initialization will fail.
    pub fn new(
        tax_ids: Vec<String>,
        parent_ids: Vec<usize>,
        names: Option<Vec<String>>,
        ranks: Option<Vec<Option<TaxRank>>>,
        dists: Option<Vec<f32>>,
    ) -> Self {
        let adj_names = names.unwrap_or_else(|| vec!["".to_string(); tax_ids.len()]);
        let adj_ranks = ranks.unwrap_or_else(|| vec![None; tax_ids.len()]);
        let adj_dists = dists.unwrap_or_else(|| vec![1.; tax_ids.len()]);
        assert_eq!(
            tax_ids.len(),
            parent_ids.len(),
            "each taxa must have only one parent"
        );
        assert_eq!(
            tax_ids.len(),
            adj_names.len(),
            "each taxa must have only one name"
        );
        assert_eq!(
            tax_ids.len(),
            adj_ranks.len(),
            "each taxa must have only one rank"
        );
        assert_eq!(
            tax_ids.len(),
            adj_dists.len(),
            "each taxa must have only one parent distance"
        );

        GeneralTaxonomy {
            tax_ids,
            parent_ids,
            parent_dists: adj_dists,
            ranks: adj_ranks,
            names: adj_names,
            tax_id_lookup: None,
            children_lookup: None,
        }
        .index(true, true)
    }

    pub fn from_taxonomy<'t, T: 't, D: 't>(taxonomy: &'t impl Taxonomy<'t, T, D>) -> Result<Self>
    where
        T: Clone + Debug + Display + PartialEq + Ord,
        D: Clone + Debug + Ord + Into<f64> + PartialOrd + PartialEq + Sum,
    {
        let tax_iter = taxonomy.traverse(taxonomy.root())?;
        let mut stack = Vec::new();

        let mut tax_ids = Vec::new();
        let mut parent_ids = Vec::new();
        let mut parent_dists = Vec::new();
        let mut ranks = Vec::new();
        let mut names = Vec::new();

        let mut ix = 0;
        for (tax_id, pre) in tax_iter {
            if !pre {
                stack.pop();
                continue;
            }
            tax_ids.push(tax_id.to_string());

            names.push(taxonomy.name(tax_id.clone())?.to_string());
            ranks.push(taxonomy.rank(tax_id.clone())?);

            if let Some((_, dist)) = taxonomy.parent(tax_id.clone())? {
                parent_dists.push(dist.into() as f32);
                parent_ids.push(*stack.last().unwrap());
            } else if tax_id == taxonomy.root() {
                parent_dists.push(0.);
                parent_ids.push(0);
            } else {
                return Err(TaxonomyError::MalformedTree {
                    tax_id: tax_id.to_string(),
                }
                .into());
            }
            stack.push(ix);
            ix += 1;
        }
        Ok(GeneralTaxonomy {
            tax_ids,
            parent_ids,
            parent_dists,
            ranks,
            names,
            tax_id_lookup: None,
            children_lookup: None,
        }
        .index(true, true))
    }

    /// Given an internal tax_id (an array position) return the
    /// corresponding external tax_id (e.g. a NCBI taxonomy ID).
    ///
    /// Because the `Taxonomy` train is implemented for internal IDs also,
    /// the can be used to speed up some operations by avoiding a string lookup.
    #[allow(clippy::wrong_self_convention)]
    #[inline]
    pub fn from_internal_id(&self, tax_id: IntTaxID) -> Result<&str> {
        // TODO: index so we'll get an Err
        Ok(&self.tax_ids[tax_id as usize])
    }

    /// Given an external tax_id (e.g. a NCBI taxonomy ID) return the
    /// corresponding internal tax_id (the position of that tax node in the
    /// internal array).
    ///
    /// Because the `Taxonomy` train is implemented for internal IDs also,
    /// the can be used to speed up some operations by avoiding a string lookup.
    #[inline]
    pub fn to_internal_id(&self, tax_id: &str) -> Result<IntTaxID> {
        match self.tax_id_lookup {
            Some(ref lookup) => lookup.get(tax_id).map_or_else(
                || {
                    Err(TaxonomyError::NoSuchKey {
                        key: tax_id.to_string(),
                    }
                    .into())
                },
                |t| Ok(*t),
            ),
            None => {
                for (ix, k) in self.tax_ids.iter().enumerate() {
                    if k == tax_id {
                        return Ok(ix as IntTaxID);
                    }
                }
                Err(TaxonomyError::NoSuchKey {
                    key: tax_id.to_string(),
                }
                .into())
            }
        }
    }

    pub fn drop(&mut self, tax_id: IntTaxID) -> Result<()> {
        // don't allow deleting root because that messes up the tree structure
        if tax_id == 0 {
            return Err(TaxonomyError::MalformedTree { tax_id: self.tax_ids[0].clone() }.into());
        }

        // reattach all the child nodes to the parent of the deleted node
        let node_parent = self.parent_ids[tax_id as usize];
        let parent_dist = self.parent_dists[tax_id as usize];
        let mut children = Vec::new();
        for (ix, (parent, dist)) in self
            .parent_ids
            .iter_mut()
            .zip(self.parent_dists.iter_mut())
            .enumerate()
        {
            if *parent == tax_id {
                *parent = node_parent;
                *dist += parent_dist;
                children.push(ix);
            }
        }

        // update the cached lookup tables
        if let Some(ref mut tax_id_lookup) = self.tax_id_lookup {
            // note: we do this before changing the table itself just so we
            // can do this lookup here
            tax_id_lookup.remove(&self.tax_ids[tax_id as usize]);
            for value in tax_id_lookup.values_mut() {
                if *value > tax_id {
                    *value -= 1;
                }
            }
        };
        if let Some(ref mut children_lookup) = self.children_lookup {
            children_lookup[node_parent].extend(children);
            children_lookup.remove(tax_id);
        };

        // and delete the node from all the other tables
        // (note we do this last so we still have the tax id above)
        self.tax_ids.remove(tax_id);
        self.parent_ids.remove(tax_id);
        self.parent_dists.remove(tax_id);
        self.ranks.remove(tax_id);
        self.names.remove(tax_id);

        // everything after `tax_id` in parents needs to get decremented by 1
        // because we've changed the actual array size
        for parent in self.parent_ids.iter_mut().skip(tax_id + 1) {
            *parent -= 1;
        }

        Ok(())
    }
}

/// This is the implementation for "internal" tax ID lookup; these IDs are
/// arbitrary (they're positions of the tax nodes in the internal array) and
/// not linked at all to the "external" (e.g. NCBI) IDs. Using these IDs
/// directly can lead to a decent speed up without having to build indicies.
impl<'s> Taxonomy<'s, IntTaxID, f32> for GeneralTaxonomy {
    fn root(&self) -> IntTaxID {
        0
    }

    fn children(&self, tax_id: IntTaxID) -> Result<Vec<IntTaxID>> {
        if let Some(ref child_lookup) = self.children_lookup {
            // O(1) implementation
            if tax_id as usize >= self.parent_ids.len() {
                return Err(TaxonomyError::NoSuchKey {
                    key: tax_id.to_string(),
                }
                .into());
            }
            Ok(child_lookup[tax_id].to_vec())
        } else {
            // O(n) implementation (slow!) -> we don't actually need
            // to use this for any of the operations we do though!
            let usize_tax_id = tax_id as usize;
            let mut children = Vec::new();
            for (i, t) in self.parent_ids.iter().enumerate() {
                if t == &usize_tax_id {
                    children.push(i as IntTaxID)
                }
            }
            Ok(children)
        }
    }

    fn parent(&self, tax_id: IntTaxID) -> Result<Option<(IntTaxID, f32)>> {
        // O(1) implementation
        if tax_id as usize >= self.parent_ids.len() {
            return Err(TaxonomyError::NoSuchKey {
                key: tax_id.to_string(),
            }
            .into());
        } else if tax_id == 0 {
            return Ok(None);
        }
        Ok(Some((
            self.parent_ids[tax_id as usize] as IntTaxID,
            self.parent_dists[tax_id as usize],
        )))
    }

    fn len(&self) -> usize {
        self.tax_ids.len()
    }

    fn name(&self, tax_id: usize) -> Result<&str> {
        if tax_id >= self.parent_ids.len() {
            return Err(TaxonomyError::NoSuchKey {
                key: tax_id.to_string(),
            }
            .into());
        }
        Ok(&self.names[tax_id])
    }

    fn rank(&self, tax_id: usize) -> Result<Option<TaxRank>> {
        if tax_id >= self.parent_ids.len() {
            return Err(TaxonomyError::NoSuchKey {
                key: tax_id.to_string(),
            }
            .into());
        }
        Ok(self.ranks[tax_id])
    }
}

/// This is the implementation for &str taxonomy access for a more
/// end-user understandable (but slightly slower) workflow.
impl<'s> Taxonomy<'s, &'s str, f32> for GeneralTaxonomy {
    fn root(&'s self) -> &'s str {
        self.from_internal_id(0).unwrap()
    }

    fn children(&self, tax_id: &str) -> Result<Vec<&str>> {
        if let Some(ref child_lookup) = self.children_lookup {
            // O(1) implementation
            let usize_tax_id = self.to_internal_id(&tax_id)?;
            child_lookup[usize_tax_id]
                .iter()
                .map(|x| self.from_internal_id(*x))
                .collect()
        } else {
            // O(n) implementation (slow!)
            let mut children = Vec::new();
            let usize_tax_id = self.to_internal_id(&tax_id)?;
            for (i, t) in self.parent_ids.iter().enumerate() {
                if t == &usize_tax_id {
                    children.push(self.from_internal_id(i as IntTaxID)?)
                }
            }
            Ok(children)
        }
    }

    fn parent(&self, tax_id: &str) -> Result<Option<(&str, f32)>> {
        // O(1) implementation
        let usize_tax_id = self.to_internal_id(&tax_id)?;
        if usize_tax_id == 0 {
            return Ok(None);
        }
        Ok(Some((
            self.from_internal_id(self.parent_ids[usize_tax_id])?,
            self.parent_dists[usize_tax_id],
        )))
    }

    fn len(&self) -> usize {
        self.tax_ids.len()
    }

    fn name(&self, tax_id: &str) -> Result<&str> {
        let tax_id = self.to_internal_id(&tax_id)?;
        Ok(&self.names[tax_id])
    }

    fn rank(&self, tax_id: &str) -> Result<Option<TaxRank>> {
        let tax_id = self.to_internal_id(&tax_id)?;
        Ok(self.ranks[tax_id])
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::taxonomy::test::MockTax;

    fn create_example() -> GeneralTaxonomy {
        GeneralTaxonomy::from_taxonomy(&MockTax).unwrap()
    }

    #[test]
    fn test_taxonomy_conversion() {
        let mock_tax = MockTax;
        let general_tax = create_example();
        assert_eq!(mock_tax.len(), Taxonomy::<&str, f32>::len(&general_tax));
    }

    #[test]
    fn test_drop() {
        let mut tax = create_example();
        assert_eq!(tax.parent("1224").unwrap(), Some(("2", 1.)));
        tax.drop(tax.to_internal_id("2").unwrap()).unwrap();
        assert_eq!(tax.parent("1224").unwrap(), Some(("131567", 2.)));
        assert!(tax.drop(0).is_err());
    }
}
