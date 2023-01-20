use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;

use serde_json::Value;

use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy;

pub type InternalIndex = usize;

/// The type that is returned when loading any taxonomies through that library.
/// It include 2 implementations of the [Taxonomy] trait: one using strings as ids
/// (easier to use but slower) and one using internal indices (harder to use but faster).
/// Most fields are public for flexibility reason.
#[derive(Clone, Debug, Deserialize, Serialize, PartialEq)]
pub struct GeneralTaxonomy {
    pub tax_ids: Vec<String>,
    pub parent_ids: Vec<InternalIndex>,
    pub parent_distances: Vec<f32>,
    pub names: Vec<String>,
    pub ranks: Vec<TaxRank>,
    // Only used by the JSON format
    pub data: Vec<HashMap<String, Value>>,

    // Lookup tables that can dramatically speed up some operations
    pub(crate) tax_id_lookup: HashMap<String, InternalIndex>,
    pub(crate) children_lookup: Vec<Vec<InternalIndex>>,
}

impl Default for GeneralTaxonomy {
    /// Create a new GeneralTaxonomy with only a root node.
    fn default() -> Self {
        let mut tax = GeneralTaxonomy {
            tax_ids: vec!["1".to_string()],
            parent_ids: vec![0],
            parent_distances: vec![1.],
            ranks: vec![TaxRank::Unspecified],
            names: vec!["root".to_string()],
            data: vec![HashMap::new()],

            tax_id_lookup: HashMap::new(),
            children_lookup: Vec::new(),
        };

        tax.index();
        tax
    }
}

impl GeneralTaxonomy {
    /// Generates indices to speed up searching by tax_id and child
    /// retrieval. We almost always want these (so `new` always
    /// calls this), but it's possible that a user might want to save
    /// memory/start-up time so it's theoretically possible to use this
    /// struct without running this.
    fn index(&mut self) {
        self.tax_id_lookup.clear();
        for (ix, tax_id) in self.tax_ids.iter().enumerate() {
            self.tax_id_lookup.insert(tax_id.clone(), ix);
        }

        for v in self.children_lookup.iter_mut() {
            v.clear();
        }

        if self.children_lookup.len() != self.tax_ids.len() {
            self.children_lookup.resize(self.tax_ids.len(), Vec::new());
        }

        for (ix, parent_ix) in self.parent_ids.iter().enumerate() {
            if ix != 0 {
                self.children_lookup[*parent_ix].push(ix);
            }
        }
    }

    /// Ensures that all nodes go back to a root and that tax ids are unique and raise an error otherwise
    fn validate(&self) -> TaxonomyResult<()> {
        if self.parent_ids.is_empty() {
            return Ok(());
        }

        let mut nodes_ok = HashSet::with_capacity(self.tax_ids.len());
        nodes_ok.insert(self.parent_ids[0]);

        for i in &self.parent_ids {
            let mut evaluated = HashSet::new();
            evaluated.insert(*i);
            let mut current = *i;

            loop {
                let parent = self.parent_ids[current];
                if nodes_ok.contains(&parent) {
                    break;
                }
                if evaluated.contains(&parent) {
                    return Err(Error::new(ErrorKind::InvalidTaxonomy(
                        "Cycle detected in taxonomy".to_owned(),
                    )));
                }

                evaluated.insert(parent);
                current = parent;
            }

            nodes_ok.extend(evaluated);
        }

        Ok(())
    }

    /// Ensures that tax ids are unique.
    /// This is not useful for all kind of taxonomies since some format have optional ids.
    pub fn validate_uniqueness(&self) -> TaxonomyResult<()> {
        if self.tax_ids.is_empty() {
            return Ok(());
        }

        let mut tax_ids = self.tax_ids.clone();
        tax_ids.sort_unstable();
        let mut dupes = HashSet::new();
        let mut last = tax_ids.get(0).unwrap();

        for i in 1..tax_ids.len() {
            let cur = tax_ids.get(i).unwrap();
            if cur == last {
                dupes.insert(cur);
                continue;
            }
            last = cur;
        }

        if !dupes.is_empty() {
            return Err(Error::new(ErrorKind::InvalidTaxonomy(format!(
                "Some tax ids are duplicated: {:?}.",
                dupes.iter().collect::<Vec<_>>()
            ))));
        }

        Ok(())
    }

    #[inline]
    pub fn to_internal_index(&self, tax_id: &str) -> TaxonomyResult<InternalIndex> {
        self.tax_id_lookup.get(tax_id).map_or_else(
            || Err(Error::new(ErrorKind::NoSuchTaxId(tax_id.to_owned()))),
            |t| Ok(*t),
        )
    }

    #[inline]
    pub fn from_internal_index(&self, tax_id: InternalIndex) -> TaxonomyResult<&str> {
        if tax_id >= self.tax_ids.len() {
            return Err(Error::new(ErrorKind::NoSuchTaxId(format!(
                "Internal ID: {}",
                tax_id
            ))));
        }
        Ok(&self.tax_ids[tax_id])
    }

    pub fn from_arrays(
        tax_ids: Vec<String>,
        parent_ids: Vec<InternalIndex>,
        names: Option<Vec<String>>,
        ranks: Option<Vec<TaxRank>>,
        distances: Option<Vec<f32>>,
        data: Option<Vec<HashMap<String, Value>>>,
    ) -> TaxonomyResult<Self> {
        let size = tax_ids.len();
        let adj_names = names.unwrap_or_else(|| vec![String::new(); tax_ids.len()]);
        let adj_ranks = ranks.unwrap_or_else(|| vec![TaxRank::Unspecified; tax_ids.len()]);
        let adj_distances = distances.unwrap_or_else(|| vec![1.0; tax_ids.len()]);
        let adj_data = data.unwrap_or_else(|| vec![HashMap::new(); tax_ids.len()]);

        if size != parent_ids.len() {
            return Err(Error::new(ErrorKind::InvalidTaxonomy(
                "Mismatched number of tax ids and parent ids".to_owned(),
            )));
        }
        if size != adj_names.len() {
            return Err(Error::new(ErrorKind::InvalidTaxonomy(
                "Mismatched number of tax ids and names".to_owned(),
            )));
        }
        if size != adj_ranks.len() {
            return Err(Error::new(ErrorKind::InvalidTaxonomy(
                "Mismatched number of tax ids and names".to_owned(),
            )));
        }
        if size != adj_distances.len() {
            return Err(Error::new(ErrorKind::InvalidTaxonomy(
                "Mismatched number of tax ids and distances".to_owned(),
            )));
        }
        if size != adj_data.len() {
            return Err(Error::new(ErrorKind::InvalidTaxonomy(
                "Mismatched number of tax ids and extra data".to_owned(),
            )));
        }

        let mut tax = GeneralTaxonomy {
            tax_ids,
            parent_ids,
            parent_distances: adj_distances,
            names: adj_names,
            ranks: adj_ranks,
            data: adj_data,

            tax_id_lookup: HashMap::with_capacity(size),
            children_lookup: vec![Vec::new(); size],
        };
        tax.index();
        tax.validate()?;
        Ok(tax)
    }

    /// Retrieves all external IDs given a name
    pub fn find_all_by_name(&self, name: &str) -> Vec<&str> {
        let name_indices = self
            .names
            .iter()
            .enumerate()
            .filter(|(_, val)| val == &name)
            .map(|(pos, _)| pos)
            .collect::<Vec<usize>>();

        name_indices
            .iter()
            .map(|&pos| &self.tax_ids[pos] as &str)
            .collect()
    }

    /// Add a new node to the taxonomy.
    pub fn add(&mut self, parent_id: &str, tax_id: &str) -> TaxonomyResult<()> {
        let parent_idx = self.to_internal_index(parent_id)?;
        let new_idx = self.tax_ids.len();

        self.tax_ids.push(tax_id.to_string());
        self.parent_ids.push(parent_idx);
        self.parent_distances.push(1.0);
        self.ranks.push(TaxRank::Unspecified);
        self.names.push(String::new());
        self.data.push(HashMap::new());

        // update the cached lookup tables
        self.tax_id_lookup.insert(tax_id.to_string(), new_idx);

        if self.children_lookup.len() != self.tax_ids.len() {
            self.children_lookup.resize(self.tax_ids.len(), Vec::new());
        }

        self.children_lookup[parent_idx].push(new_idx);

        Ok(())
    }

    /// Remove a single node from the taxonomy.
    ///
    /// Unlike pruning the children of the node are kept, but are rejoined
    /// onto the node's parent. The root node can not be removed
    pub fn remove(&mut self, tax_id: &str) -> TaxonomyResult<()> {
        let idx = self.to_internal_index(tax_id)?;
        // don't allow deleting root because that messes up the tree structure
        if idx == 0 {
            return Err(Error::new(ErrorKind::OperationNotAllowed(
                "Cannot delete root of taxonomy".to_owned(),
            )));
        }

        // reattach all the child nodes to the parent of the deleted node
        let node_parent = self.parent_ids[idx];
        let parent_dist = self.parent_distances[idx];
        for (parent, dist) in self
            .parent_ids
            .iter_mut()
            .zip(self.parent_distances.iter_mut())
        {
            if *parent == idx {
                *parent = node_parent;
                *dist += parent_dist;
            }
        }

        // and delete the node from all the other tables
        // (note we do this last so we still have the tax id above)
        self.tax_ids.remove(idx);
        self.parent_ids.remove(idx);
        self.parent_distances.remove(idx);
        self.ranks.remove(idx);
        self.names.remove(idx);

        // everything after `tax_id` in parents needs to get decremented by 1
        // because we've changed the actual array size
        for parent in self.parent_ids.iter_mut().skip(idx + 1) {
            if *parent > 0 {
                *parent -= 1;
            }
        }

        // we could try to update the lookups, but all of the values in
        // tax_id_lookup need to be scanned and all of the keys and values
        // in children_lookup have changed so it's easier to just rebuild
        self.index();
        Ok(())
    }
}

/// This is the implementation for &str taxonomy access for a more
/// end-user understandable (but slightly slower) workflow.
impl<'t> Taxonomy<'t, &'t str> for GeneralTaxonomy {
    fn root(&'t self) -> &'t str {
        &self.tax_ids[0]
    }

    fn children(&'t self, tax_id: &str) -> TaxonomyResult<Vec<&'t str>> {
        let idx = self.to_internal_index(tax_id)?;
        self.children_lookup[idx]
            .iter()
            .map(|x| self.from_internal_index(*x))
            .collect()
    }

    fn parent(&'t self, tax_id: &str) -> TaxonomyResult<Option<(&'t str, f32)>> {
        let idx = self.to_internal_index(tax_id)?;
        if idx == 0 {
            return Ok(None);
        }

        Ok(Some((
            self.from_internal_index(self.parent_ids[idx])?,
            self.parent_distances[idx],
        )))
    }

    fn name(&'t self, tax_id: &str) -> TaxonomyResult<&str> {
        let idx = self.to_internal_index(tax_id)?;
        Ok(&self.names[idx])
    }

    fn data(&'t self, tax_id: &str) -> TaxonomyResult<Cow<'t, HashMap<String, Value>>> {
        let idx = self.to_internal_index(tax_id)?;
        Ok(Cow::Borrowed(&self.data[idx]))
    }

    fn rank(&'t self, tax_id: &str) -> TaxonomyResult<TaxRank> {
        let idx = self.to_internal_index(tax_id)?;
        Ok(self.ranks[idx])
    }

    fn len(&'t self) -> usize
    where
        Self: Sized,
    {
        self.tax_ids.len()
    }
}

/// This is the implementation for "internal" tax ID lookup; these IDs are
/// arbitrary (they're positions of the tax nodes in the internal array) and
/// not linked at all to the "external" (e.g. NCBI) IDs. Using these IDs
/// directly can lead to a decent speed up without having to build indices.
/// This is about 3-8x faster than the &str impl, depending on the usage.
impl<'t> Taxonomy<'t, InternalIndex> for GeneralTaxonomy {
    fn root(&'t self) -> InternalIndex {
        0
    }

    fn children(&'t self, tax_id: InternalIndex) -> TaxonomyResult<Vec<InternalIndex>> {
        if let Some(children) = self.children_lookup.get(tax_id) {
            Ok(children.to_vec())
        } else {
            Err(Error::new(ErrorKind::NoSuchInternalIndex(tax_id)))
        }
    }

    fn parent(&'t self, idx: InternalIndex) -> TaxonomyResult<Option<(InternalIndex, f32)>> {
        if idx == 0 {
            return Ok(None);
        }
        if idx >= self.parent_ids.len() {
            return Err(Error::new(ErrorKind::NoSuchInternalIndex(idx)));
        }
        Ok(Some((self.parent_ids[idx], self.parent_distances[idx])))
    }

    fn name(&'t self, idx: InternalIndex) -> TaxonomyResult<&str> {
        if let Some(name) = self.names.get(idx) {
            Ok(name)
        } else {
            Err(Error::new(ErrorKind::NoSuchInternalIndex(idx)))
        }
    }

    fn data(&'t self, idx: InternalIndex) -> TaxonomyResult<Cow<'t, HashMap<String, Value>>> {
        if let Some(data) = self.data.get(idx) {
            Ok(Cow::Borrowed(data))
        } else {
            Err(Error::new(ErrorKind::NoSuchInternalIndex(idx)))
        }
    }

    fn rank(&'t self, idx: InternalIndex) -> TaxonomyResult<TaxRank> {
        if let Some(rank) = self.ranks.get(idx) {
            Ok(*rank)
        } else {
            Err(Error::new(ErrorKind::NoSuchInternalIndex(idx)))
        }
    }

    fn len(&'t self) -> usize
    where
        Self: Sized,
    {
        self.tax_ids.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::formats::json::load;
    use std::io::Cursor;

    fn create_test_taxonomy() -> GeneralTaxonomy {
        let example = r#"{
            "nodes": [
                {"id": "1", "name": "root", "readcount": 1000},
                {"id": "2", "name": "Bacteria", "rank": "no rank", "readcount": 1000},
                {"id": "562", "name": "Escherichia coli", "rank": "species", "readcount": 1000},
                {"id": "1000", "name": "Covid", "rank": "", "readcount": 1000},
                {"id": "101", "name": "Proteus", "rank": "", "readcount": 1000},
                {"id": "102", "name": "Proteus", "rank": "", "readcount": 1000}
            ],
            "links": [
                {"source": 1, "target": 0},
                {"source": 2, "target": 1},
                {"source": 3, "target": 0},
                {"source": 4, "target": 0},
                {"source": 5, "target": 0}
            ]
        }"#;

        load(Cursor::new(example), None).unwrap()
    }

    #[test]
    fn implements_taxonomy_correctly() {
        let tax = create_test_taxonomy();
        assert_eq!(Taxonomy::<&str>::len(&tax), 6);
        assert_eq!(tax.children("1").unwrap(), vec!["2", "1000", "101", "102"]);
        assert_eq!(tax.name("562").unwrap(), "Escherichia coli");
        assert_eq!(tax.rank("562").unwrap(), TaxRank::Species);
        assert_eq!(tax.parent("562").unwrap(), Some(("2", 1.0)));
    }

    #[test]
    fn can_find_all_by_name() {
        let tax = create_test_taxonomy();
        let res = tax.find_all_by_name("Bacteria");
        assert_eq!(res, vec!["2"]);
        let res = tax.find_all_by_name("Proteus");
        assert_eq!(res, vec!["101", "102"]);
    }

    #[test]
    fn can_add_node() {
        let mut tax = create_test_taxonomy();
        let tax_size = Taxonomy::<&str>::len(&tax);
        tax.add("2", "200").unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), tax_size + 1);
        assert_eq!(tax.parent("200").unwrap(), Some(("2", 1.0)));
        assert_eq!(tax.lineage("200").unwrap(), vec!["200", "2", "1"]);

        // parent_id, tax_id
        //
        // make sure we can continue to add nodes aftewards
        //
        let tax_size = Taxonomy::<&str>::len(&tax);
        tax.add("200", "1024").unwrap();

        assert_eq!(Taxonomy::<&str>::len(&tax), tax_size + 1);
        assert_eq!(tax.parent("1024").unwrap(), Some(("200", 1.0)));
        assert_eq!(tax.lineage("1024").unwrap(), vec!["1024", "200", "2", "1"]);
    }

    #[test]
    fn can_remove_node() {
        let mut tax = create_test_taxonomy();
        let tax_size = Taxonomy::<&str>::len(&tax);
        tax.remove("2").unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), tax_size - 1);
        assert_eq!(tax.parent("562").unwrap(), Some(("1", 2.0)));
        assert_eq!(tax.lineage("562").unwrap(), vec!["562", "1"]);
        // can't remove root
        assert!(tax.remove("1").is_err());
    }

    #[test]
    fn errors_on_taxonomy_with_cycle() {
        let example = r#"{
            "nodes": [
                {"id": "1", "name": "root", "readcount": 1000},
                {"id": "2", "name": "Bacteria", "rank": "no rank", "readcount": 1000},
                {"id": "562", "name": "Escherichia coli", "rank": "species", "readcount": 1000},
                {"id": "1000", "name": "Covid", "rank": "", "readcount": 1000},
                {"id": "2000", "name": "Covid-2", "rank": "", "readcount": 1000}
            ],
            "links": [
                {"source": 1, "target": 0},
                {"source": 2, "target": 1},
                {"source": 3, "target": 4},
                {"source": 4, "target": 3}
            ]
        }"#;

        let res = load(Cursor::new(example), None);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().kind,
            ErrorKind::InvalidTaxonomy("Cycle detected in taxonomy".to_owned())
        );
    }

    #[test]
    fn can_validate_uniqueness() {
        let example = r#"{
            "nodes": [
                {"id": "1", "name": "root", "readcount": 1000},
                {"id": "1", "name": "Bacteria", "rank": "no rank", "readcount": 1000}
            ],
            "links": [
                {"source": 1, "target": 0}
            ]
        }"#;

        let res = load(Cursor::new(example), None);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().kind,
            ErrorKind::InvalidTaxonomy("Some tax ids are duplicated: [\"1\"].".to_owned())
        );
    }
}
