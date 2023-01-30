use crate::errors::TaxonomyResult;
use crate::rank::TaxRank;
use serde_json::Value;
use std::borrow::Cow;
use std::collections::{HashMap, VecDeque};
use std::fmt::{Debug, Display};

/// The trait to implement for each form we implement.
/// It does include some default implementations for a few functions as well as a way
/// to iterate through it, starting with the specified node.
pub trait Taxonomy<'t, T: 't>
where
    T: Clone + Debug + Display + PartialEq,
{
    /// Returns the root node of the entire tree.
    fn root(&'t self) -> T;

    /// Returns a [Vec] of all the direct children IDs of the given tax_id.
    fn children(&'t self, tax_id: T) -> TaxonomyResult<Vec<T>>;

    /// Returns a [Vec] of all the children IDs of the given tax_id.
    fn descendants(&'t self, tax_id: T) -> TaxonomyResult<Vec<T>>;

    /// Returns the parent of the given taxonomic node and the distance to said parent.
    /// The parent of the root node will return [None].
    fn parent(&'t self, tax_id: T) -> TaxonomyResult<Option<(T, f32)>>;

    /// Returns a [Vec] of taxonomy nodes from the one provided back to root.
    /// This method must return the node itself as the first entry in the list
    /// and the root node as the last entry in the list.
    fn lineage(&'t self, tax_id: T) -> TaxonomyResult<Vec<T>> {
        let mut parents = Vec::new();
        let mut last_parent = tax_id.clone();
        parents.push(tax_id);
        while let Some(p) = self.parent(last_parent)? {
            last_parent = p.0.clone();
            parents.push(p.0);
        }
        Ok(parents)
    }

    /// Returns the parent at a given taxonomic rank (note this may not be
    /// the immediate parent). This also returns the distance *to* that parent.
    fn parent_at_rank(&'t self, tax_id: T, rank: TaxRank) -> TaxonomyResult<Option<(T, f32)>> {
        // if this node is at the rank, just return it
        if self.rank(tax_id.clone())? == rank {
            return Ok(Some((tax_id, 0.0)));
        }

        // traverse up the tree looking for it
        let mut cur_id = tax_id;
        let mut dists = Vec::new();
        while let Some(p) = self.parent(cur_id)? {
            dists.push(p.1);
            if self.rank(p.0.clone())? == rank {
                return Ok(Some((p.0, dists.into_iter().sum())));
            }
            cur_id = p.0.clone();
        }
        Ok(None)
    }

    /// Returns the first common parent between two nodes. E.g. for the tree:
    ///
    /// #          /-- C -- D
    /// # A -- B --
    /// #          \-- E
    ///
    /// The LCA (lowest common ancestor) of nodes D and E is node B and the
    /// LCA of nodes D and C is node C itself.
    fn lca(&'t self, id1: T, id2: T) -> TaxonomyResult<T> {
        // make a vec of parents of id1
        let mut id1_parents = VecDeque::new();
        id1_parents.push_front(id1);
        while let Some(p) = self.parent(id1_parents.front().unwrap().clone())? {
            id1_parents.push_front(p.0);
        }

        // make a vec of parents of id2
        let mut id2_parents = VecDeque::new();
        id2_parents.push_front(id2);
        while let Some(p) = self.parent(id2_parents.front().unwrap().clone())? {
            id2_parents.push_front(p.0);
        }

        // find the lowest common ancestor
        let mut common = self.root();
        for (pid1, pid2) in id1_parents.into_iter().zip(id2_parents.into_iter()) {
            if pid1 != pid2 {
                break;
            }
            common = pid1;
        }
        Ok(common)
    }

    /// Returns the name of the tax_id provided.
    fn name(&'t self, tax_id: T) -> TaxonomyResult<&str>;

    /// Returns the additional data at the given tax id node
    /// This is only used by the json taxonomy so the value is serde_json::Value
    /// By default it just returns an empty hashmap
    fn data(&'t self, _tax_id: T) -> TaxonomyResult<Cow<'t, HashMap<String, Value>>> {
        Ok(Cow::Owned(HashMap::new()))
    }

    /// Returns the taxonomic rank of the tax_id provided.
    fn rank(&'t self, tax_id: T) -> TaxonomyResult<TaxRank>;

    /// Generates an iterator that traces over the entire taxonomic tree from the given node. During
    /// preorder traversal, it returns `(T, true)` and during postorder traversal
    /// it returns `(T, false)`
    fn traverse(&'t self, node: T) -> TaxonomyResult<TaxonomyIterator<'t, T>>
    where
        Self: Sized,
    {
        Ok(TaxonomyIterator::new(self, node))
    }

    /// Returns the number of nodes in the taxonomy
    fn len(&'t self) -> usize
    where
        Self: Sized,
    {
        self.traverse(self.root()).unwrap().count() / 2
    }

    /// Determines if there are any nodes at all in the taxonomy. This should
    /// almost always be implemented for performance reasons.
    fn is_empty(&'t self) -> bool
    where
        Self: Sized,
    {
        self.len() == 0
    }
}

pub struct TaxonomyIterator<'t, T: 't> {
    nodes_left: Vec<T>,
    visited_nodes: Vec<T>,
    tax: &'t dyn Taxonomy<'t, T>,
}

impl<'t, T> TaxonomyIterator<'t, T> {
    pub fn new(tax: &'t dyn Taxonomy<'t, T>, root_node: T) -> Self {
        TaxonomyIterator {
            nodes_left: vec![root_node],
            visited_nodes: vec![],
            tax,
        }
    }
}

impl<'t, T> Iterator for TaxonomyIterator<'t, T>
where
    T: Clone + Debug + Display + PartialEq,
{
    type Item = (T, bool);

    fn next(&mut self) -> Option<Self::Item> {
        if self.nodes_left.is_empty() {
            return None;
        }

        let cur_node = self.nodes_left.last().unwrap().clone();
        let node_visited = {
            let last_visited = self.visited_nodes.last();
            Some(&cur_node) == last_visited
        };
        if node_visited {
            self.visited_nodes.pop();
            Some((self.nodes_left.pop().unwrap(), false))
        } else {
            self.visited_nodes.push(cur_node.clone());
            let children = self.tax.children(cur_node.clone()).unwrap();
            if !children.is_empty() {
                self.nodes_left.extend(children);
            }
            Some((cur_node, true))
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use std::collections::HashSet;

    pub(crate) struct MockTax;

    impl<'t> Taxonomy<'t, u32> for MockTax {
        fn root(&self) -> u32 {
            1
        }

        fn children(&self, tax_id: u32) -> TaxonomyResult<Vec<u32>> {
            Ok(match tax_id {
                1 => vec![131567],
                131567 => vec![2],
                2 => vec![1224],
                1224 => vec![1236],
                1236 => vec![135622, 135613], // Gammaproteobacteria
                135622 => vec![22],
                22 => vec![62322],
                62322 => vec![56812], // Shewanella
                135613 => vec![1046],
                1046 => vec![53452],
                53452 => vec![61598], // Lamprocystis
                61598 => vec![765909],
                _ => vec![],
            })
        }

        fn descendants(&'t self, tax_id: u32) -> TaxonomyResult<Vec<u32>> {
            let children: HashSet<u32> = self
                .traverse(tax_id)?
                .map(|(n, _)| n)
                .filter(|n| *n != tax_id)
                .collect();
            let mut children: Vec<u32> = children.into_iter().collect();
            children.sort_unstable();
            Ok(children)
        }

        fn parent(&self, tax_id: u32) -> TaxonomyResult<Option<(u32, f32)>> {
            Ok(match tax_id {
                131567 => Some((1, 1.)),
                2 => Some((131567, 1.)),
                1224 => Some((2, 1.)),
                1236 => Some((1224, 1.)), // Gammaproteobacteria
                135622 => Some((1236, 1.)),
                22 => Some((135622, 1.)),
                62322 => Some((22, 1.)),
                56812 => Some((62322, 1.)), // Shewanella frigidimarina
                135613 => Some((1236, 1.)),
                1046 => Some((135613, 1.)),
                53452 => Some((1046, 1.)),
                61598 => Some((53452, 1.)),  // Lamprocystis purpurea
                765909 => Some((61598, 1.)), // Lamprocystis purpurea
                _ => None,
            })
        }

        fn name(&self, tax_id: u32) -> TaxonomyResult<&str> {
            Ok(match tax_id {
                1 => "root",
                131567 => "cellular organisms",
                2 => "Bacteria",
                1224 => "Proteobacteria",
                1236 => "Gammaproteobacteria",
                135613 => "Chromatiales",
                1046 => "Chromatiaceae",
                53452 => "Lamprocystis",
                61598 => "Lamprocystis purpurea",
                765909 => "Lamprocystis purpurea DSM 4197",
                _ => "",
            })
        }

        fn rank(&self, tax_id: u32) -> TaxonomyResult<TaxRank> {
            Ok(match tax_id {
                2 => TaxRank::Superkingdom,
                1224 => TaxRank::Phylum,
                1236 => TaxRank::Class,
                135613 => TaxRank::Order,
                1046 => TaxRank::Family,
                53452 => TaxRank::Genus,
                61598 => TaxRank::Species,
                _ => TaxRank::Unspecified,
            })
        }
    }

    #[test]
    fn test_len() {
        let tax = MockTax;
        assert_eq!(tax.root(), 1);
        assert_eq!(tax.len(), 14);
        assert_eq!(tax.is_empty(), false);
    }

    #[test]
    fn test_descendants() {
        let tax = MockTax;
        assert_eq!(
            tax.descendants(2).unwrap(),
            vec![22, 1046, 1224, 1236, 53452, 56812, 61598, 62322, 135613, 135622, 765909]
        );
    }

    #[test]
    fn test_lca() {
        let tax = MockTax;
        assert_eq!(tax.lca(56812, 22).unwrap(), 22);
        assert_eq!(tax.lca(56812, 765909).unwrap(), 1236);
    }

    #[test]
    fn test_lineage() {
        let tax = MockTax;
        assert_eq!(tax.lineage(1).unwrap(), vec![1]);
        assert_eq!(
            tax.lineage(61598).unwrap(),
            vec![61598, 53452, 1046, 135613, 1236, 1224, 2, 131567, 1]
        );
    }

    #[test]
    fn test_parent_at_rank() {
        let tax = MockTax;

        assert_eq!(
            tax.parent_at_rank(765909, TaxRank::Genus).unwrap(),
            Some((53452, 2.))
        );
        assert_eq!(
            tax.parent_at_rank(765909, TaxRank::Class).unwrap(),
            Some((1236, 5.))
        );
        assert_eq!(
            tax.parent_at_rank(1224, TaxRank::Phylum).unwrap(),
            Some((1224, 0.))
        );
        assert_eq!(tax.parent_at_rank(1224, TaxRank::Genus).unwrap(), None,);
    }

    #[test]
    fn test_traversal() {
        let tax = MockTax;
        let mut visited = HashSet::new();
        let n_nodes = tax
            .traverse(tax.root())
            .unwrap()
            .enumerate()
            .map(|(ix, (tid, pre))| match ix {
                0 => {
                    assert_eq!(tid, 1, "root is first");
                    assert_eq!(pre, true, "preorder visits happen first");
                }
                27 => {
                    assert_eq!(tid, 1, "root is last too");
                    assert_eq!(pre, false, "postorder visits happen last");
                }
                _ => {
                    if pre {
                        visited.insert(tid);
                    } else {
                        assert!(visited.contains(&tid));
                    }
                }
            })
            .count();

        assert_eq!(n_nodes, 28, "Each node appears twice");
    }
}
