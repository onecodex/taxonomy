//! The Taxonomy trait defines a large suite of methods that can be used by
//! implementing a limited subset of required methods.
use std::collections::VecDeque;
use std::fmt::{Debug, Display};
use std::iter::Sum;

use crate::rank::TaxRank;
use crate::Result;

/// Taxonomy is a trait exposing a number of traversal and informational
/// methods given a smaller number of required methods. This allows
/// different taxonomic data structures with different memory/speed
/// requirements to have a common set of methods.
/// `T` is the type for the taxonomy ID and `D` for the distance
pub trait Taxonomy<'t, T: 't, D: 't>
where
    T: Clone + Debug + Display + PartialEq,
    D: Debug + PartialOrd + PartialEq + Sum,
{
    /// Returns the root node of the entire tree.
    fn root(&'t self) -> T;

    /// Returns a [Vec] of all the child IDs of the given tax_id.
    fn children(&'t self, tax_id: T) -> Result<Vec<T>>;

    /// Returns the parent of the given taxonomic node and the distance to said parent.
    /// The parent of the root node should always be [None].
    fn parent(&'t self, tax_id: T) -> Result<Option<(T, D)>>; // None if root or not in taxonomy

    /// Returns the name of the tax_id provided.
    ///
    /// Although not strictly required for many taxonomic operations,
    /// this method allows taxonomic exports to conveniently have access
    /// to names in a standardized fashion.
    fn name(&self, tax_id: T) -> Result<&str>;

    /// Returns the taxonomic rank of the tax_id provided.
    ///
    /// Although not strictly required for many of the operations we implement
    /// here, this is primarily here to allow taxonomic exports to conveniently
    /// have access to ranks in a standardized fashion.
    fn rank(&self, tax_id: T) -> Result<Option<TaxRank>>;

    /// Returns a [Vec] of taxonomy nodes from the one provided back to root.
    /// This method must return the node itself as the first entry in the list
    /// and the root node as the last entry in the list.
    fn lineage(&'t self, tax_id: T) -> Result<Vec<T>> {
        // make a vec of parents of id1
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
    fn parent_at_rank(&'t self, tax_id: T, rank: TaxRank) -> Result<Option<(T, D)>> {
        // if this node is at the rank, just return it
        if let Some(cur_rank) = self.rank(tax_id.clone())? {
            if cur_rank == rank {
                return Ok(Some((tax_id, vec![].into_iter().sum())));
            }
        }

        // traverse up the tree looking for it
        let mut cur_id = tax_id;
        let mut dists = Vec::new();
        while let Some(p) = self.parent(cur_id)? {
            dists.push(p.1);
            if let Some(cur_rank) = self.rank(p.0.clone())? {
                if cur_rank == rank {
                    return Ok(Some((p.0, dists.into_iter().sum())));
                }
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
    fn lca(&'t self, id1: T, id2: T) -> Result<T> {
        // make a vec of parents of id1
        let mut id1_parents: VecDeque<T> = VecDeque::new();
        id1_parents.push_front(id1);
        while let Some(p) = self.parent(id1_parents.front().unwrap().clone())? {
            id1_parents.push_front(p.0);
        }

        // make a vec of parents of id2
        let mut id2_parents: VecDeque<T> = VecDeque::new();
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

    /// Generates an iterator that traces over the entire taxonomic tree. During
    /// preorder traversal, it returns `(T, true)` and during postorder traversal
    /// it returns `(T, false)`
    fn traverse(&'t self, node: T) -> Result<TaxonomyIterator<'t, T, D>>
    where
        Self: Sized,
    {
        Ok(TaxonomyIterator::new(self, node))
    }

    /// Returns the number of nodes in the taxonomy
    ///
    /// A default implementation is specified, but you almost always want
    /// to provide your own for speed reasons.
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

pub struct TaxonomyIterator<'t, T: 't, D: 't> {
    nodes_left: Vec<T>,
    visited_nodes: Vec<T>,
    tax: &'t dyn Taxonomy<'t, T, D>,
}

impl<'t, T, D> TaxonomyIterator<'t, T, D> {
    pub fn new(tax: &'t dyn Taxonomy<'t, T, D>, root_node: T) -> Self {
        TaxonomyIterator {
            nodes_left: vec![root_node],
            visited_nodes: vec![],
            tax,
        }
    }
}

impl<'t, T, D> Iterator for TaxonomyIterator<'t, T, D>
where
    T: Clone + Debug + Display + PartialEq,
    D: Debug + PartialOrd + PartialEq + Sum,
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
pub(crate) mod test {
    use std::collections::HashSet;

    use crate::rank::TaxRank;
    use crate::taxonomy::Taxonomy;
    use crate::Result;

    pub(crate) struct MockTax;

    impl<'t> Taxonomy<'t, u32, u16> for MockTax {
        fn root(&self) -> u32 {
            1
        }

        fn children(&self, tax_id: u32) -> Result<Vec<u32>> {
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

        fn parent(&self, tax_id: u32) -> Result<Option<(u32, u16)>> {
            Ok(match tax_id {
                131567 => Some((1, 1)),
                2 => Some((131567, 1)),
                1224 => Some((2, 1)),
                1236 => Some((1224, 1)), // Gammaproteobacteria
                135622 => Some((1236, 1)),
                22 => Some((135622, 1)),
                62322 => Some((22, 1)),
                56812 => Some((62322, 1)), // Shewanella frigidimarina
                135613 => Some((1236, 1)),
                1046 => Some((135613, 1)),
                53452 => Some((1046, 1)),
                61598 => Some((53452, 1)),  // Lamprocystis purpurea
                765909 => Some((61598, 1)), // Lamprocystis purpurea
                _ => None,
            })
        }

        fn name(&self, tax_id: u32) -> Result<&str> {
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

        fn rank(&self, tax_id: u32) -> Result<Option<TaxRank>> {
            Ok(match tax_id {
                2 => Some(TaxRank::Superkingdom),
                1224 => Some(TaxRank::Phylum),
                1236 => Some(TaxRank::Class),
                135613 => Some(TaxRank::Order),
                1046 => Some(TaxRank::Family),
                53452 => Some(TaxRank::Genus),
                61598 => Some(TaxRank::Species),
                _ => None,
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
    fn test_lca() -> Result<()> {
        let tax = MockTax;
        assert_eq!(tax.lca(56812, 22)?, 22);
        assert_eq!(tax.lca(56812, 765909)?, 1236);
        Ok(())
    }

    #[test]
    fn test_lineage() -> Result<()> {
        let tax = MockTax;
        assert_eq!(tax.lineage(1)?, vec![1]);
        assert_eq!(
            tax.lineage(61598)?,
            vec![61598, 53452, 1046, 135613, 1236, 1224, 2, 131567, 1]
        );
        Ok(())
    }

    #[test]
    fn test_parent_at_rank() -> Result<()> {
        let tax = MockTax;

        assert_eq!(
            tax.parent_at_rank(765909, TaxRank::Genus)?,
            Some((53452, 2))
        );
        assert_eq!(tax.parent_at_rank(765909, TaxRank::Class)?, Some((1236, 5)));
        assert_eq!(tax.parent_at_rank(1224, TaxRank::Phylum)?, Some((1224, 0)));
        assert_eq!(tax.parent_at_rank(1224, TaxRank::Genus)?, None,);
        Ok(())
    }

    #[test]
    fn test_traversal() -> Result<()> {
        let tax = MockTax;
        let mut visited = HashSet::new();
        let n_nodes = tax
            .traverse(tax.root())?
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
        Ok(())
    }
}
