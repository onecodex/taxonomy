use crate::errors::TaxonomyResult;
use crate::rank::TaxRank;
use std::collections::VecDeque;

pub trait Taxonomy<'t> {
    /// Returns the root node of the entire tree.
    fn root(&'t self) -> &'t str;

    /// Returns a [Vec] of all the child IDs of the given tax_id.
    fn children(&'t self, tax_id: &str) -> TaxonomyResult<Vec<&'t str>>;

    /// Returns the parent of the given taxonomic node and the distance to said parent.
    /// The parent of the root node will return [None].
    fn parent(&'t self, tax_id: &str) -> TaxonomyResult<Option<(&'t str, f32)>>;

    /// Returns a [Vec] of taxonomy nodes from the one provided back to root.
    /// This method must return the node itself as the first entry in the list
    /// and the root node as the last entry in the list.
    fn lineage(&'t self, tax_id: &'t str) -> TaxonomyResult<Vec<&'t str>> {
        let mut parents = Vec::new();
        let mut last_parent = tax_id;
        parents.push(tax_id);
        while let Some(p) = self.parent(last_parent)? {
            last_parent = p.0;
            parents.push(p.0);
        }
        Ok(parents)
    }

    /// Returns the parent at a given taxonomic rank (note this may not be
    /// the immediate parent). This also returns the distance *to* that parent.
    fn parent_at_rank(
        &'t self,
        tax_id: &'t str,
        rank: TaxRank,
    ) -> TaxonomyResult<Option<(&'t str, f32)>> {
        // if this node is at the rank, just return it
        if self.rank(tax_id)? == rank {
            return Ok(Some((tax_id, 0.0)));
        }

        // traverse up the tree looking for it
        let mut cur_id = tax_id;
        let mut dists = Vec::new();
        while let Some(p) = self.parent(cur_id)? {
            dists.push(p.1);
            if self.rank(p.0)? == rank {
                return Ok(Some((p.0, dists.into_iter().sum())));
            }
            cur_id = p.0;
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
    fn lca(&'t self, id1: &'t str, id2: &'t str) -> TaxonomyResult<&'t str> {
        // make a vec of parents of id1
        let mut id1_parents = VecDeque::new();
        id1_parents.push_front(id1);
        while let Some(p) = self.parent(*id1_parents.front().unwrap())? {
            id1_parents.push_front(p.0);
        }

        // make a vec of parents of id2
        let mut id2_parents = VecDeque::new();
        id2_parents.push_front(id2);
        while let Some(p) = self.parent(*id2_parents.front().unwrap())? {
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
    fn name(&'t self, tax_id: &str) -> TaxonomyResult<&str>;

    /// Returns the taxonomic rank of the tax_id provided.
    fn rank(&'t self, tax_id: &str) -> TaxonomyResult<TaxRank>;

    /// Used by traverse, to grab the internal tax id for the given tax id. Only there
    /// for lifetime reasons. Raise an error if not found.
    fn get_internal_tax_id(&'t self, node: &str) -> TaxonomyResult<&'t str>;

    /// Generates an iterator that traces over the entire taxonomic tree from the given node. During
    /// preorder traversal, it returns `(&str, true)` and during postorder traversal
    /// it returns `(&str, false)`
    fn traverse(&'t self, node: &str) -> TaxonomyResult<TaxonomyIterator<'t>>
    where
        Self: Sized,
    {
        let internal = self.get_internal_tax_id(node)?;
        Ok(TaxonomyIterator::new(self, internal))
    }

    /// Returns the number of nodes in the taxonomy
    fn len(&'t self) -> usize;

    /// Determines if there are any nodes at all in the taxonomy. This should
    /// almost always be implemented for performance reasons.
    fn is_empty(&'t self) -> bool {
        self.len() == 0
    }
}

pub struct TaxonomyIterator<'t> {
    nodes_left: Vec<&'t str>,
    visited_nodes: Vec<&'t str>,
    tax: &'t dyn Taxonomy<'t>,
}

impl<'t> TaxonomyIterator<'t> {
    pub fn new(tax: &'t dyn Taxonomy<'t>, root_node: &'t str) -> Self {
        TaxonomyIterator {
            nodes_left: vec![root_node],
            visited_nodes: vec![],
            tax,
        }
    }
}

impl<'t> Iterator for TaxonomyIterator<'t> {
    type Item = (&'t str, bool);

    fn next(&mut self) -> Option<Self::Item> {
        if self.nodes_left.is_empty() {
            return None;
        }

        let cur_node = *self.nodes_left.last().unwrap();
        let node_visited = {
            let last_visited = self.visited_nodes.last();
            Some(&cur_node) == last_visited
        };
        if node_visited {
            self.visited_nodes.pop();
            Some((self.nodes_left.pop().unwrap(), false))
        } else {
            self.visited_nodes.push(cur_node);
            let children = self.tax.children(cur_node).unwrap();
            if !children.is_empty() {
                self.nodes_left.extend(children);
            }
            Some((cur_node, true))
        }
    }
}
