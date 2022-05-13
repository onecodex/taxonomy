use std::collections::HashSet;
use std::fmt::{Debug, Display};
use std::hash::Hash;

use crate::errors::TaxonomyResult;
use crate::{GeneralTaxonomy, Taxonomy};

/// Return a tree with these tax_ids and their children removed.
///
/// Note this uses internal indices so you'll have to convert "string"
/// keys using `to_internal_index` if you have those.
pub fn prune_away<'t, T: 't + Clone + Copy + Debug + Display + Eq + Hash + PartialEq>(
    tax: &'t impl Taxonomy<'t, T>,
    tax_ids: &[T],
) -> TaxonomyResult<GeneralTaxonomy> {
    let mut new_ids = Vec::new();
    let mut parent_ids = Vec::new();
    let mut dists = Vec::new();
    let mut names = Vec::new();
    let mut ranks = Vec::new();
    let mut data = Vec::new();

    let tax_set: HashSet<T> = tax_ids.iter().cloned().collect();
    let mut dropping: u8 = 0;
    let mut cur_lineage = Vec::new();
    for (node, pre) in tax.traverse(tax.root())? {
        if tax_set.contains(&node) {
            if pre {
                dropping += 1;
            } else {
                dropping -= 1;
            }
        }
        if dropping == 0 {
            if pre {
                new_ids.push(node.to_string());
                parent_ids.push(cur_lineage.last().map(|x| x - 1).unwrap_or(0));
                dists.push(tax.parent(node)?.map(|x| x.1).unwrap_or(0.));
                names.push(tax.name(node)?.to_string());
                ranks.push(tax.rank(node)?);
                data.push((*tax.data(node)?).to_owned());

                cur_lineage.push(new_ids.len());
            } else {
                cur_lineage.pop();
            }
        }
    }
    GeneralTaxonomy::from_arrays(
        new_ids,
        parent_ids,
        Some(names),
        Some(ranks),
        Some(dists),
        Some(data),
    )
}

/// Return a tree containing only the given tax_ids and their parents.
///
/// Note this uses internal indices so you'll have to convert "string"
/// keys using `to_internal_index` if you have those.
pub fn prune_to<'t, T: 't>(
    tax: &'t impl Taxonomy<'t, T>,
    tax_ids: &[T],
    include_children: bool,
) -> TaxonomyResult<GeneralTaxonomy>
where
    T: Clone + Copy + Debug + Display + Eq + Hash + PartialEq,
{
    let mut good_ids: HashSet<T> = tax_ids.iter().cloned().collect();
    for (node, pre) in tax.traverse(tax.root())? {
        if let Some((parent_node, _)) = tax.parent(node)? {
            if pre && include_children && good_ids.contains(&parent_node) {
                // insert child nodes on the traverse down (add node if parent is in)
                good_ids.insert(node);
            } else if !pre && good_ids.contains(&node) {
                // insert parent nodes of stuff we've seen on the traverse back
                // (note this will add some duplicates of the "children" nodes
                // but since this is a set that still works)
                good_ids.insert(parent_node);
            }
        }
    }

    // create a new taxonomy
    let mut new_ids = Vec::new();
    let mut parent_ids = Vec::new();
    let mut dists = Vec::new();
    let mut names = Vec::new();
    let mut ranks = Vec::new();
    let mut data = Vec::new();

    let mut cur_lineage = Vec::new();
    for (node, pre) in tax.traverse(tax.root())? {
        if pre {
            if good_ids.contains(&node) {
                new_ids.push(node.to_string());
                parent_ids.push(cur_lineage.last().map(|x| x - 1).unwrap_or(0));
                dists.push(tax.parent(node)?.map(|x| x.1).unwrap_or(0.));
                names.push(tax.name(node)?.to_string());
                ranks.push(tax.rank(node)?);
                data.push((*tax.data(node)?).to_owned());
            }
            cur_lineage.push(new_ids.len());
        } else {
            cur_lineage.pop();
        }
    }

    GeneralTaxonomy::from_arrays(
        new_ids,
        parent_ids,
        Some(names),
        Some(ranks),
        Some(dists),
        Some(data),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::taxonomy::tests::MockTax;

    #[test]
    fn test_prune_away() {
        let tax = MockTax;
        assert_eq!(Taxonomy::len(&tax), 14);

        let pruned = prune_away(&tax, &[135622]).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&pruned), 10);
        assert_eq!(pruned.parent("131567").unwrap().unwrap().0, "1");
    }

    #[test]
    fn test_prune_to() {
        let tax = MockTax;
        assert_eq!(Taxonomy::len(&tax), 14);

        let pruned = prune_to(&tax, &[62322], true).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&pruned), 9);
        assert_eq!(pruned.parent("131567").unwrap().unwrap().0, "1");
    }
}
