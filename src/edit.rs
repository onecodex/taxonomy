//! Functions for changing the topology or contents of a taxonomy
use std::collections::HashSet;

use crate::base::{GeneralTaxonomy, IntTaxID};
use crate::taxonomy::Taxonomy;
use crate::Result;

fn filter_to_nodes(tax: &GeneralTaxonomy, tax_nodes: &[IntTaxID]) -> GeneralTaxonomy {
    let mut tax_ids = Vec::new();
    let mut parent_ids = Vec::new();
    let mut dists = Vec::new();
    let mut names = Vec::new();
    let mut ranks = Vec::new();
    for node in tax_nodes {
        tax_ids.push(tax.tax_ids[*node].clone());
        parent_ids.push(tax.parent_ids[*node]);
        dists.push(tax.parent_dists[*node]);
        names.push(tax.names[*node].clone());
        ranks.push(tax.ranks[*node].clone());
    }

    GeneralTaxonomy {
        tax_ids,
        parent_ids,
        parent_dists: dists,
        ranks,
        names,
        tax_id_lookup: None,
        children_lookup: None,
    }
}

/// Return a tree with these tax_ids and their children removed.
///
/// Note this uses internal indices so you'll have to convert "string"
/// keys using `to_internal_id` if you have those.
pub fn prune_away(tax: &GeneralTaxonomy, tax_ids: &[IntTaxID]) -> Result<GeneralTaxonomy> {
    let tax_set: HashSet<IntTaxID> = tax_ids.iter().cloned().collect();
    let mut good_ids: Vec<IntTaxID> = Vec::new();
    let mut dropping: u8 = 0;
    for (node, pre) in tax.traverse(tax.root())? {
        if tax_set.contains(&node) {
            if pre {
                dropping += 1;
            } else {
                dropping -= 1;
            }
        }
        if dropping == 0 {
            good_ids.push(node);
        }
    }
    Ok(filter_to_nodes(tax, &good_ids))
}

/// Return a tree containing only the given tax_ids and their parents.
///
/// Note this uses internal indices so you'll have to convert "string"
/// keys using `to_internal_id` if you have those.
pub fn prune_to(
    tax: &GeneralTaxonomy,
    tax_ids: &[IntTaxID],
    include_children: bool,
) -> Result<GeneralTaxonomy> {
    let mut good_ids: HashSet<IntTaxID> = tax_ids.iter().cloned().collect();
    for (node, pre) in tax.traverse(tax.root())? {
        if let Some((parent_node, _)) = tax.parent(node)? {
            // insert child nodes on the traverse down (add node if parent is in)
            if pre && include_children && good_ids.contains(&parent_node) {
                good_ids.insert(node);
            }

            // insert parent nodes of stuff we've seen on the traverse back
            // (note this will add some duplicates of the "children" nodes
            // but since this is a set that still works)
            if !pre && good_ids.contains(&node) {
                good_ids.insert(parent_node);
            }
        }
    }
    let good_ids: Vec<IntTaxID> = good_ids.into_iter().collect();
    Ok(filter_to_nodes(tax, &good_ids))
}

// FIXME: hammer this out
// pub fn filter(
//     tax: &GeneralTaxonomy,
//     tax_ids: &[IntTaxID],
//     remove_ids: bool,
// ) -> GeneralTaxonomy {
//     let good_ids: HashSet<IntTaxID> = tax_ids.iter().cloned().collect();
//     let mut cur_lineage: Vec<IntTaxID> = Vec::new();
//     // TODO: need to track the names, ranks, parents, dists
//     for (node, pre) in tax.traverse(tax.root()) {
//         if !remove_ids && good_ids.contains(&node) {
//             if pre {
//                 // TODO: add node to new taxonomy?
//                 // TODO: add node to cur_lineage
//             } else {
//                 // TODO: pop last node off cur_lineage (and check it matches?)
//             }
//         }
//     }
// }
