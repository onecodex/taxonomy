//! Code for working with distance relationships in taxonomic trees
//! and tracking weights assigned to different nodes in the tree.
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fmt::{Debug, Display};
use std::hash::{BuildHasher, Hash};
use std::iter::Sum;

use crate::taxonomy::Taxonomy;
use crate::Result;

/// Calculates the summed weight of every path to root for a set of node
/// weights given a corresponding taxonomy.
///
/// This is like `maximum_weighted_path`, but rather than just returning the
/// _best_ path, it returns all of the paths.
pub fn all_weighted_paths<'t, D: 't, T: 't, W: 't, S: BuildHasher>(
    taxonomy: &'t impl Taxonomy<'t, T, D>,
    weights: &HashMap<T, W, S>,
) -> Result<Vec<(T, W)>>
where
    D: Debug + PartialOrd + Sum,
    T: Clone + Debug + Display + Hash + PartialEq + Ord,
    W: Clone + Debug + Ord + PartialOrd + PartialEq + Sum,
{
    let mut taxs_by_score: BinaryHeap<(W, T)> = BinaryHeap::new();
    let mut stem_ids: HashSet<T> = HashSet::new();
    for tax_id in weights.keys() {
        if stem_ids.contains(tax_id) {
            continue;
        }
        let score = taxonomy
            .lineage((*tax_id).clone())?
            .iter()
            .filter_map(|t| {
                if t != tax_id {
                    stem_ids.insert((*t).clone());
                }
                weights.get(t).cloned()
            })
            .sum();

        taxs_by_score.push((score, tax_id.clone()));
    }
    // nodes with higher scores should be first so we need to reverse
    let mut sorted_taxs = taxs_by_score.into_sorted_vec();
    sorted_taxs.reverse();
    // switch the weight/tax_id order and remove remaining stems
    let taxs = sorted_taxs
        .into_iter()
        .filter_map(|(weight, tax_id)| {
            if stem_ids.contains(&tax_id) {
                None
            } else {
                Some((tax_id, weight))
            }
        })
        .collect();
    Ok(taxs)
}

/// Find the lineage that has the greatest summed weight from all of the
/// weights and a corresponding taxonomy.
///
/// Note that this implementation doesn't use the "distances" in the taxonomy
/// itself, but only uses the weights provided. This can greatly speed up the
/// calculation because it reduces the number of possible paths that need to
/// be checked.
///
/// If `take_first_in_tie` is set the tax_id returned will be the most
/// specific possible (e.g. a strain), but may not be completely correct if
/// there are multiple specific paths (e.g. multiple strains with the same
/// path weight). Setting to false will return the common ancestor of those
/// strains.
pub fn maximum_weighted_path<'t, D: 't, T: 't, W: 't, S: BuildHasher>(
    taxonomy: &'t impl Taxonomy<'t, T, D>,
    weights: &HashMap<T, W, S>,
    take_first_in_tie: bool,
) -> Result<Option<(T, W)>>
where
    D: Debug + PartialOrd + Sum,
    T: Clone + Debug + Display + Hash + PartialEq + Ord,
    W: Clone + Debug + PartialOrd + PartialEq + Sum,
{
    let mut max_taxes: Vec<T> = Vec::new();
    // this is gross, but there's no "Zero" trait we can define to init this
    let mut max_score: W = Vec::new().into_iter().sum();
    for tax_id in weights.keys() {
        let mut tax_node = tax_id.clone();
        let mut scores = Vec::new();
        loop {
            if let Some(score) = weights.get(&tax_node) {
                scores.push((*score).clone());
            }
            match taxonomy.parent(tax_node)? {
                Some(p) => tax_node = p.0,
                None => break,
            }
        }

        let score: W = scores.into_iter().sum();
        if score > max_score {
            max_score = score.clone();
            max_taxes.clear();
        }

        if score >= max_score {
            max_taxes.push(tax_id.clone());
        }
    }

    if take_first_in_tie {
        return Ok(max_taxes.into_iter().next().map(|a| (a, max_score)));
    }

    let first_child = max_taxes.pop();
    let ancestor =
        max_taxes
            .into_iter()
            .try_fold(first_child, |ancestor, child| match ancestor {
                None => Ok(None),
                Some(a) => taxonomy.lca(a, child).map(Some),
            })?;

    Ok(ancestor.map(|a| (a, max_score)))
    // is max_score going to be a little low here because we're not counting
    // all the leaf nodes?
}

/// Coverts a set of weights into the set of weights including their children
/// using the corresponding taxonomy.
///
/// For example, for a taxonomic classification of a set of sequencing reads
/// where a classification may be non-specific to a leaf node, this would
/// turn the raw read counts at each taxonomic node into the set of read
/// counts including all of the children of that node (which is probably more
/// useful for most use cases).
pub fn rollup_weights<'t, D: 't, T: 't, W: 't, S: BuildHasher>(
    taxonomy: &'t impl Taxonomy<'t, T, D>,
    weights: &HashMap<T, W, S>,
) -> Result<Vec<(T, W)>>
where
    D: Debug + PartialOrd + Sum,
    T: Clone + Debug + Display + Hash + PartialEq + Ord,
    W: Clone + Debug + PartialOrd + PartialEq + Sum,
{
    let mut all_weights: HashMap<T, Vec<W>> = HashMap::new();
    for (leaf_id, weight) in weights {
        for tax_id in taxonomy.lineage(leaf_id.clone())? {
            let tax_weights = all_weights.entry(tax_id).or_insert_with(Vec::new);
            tax_weights.push(weight.clone());
        }
    }
    Ok(all_weights
        .into_iter()
        .map(|(tax_id, all_weights)| (tax_id, all_weights.into_iter().sum()))
        .collect())
}

#[test]
fn test_all_weighted_path() -> Result<()> {
    use crate::taxonomy::test::MockTax;
    let tax = MockTax;
    let mut hits: HashMap<u32, u16> = HashMap::new();
    hits.insert(765909, 41);
    hits.insert(1, 25);
    hits.insert(131567, 233);
    hits.insert(2, 512);
    hits.insert(1224, 33);
    hits.insert(1236, 275);
    hits.insert(135622, 59);
    hits.insert(22, 270);
    hits.insert(62322, 49);
    hits.insert(56812, 1);
    let weights = all_weighted_paths(&tax, &hits)?;
    assert_eq!(weights, vec![(56812, 1457), (765909, 1119)]);
    Ok(())
}

#[test]
fn test_maximum_weighted_path() -> Result<()> {
    use crate::taxonomy::test::MockTax;
    let tax = MockTax;
    let mut hits: HashMap<u32, u16> = HashMap::new();
    hits.insert(765909, 41);
    hits.insert(1, 25);
    hits.insert(131567, 233);
    hits.insert(2, 512);
    hits.insert(1224, 33);
    hits.insert(1236, 275);
    hits.insert(135622, 59);
    hits.insert(22, 270);
    hits.insert(62322, 49);
    hits.insert(56812, 1);
    let (node, weight) = maximum_weighted_path(&tax, &hits, false)?.unwrap();
    assert_eq!(node, 56812);
    assert_eq!(weight, 25 + 233 + 512 + 33 + 275 + 59 + 270 + 49 + 1);

    let mut eq_hits: HashMap<u32, u16> = HashMap::new();
    eq_hits.insert(2, 100);
    eq_hits.insert(56812, 10);
    eq_hits.insert(765909, 10);
    let (node, weight) = maximum_weighted_path(&tax, &eq_hits, true)?.unwrap();
    assert!(node == 56812 || node == 765909);
    assert_eq!(weight, 110);

    Ok(())
}

#[test]
fn test_rollup() -> Result<()> {
    use crate::taxonomy::test::MockTax;
    let tax = MockTax;
    let mut hits: HashMap<u32, u16> = HashMap::new();
    hits.insert(1, 25);
    hits.insert(2, 512);
    hits.insert(1224, 33);
    hits.insert(56812, 1);
    hits.insert(765909, 41);
    let mut rolled_hits: Vec<(u32, u16)> = rollup_weights(&tax, &hits)?.into_iter().collect();
    rolled_hits.sort();
    assert_eq!(
        rolled_hits,
        vec![
            (1, 25 + 512 + 33 + 1 + 41),
            (2, 512 + 33 + 1 + 41),
            (22, 1),
            (1046, 41),
            (1224, 33 + 1 + 41),
            (1236, 1 + 41),
            (53452, 41),
            (56812, 1),
            (61598, 41),
            (62322, 1),
            (131567, 587),
            (135613, 41),
            (135622, 1),
            (765909, 41),
        ]
    );
    Ok(())
}
