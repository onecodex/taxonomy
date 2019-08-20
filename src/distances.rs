//! Code for working with distance relationships in taxonomic trees
//! and tracking weights assigned to different nodes in the tree.
use std::collections::{BinaryHeap, HashMap};
use std::fmt::{Debug, Display};
use std::hash::{BuildHasher, Hash};
use std::iter::Sum;

use crate::taxonomy::Taxonomy;
use crate::Result;

/// Given a taxonomy and a set of "weights," calculate the summed weight for
/// every path to root. This is like `maximum_weighted_path`, but rather than
/// just returning the _best_ path, it returns all of the paths.
pub fn all_weighted_paths<'t, T: 't, D: 't, S: BuildHasher>(
    taxonomy: &'t impl Taxonomy<'t, T, D>,
    weights: &HashMap<T, D, S>,
) -> Result<Vec<(D, T)>>
where
    T: Clone + Debug + Display + Hash + PartialEq + Ord,
    D: Clone + Debug + Ord + PartialOrd + PartialEq + Sum,
{
    let mut taxs_by_score: BinaryHeap<(D, T)> = BinaryHeap::new();

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
        taxs_by_score.push((scores.into_iter().sum(), tax_id.clone()));
    }
    // nodes with higher scores should be first so we need to reverse
    let mut sorted_taxs = taxs_by_score.into_sorted_vec();
    sorted_taxs.reverse();
    Ok(sorted_taxs)
}

/// Given a taxonomy and a set of "weights," find the lineage that has the
/// greatest summed weight.
///
/// Note that this implementation doesn't use the "distances" in the taxonomy
/// itself, but only uses the weights provided. This can greatly speed up the
/// calculation because it reduces the number of possible paths that need to
/// be checked.
pub fn maximum_weighted_path<'t, T: 't, D: 't, S: BuildHasher>(
    taxonomy: &'t impl Taxonomy<'t, T, D>,
    weights: &HashMap<T, D, S>,
    take_first_in_tie: bool,
) -> Result<(Option<T>, D)>
where
    T: Clone + Debug + Display + Hash + PartialEq + Ord,
    D: Clone + Debug + Ord + PartialOrd + PartialEq + Sum,
{
    let mut max_taxes: Vec<T> = Vec::new();
    // this is gross, but there's no "Zero" trait we can define to init this
    let mut max_score: D = Vec::new().into_iter().sum();
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

        let score: D = scores.into_iter().sum();
        if score > max_score {
            max_score = score.clone();
            max_taxes.clear();
        }

        if score >= max_score {
            max_taxes.push(tax_id.clone());
        }
    }

    if take_first_in_tie {
        return Ok((max_taxes.into_iter().next(), max_score));
    }

    let first_child = max_taxes.pop();
    let ancestor =
        max_taxes
            .into_iter()
            .try_fold(first_child, |ancestor, child| match ancestor {
                None => Ok(None),
                Some(a) => taxonomy.lca(a, child).map(Some),
            })?;
    Ok((ancestor, max_score))
    // is max_score going to be a little low here because we're not counting
    // all the leaf nodes?
}

#[test]
fn test_maximum_weighted_path() {
    use crate::taxonomy::test::MockTax;
    let tax = MockTax;
    let mut hits: HashMap<u32, u16> = HashMap::new();
    hits.insert(0, 41);
    hits.insert(1, 25);
    hits.insert(131567, 233);
    hits.insert(2, 512);
    hits.insert(1224, 33);
    hits.insert(1236, 275);
    hits.insert(135622, 59);
    hits.insert(22, 270);
    hits.insert(62322, 49);
    hits.insert(56812, 1);
    assert_eq!(
        maximum_weighted_path(&tax, &hits, false).unwrap().0,
        Some(56812)
    );
}
