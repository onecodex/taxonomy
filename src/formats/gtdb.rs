use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

use crate::base::{GeneralTaxonomy, InternalIndex};
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;

/// Read GTDB format into a Taxonomy object out of a `reader`.
///
/// Still somewhat experimental and may not support all GTDB features.
pub fn load<R: Read>(reader: &mut R) -> TaxonomyResult<GeneralTaxonomy> {
    let mut tax_names = Vec::new();
    // Maps tax name to (index, lineage prefix)
    let mut tax_names_idx: HashMap<String, (InternalIndex, String)> = HashMap::new();
    let mut parent_ids = Vec::new();
    let mut tax_ranks = Vec::new();

    for (row_idx, row_result) in BufReader::new(reader).lines().enumerate() {
        let line_num = row_idx + 1;
        let row = row_result?;
        // parts[0] -> accession number
        // parts[1] -> lineage, `;` separated
        let parts: Vec<_> = row.split('\t').collect();
        if parts.len() != 2 {
            return Err(Error::new(ErrorKind::ImportError {
                line: line_num,
                msg: "Expected tab-delimited line with exactly two parts (accession and lineage)"
                    .to_owned(),
            }));
        }
        let lineage: Vec<_> = parts[1].split(';').collect();

        for (i, level) in lineage.iter().enumerate() {
            let curr_prefix = lineage[..i].join(";");
            match tax_names_idx.get(*level) {
                Some((_, prev_prefix)) => {
                    if curr_prefix != *prev_prefix {
                        return Err(Error::new(ErrorKind::ImportError {
                            line: line_num,
                            msg: format!(
                                "Inconsistent lineages for taxon {}: {} != {}",
                                level, curr_prefix, prev_prefix
                            ),
                        }));
                    }
                    continue;
                }
                None => {
                    let tax_rank = match &level[..3] {
                        "d__" => TaxRank::Domain,
                        "p__" => TaxRank::Phylum,
                        "c__" => TaxRank::Class,
                        "o__" => TaxRank::Order,
                        "f__" => TaxRank::Family,
                        "g__" => TaxRank::Genus,
                        "s__" => TaxRank::Species,
                        _ => TaxRank::Unspecified,
                    };
                    tax_names.push(level.to_string());
                    tax_ranks.push(tax_rank);
                    let idx = tax_names.len() - 1;
                    tax_names_idx.insert(level.to_string(), (idx, curr_prefix));
                    if i > 0 {
                        let parent_id = tax_names_idx.get(lineage[i - 1]).unwrap().0;
                        parent_ids.push(parent_id);
                    } else {
                        parent_ids.push(0);
                    }
                }
            }
        }
    }

    GeneralTaxonomy::from_arrays(
        tax_names.clone(),
        parent_ids,
        Some(tax_names),
        Some(tax_ranks),
        None,
        None,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::taxonomy::Taxonomy;
    use std::fs::File;

    #[test]
    fn can_load_gtdb_format() {
        let mut file = File::open("tests/data/gtdb_sample.tsv").unwrap();
        let tax = load(&mut file).unwrap();

        let mut tax_id = "d__Bacteria";
        assert_eq!(tax.rank(tax_id).unwrap(), TaxRank::Domain);
        assert_eq!(tax.parent(tax_id).unwrap(), None);
        assert_eq!(tax.lineage(tax_id).unwrap(), ["d__Bacteria"]);

        tax_id = "c__Bacilli";
        assert_eq!(tax.rank(tax_id).unwrap(), TaxRank::Class);
        assert_eq!(
            tax.lineage(tax_id).unwrap(),
            ["c__Bacilli", "p__Firmicutes", "d__Bacteria"]
        );

        tax_id = "s__Escherichia coli";
        assert_eq!(tax.rank(tax_id).unwrap(), TaxRank::Species);
        assert_eq!(
            tax.lineage(tax_id).unwrap(),
            [
                "s__Escherichia coli",
                "g__Escherichia",
                "f__Enterobacteriaceae",
                "o__Enterobacterales",
                "c__Gammaproteobacteria",
                "p__Proteobacteria",
                "d__Bacteria"
            ]
        );
    }

    #[test]
    fn invalid_gtdb_format() {
        let mut file = File::open("tests/data/gtdb_invalid.tsv").unwrap();
        assert!(load(&mut file).is_err());
    }

    #[test]
    fn inconsistent_gtdb_taxonomy() {
        let mut file = File::open("tests/data/gtdb_inconsistent.tsv").unwrap();
        assert!(load(&mut file).is_err());
    }
}
