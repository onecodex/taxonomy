use std::collections::HashMap;
use std::io::Read;

use crate::base::GeneralTaxonomy;
use crate::errors::TaxonomyResult;
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy;

// TODO: load from file or str? or both?
pub fn load(data: &str) -> TaxonomyResult<GeneralTaxonomy> {
    // no distance, no tax id

    // Should use a linked hashset for more efficient contains
    let mut tax_ids = Vec::new();
    let mut tax_names = Vec::new();
    let mut tax_names_idx = HashMap::new();
    let mut parent_ids = Vec::new();
    let mut tax_ranks = Vec::new();

    for row in data.lines() {
        // parts[0] -> accession number
        // parts[1] -> lineage, `;` separated
        let parts: Vec<_> = row.split('\t').collect();
        assert_eq!(parts.len(), 2);
        let accession_number = parts[0];
        tax_ids.push(accession_number.to_string());
        let lineage: Vec<_> = parts[1].split(';').collect();

        for (i, level) in lineage.iter().enumerate() {
            if tax_names_idx.contains_key(*level) {
                continue;
            }
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
            tax_names_idx.insert(level.to_string(), idx);
            if i > 0 {
                let parent_id = tax_names_idx.get(lineage[i - 1]).unwrap();
                parent_ids.push(*parent_id);
            } else {
                parent_ids.push(0);
            }
        }
    }
    // println!("{:?}", tax_names);

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

    #[test]
    fn can_load_gtdb_format() {
        let data = std::fs::read_to_string("tests/data/gtdb_sample.tsv").unwrap();
        let tax = load(&data).unwrap();

        let lineage = tax.lineage("s__Escherichia coli").unwrap();
        println!("{:?}", lineage);
        assert!(false);
    }
}
