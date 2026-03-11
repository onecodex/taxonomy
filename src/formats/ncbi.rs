use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

use crate::base::{GeneralTaxonomy, TaxonomyValue};
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy;

const NODES_FILENAME: &str = "nodes.dmp";
const NAMES_FILENAME: &str = "names.dmp";

/// Represents a taxonomic name with its class and optional unique name
#[derive(Debug, Clone)]
struct TaxonomyName {
    name_class: String,
    name_text: String,
    unique_name: Option<String>,
}

impl TaxonomyName {
    fn to_taxonomy_value(&self) -> TaxonomyValue {
        let mut map = HashMap::new();
        map.insert(
            "name_class".to_string(),
            TaxonomyValue::String(self.name_class.clone()),
        );
        map.insert(
            "name_text".to_string(),
            TaxonomyValue::String(self.name_text.clone()),
        );
        if let Some(ref unique) = self.unique_name {
            map.insert(
                "unique_name".to_string(),
                TaxonomyValue::String(unique.clone()),
            );
        }
        TaxonomyValue::Object(map)
    }

    fn from_taxonomy_value(value: &TaxonomyValue) -> Option<Self> {
        if let TaxonomyValue::Object(map) = value {
            let name_class = map.get("name_class")?.as_str()?.to_string();
            let name_text = map.get("name_text")?.as_str()?.to_string();
            let unique_name = map
                .get("unique_name")
                .and_then(|v| v.as_str())
                .map(|s| s.to_string());
            Some(TaxonomyName {
                name_class,
                name_text,
                unique_name,
            })
        } else {
            None
        }
    }
}

/// Data for a single row in nodes.dmp
struct NodesRow<'a> {
    tax_id: &'a str,
    parent: &'a str,
    rank: &'a str,
    embl_code: &'a str,
    division_id: &'a str,
    inherited_div: &'a str,
    genetic_code: &'a str,
    inherited_gc: &'a str,
    mito_gc: &'a str,
    inherited_mgc: &'a str,
    genbank_hidden: &'a str,
    subtree_hidden: &'a str,
    comments: &'a str,
}

/// Loads a NCBI taxonomy from the given directory.
/// The directory should contain at least two files: `nodes.dmp` and `names.dmp`.
pub fn load<P: AsRef<Path>>(ncbi_directory: P) -> TaxonomyResult<GeneralTaxonomy> {
    let dir = ncbi_directory.as_ref();
    let nodes_file = std::fs::File::open(dir.join(NODES_FILENAME))?;
    let names_file = std::fs::File::open(dir.join(NAMES_FILENAME))?;

    // First we go through the nodes
    let mut tax_ids: Vec<String> = Vec::new();
    let mut parents: Vec<String> = Vec::new();
    let mut ranks: Vec<TaxRank> = Vec::new();
    let mut data: Vec<HashMap<String, TaxonomyValue>> = Vec::new();
    let mut tax_to_idx: HashMap<String, usize> = HashMap::new();

    for (ix, line) in BufReader::new(nodes_file).lines().enumerate() {
        let mut fields: Vec<String> = line?.split("\t|\t").map(|x| x.to_string()).collect();
        if fields.len() < 12 {
            return Err(Error::new(ErrorKind::ImportError {
                line: ix,
                msg: "Not enough fields in nodes.dmp; bad line?".to_owned(),
            }));
        }
        let tax_id = fields.remove(0).trim().to_string();
        let parent_tax_id = fields.remove(0).trim().to_string();
        let rank = TaxRank::from_str(fields.remove(0).trim())?;
        let embl_code = fields.remove(0).trim().to_string();
        let division_id = fields.remove(0).trim().parse::<i64>().ok();
        let inherited_div_flag = fields.remove(0).trim() == "1";
        let genetic_code_id = fields.remove(0).trim().parse::<i64>().ok();
        let inherited_gc_flag = fields.remove(0).trim() == "1";
        let mitochondrial_genetic_code_id = fields.remove(0).trim().parse::<i64>().ok();
        let inherited_mgc_flag = fields.remove(0).trim() == "1";
        let genbank_hidden_flag = fields.remove(0).trim() == "1";
        let hidden_subtree_root_flag = fields.remove(0).trim() == "1";
        let comments = if !fields.is_empty() {
            fields
                .remove(0)
                .trim()
                .trim_end_matches('|')
                .trim()
                .to_string()
        } else {
            String::new()
        };

        tax_ids.push(tax_id.clone());
        parents.push(parent_tax_id.to_string());
        ranks.push(rank);

        // Store all node fields in data HashMap with proper types
        let mut node_data = HashMap::new();
        if !embl_code.is_empty() {
            node_data.insert("embl_code".to_string(), TaxonomyValue::String(embl_code));
        }
        if let Some(div_id) = division_id {
            node_data.insert("division_id".to_string(), TaxonomyValue::Integer(div_id));
        }
        node_data.insert(
            "inherited_div_flag".to_string(),
            TaxonomyValue::Bool(inherited_div_flag),
        );
        if let Some(gc_id) = genetic_code_id {
            node_data.insert("genetic_code_id".to_string(), TaxonomyValue::Integer(gc_id));
        }
        node_data.insert(
            "inherited_gc_flag".to_string(),
            TaxonomyValue::Bool(inherited_gc_flag),
        );
        if let Some(mgc_id) = mitochondrial_genetic_code_id {
            node_data.insert(
                "mitochondrial_genetic_code_id".to_string(),
                TaxonomyValue::Integer(mgc_id),
            );
        }
        node_data.insert(
            "inherited_mgc_flag".to_string(),
            TaxonomyValue::Bool(inherited_mgc_flag),
        );
        node_data.insert(
            "genbank_hidden_flag".to_string(),
            TaxonomyValue::Bool(genbank_hidden_flag),
        );
        node_data.insert(
            "hidden_subtree_root_flag".to_string(),
            TaxonomyValue::Bool(hidden_subtree_root_flag),
        );
        if !comments.is_empty() {
            node_data.insert("comments".to_string(), TaxonomyValue::String(comments));
        }

        data.push(node_data);
        tax_to_idx.insert(tax_id, ix);
    }

    // TODO: fixme? this fails if we have unmapped parent nodes (i.e. file is truncated?)
    let mut parent_ids = Vec::with_capacity(parents.len());
    for (i, parent) in parents.into_iter().enumerate() {
        if let Some(idx) = tax_to_idx.get(&parent) {
            parent_ids.push(*idx);
        } else {
            return Err(Error::new(ErrorKind::ImportError {
                line: i + 1,
                msg: format!("Parent ID {} could not be found in nodes.dmp", parent),
            }));
        }
    }

    // Grab scientific names and store all other name types in data
    let mut names: Vec<String> = vec![String::new(); tax_ids.len()];
    for (ix, line) in BufReader::new(names_file).lines().enumerate() {
        let mut fields: Vec<String> = line?.split("\t|\t").map(|x| x.to_string()).collect();
        if fields.len() < 4 {
            return Err(Error::new(ErrorKind::ImportError {
                line: ix,
                msg: "Not enough fields in names.dmp".to_owned(),
            }));
        }
        let tax_id = fields.remove(0).trim().to_string();
        let name_txt = fields.remove(0).trim().to_string();
        let unique_name = fields.remove(0).trim().to_string();
        let name_class = fields
            .remove(0)
            .trim()
            .trim_end_matches('|')
            .trim()
            .to_string();

        if let Some(&idx) = tax_to_idx.get(&tax_id) {
            if name_class == "scientific name" {
                names[idx] = name_txt.clone();
            } else {
                // Store non-scientific name types in data to avoid duplication
                // (scientific name is already in the main names vector)
                let taxonomy_name = TaxonomyName {
                    name_class: name_class.clone(),
                    name_text: name_txt,
                    unique_name: if unique_name.is_empty() {
                        None
                    } else {
                        Some(unique_name)
                    },
                };

                // Get or create the names vector for this taxon
                let names_key = "names";
                let names_array = data[idx]
                    .entry(names_key.to_string())
                    .or_insert_with(|| TaxonomyValue::Array(Vec::new()));

                if let Some(arr) = names_array.as_array_mut() {
                    arr.push(taxonomy_name.to_taxonomy_value());
                }
            }
        }
    }

    let gt = GeneralTaxonomy::from_arrays(
        tax_ids,
        parent_ids,
        Some(names),
        Some(ranks),
        None,
        Some(data),
    )?;
    gt.validate_uniqueness()?;
    Ok(gt)
}

/// Helper function to write a single row to nodes.dmp
fn write_nodes_row(writer: &mut BufWriter<std::fs::File>, row: &NodesRow) -> std::io::Result<()> {
    writeln!(
        writer,
        "{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|",
        row.tax_id,
        row.parent,
        row.rank,
        row.embl_code,
        row.division_id,
        row.inherited_div,
        row.genetic_code,
        row.inherited_gc,
        row.mito_gc,
        row.inherited_mgc,
        row.genbank_hidden,
        row.subtree_hidden,
        row.comments
    )
}

/// Helper function to write a single row to names.dmp
fn write_names_row(
    writer: &mut BufWriter<std::fs::File>,
    tax_id: &str,
    name: &str,
    unique_name: &str,
    name_class: &str,
) -> std::io::Result<()> {
    writeln!(
        writer,
        "{}\t|\t{}\t|\t{}\t|\t{}\t|",
        tax_id, name, unique_name, name_class
    )
}

pub fn save<'t, T, P: AsRef<Path>, X: Taxonomy<'t, T>>(
    tax: &'t X,
    out_dir: P,
    include_extended: bool,
) -> TaxonomyResult<()>
where
    T: 't + Clone + Debug + Display + PartialEq,
{
    let dir = out_dir.as_ref();
    std::fs::create_dir_all(dir)?;
    let mut node_writer = BufWriter::new(std::fs::File::create(dir.join(NODES_FILENAME))?);
    let mut name_writer = BufWriter::new(std::fs::File::create(dir.join(NAMES_FILENAME))?);

    let root = tax.root();
    for key in tax.traverse(root.clone())?.filter(|x| x.1).map(|x| x.0) {
        let name = tax.name(key.clone())?;
        let rank = tax.rank(key.clone())?;
        let parent = if key == root {
            format!("{}", key)
        } else {
            tax.parent(key.clone())?
                .map(|(x, _)| format!("{}", x))
                .unwrap_or_default()
        };

        // Extract data fields with proper types
        let node_data = tax.data(key.clone())?;

        let (
            embl_code,
            division_id,
            inherited_div_flag,
            genetic_code_id,
            inherited_gc_flag,
            mitochondrial_genetic_code_id,
            inherited_mgc_flag,
            genbank_hidden_flag,
            hidden_subtree_root_flag,
            comments,
        ) = if include_extended {
            // Extract extended fields from data
            let embl_code = node_data
                .get("embl_code")
                .and_then(|v| v.as_str())
                .unwrap_or("");
            let division_id = node_data
                .get("division_id")
                .and_then(|v| v.as_i64())
                .map(|i| i.to_string())
                .unwrap_or_else(|| "".to_string());
            let inherited_div_flag = node_data
                .get("inherited_div_flag")
                .and_then(|v| v.as_bool())
                .map(|b| if b { "1" } else { "0" })
                .unwrap_or("");
            let genetic_code_id = node_data
                .get("genetic_code_id")
                .and_then(|v| v.as_i64())
                .map(|i| i.to_string())
                .unwrap_or_else(|| "".to_string());
            let inherited_gc_flag = node_data
                .get("inherited_gc_flag")
                .and_then(|v| v.as_bool())
                .map(|b| if b { "1" } else { "0" })
                .unwrap_or("");
            let mitochondrial_genetic_code_id = node_data
                .get("mitochondrial_genetic_code_id")
                .and_then(|v| v.as_i64())
                .map(|i| i.to_string())
                .unwrap_or_else(|| "".to_string());
            let inherited_mgc_flag = node_data
                .get("inherited_mgc_flag")
                .and_then(|v| v.as_bool())
                .map(|b| if b { "1" } else { "0" })
                .unwrap_or("");
            let genbank_hidden_flag = node_data
                .get("genbank_hidden_flag")
                .and_then(|v| v.as_bool())
                .map(|b| if b { "1" } else { "0" })
                .unwrap_or("");
            let hidden_subtree_root_flag = node_data
                .get("hidden_subtree_root_flag")
                .and_then(|v| v.as_bool())
                .map(|b| if b { "1" } else { "0" })
                .unwrap_or("");
            let comments = node_data
                .get("comments")
                .and_then(|v| v.as_str())
                .unwrap_or("");

            (
                embl_code,
                division_id,
                inherited_div_flag,
                genetic_code_id,
                inherited_gc_flag,
                mitochondrial_genetic_code_id,
                inherited_mgc_flag,
                genbank_hidden_flag,
                hidden_subtree_root_flag,
                comments,
            )
        } else {
            // Use empty/default values for backwards compatibility
            (
                "",
                "".to_string(),
                "",
                "".to_string(),
                "",
                "".to_string(),
                "",
                "",
                "",
                "",
            )
        };

        // Write scientific name
        write_names_row(
            &mut name_writer,
            &format!("{}", &key),
            &name,
            "",
            "scientific name",
        )?;

        // Write all alternative names from data (only if extended fields are included)
        if include_extended {
            if let Some(names_value) = node_data.get("names") {
                if let Some(names_array) = names_value.as_array() {
                    for name_value in names_array {
                        if let Some(tax_name) = TaxonomyName::from_taxonomy_value(name_value) {
                            // Skip scientific name as it's already written above
                            if tax_name.name_class != "scientific name" {
                                write_names_row(
                                    &mut name_writer,
                                    &format!("{}", &key),
                                    &tax_name.name_text,
                                    tax_name.unique_name.as_deref().unwrap_or(""),
                                    &tax_name.name_class,
                                )?;
                            }
                        }
                    }
                }
            }
        }

        // Write nodes.dmp entry with all fields
        let tax_id_str = format!("{}", &key);
        write_nodes_row(
            &mut node_writer,
            &NodesRow {
                tax_id: &tax_id_str,
                parent: &parent,
                rank: rank.to_ncbi_rank(),
                embl_code,
                division_id: &division_id,
                inherited_div: inherited_div_flag,
                genetic_code: &genetic_code_id,
                inherited_gc: inherited_gc_flag,
                mito_gc: &mitochondrial_genetic_code_id,
                inherited_mgc: inherited_mgc_flag,
                genbank_hidden: genbank_hidden_flag,
                subtree_hidden: hidden_subtree_root_flag,
                comments,
            },
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::taxonomy::Taxonomy;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn can_import_ncbi() {
        let nodes =
            "1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
    10239\t|\t1\t|\tno_rank\t|\t\t|\t9\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
    2\t|\t131567\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
    543\t|\t91347\t|\tfamily\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
    561\t|\t543\t|\tgenus\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
    562\t|\t561\t|\tspecies\t|\tEC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\t\t|
    1224\t|\t2\t|\tphylum\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
    1236\t|\t1224\t|\tclass\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
    91347\t|\t1236\t|\torder\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|
    131567\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t1\t|\t1\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\t\t|";
        let names = "1\t|\tall\t|\t\t|\tsynonym\t|
    1\t|\troot\t|\t\t|\tscientific name\t|
    10239\t|\tViruses\t|\t\t|\tscientific name\t|
    2\t|\tBacteria\t|\tBacteria <prokaryotes>\t|\tscientific name\t|
    543\t|\tEnterobacteriaceae\t|\t\t|\tscientific name\t|
    561\t|\tEscherchia\t|\t\t|\tmisspelling\t|
    561\t|\tEscherichia\t|\t\t|\tscientific name\t|
    561\t|\tEscherichia Castellani and Chalmers 1919\t|\t\t|\tauthority\t|
    562\t|\t\"Bacillus coli\" Migula 1895\t|\t\t|\tauthority\t|
    562\t|\t\"Bacterium coli commune\" Escherich 1885\t|\t\t|\tauthority\t|
    562\t|\t\"Bacterium coli\" (Migula 1895) Lehmann and Neumann 1896\t|\t\t|\tauthority\t|
    562\t|\tATCC 11775\t|\t\t|\ttype material\t|
    562\t|\tBacillus coli\t|\t\t|\tsynonym\t|
    562\t|\tBacterium coli\t|\t\t|\tsynonym\t|
    562\t|\tBacterium coli commune\t|\t\t|\tsynonym\t|
    562\t|\tCCUG 24\t|\t\t|\ttype material\t|
    562\t|\tCCUG 29300\t|\t\t|\ttype material\t|
    562\t|\tCIP 54.8\t|\t\t|\ttype material\t|
    562\t|\tDSM 30083\t|\t\t|\ttype material\t|
    562\t|\tE. coli\t|\t\t|\tcommon name\t|
    562\t|\tEnterococcus coli\t|\t\t|\tsynonym\t|
    562\t|\tEscherchia coli\t|\t\t|\tmisspelling\t|
    562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|
    562\t|\tEscherichia coli (Migula 1895) Castellani and Chalmers 1919\t|\t\t|\tauthority\t|
    562\t|\tEscherichia sp. 3_2_53FAA\t|\t\t|\tincludes\t|
    562\t|\tEscherichia sp. MAR\t|\t\t|\tincludes\t|
    562\t|\tEscherichia/Shigella coli\t|\t\t|\tequivalent name\t|
    562\t|\tEschericia coli\t|\t\t|\tmisspelling\t|
    562\t|\tJCM 1649\t|\t\t|\ttype material\t|
    562\t|\tLMG 2092\t|\t\t|\ttype material\t|
    562\t|\tNBRC 102203\t|\t\t|\ttype material\t|
    562\t|\tNCCB 54008\t|\t\t|\ttype material\t|
    562\t|\tNCTC 9001\t|\t\t|\ttype material\t|
    562\t|\tbacterium 10a\t|\t\t|\tincludes\t|
    562\t|\tbacterium E3\t|\t\t|\tincludes\t|
    1224\t|\tAlphaproteobacteraeota\t|\t\t|\tsynonym\t|
    1224\t|\tAlphaproteobacteraeota Oren et al. 2015\t|\t\t|\tauthority\t|
    1224\t|\tProteobacteria\t|\t\t|\tscientific name\t|
    1224\t|\tProteobacteria Garrity et al. 2005\t|\t\t|\tauthority\t|
    1224\t|\tProteobacteria [class] Stackebrandt et al. 1988\t|\t\t|\tauthority\t|
    1224\t|\tnot Proteobacteria Cavalier-Smith 2002\t|\t\t|\tauthority\t|
    1224\t|\tproteobacteria\t|\tproteobacteria<blast1224>\t|\tblast name\t|
    1224\t|\tpurple bacteria\t|\t\t|\tcommon name\t|
    1224\t|\tpurple bacteria and relatives\t|\t\t|\tcommon name\t|
    1224\t|\tpurple non-sulfur bacteria\t|\t\t|\tcommon name\t|
    1224\t|\tpurple photosynthetic bacteria\t|\t\t|\tcommon name\t|
    1224\t|\tpurple photosynthetic bacteria and relatives\t|\t\t|\tcommon name\t|
    1236\t|\tGammaproteobacteria\t|\t\t|\tscientific name\t|
    91347\t|\tEnterobacterales\t|\t\t|\tscientific name\t|
    131567\t|\tbiota\t|\t\t|\tsynonym\t|
    131567\t|\tcellular organisms\t|\t\t|\tscientific name\t|";

        let dir = tempdir().unwrap();
        let path = dir.path();
        let mut nodes_file = std::fs::File::create(path.join(NODES_FILENAME)).unwrap();
        writeln!(nodes_file, "{}", nodes).unwrap();
        let mut names_file = std::fs::File::create(path.join(NAMES_FILENAME)).unwrap();
        writeln!(names_file, "{}", names).unwrap();

        let tax = load(path).unwrap();
        assert_eq!(
            Taxonomy::<&str>::name(&tax, "562").unwrap(),
            "Escherichia coli"
        );
        assert_eq!(
            Taxonomy::<&str>::rank(&tax, "562").unwrap(),
            TaxRank::Species
        );
        assert_eq!(
            Taxonomy::<&str>::children(&tax, "561").unwrap(),
            vec!["562"]
        );
        assert_eq!(
            Taxonomy::<&str>::parent(&tax, "562").unwrap(),
            Some(("561", 1.))
        );

        // Check that node data fields are loaded with proper types
        let data_562 = Taxonomy::<&str>::data(&tax, "562").unwrap();
        assert_eq!(
            data_562.get("genetic_code_id").and_then(|v| v.as_i64()),
            Some(11)
        );
        // Check boolean flags are parsed correctly
        assert_eq!(
            data_562
                .get("genbank_hidden_flag")
                .and_then(|v| v.as_bool()),
            Some(true)
        );

        // Check that alternative names are stored
        assert!(data_562.contains_key("names"));
        let names_array = data_562.get("names").and_then(|v| v.as_array()).unwrap();
        let common_name = names_array
            .iter()
            .find_map(|v| {
                let name = TaxonomyName::from_taxonomy_value(v)?;
                if name.name_class == "common name" {
                    Some(name.name_text)
                } else {
                    None
                }
            })
            .unwrap();
        assert_eq!(common_name, "E. coli");

        let out = path.join("out");
        save::<&str, _, _>(&tax, &out, true).unwrap();

        // now load again and validate a few taxids
        let tax2 = load(&out).unwrap();

        // Check E. coli (562)
        assert_eq!(
            Taxonomy::<&str>::name(&tax2, "562").unwrap(),
            "Escherichia coli"
        );
        assert_eq!(
            Taxonomy::<&str>::rank(&tax2, "562").unwrap(),
            TaxRank::Species
        );
        assert_eq!(
            Taxonomy::<&str>::parent(&tax2, "562").unwrap(),
            Some(("561", 1.))
        );

        // Check Escherichia (561)
        assert_eq!(Taxonomy::<&str>::name(&tax2, "561").unwrap(), "Escherichia");
        assert_eq!(
            Taxonomy::<&str>::rank(&tax2, "561").unwrap(),
            TaxRank::Genus
        );
        assert_eq!(
            Taxonomy::<&str>::parent(&tax2, "561").unwrap(),
            Some(("543", 1.))
        );

        // Check root (1)
        assert_eq!(Taxonomy::<&str>::name(&tax2, "1").unwrap(), "root");

        // Check children relationship preserved
        assert_eq!(
            Taxonomy::<&str>::children(&tax2, "561").unwrap(),
            vec!["562"]
        );

        // Check that data fields are preserved through save/load cycle with proper types
        let data_562_after = Taxonomy::<&str>::data(&tax2, "562").unwrap();
        assert_eq!(
            data_562_after
                .get("genetic_code_id")
                .and_then(|v| v.as_i64()),
            Some(11)
        );

        // Check that alternative names are preserved
        assert!(data_562_after.contains_key("names"));
        let names_array_after = data_562_after
            .get("names")
            .and_then(|v| v.as_array())
            .unwrap();
        let common_name_after = names_array_after
            .iter()
            .find_map(|v| {
                let name = TaxonomyName::from_taxonomy_value(v)?;
                if name.name_class == "common name" {
                    Some(name.name_text)
                } else {
                    None
                }
            })
            .unwrap();
        assert_eq!(common_name_after, "E. coli");

        // Test that alternative names can be found with find_all_by_name
        let found_by_common = tax.find_all_by_name("E. coli");
        assert_eq!(found_by_common.len(), 1);
        assert_eq!(found_by_common[0], "562");

        let found_by_scientific = tax.find_all_by_name("Escherichia coli");
        assert_eq!(found_by_scientific.len(), 1);
        assert_eq!(found_by_scientific[0], "562");

        // Test finding by synonym
        let found_by_synonym = tax.find_all_by_name("Bacillus coli");
        assert_eq!(found_by_synonym.len(), 1);
        assert_eq!(found_by_synonym[0], "562");
    }
}
