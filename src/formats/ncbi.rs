use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

use crate::base::GeneralTaxonomy;
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy;

const NODES_FILENAME: &str = "nodes.dmp";
const NAMES_FILENAME: &str = "names.dmp";

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
    let mut data: Vec<HashMap<String, serde_json::Value>> = Vec::new();
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
        let rank = fields.remove(0).trim().to_string();
        let embl_code = fields.remove(0).trim().to_string();
        let division_id = fields.remove(0).trim().to_string();
        let inherited_div_flag = fields.remove(0).trim().to_string();
        let genetic_code_id = fields.remove(0).trim().to_string();
        let inherited_GC_flag = fields.remove(0).trim().to_string();
        let mitochondrial_genetic_code_id = fields.remove(0).trim().to_string();
        let inherited_MGC_flag = fields.remove(0).trim().to_string();
        let GenBank_hidden_flag = fields.remove(0).trim().to_string();
        let hidden_subtree_root_flag = fields.remove(0).trim().to_string();
        let comments = if !fields.is_empty() {
            fields.remove(0).trim().trim_end_matches('|').trim().to_string()
        } else {
            String::new()
        };

        tax_ids.push(tax_id.clone());
        parents.push(parent_tax_id.to_string());
        ranks.push(TaxRank::from_str(&rank)?);

        // Store all node fields in data HashMap
        let mut node_data = HashMap::new();
        if !embl_code.is_empty() {
            node_data.insert("embl_code".to_string(), serde_json::Value::String(embl_code));
        }
        if !division_id.is_empty() {
            node_data.insert("division_id".to_string(), serde_json::Value::String(division_id));
        }
        if !inherited_div_flag.is_empty() {
            node_data.insert("inherited_div_flag".to_string(), serde_json::Value::String(inherited_div_flag));
        }
        if !genetic_code_id.is_empty() {
            node_data.insert("genetic_code_id".to_string(), serde_json::Value::String(genetic_code_id));
        }
        if !inherited_GC_flag.is_empty() {
            node_data.insert("inherited_GC_flag".to_string(), serde_json::Value::String(inherited_GC_flag));
        }
        if !mitochondrial_genetic_code_id.is_empty() {
            node_data.insert("mitochondrial_genetic_code_id".to_string(), serde_json::Value::String(mitochondrial_genetic_code_id));
        }
        if !inherited_MGC_flag.is_empty() {
            node_data.insert("inherited_MGC_flag".to_string(), serde_json::Value::String(inherited_MGC_flag));
        }
        if !GenBank_hidden_flag.is_empty() {
            node_data.insert("GenBank_hidden_flag".to_string(), serde_json::Value::String(GenBank_hidden_flag));
        }
        if !hidden_subtree_root_flag.is_empty() {
            node_data.insert("hidden_subtree_root_flag".to_string(), serde_json::Value::String(hidden_subtree_root_flag));
        }
        if !comments.is_empty() {
            node_data.insert("comments".to_string(), serde_json::Value::String(comments));
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
        let name_class = fields.remove(0).trim().trim_end_matches('|').trim().to_string();

        if let Some(&idx) = tax_to_idx.get(&tax_id) {
            if name_class == "scientific name" {
                names[idx] = name_txt.clone();
            }

            // Store all name types in data
            let name_key = format!("name_{}", name_class.replace(" ", "_"));
            data[idx].insert(name_key, serde_json::Value::String(name_txt));

            if !unique_name.is_empty() {
                let unique_key = format!("unique_name_{}", name_class.replace(" ", "_"));
                data[idx].insert(unique_key, serde_json::Value::String(unique_name));
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
fn write_nodes_row(
    writer: &mut BufWriter<std::fs::File>,
    tax_id: &str,
    parent: &str,
    rank: &str,
    embl_code: &str,
    division_id: &str,
    inherited_div: &str,
    genetic_code: &str,
    inherited_gc: &str,
    mito_gc: &str,
    inherited_mgc: &str,
    genbank_hidden: &str,
    subtree_hidden: &str,
    comments: &str,
) -> std::io::Result<()> {
    write!(
        writer,
        "{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\n",
        tax_id, parent, rank, embl_code, division_id, inherited_div,
        genetic_code, inherited_gc, mito_gc, inherited_mgc,
        genbank_hidden, subtree_hidden, comments
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
    write!(
        writer,
        "{}\t|\t{}\t|\t{}\t|\t{}\t|\n",
        tax_id, name, unique_name, name_class
    )
}

pub fn save<'t, T: 't, P: AsRef<Path>, X: Taxonomy<'t, T>>(
    tax: &'t X,
    out_dir: P,
) -> TaxonomyResult<()>
where
    T: Clone + Debug + Display + PartialEq,
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

        // Extract data fields
        let node_data = tax.data(key.clone())?;
        let embl_code = node_data
            .get("embl_code")
            .and_then(|v| v.as_str())
            .unwrap_or("");
        let division_id = node_data
            .get("division_id")
            .and_then(|v| v.as_str())
            .unwrap_or("0");
        let inherited_div_flag = node_data
            .get("inherited_div_flag")
            .and_then(|v| v.as_str())
            .unwrap_or("0");
        let genetic_code_id = node_data
            .get("genetic_code_id")
            .and_then(|v| v.as_str())
            .unwrap_or("1");
        let inherited_GC_flag = node_data
            .get("inherited_GC_flag")
            .and_then(|v| v.as_str())
            .unwrap_or("0");
        let mitochondrial_genetic_code_id = node_data
            .get("mitochondrial_genetic_code_id")
            .and_then(|v| v.as_str())
            .unwrap_or("0");
        let inherited_MGC_flag = node_data
            .get("inherited_MGC_flag")
            .and_then(|v| v.as_str())
            .unwrap_or("0");
        let GenBank_hidden_flag = node_data
            .get("GenBank_hidden_flag")
            .and_then(|v| v.as_str())
            .unwrap_or("0");
        let hidden_subtree_root_flag = node_data
            .get("hidden_subtree_root_flag")
            .and_then(|v| v.as_str())
            .unwrap_or("0");
        let comments = node_data
            .get("comments")
            .and_then(|v| v.as_str())
            .unwrap_or("");

        // Write scientific name
        write_names_row(&mut name_writer, &format!("{}", &key), &name, "", "scientific name")?;

        // Write all alternative names from data
        for (data_key, value) in node_data.iter() {
            if data_key.starts_with("name_") && data_key != "name_scientific_name" {
                let name_class = data_key.strip_prefix("name_").unwrap().replace("_", " ");
                let name_txt = value.as_str().unwrap_or("");
                let unique_name_key = format!("unique_name_{}", data_key.strip_prefix("name_").unwrap());
                let unique_name = node_data
                    .get(&unique_name_key)
                    .and_then(|v| v.as_str())
                    .unwrap_or("");
                write_names_row(&mut name_writer, &format!("{}", &key), name_txt, unique_name, &name_class)?;
            }
        }

        // Write nodes.dmp entry with all fields
        write_nodes_row(
            &mut node_writer,
            &format!("{}", &key),
            &parent,
            rank.to_ncbi_rank(),
            embl_code,
            division_id,
            inherited_div_flag,
            genetic_code_id,
            inherited_GC_flag,
            mitochondrial_genetic_code_id,
            inherited_MGC_flag,
            GenBank_hidden_flag,
            hidden_subtree_root_flag,
            comments,
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

        // Check that node data fields are loaded
        let data_562 = Taxonomy::<&str>::data(&tax, "562").unwrap();
        assert_eq!(
            data_562.get("genetic_code_id").and_then(|v| v.as_str()),
            Some("11")
        );

        // Check that alternative names are stored
        assert!(data_562.contains_key("name_common_name"));
        assert_eq!(
            data_562.get("name_common_name").and_then(|v| v.as_str()),
            Some("E. coli")
        );

        let out = path.join("out");
        save::<&str, _, _>(&tax, &out).unwrap();

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

        // Check that data fields are preserved through save/load cycle
        let data_562_after = Taxonomy::<&str>::data(&tax2, "562").unwrap();
        assert_eq!(
            data_562_after.get("genetic_code_id").and_then(|v| v.as_str()),
            Some("11")
        );

        // Check that alternative names are preserved
        assert!(data_562_after.contains_key("name_common_name"));
        assert_eq!(
            data_562_after.get("name_common_name").and_then(|v| v.as_str()),
            Some("E. coli")
        );
    }
}
