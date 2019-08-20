use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::iter::Sum;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use crate::base::GeneralTaxonomy;
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy;
use crate::{Result, TaxonomyError};

// TODO: we should allow read and write methods that take a nodes.dmp file
// followed by a names.dmp file (separated by a `--` line?) to allow piping
// NCBI taxonomies in and out

// TODO: add save_ncbi_writers fn?

pub fn load_ncbi_files<P>(ncbi_directory: P) -> Result<GeneralTaxonomy>
where
    P: AsRef<Path>,
{
    let mut nodes_path = PathBuf::from(ncbi_directory.as_ref());
    let mut names_path = nodes_path.clone();
    nodes_path.push("nodes.dmp");
    names_path.push("names.dmp");

    let nodes_file = File::open(nodes_path)?;
    let names_file = File::open(names_path)?;

    load_ncbi(nodes_file, names_file)
}

pub fn load_ncbi<R>(node_reader: R, name_reader: R) -> Result<GeneralTaxonomy>
where
    R: Read,
{
    let nodes_buf = BufReader::new(node_reader);
    let mut tax_ids: Vec<String> = Vec::new();
    let mut parents: Vec<String> = Vec::new();
    let mut ranks: Vec<Option<TaxRank>> = Vec::new();
    let mut tax_to_idx: HashMap<String, usize> = HashMap::new();

    for (ix, line) in nodes_buf.lines().enumerate() {
        let mut fields: Vec<String> = line?.split("\t|\t").map(|x| x.to_string()).collect();
        if fields.len() < 10 {
            // should be at least 14
            let msg = if ix == 0 {
                "Not enough fields; perhaps names and nodes files are switched?"
            } else {
                "Not enough fields; nodes file is bad?"
            };
            return Err(TaxonomyError::ImportError {
                file: "nodes.dmp".to_string(),
                line: ix,
                msg: msg.to_string(),
            }
            .into());
        }
        let tax_id = fields.remove(0).to_string();
        let parent_tax_id = fields.remove(0);
        let rank = fields.remove(0);
        // let embl_code = fields[3];
        // let gencode = fields[6];
        // let mito_gencode = fields[8];

        tax_ids.push(tax_id.clone());
        parents.push(parent_tax_id.to_string());
        ranks.push(TaxRank::from_str(&rank).ok());
        tax_to_idx.insert(tax_id, ix);
    }

    // TODO: fixme? this fails if we have unmapped parent nodes (i.e. file is truncated?)
    let parent_ids: Result<Vec<usize>> = parents
        .iter()
        .enumerate()
        .map(|(ix, x)| {
            tax_to_idx.get(&*x.as_str()).copied().ok_or_else(|| {
                TaxonomyError::ImportError {
                    file: "nodes.dmp".to_string(),
                    line: ix + 1,
                    msg: format!("Parent ID {} could not be found", x),
                }
                .into()
            })
        })
        .collect();

    let parent_ids = parent_ids?;

    let names_buf = BufReader::new(name_reader);
    let mut names: Vec<String> = vec!["".to_string(); tax_ids.len()];
    for (ix, line) in names_buf.lines().enumerate() {
        let mut fields: Vec<String> = line?.split("\t|\t").map(|x| x.to_string()).collect();
        if fields.len() > 10 {
            // should only be 5
            return Err(TaxonomyError::ImportError {
                file: "names.dmp".to_string(),
                line: ix,
                msg: "Too many fields?".to_string(),
            }
            .into());
        }
        let tax_id = fields.remove(0);
        let name = fields.remove(0);
        let name_class = fields.remove(1);
        if name_class.starts_with("scientific name") {
            let name = name.to_string();
            names[tax_to_idx[&*tax_id]] = name;
        }
    }

    Ok(GeneralTaxonomy::new(
        tax_ids,
        parent_ids,
        Some(names),
        Some(ranks),
        None,
    ))
}

// TODO: add root_node
pub fn save_ncbi_files<'t, P, T: 't, D: 't, X>(tax: &'t X, name_file: P, node_file: P) -> Result<()>
where
    P: AsRef<Path>,
    T: Clone + Debug + Display + PartialEq,
    D: Debug + PartialOrd + Sum,
    X: Taxonomy<'t, T, D>,
{
    let nodef = File::create(node_file)?;
    let mut node_writer = BufWriter::new(nodef);
    let namef = File::create(name_file)?;
    let mut name_writer = BufWriter::new(namef);
    for key in tax.traverse(tax.root())?.filter(|x| x.1).map(|x| x.0) {
        let name = tax.name(key.clone())?;
        let rank = tax.rank(key.clone())?;
        name_writer
            .write_all(format!("{}\t|\t{}\t|\tscientific name\t|", &key, name).as_bytes())?;
        node_writer.write_all(
            format!(
                "{}\t|\t{}\t|\t{}\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|",
                &key,
                tax.parent(key.clone())?
                    .map(|x| format!("{}", x.0))
                    .unwrap_or_else(|| "".to_string()),
                rank.map(|x| x.to_ncbi_rank()).unwrap_or(""),
            )
            .as_bytes(),
        )?;
    }

    Ok(())
}

#[test]
fn test_ncbi_importer() {
    use std::io::Cursor;

    let nodes = "1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|
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
    let tax = load_ncbi(Cursor::new(nodes), Cursor::new(names)).unwrap();
    assert_eq!(tax.name("562").unwrap(), "Escherichia coli");
    assert_eq!(tax.rank("562").unwrap(), Some(TaxRank::Species));
    assert_eq!(tax.parent("562").unwrap(), Some(("561", 1.)));
    assert_eq!(tax.children("561").unwrap(), vec!["562"]);
}
