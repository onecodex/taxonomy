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

pub fn load<P: AsRef<Path>>(ncbi_directory: P) -> TaxonomyResult<GeneralTaxonomy> {
    let dir = ncbi_directory.as_ref();
    let nodes_file = std::fs::File::open(dir.join(NODES_FILENAME))?;
    let names_file = std::fs::File::open(dir.join(NAMES_FILENAME))?;

    // First we go through the nodes
    let mut tax_ids: Vec<String> = Vec::new();
    let mut parents: Vec<String> = Vec::new();
    let mut ranks: Vec<TaxRank> = Vec::new();
    let mut tax_to_idx: HashMap<String, usize> = HashMap::new();

    for (ix, line) in BufReader::new(nodes_file).lines().enumerate() {
        let mut fields: Vec<String> = line?.split("\t|\t").map(|x| x.to_string()).collect();
        if fields.len() < 10 {
            // should be at least 14
            return Err(Error::new(ErrorKind::ImportError {
                line: ix,
                msg: "Not enough fields in nodes.dmp; bad line?".to_owned(),
            }));
        }
        let tax_id = fields.remove(0).trim().to_string();
        let parent_tax_id = fields.remove(0).trim().to_string();
        let rank = fields.remove(0);

        tax_ids.push(tax_id.clone());
        parents.push(parent_tax_id.to_string());
        ranks.push(TaxRank::from_str(&rank)?);
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

    // And then grab their names by their idx
    let mut names: Vec<String> = vec![String::new(); tax_ids.len()];
    for (ix, line) in BufReader::new(names_file).lines().enumerate() {
        let mut fields: Vec<String> = line?.split("\t|\t").map(|x| x.to_string()).collect();
        if fields.len() > 10 {
            // should only be 5
            return Err(Error::new(ErrorKind::ImportError {
                line: ix,
                msg: "Too many fields in names.dmp".to_owned(),
            }));
        }
        let tax_id = fields.remove(0).trim().to_string();
        let name = fields.remove(0).trim().to_string();
        let name_class = fields.remove(1);
        if name_class.starts_with("scientific name") {
            let name = name.to_string();
            names[tax_to_idx[&*tax_id]] = name;
        }
    }

    GeneralTaxonomy::from_arrays(tax_ids, parent_ids, Some(names), Some(ranks), None, None)
}

/// Only supports dumping with the &str impl as there seems to be no way to specify `T` other
/// than passing a dummy argument with `impl trait`?
pub fn save<'t, T: 't, P: AsRef<Path>, X: Taxonomy<'t, T>>(
    tax: &'t X,
    out_dir: P,
) -> TaxonomyResult<()>
where
    T: Clone + Debug + Display + PartialEq,
{
    let dir = out_dir.as_ref();
    std::fs::create_dir_all(&dir)?;
    let mut node_writer = BufWriter::new(std::fs::File::create(dir.join(NODES_FILENAME))?);
    let mut name_writer = BufWriter::new(std::fs::File::create(dir.join(NAMES_FILENAME))?);

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
                    .map(|(x, _)| format!("{}", x))
                    .unwrap_or_default(),
                rank.to_ncbi_rank(),
            )
            .as_bytes(),
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

        let out = path.join("out");
        save::<&str, _, _>(&tax, &out).unwrap();
    }
}
