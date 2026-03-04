use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::path::Path;
use std::str::FromStr;
use std::sync::Arc;

use arrow::array::{Array, Float32Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;

use crate::base::GeneralTaxonomy;
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy;

fn taxonomy_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("tax_id", DataType::Utf8, false),
        Field::new("parent_id", DataType::Utf8, false),
        Field::new("name", DataType::Utf8, false),
        Field::new("rank", DataType::Utf8, false),
        Field::new("parent_distance", DataType::Float32, false),
    ]))
}

fn missing_column_err(name: &str) -> Error {
    Error::new(ErrorKind::ImportError {
        line: 0,
        msg: format!("Missing '{}' column", name),
    })
}

fn get_string_column<'a>(batch: &'a RecordBatch, name: &str) -> TaxonomyResult<&'a StringArray> {
    let col = batch
        .column_by_name(name)
        .ok_or_else(|| missing_column_err(name))?;
    col.as_any().downcast_ref::<StringArray>().ok_or_else(|| {
        Error::new(ErrorKind::ImportError {
            line: 0,
            msg: format!("'{}' column must be a string type", name),
        })
    })
}

/// Loads a taxonomy from a Parquet file.
///
/// The file must contain the following columns:
/// - `tax_id`: string — the unique identifier for each node
/// - `parent_id`: string — the `tax_id` of the parent (same as `tax_id` for the root)
/// - `name`: string — the scientific name
/// - `rank`: string — the taxonomic rank (e.g. "species", "genus")
/// - `parent_distance`: float32 — branch length to the parent
pub fn load<P: AsRef<Path>>(path: P) -> TaxonomyResult<GeneralTaxonomy> {
    let file = std::fs::File::open(path.as_ref())?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let reader = builder.build()?;

    let mut tax_ids: Vec<String> = Vec::new();
    let mut parent_id_strings: Vec<String> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut ranks_str: Vec<String> = Vec::new();
    let mut parent_distances: Vec<f32> = Vec::new();

    for result in reader {
        let batch = result?;
        let n = batch.num_rows();

        let tax_id_col = get_string_column(&batch, "tax_id")?;
        let parent_id_col = get_string_column(&batch, "parent_id")?;
        let name_col = get_string_column(&batch, "name")?;
        let rank_col = get_string_column(&batch, "rank")?;
        let distance_col = batch
            .column_by_name("parent_distance")
            .ok_or_else(|| missing_column_err("parent_distance"))?
            .as_any()
            .downcast_ref::<Float32Array>()
            .ok_or_else(|| {
                Error::new(ErrorKind::ImportError {
                    line: 0,
                    msg: "'parent_distance' column must be Float32 type".to_owned(),
                })
            })?;

        for i in 0..n {
            tax_ids.push(tax_id_col.value(i).to_string());
            parent_id_strings.push(parent_id_col.value(i).to_string());
            names.push(name_col.value(i).to_string());
            ranks_str.push(rank_col.value(i).to_string());
            parent_distances.push(distance_col.value(i));
        }
    }

    let mut tax_to_idx: HashMap<String, usize> = HashMap::with_capacity(tax_ids.len());
    for (ix, tax_id) in tax_ids.iter().enumerate() {
        tax_to_idx.insert(tax_id.clone(), ix);
    }

    let mut parent_ids: Vec<usize> = Vec::with_capacity(parent_id_strings.len());
    for (i, parent) in parent_id_strings.iter().enumerate() {
        if let Some(idx) = tax_to_idx.get(parent) {
            parent_ids.push(*idx);
        } else {
            return Err(Error::new(ErrorKind::ImportError {
                line: i,
                msg: format!("Parent ID '{}' not found in tax_id column", parent),
            }));
        }
    }

    let ranks: Vec<TaxRank> = ranks_str
        .iter()
        .map(|r| TaxRank::from_str(r))
        .collect::<Result<Vec<_>, _>>()?;

    let gt = GeneralTaxonomy::from_arrays(
        tax_ids,
        parent_ids,
        Some(names),
        Some(ranks),
        Some(parent_distances),
        None,
    )?;
    gt.validate_uniqueness()?;
    Ok(gt)
}

/// Saves a taxonomy to a Parquet file.
///
/// The output file will contain the following columns:
/// - `tax_id`: string
/// - `parent_id`: string (same as `tax_id` for the root node)
/// - `name`: string
/// - `rank`: string
/// - `parent_distance`: float32
pub fn save<'t, T: 't, P: AsRef<Path>, X: Taxonomy<'t, T>>(
    tax: &'t X,
    path: P,
) -> TaxonomyResult<()>
where
    T: Clone + Debug + Display + Eq + Hash,
{
    let schema = taxonomy_schema();

    let mut tax_ids: Vec<String> = Vec::new();
    let mut parent_ids: Vec<String> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut ranks: Vec<String> = Vec::new();
    let mut distances: Vec<f32> = Vec::new();

    let root = tax.root();
    for key in tax.traverse(root.clone())?.filter(|x| x.1).map(|x| x.0) {
        let name = tax.name(key.clone())?;
        let rank = tax.rank(key.clone())?;
        let (parent_id, distance) = if key == root {
            (format!("{}", key), 1.0f32)
        } else {
            tax.parent(key.clone())?
                .map(|(p, d)| (format!("{}", p), d))
                .unwrap_or_else(|| (format!("{}", key), 1.0f32))
        };

        tax_ids.push(format!("{}", key));
        parent_ids.push(parent_id);
        names.push(name.to_string());
        ranks.push(rank.to_string());
        distances.push(distance);
    }

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(tax_ids)),
            Arc::new(StringArray::from(parent_ids)),
            Arc::new(StringArray::from(names)),
            Arc::new(StringArray::from(ranks)),
            Arc::new(Float32Array::from(distances)),
        ],
    )?;

    let file = std::fs::File::create(path.as_ref())?;
    let mut writer = ArrowWriter::try_new(file, schema, None)?;
    writer.write(&batch)?;
    writer.close()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::taxonomy::Taxonomy;
    use tempfile::tempdir;

    fn make_test_taxonomy() -> GeneralTaxonomy {
        GeneralTaxonomy::from_arrays(
            vec!["1".to_string(), "2".to_string(), "3".to_string()],
            vec![0, 0, 1],
            Some(vec![
                "root".to_string(),
                "Bacteria".to_string(),
                "E. coli".to_string(),
            ]),
            Some(vec![
                TaxRank::Unspecified,
                TaxRank::Superkingdom,
                TaxRank::Species,
            ]),
            Some(vec![1.0, 1.0, 1.0]),
            None,
        )
        .unwrap()
    }

    #[test]
    fn test_roundtrip() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("taxonomy.parquet");

        let original = make_test_taxonomy();
        save::<&str, _, _>(&original, &path).unwrap();

        let loaded = load(&path).unwrap();

        assert_eq!(Taxonomy::<&str>::name(&loaded, "1").unwrap(), "root");
        assert_eq!(Taxonomy::<&str>::name(&loaded, "2").unwrap(), "Bacteria");
        assert_eq!(Taxonomy::<&str>::name(&loaded, "3").unwrap(), "E. coli");
        assert_eq!(
            Taxonomy::<&str>::rank(&loaded, "3").unwrap(),
            TaxRank::Species
        );
        assert_eq!(Taxonomy::<&str>::children(&loaded, "1").unwrap(), vec!["2"]);
        assert_eq!(
            Taxonomy::<&str>::parent(&loaded, "3").unwrap(),
            Some(("2", 1.0))
        );
    }
}
