use std::collections::HashMap;
use std::fmt;
use std::io::{Read, Write};
use std::str::FromStr;

use serde::{de, Deserialize, Deserializer, Serialize};
use serde_json::{json, to_value, to_writer, Value};

use crate::base::GeneralTaxonomy;
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;
use crate::Taxonomy;

#[derive(Eq, PartialEq)]
pub enum JsonFormat {
    NodeLink,
    Tree,
}

// Sometimes the tax ID might be an integer but we only want strings
struct DeserializeU64OrStringVisitor;

impl<'de> de::Visitor<'de> for DeserializeU64OrStringVisitor {
    type Value = String;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("an integer or a string")
    }

    fn visit_u64<E>(self, v: u64) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(format!("{}", v))
    }

    fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(v.to_string())
    }
}

fn deserialize_u64_or_string<'de, D>(deserializer: D) -> Result<String, D::Error>
where
    D: Deserializer<'de>,
{
    deserializer.deserialize_any(DeserializeU64OrStringVisitor)
}

fn deserialize_tax_rank<'de, D>(deserializer: D) -> Result<TaxRank, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = Deserialize::deserialize(deserializer)?;
    if s.is_empty() {
        return Ok(TaxRank::Unspecified);
    }
    TaxRank::from_str(&s).map_err(de::Error::custom)
}

fn default_tax_rank() -> TaxRank {
    TaxRank::Unspecified
}

#[derive(Debug, Serialize, Deserialize)]
struct TaxNode {
    #[serde(deserialize_with = "deserialize_u64_or_string")]
    id: String,
    name: String,
    #[serde(deserialize_with = "deserialize_tax_rank")]
    #[serde(default = "default_tax_rank")]
    rank: TaxRank,

    #[serde(flatten)]
    extra: HashMap<String, Value>,
}

#[derive(Debug, Serialize, Deserialize)]
struct Link {
    source: usize,
    target: usize,
}

/// Requires `nodes` to be an array of object with at least {rank, name, id}
fn load_node_link_json(tax_json: &Value) -> TaxonomyResult<GeneralTaxonomy> {
    let json_tax_nodes = tax_json["nodes"]
        .as_array()
        .ok_or_else(|| {
            Error::new(ErrorKind::ImportError {
                line: 0,
                msg: "'nodes' not in JSON".to_owned(),
            })
        })?
        .clone();

    let mut tax_nodes = Vec::with_capacity(json_tax_nodes.len());
    for n in json_tax_nodes {
        let node: TaxNode = serde_json::from_value(n)?;
        if node.id.parse::<usize>().is_err() {
            return Err(Error::new(ErrorKind::ImportError {
                line: 0,
                msg: format!("Tax ID {} cannot be converted to an integer", node.id),
            }));
        }
        tax_nodes.push(node);
    }
    // TODO: is that needed?
    tax_nodes.sort_unstable_by(|a, b| {
        let tax_id_a = a.id.parse::<usize>().unwrap();
        let tax_id_b = b.id.parse::<usize>().unwrap();
        tax_id_a.cmp(&tax_id_b)
    });

    let tax_links = tax_json["links"]
        .as_array()
        .ok_or_else(|| {
            Error::new(ErrorKind::ImportError {
                line: 0,
                msg: "'links' not in JSON".to_owned(),
            })
        })?
        .clone();

    let num_nodes = tax_nodes.len();
    let mut tax_ids: Vec<String> = Vec::with_capacity(tax_nodes.len());
    let mut names: Vec<String> = Vec::with_capacity(tax_nodes.len());
    let mut ranks: Vec<TaxRank> = Vec::with_capacity(tax_nodes.len());
    let mut data: Vec<HashMap<String, Value>> = Vec::with_capacity(tax_nodes.len());
    let mut parent_ids = vec![0; tax_nodes.len()];

    for node in tax_nodes {
        tax_ids.push(node.id);
        names.push(node.name);
        ranks.push(node.rank);
        data.push(node.extra);
    }

    for l in tax_links {
        let link: Link = serde_json::from_value(l)?;
        if link.source >= num_nodes {
            return Err(Error::new(ErrorKind::ImportError {
                line: 0,
                msg: format!(
                    "JSON source index {} specified, but there are only {} nodes",
                    link.source, num_nodes
                ),
            }));
        }
        if link.target >= num_nodes {
            return Err(Error::new(ErrorKind::ImportError {
                line: 0,
                msg: format!(
                    "JSON target index {} specified, but there are only {} nodes",
                    link.target, num_nodes
                ),
            }));
        }
        parent_ids[link.source] = link.target;
    }

    GeneralTaxonomy::from_arrays(
        tax_ids,
        parent_ids,
        Some(names),
        Some(ranks),
        None,
        Some(data),
    )
}

#[derive(Debug, Serialize, Deserialize)]
struct TaxNodeTree {
    #[serde(deserialize_with = "deserialize_u64_or_string")]
    id: String,
    name: String,
    #[serde(deserialize_with = "deserialize_tax_rank")]
    #[serde(default = "default_tax_rank")]
    rank: TaxRank,
    children: Vec<TaxNodeTree>,

    #[serde(flatten)]
    extra: HashMap<String, Value>,
}

fn load_tree_json(tax_json: &Value) -> TaxonomyResult<GeneralTaxonomy> {
    let root_node: TaxNodeTree = serde_json::from_value(tax_json.clone())?;

    fn add_node(
        parent_loc: usize,
        node: TaxNodeTree,
        tax_ids: &mut Vec<String>,
        parent_ids: &mut Vec<usize>,
        names: &mut Vec<String>,
        ranks: &mut Vec<TaxRank>,
        data: &mut Vec<HashMap<String, Value>>,
    ) -> TaxonomyResult<()> {
        tax_ids.push(node.id);
        parent_ids.push(parent_loc);
        names.push(node.name);
        ranks.push(node.rank);
        data.push(node.extra);
        let loc = tax_ids.len() - 1;
        for child in node.children {
            add_node(loc, child, tax_ids, parent_ids, names, ranks, data)?;
        }

        Ok(())
    }

    let mut tax_ids = Vec::new();
    let mut parent_ids = Vec::new();
    let mut names = Vec::new();
    let mut ranks = Vec::new();
    let mut data = Vec::new();

    add_node(
        0,
        root_node,
        &mut tax_ids,
        &mut parent_ids,
        &mut names,
        &mut ranks,
        &mut data,
    )?;

    GeneralTaxonomy::from_arrays(
        tax_ids,
        parent_ids,
        Some(names),
        Some(ranks),
        None,
        Some(data),
    )
}

pub fn load<R: Read>(reader: R, json_pointer: Option<&str>) -> TaxonomyResult<GeneralTaxonomy> {
    let tax_json: Value = serde_json::from_reader(reader)?;
    let actual_tax_json = if let Some(p) = json_pointer {
        tax_json.pointer(p).ok_or_else(|| {
            Error::new(ErrorKind::ImportError {
                line: 0,
                msg: format!("JSON path {} does not correspond to a value", p),
            })
        })?
    } else {
        &tax_json
    };

    // determine the JSON type
    if actual_tax_json.get("nodes").is_some() {
        load_node_link_json(actual_tax_json)
    } else {
        load_tree_json(actual_tax_json)
    }
}

pub fn save<W: Write>(
    writer: W,
    taxonomy: &GeneralTaxonomy,
    format: JsonFormat,
    root_node: Option<&str>,
) -> TaxonomyResult<()> {
    let root_node = root_node.unwrap_or_else(|| taxonomy.root());
    let json_data = match format {
        JsonFormat::NodeLink => serialize_as_node_links(taxonomy, root_node),
        JsonFormat::Tree => serialize_as_tree(taxonomy, root_node),
    }?;
    to_writer(writer, &json_data)?;

    Ok(())
}

fn serialize_as_tree(taxonomy: &GeneralTaxonomy, tax_id: &str) -> TaxonomyResult<Value> {
    fn inner(tax: &GeneralTaxonomy, tax_id: &str) -> TaxonomyResult<TaxNodeTree> {
        let mut children = Vec::new();
        for child in tax.children(tax_id)? {
            children.push(inner(tax, child)?);
        }
        let node = TaxNodeTree {
            id: tax_id.to_owned(),
            name: tax.name(tax_id)?.to_owned(),
            rank: tax.rank(tax_id)?,
            children,
            extra: tax.data(tax_id)?.clone(),
        };
        Ok(node)
    }

    Ok(to_value(inner(taxonomy, tax_id)?)?)
}

fn serialize_as_node_links(taxonomy: &GeneralTaxonomy, root_id: &str) -> TaxonomyResult<Value> {
    let mut nodes = Vec::new();
    let mut links = Vec::new();
    let mut id_to_idx = HashMap::new();

    for (ix, (tid, _pre)) in taxonomy.traverse(root_id)?.filter(|x| x.1).enumerate() {
        let node = TaxNode {
            id: tid.to_owned(),
            name: taxonomy.name(tid)?.to_owned(),
            rank: taxonomy.rank(tid)?,
            extra: taxonomy.data(tid)?.clone(),
        };
        nodes.push(to_value(&node).unwrap());
        id_to_idx.insert(tid, ix);
        if let Some(parent_id) = taxonomy.parent(tid)? {
            links.push(json!({
                "source": ix,
                "target": id_to_idx[&parent_id.0],
            }));
        }
    }

    let tax_json = json!({
        "nodes": nodes,
        "links": links,
        "directed": true,
        "multigraph": false,
        "graph": [],
    });

    Ok(tax_json)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::taxonomy::Taxonomy;
    use serde_json::from_str;
    use std::io::Cursor;

    #[test]
    fn can_load_empty_node_link_format() {
        let example = r#"{"nodes": [], "links": []}"#;
        let tax = load(Cursor::new(example), None).unwrap();
        assert_eq!(tax.len(), 0);
    }

    #[derive(Debug, Serialize, Deserialize)]
    struct NodeLinkFormat {
        nodes: Vec<TaxNode>,
        links: Vec<Link>,
    }

    #[test]
    fn can_load_small_node_link_format() {
        let example = r#"{
            "nodes": [
                {"id": "1", "name": "root", "readcount": 1000},
                {"id": "2", "name": "Bacteria", "rank": "no rank", "readcount": 1000},
                {"id": "562", "name": "Escherichia coli", "rank": "species", "readcount": 1000}
            ],
            "links": [
                {"source": 1, "target": 0},
                {"source": 2, "target": 1}
            ]
        }"#;
        let tax = load(Cursor::new(example), None).unwrap();
        assert_eq!(tax.len(), 3);
        assert_eq!(tax.root(), "1");
        assert_eq!(tax.children("1").unwrap(), vec!["2"]);
        assert_eq!(tax.lineage("562").unwrap(), vec!["562", "2", "1"]);
        assert_eq!(tax.data("1").unwrap()["readcount"].as_i64().unwrap(), 1000);
        // Converting to Value should give us the same data
        let serialized = serialize_as_node_links(&tax, tax.root()).unwrap();
        let serialized_input: NodeLinkFormat = from_str(example).unwrap();
        assert_eq!(
            serialized.as_object().unwrap()["nodes"],
            to_value(serialized_input.nodes).unwrap()
        );
        assert_eq!(
            serialized.as_object().unwrap()["links"],
            to_value(serialized_input.links).unwrap()
        );
    }

    #[test]
    fn can_load_tree_format() {
        let example = r#"{
            "id": "1",
            "name": "root",
            "rank": "no rank",
            "readcount": 1000,
            "children": [
                {
                    "id": "2",
                    "name": "Bacteria",
                    "rank": "no rank",
                    "readcount": 2000,
                    "children": [
                        {
                            "id": "562",
                            "name": "Escherichia coli",
                            "rank": "species",
                            "readcount": 3000,
                            "children": []
                        }
                    ]
                }
            ]
        }"#;

        let tax = load(Cursor::new(example), None).unwrap();
        assert_eq!(tax.len(), 3);
        assert_eq!(tax.root(), "1");
        assert_eq!(tax.children("1").unwrap(), vec!["2"]);
        assert_eq!(tax.lineage("562").unwrap(), vec!["562", "2", "1"]);
        assert_eq!(tax.data("1").unwrap()["readcount"].as_i64().unwrap(), 1000);
        // Converting to Value should give us the same data
        let serialized = serialize_as_tree(&tax, tax.root()).unwrap();
        let input_val: TaxNodeTree = from_str(&example).unwrap();
        assert_eq!(serialized, to_value(&input_val).unwrap());
    }

    #[test]
    fn errors_on_invalid_tree_format() {
        // no id at root
        let example = r#"{}"#;
        assert!(load(Cursor::new(example), None).is_err());

        // rank as a number
        let example = r#"{"id": "1", "rank": 5}"#;
        assert!(load(Cursor::new(example), None).is_err());
    }

    #[test]
    fn can_load_node_link_and_fix_order_matter() {
        let example = r#"
{"multigraph": false, "directed": true, "graph": [], "nodes": [{"name": "Bacillus subtilis subsp. subtilis", "rank": "subspecies", "id": 135461}, {"name": "Staphylococcus aureus", "rank": "species", "id": 1280}, {"name": "Mycobacterium", "rank": "genus", "id": 1763}, {"name": "Klebsiella", "rank": "genus", "id": 570}, {"name": "Bacillus anthracis str. 'Ames Ancestor'", "rank": "strain", "id": 261594}, {"name": "Deinococcaceae", "rank": "family", "id": 183710}, {"name": "Escherichia coli str. K-12 substr. MG1655", "rank": "no rank", "id": 511145}, {"name": "Escherichia coli K-12", "rank": "strain", "id": 83333}, {"name": "Mycobacterium tuberculosis H37Rv", "rank": "strain", "id": 83332}, {"name": "Salmonella enterica subsp. enterica serovar Typhimurium str. LT2", "rank": "strain", "id": 99287}, {"name": "Bacillus subtilis group", "rank": "species group", "id": 653685}, {"name": "root", "rank": "no rank", "id": 1}, {"name": "cellular organisms", "rank": "no rank", "id": 131567}, {"name": "Klebsiella pneumoniae subsp. pneumoniae", "rank": "subspecies", "id": 72407}, {"name": "Deinococci", "rank": "class", "id": 188787}, {"name": "Terrabacteria group", "rank": "clade", "id": 1783272}, {"name": "Klebsiella pneumoniae subsp. pneumoniae HS11286", "rank": "strain", "id": 1125630}, {"name": "Actinomycetia", "rank": "class", "id": 1760}, {"name": "Salmonella enterica subsp. enterica", "rank": "subspecies", "id": 59201}, {"name": "Bacillus anthracis", "rank": "species", "id": 1392}, {"name": "Bacillus", "rank": "genus", "id": 1386}, {"name": "Wolbachia", "rank": "genus", "id": 953}, {"name": "Deinococcus", "rank": "genus", "id": 1298}, {"name": "Firmicutes", "rank": "phylum", "id": 1239}, {"name": "Escherichia", "rank": "genus", "id": 561}, {"name": "Proteobacteria", "rank": "phylum", "id": 1224}, {"name": "Vibrio cholerae MS6", "rank": "strain", "id": 1420885}, {"name": "Escherichia coli", "rank": "species", "id": 562}, {"name": "Staphylococcus", "rank": "genus", "id": 1279}, {"name": "Wolbachieae", "rank": "tribe", "id": 952}, {"name": "Vibrionales", "rank": "order", "id": 135623}, {"name": "Salmonella enterica", "rank": "species", "id": 28901}, {"name": "Vibrio", "rank": "genus", "id": 662}, {"name": "Bacillus cereus group", "rank": "species group", "id": 86661}, {"name": "Mycobacteriaceae", "rank": "family", "id": 1762}, {"name": "Enterobacteriaceae", "rank": "family", "id": 543}, {"name": "Bacillaceae", "rank": "family", "id": 186817}, {"name": "Bacillus subtilis", "rank": "species", "id": 1423}, {"name": "Vibrionaceae", "rank": "family", "id": 641}, {"name": "Enterobacterales", "rank": "order", "id": 91347}, {"name": "Salmonella", "rank": "genus", "id": 590}, {"name": "Wolbachia pipientis", "rank": "species", "id": 955}, {"name": "Vibrio cholerae", "rank": "species", "id": 666}, {"name": "Anaplasmataceae", "rank": "family", "id": 942}, {"name": "Bacillus subtilis subsp. subtilis str. 168", "rank": "strain", "id": 224308}, {"name": "Deinococcus-Thermus", "rank": "phylum", "id": 1297}, {"name": "Klebsiella/Raoultella group", "rank": "no rank", "id": 2890311}, {"name": "Gammaproteobacteria", "rank": "class", "id": 1236}, {"name": "Deinococcus radiodurans", "rank": "species", "id": 1299}, {"name": "Actinobacteria", "rank": "phylum", "id": 201174}, {"name": "Staphylococcaceae", "rank": "family", "id": 90964}, {"name": "Klebsiella pneumoniae", "rank": "species", "id": 573}, {"name": "Deinococcales", "rank": "order", "id": 118964}, {"name": "Mycobacterium tuberculosis", "rank": "species", "id": 1773}, {"name": "Salmonella enterica subsp. enterica serovar Typhimurium", "rank": "no rank", "id": 90371}, {"name": "Bacillales", "rank": "order", "id": 1385}, {"name": "Alphaproteobacteria", "rank": "class", "id": 28211}, {"name": "Corynebacteriales", "rank": "order", "id": 85007}, {"name": "Mycobacterium tuberculosis complex", "rank": "species group", "id": 77643}, {"name": "Bacteria", "rank": "superkingdom", "id": 2}, {"name": "Staphylococcus aureus subsp. aureus NCTC 8325", "rank": "strain", "id": 93061}, {"name": "Bacilli", "rank": "class", "id": 91061}, {"name": "Rickettsiales", "rank": "order", "id": 766}], "links": [{"source": 59, "target": 12}, {"source": 35, "target": 39}, {"source": 24, "target": 35}, {"source": 27, "target": 24}, {"source": 3, "target": 46}, {"source": 51, "target": 3}, {"source": 40, "target": 35}, {"source": 38, "target": 30}, {"source": 32, "target": 38}, {"source": 42, "target": 32}, {"source": 62, "target": 56}, {"source": 43, "target": 62}, {"source": 29, "target": 43}, {"source": 21, "target": 29}, {"source": 41, "target": 21}, {"source": 25, "target": 59}, {"source": 47, "target": 25}, {"source": 23, "target": 15}, {"source": 28, "target": 50}, {"source": 1, "target": 28}, {"source": 45, "target": 15}, {"source": 22, "target": 5}, {"source": 48, "target": 22}, {"source": 55, "target": 61}, {"source": 20, "target": 36}, {"source": 19, "target": 33}, {"source": 37, "target": 10}, {"source": 17, "target": 49}, {"source": 34, "target": 57}, {"source": 2, "target": 34}, {"source": 53, "target": 58}, {"source": 56, "target": 25}, {"source": 31, "target": 40}, {"source": 18, "target": 31}, {"source": 13, "target": 51}, {"source": 58, "target": 2}, {"source": 8, "target": 53}, {"source": 7, "target": 27}, {"source": 57, "target": 17}, {"source": 33, "target": 20}, {"source": 54, "target": 18}, {"source": 50, "target": 55}, {"source": 61, "target": 23}, {"source": 39, "target": 47}, {"source": 60, "target": 1}, {"source": 9, "target": 54}, {"source": 52, "target": 14}, {"source": 12, "target": 11}, {"source": 0, "target": 37}, {"source": 30, "target": 47}, {"source": 5, "target": 52}, {"source": 36, "target": 55}, {"source": 14, "target": 45}, {"source": 49, "target": 15}, {"source": 44, "target": 0}, {"source": 4, "target": 19}, {"source": 6, "target": 7}, {"source": 10, "target": 20}, {"source": 16, "target": 13}, {"source": 26, "target": 42}, {"source": 15, "target": 59}, {"source": 46, "target": 35}]}
        "#;
        let tax = load(Cursor::new(example), None).unwrap();
        assert_eq!(tax.len(), 63);
        assert_eq!(tax.root(), "1");
        assert_eq!(tax.lineage("562").unwrap(), vec!["562"]);
    }

    #[test]
    fn can_load_from_json_path() {
        let example = r#"{"test": {"sub": {"nodes": [], "links": []}}}"#;
        let tax = load(Cursor::new(example), Some("/test/sub")).unwrap();
        assert_eq!(tax.len(), 0);
    }
}
