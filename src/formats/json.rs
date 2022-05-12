use std::collections::HashMap;
use std::fmt;
use std::io::{Read, Write};
use std::str::FromStr;

use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
use serde_json::{json, to_value, to_writer, Value};

use crate::base::GeneralTaxonomy;
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;
use crate::Taxonomy;

#[derive(Eq, PartialEq)]
pub enum JsonFormat {
    /// The node link format is made of 2 arrays:
    /// 1. the nodes with their taxonomic info, in theory sorted by tax id
    /// 2. the links between each node: a `source` node has a `parent` node.
    /// The nodes are represented by indices in the nodes array
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
    let opt: Option<String> = Option::deserialize(deserializer)?;
    if let Some(s) = opt {
        if s.is_empty() {
            return Ok(TaxRank::Unspecified);
        }
        TaxRank::from_str(&s).map_err(de::Error::custom)
    } else {
        Ok(TaxRank::Unspecified)
    }

}

fn serialize_tax_rank<S>(x: &TaxRank, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    s.serialize_str(x.to_ncbi_rank())
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
    #[serde(serialize_with = "serialize_tax_rank")]
    #[serde(default = "default_tax_rank")]
    rank: TaxRank,

    /// We want to keep any other fields that was in JSON, they all get put in this hashmap
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

    // The JSON import might be messed up so we sort the nodes by ids
    // but since the links refer to their position, we need to first record them
    #[allow(clippy::needless_collect)]
    let initial_positions: Vec<_> = tax_nodes.iter().map(|n| n.id.clone()).collect();
    tax_nodes.sort_unstable_by(|a, b| {
        let tax_id_a = a.id.parse::<usize>().unwrap();
        let tax_id_b = b.id.parse::<usize>().unwrap();
        tax_id_a.cmp(&tax_id_b)
    });
    let after_positions: Vec<_> = tax_nodes.iter().map(|n| &n.id).collect();

    let mut idx_mapping_after_sorting = HashMap::new();
    for (i, id) in initial_positions.into_iter().enumerate() {
        idx_mapping_after_sorting
            .insert(i, after_positions.iter().position(|p| *p == &id).unwrap());
    }

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

        // Then we need to fix up the parents if needed since we might have re-ordered
        // content
        let fixed_source = idx_mapping_after_sorting[&link.source];
        let fixed_target = idx_mapping_after_sorting[&link.target];

        parent_ids[fixed_source] = fixed_target;
        // parent should always be higher than child in the tree and since we might
        // have moved things around we need to ensure it stays the case
        if fixed_target > fixed_source {
            parent_ids[fixed_target] = fixed_source;
        } else {
            parent_ids[fixed_source] = fixed_target;
        }
        parent_ids[idx_mapping_after_sorting[&link.source]] =
            idx_mapping_after_sorting[&link.target];
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
    #[serde(serialize_with = "serialize_tax_rank")]
    #[serde(default = "default_tax_rank")]
    rank: TaxRank,
    #[serde(default)]
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

        let tax2 = load_node_link_json(&serialized).unwrap();
        assert_eq!(tax2.len(), 3);
        assert_eq!(tax2.root(), "1");
        assert_eq!(tax2.children("1").unwrap(), vec!["2"]);
        assert_eq!(tax2.lineage("562").unwrap(), vec!["562", "2", "1"]);
        assert_eq!(tax2.data("1").unwrap()["readcount"].as_i64().unwrap(), 1000);
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
                            "readcount": 3000
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
        let example = r#"{"id": "1", "rank": 5, "name": ""}"#;
        assert!(load(Cursor::new(example), None).is_err());
    }

    #[test]
    fn can_load_node_link_and_fix_order_them() {
        let example = r#"
{
  "multigraph": false,
  "directed": true,
  "graph": [],
  "nodes": [
    {
      "name": "root",
      "rank": "no rank",
      "id": 1
    },
    {
      "name": "genus 2",
      "rank": "genus",
      "id": 9
    },
    {
      "name": "superkingdom 1",
      "rank": "superkingdom",
      "id": 2
    },
    {
      "name": "species 2.1",
      "rank": "species",
      "id": 11
    },
    {
      "name": "genus 1",
      "rank": "genus",
      "id": 8
    },
    {
      "name": "class 1",
      "rank": "class",
      "id": 5
    },
    {
      "name": "kingdom 1",
      "rank": "kingdom",
      "id": 3
    },
    {
      "name": "phylum 1",
      "rank": "phylum",
      "id": 4
    },
    {
      "name": "order 1",
      "rank": "order",
      "id": 6
    },
    {
      "name": "family 1",
      "rank": "family",
      "id": 7
    },
    {
      "name": "species 1.1",
      "rank": "species",
      "id": 10
    }
  ],
  "links": [
    {
      "source": 3,
      "target": 1
    },
    {
      "source": 10,
      "target": 4
    },
    {
      "source": 1,
      "target": 9
    },
    {
      "source": 4,
      "target": 9
    },
    {
      "source": 9,
      "target": 8
    },
    {
      "source": 8,
      "target": 5
    },
    {
      "source": 5,
      "target": 7
    },
    {
      "source": 7,
      "target": 6
    },
    {
      "source": 6,
      "target": 2
    },
    {
      "source": 2,
      "target": 0
    }
  ]
}        "#;
        let tax = load(Cursor::new(example), None).unwrap();
        assert_eq!(tax.len(), 11);
        assert_eq!(tax.root(), "1");
        assert_eq!(
            tax.lineage("10").unwrap(),
            vec!["10", "8", "7", "6", "5", "4", "3", "2", "1"]
        );
    }

    #[test]
    fn can_load_from_json_path() {
        let example = r#"{"test": {"sub": {"nodes": [], "links": []}}}"#;
        let tax = load(Cursor::new(example), Some("/test/sub")).unwrap();
        assert_eq!(tax.len(), 0);
    }

    #[test]
    fn can_handle_null_ranks() {
        let example = r#"{"id": "1", "rank": null, "name": ""}"#;
        assert!(load(Cursor::new(example), None).is_ok());
    }
}
