use std::collections::HashMap;
use std::fmt;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::io::{Read, Write};
use std::str::FromStr;

use crate::base::GeneralTaxonomy;
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;
use crate::Taxonomy;
use serde::{de, Deserialize, Deserializer, Serialize, Serializer};
use serde_json::{json, to_value, to_writer, Value};

/// We can handle 2 kinds of JSON formats:
/// 1. The Node Link format
/// 2. The Tree format
///
/// The Node Link Format:
///
/// The taxonomy is represented as a directed acyclic graph (DAG) in JSON using a node-link
/// structure:
/// ```json
/// {
///   "nodes": [
///     { "id": 0, "name": "root", "rank": "no rank" },
///     { "id": 1, "name": "superkingdom A", "rank": "superkingdom" },
///     { "id": 2, "name": "kingdom A", "rank": "kingdom" },
///     { "id": 3, "name": "phylum A", "rank": "phylum" },
///     { "id": 4, "name": "genus A", "rank": "genus" },
///     { "id": 5, "name": "species A1", "rank": "species" }
///   ],
///   "links": [
///     { "source": 5, "target": 4 },
///     { "source": 4, "target": 3 },
///     { "source": 3, "target": 2 },
///     { "source": 2, "target": 1 },
///     { "source": 1, "target": 0 }
///   ]
/// }
/// ```
///
/// nodes:
/// A list of taxonomic units. Each node has:
///   id – unique integer identifier
///   name – the scientific name or placeholder label
///   rank – the taxonomic rank (e.g., "species", "genus", "family")
///
/// links:
/// A list of directed edges between nodes. Each link has:
///   source – the child node’s index in the nodes array
///   target – the parent node’s index in the nodes array
///
/// The Tree Format
///
/// In the tree format, child nodes are nested within their parent node. This format is more
/// natural as the taxonomic structure is immediately apparent from reading it.
///
/// ```json
/// {
///     "id": "1",
///     "name": "root",
///     "rank": "no rank",
///     "children": [
///         {
///             "id": "2",
///             "name": "Bacteria",
///             "rank": "no rank",
///             "children": [
///                 {
///                     "id": "562",
///                     "name": "Escherichia coli",
///                     "rank": "species",
///                 }
///             ]
///         }
///     ]
/// }
/// ```
//
/// For both formats, you can add more data on each node object and these will be available after loading.
/// If a `rank` propery is present, it will be parsed as a NCBI rank.
#[derive(Eq, PartialEq)]
pub enum JsonFormat {
    /// The node link format is made of 2 arrays:
    /// 1. the nodes with their taxonomic info. Internal (integer) IDs are the node's position in
    ///    the nodes array
    /// 2. the links between each node: a `source` node has a `parent` node.
    /// The nodes are represented by indices in the nodes array
    /// Only use that format if you have existing taxonomies in that format.
    NodeLink,
    /// The preferred format
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
    let gt = if actual_tax_json.get("nodes").is_some() {
        load_node_link_json(actual_tax_json)
    } else {
        load_tree_json(actual_tax_json)
    }?;
    gt.validate_uniqueness()?;
    Ok(gt)
}

pub fn save<'t, W: Write, T: 't, X: Taxonomy<'t, T>>(
    writer: W,
    taxonomy: &'t X,
    format: JsonFormat,
    root_node: Option<T>,
) -> TaxonomyResult<()>
where
    T: Clone + Debug + Display + Eq + Hash + PartialEq,
{
    // Defaults to root but root exists only if taxonomy is not empty.
    let root_node = root_node.or(if taxonomy.is_empty() {
        None
    } else {
        Some(taxonomy.root())
    });

    let json_data = match format {
        JsonFormat::NodeLink => serialize_as_node_links(taxonomy, root_node),
        JsonFormat::Tree => serialize_as_tree(taxonomy, root_node),
    }?;
    to_writer(writer, &json_data)?;

    Ok(())
}

/// Memory-efficient streaming tree export that writes JSON incrementally
/// instead of building the entire tree structure in memory first.
pub fn save_tree_streaming<'t, W: Write, T: 't, X: Taxonomy<'t, T>>(
    mut writer: W,
    taxonomy: &'t X,
    root_node: Option<T>,
) -> TaxonomyResult<()>
where
    T: Clone + Debug + Display + Eq + Hash + PartialEq,
{
    let tax_id = root_node
        .or_else(|| {
            if taxonomy.is_empty() {
                None
            } else {
                Some(taxonomy.root())
            }
        })
        .ok_or(Error::new(ErrorKind::InvalidTaxonomy(
            "Taxonomy must have a root node.".to_string(),
        )))?;

    write_node_streaming(taxonomy, tax_id, &mut writer)?;
    writer.flush()?;
    Ok(())
}

fn write_node_streaming<'t, W: Write, T: 't>(
    tax: &'t impl Taxonomy<'t, T>,
    tax_id: T,
    writer: &mut W,
) -> TaxonomyResult<()>
where
    T: Clone + Debug + Display + Eq + Hash + PartialEq,
{
    write!(writer, "{{")?;

    // Write id
    write!(writer, "\"id\":")?;
    to_writer(&mut *writer, &tax_id.to_string())?;

    // Write name
    write!(writer, ",\"name\":")?;
    to_writer(&mut *writer, &tax.name(tax_id.clone())?)?;

    // Write rank
    write!(writer, ",\"rank\":")?;
    to_writer(&mut *writer, tax.rank(tax_id.clone())?.to_ncbi_rank())?;

    // Write extra data fields
    let data = tax.data(tax_id.clone())?;
    for (key, value) in data.iter() {
        write!(writer, ",{}", serde_json::to_string(key)?)?;
        write!(writer, ":")?;
        to_writer(&mut *writer, value)?;
    }

    // Write children (always include, even if empty, to match regular format)
    let children = tax.children(tax_id)?;
    write!(writer, ",\"children\":[")?;
    for (i, child) in children.into_iter().enumerate() {
        if i > 0 {
            write!(writer, ",")?;
        }
        write_node_streaming(tax, child, writer)?;
    }
    write!(writer, "]")?;

    write!(writer, "}}")?;
    Ok(())
}

fn serialize_as_tree<'t, T: 't>(
    taxonomy: &'t impl Taxonomy<'t, T>,
    root_node: Option<T>,
) -> TaxonomyResult<Value>
where
    T: Clone + Debug + Display + Eq + Hash + PartialEq,
{
    let tax_id = root_node.ok_or(Error::new(ErrorKind::InvalidTaxonomy(
        "Taxonomy must have a root node.".to_string(),
    )))?;

    fn inner<'t, T: 't>(tax: &'t impl Taxonomy<'t, T>, tax_id: T) -> TaxonomyResult<TaxNodeTree>
    where
        T: Clone + Debug + Display + Eq + Hash + PartialEq,
    {
        let mut children = Vec::new();
        for child in tax.children(tax_id.clone())? {
            children.push(inner(tax, child)?);
        }
        let node = TaxNodeTree {
            id: tax_id.to_string(),
            name: tax.name(tax_id.clone())?.to_string(),
            rank: tax.rank(tax_id.clone())?,
            children,
            extra: (*tax.data(tax_id)?).clone(),
        };
        Ok(node)
    }

    Ok(to_value(inner(taxonomy, tax_id)?)?)
}

fn serialize_as_node_links<'t, T: 't>(
    tax: &'t impl Taxonomy<'t, T>,
    root_node: Option<T>,
) -> TaxonomyResult<Value>
where
    T: Clone + Debug + Display + Eq + Hash + PartialEq,
{
    let mut nodes = Vec::new();
    let mut links = Vec::new();
    let mut id_to_idx = HashMap::new();

    if let Some(root_id) = root_node {
        for (ix, (tid, _pre)) in tax.traverse(root_id)?.filter(|x| x.1).enumerate() {
            let node = TaxNode {
                id: tid.to_string(),
                name: tax.name(tid.clone())?.to_string(),
                rank: tax.rank(tid.clone())?,
                extra: (*tax.data(tid.clone())?).clone(),
            };
            nodes.push(to_value(&node).unwrap());
            id_to_idx.insert(tid.clone(), ix);
            if let Some(parent_id) = tax.parent(tid)? {
                links.push(json!({
                    "source": ix,
                    "target": id_to_idx[&parent_id.0],
                }));
            }
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
    fn can_load_and_save_empty_node_link_format() {
        // Savings
        let example = r#"{"nodes": [], "links": []}"#;
        let tax = load(Cursor::new(example), None).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 0);

        // Loading
        let mut out = Vec::new();
        save::<_, &str, _>(&mut out, &tax, JsonFormat::NodeLink, None).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            json!({ "nodes": [], "links": [], "directed": true, "multigraph": false, "graph": [] })
                .to_string()
        )
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
        assert_eq!(Taxonomy::<&str>::len(&tax), 3);
        assert_eq!(Taxonomy::<&str>::root(&tax), "1");
        assert_eq!(Taxonomy::<&str>::children(&tax, "1").unwrap(), vec!["2"]);
        assert_eq!(
            Taxonomy::<&str>::lineage(&tax, "562").unwrap(),
            vec!["562", "2", "1"]
        );
        assert_eq!(
            Taxonomy::<&str>::data(&tax, "1").unwrap()["readcount"]
                .as_i64()
                .unwrap(),
            1000
        );
        // Converting to Value should give us the same data
        let serialized = serialize_as_node_links(&tax, Some(Taxonomy::<&str>::root(&tax))).unwrap();

        let tax2 = load_node_link_json(&serialized).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax2), 3);
        assert_eq!(Taxonomy::<&str>::root(&tax2), "1");
        assert_eq!(Taxonomy::<&str>::children(&tax2, "1").unwrap(), vec!["2"]);
        assert_eq!(
            Taxonomy::<&str>::lineage(&tax2, "562").unwrap(),
            vec!["562", "2", "1"]
        );
        assert_eq!(
            Taxonomy::<&str>::data(&tax2, "1").unwrap()["readcount"]
                .as_i64()
                .unwrap(),
            1000
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
                            "readcount": 3000
                        }
                    ]
                }
            ]
        }"#;

        let tax = load(Cursor::new(example), None).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 3);
        assert_eq!(Taxonomy::<&str>::root(&tax), "1");
        assert_eq!(Taxonomy::<&str>::children(&tax, "1").unwrap(), vec!["2"]);
        assert_eq!(
            Taxonomy::<&str>::lineage(&tax, "562").unwrap(),
            vec!["562", "2", "1"]
        );
        assert_eq!(
            Taxonomy::<&str>::data(&tax, "1").unwrap()["readcount"]
                .as_i64()
                .unwrap(),
            1000
        );
        // Converting to Value should give us the same data
        let serialized = serialize_as_tree(&tax, Some(Taxonomy::<&str>::root(&tax))).unwrap();
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

    fn get_test_node_link_data() -> &'static str {
        r#"
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
}        "#
    }

    #[test]
    fn can_load_node_link() {
        let example = get_test_node_link_data();
        let tax = load(Cursor::new(example), None).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 11);
        assert_eq!(Taxonomy::<&str>::root(&tax), "1");
        assert_eq!(
            Taxonomy::<&str>::lineage(&tax, "10").unwrap(),
            vec!["10", "8", "7", "6", "5", "4", "3", "2", "1"]
        );
    }

    #[test]
    fn can_use_node_index_as_internal_id() {
        let example = get_test_node_link_data();
        let tax = load(Cursor::new(example), None).unwrap();
        let tax_ids = vec!["1", "9", "2", "11", "8", "5", "3", "4", "6", "7", "10"];
        assert_eq!(
            tax_ids
                .into_iter()
                .map(|i| tax.to_internal_index(i).unwrap())
                .collect::<Vec<_>>(),
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        )
    }

    #[test]
    fn can_load_from_json_path() {
        let example = r#"{"test": {"sub": {"nodes": [], "links": []}}}"#;
        let tax = load(Cursor::new(example), Some("/test/sub")).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 0);
    }

    #[test]
    fn can_handle_null_ranks() {
        let example = r#"{"id": "1", "rank": null, "name": ""}"#;
        assert!(load(Cursor::new(example), None).is_ok());
    }

    #[test]
    fn streaming_tree_matches_regular_tree() {
        // Create a taxonomy with multiple nodes and extra data
        let example = r#"{
            "id": "1",
            "name": "root",
            "rank": "no rank",
            "custom_field": "value1",
            "children": [
                {
                    "id": "2",
                    "name": "Bacteria",
                    "rank": "superkingdom",
                    "genetic_code_id": "11",
                    "children": [
                        {
                            "id": "562",
                            "name": "Escherichia coli",
                            "rank": "species",
                            "genetic_code_id": "11",
                            "embl_code": "EC",
                            "division_id": "0",
                            "name_common_name": "E. coli"
                        }
                    ]
                }
            ]
        }"#;

        let tax = load(Cursor::new(example), None).unwrap();

        // Save using regular method
        let mut regular_output = Vec::new();
        save::<_, &str, _>(&mut regular_output, &tax, JsonFormat::Tree, None).unwrap();

        // Save using streaming method
        let mut streaming_output = Vec::new();
        save_tree_streaming::<_, &str, _>(&mut streaming_output, &tax, None).unwrap();

        // Parse both outputs
        let regular_json: Value = serde_json::from_slice(&regular_output).unwrap();
        let streaming_json: Value = serde_json::from_slice(&streaming_output).unwrap();

        // They should be equivalent
        assert_eq!(regular_json, streaming_json);

        // Verify we can reload from streaming output
        let tax2 = load(Cursor::new(&streaming_output[..]), None).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax2), 3);
        assert_eq!(
            Taxonomy::<&str>::name(&tax2, "562").unwrap(),
            "Escherichia coli"
        );

        // Verify extra data is preserved
        let data = Taxonomy::<&str>::data(&tax2, "562").unwrap();
        assert_eq!(
            data.get("genetic_code_id").and_then(|v| v.as_str()),
            Some("11")
        );
        assert_eq!(data.get("embl_code").and_then(|v| v.as_str()), Some("EC"));
        assert_eq!(
            data.get("name_common_name").and_then(|v| v.as_str()),
            Some("E. coli")
        );
    }

    #[test]
    fn streaming_handles_large_data_fields() {
        // Test that streaming handles nodes with lots of extra fields (like NCBI data)
        let example = r#"{
            "id": "562",
            "name": "Escherichia coli",
            "rank": "species",
            "genetic_code_id": "11",
            "division_id": "0",
            "inherited_div_flag": "1",
            "mitochondrial_genetic_code_id": "0",
            "inherited_MGC_flag": "1",
            "GenBank_hidden_flag": "0",
            "hidden_subtree_root_flag": "0",
            "comments": "",
            "name_scientific_name": "Escherichia coli",
            "name_common_name": "E. coli",
            "name_synonym": "Bacterium coli",
            "name_authority": "Escherichia coli (Migula 1895) Castellani and Chalmers 1919"
        }"#;

        let tax = load(Cursor::new(example), None).unwrap();

        let mut streaming_output = Vec::new();
        save_tree_streaming::<_, &str, _>(&mut streaming_output, &tax, None).unwrap();

        // Reload and verify all fields preserved
        let tax2 = load(Cursor::new(&streaming_output[..]), None).unwrap();
        let data = Taxonomy::<&str>::data(&tax2, "562").unwrap();

        assert_eq!(
            data.get("genetic_code_id").and_then(|v| v.as_str()),
            Some("11")
        );
        assert_eq!(
            data.get("name_common_name").and_then(|v| v.as_str()),
            Some("E. coli")
        );
        assert_eq!(
            data.get("name_synonym").and_then(|v| v.as_str()),
            Some("Bacterium coli")
        );
        assert_eq!(data.get("comments").and_then(|v| v.as_str()), Some(""));
    }
}
