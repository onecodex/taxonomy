use std::io::{BufReader, Read};

use quick_xml::events::Event;

use crate::base::GeneralTaxonomy;
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;
use std::collections::HashMap;
use std::str::FromStr;

/// Read PhyloXML format into a Taxonomy object out of a `reader`.
///
/// Still somewhat experimental and may not support all PhyloXML features.
pub fn load<R: Read>(reader: &mut R) -> TaxonomyResult<GeneralTaxonomy> {
    let mut xml_reader = quick_xml::Reader::from_reader(BufReader::new(reader));
    xml_reader.trim_text(true);

    // find the <phylogeny> inside the XML
    let mut buf = Vec::new();
    loop {
        let evt = xml_reader.read_event(&mut buf)?;
        match evt {
            Event::Start(ref e) => {
                if e.name() == b"phylogeny" {
                    break;
                }
            }
            Event::Eof => {
                return Err(Error::new(ErrorKind::ImportError {
                    line: 0,
                    msg: "No valid phyloxml taxonomy found".to_owned(),
                }));
            }
            _ => continue,
        }
        buf.clear();
    }

    let mut tax_ids: Vec<String> = vec![];
    let mut names: Vec<Option<String>> = vec![];
    let mut parent_ids: Vec<usize> = vec![];
    let mut ranks: Vec<TaxRank> = vec![];
    let mut dists: Vec<f32> = vec![];

    let mut cur_lineage: Vec<usize> = Vec::new();

    let mut current_tag: Vec<u8> = b"".to_vec();
    loop {
        match xml_reader.read_event(&mut buf)? {
            Event::Start(ref e) => {
                match e.name() {
                    b"phylogeny" => {
                        return Err(Error::new(ErrorKind::ImportError {
                            line: 0,
                            msg: "Nested phylogeny not permitted".to_owned(),
                        }));
                    }
                    b"clade" => {
                        let attributes: HashMap<&[u8], String> = e
                            .attributes()
                            .map(|a| {
                                // TODO: sometimes this will blow up on an invalid attribute
                                // (e.g. <test a=">) but there's no way to pass the error
                                // out easily; we should probably figure out how to do that though
                                let att = a.unwrap();
                                (att.key, att.unescape_and_decode_value(&xml_reader).unwrap())
                            })
                            .collect();
                        cur_lineage.push(if tax_ids.is_empty() {
                            0 // the very first node is the root node
                        } else {
                            tax_ids.len() - 1
                        });
                        tax_ids.push("".to_string());
                        names.push(None);
                        parent_ids.push(*cur_lineage.last().unwrap_or(&0));
                        dists.push(
                            attributes
                                .get(&&b"branch_length"[..])
                                .unwrap_or(&"1".to_string())
                                .parse()
                                .map_err(|_| {
                                    Error::new(ErrorKind::ImportError {
                                        line: 0,
                                        msg: "Could not interpret branch length as a number"
                                            .to_owned(),
                                    })
                                })?,
                        );
                        ranks.push(TaxRank::Unspecified);
                    }
                    t => current_tag = t.to_vec(),
                }
            }
            Event::End(ref e) => {
                match e.name() {
                    b"phylogeny" => break,
                    b"clade" => {
                        cur_lineage.pop();
                    }
                    _ => current_tag = b"".to_vec(),
                };
            }
            Event::Text(e) => {
                if names.is_empty() {
                    // sometimes the phylogeny itself has a <name>
                    continue;
                }
                let text = e.unescape_and_decode(&xml_reader)?;
                match &current_tag[..] {
                    b"name" => *names.last_mut().unwrap() = Some(text),
                    b"id" => *tax_ids.last_mut().unwrap() = text,
                    b"branch_length" => {
                        *dists.last_mut().unwrap() = text.parse().map_err(|_| {
                            Error::new(ErrorKind::ImportError {
                                line: 0,
                                msg: "Could not interpret branch length as a number".to_string(),
                            })
                        })?;
                    }
                    b"rank" => *ranks.last_mut().unwrap() = TaxRank::from_str(&text)?,
                    // TODO: do something with confidence scores?
                    // b"confidence" => {},
                    // TODO: build up a dict of additional metadata we can scrape out?
                    // b"domain" | b"code" | b"id" | b"common_name" | b"scientific_name"
                    // | b"accession" | b"description" | b"desc" | b"symbol" | b"mol_seq" | b"property" => {},
                    _ => {} // bail!("Text outside in unexpected location"),
                }
            }
            Event::Eof => break,
            _ => (),
        }
        buf.clear();
    }

    let has_any_names = names.iter().all(|n| n == &None);

    let cleaned_names: Option<Vec<String>> = if has_any_names {
        None
    } else {
        Some(names.into_iter().map(|n| n.unwrap_or_default()).collect())
    };

    // if not everything has a tax id, should we try to use the names instead?
    // or should we should have a fallback that makes up an ID for this (and
    // for Newick trees) for anything that doesn't have any kind of identifier?

    GeneralTaxonomy::from_arrays(
        tax_ids,
        parent_ids,
        cleaned_names,
        Some(ranks),
        Some(dists),
        None,
    )
}

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use super::*;
    use crate::taxonomy::Taxonomy;

    #[test]
    fn test_load_phyloxml() {
        let text_xml = r#"
        <phylogeny rooted="true">
          <name>test taxonomy</name>
          <clade>
            <id>E</id>
            <clade>
              <id>D</id>
              <branch_length>0.3</branch_length>
              <clade>
                <name>A</name>
                <id>A</id>
                <branch_length>0.1</branch_length>
              </clade>
              <clade>
                <name>B</name>
                <id>B</id>
                <branch_length>0.2</branch_length>
              </clade>
            </clade>
            <clade>
              <name>C</name>
              <id>C</id>
              <branch_length>0.4</branch_length>
            </clade>
          </clade>
        </phylogeny>
        "#;
        let mut text_cursor = Cursor::new(text_xml);
        let tax = load(&mut text_cursor).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 5);

        let text_xml = r#"
        <phylogeny rooted="true">
           <name>test taxonomy</name>
           <clade>
              <clade branch_length="0.3">
                 <clade branch_length="0.1">
                    <name>A</name>
                 </clade>
                 <clade branch_length="0.2">
                    <name>B</name>
                 </clade>
              </clade>
              <clade branch_length="0.4">
                 <name>C</name>
              </clade>
           </clade>
        </phylogeny>
        "#;
        let mut text_cursor = Cursor::new(text_xml);
        let tax = load(&mut text_cursor).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 5);
    }

    #[test]
    fn test_no_valid_phyloxml() {
        let text_xml = r#"
        <document></document>
        "#;
        let mut text_cursor = Cursor::new(text_xml);
        assert!(load(&mut text_cursor).is_err());
    }
}
