// On PhyloXML
//
// http://etetoolkit.org/docs/latest/tutorial/tutorial_phyloxml.html
// https://biopython.org/wiki/PhyloXML
// http://www.phyloxml.org/

// On Rust XML parsing
//
// https://github.com/tafia/quick-xml
// https://docs.rs/quick-xml/0.12.4/quick_xml/
// https://docs.rs/quick-xml/0.12.4/quick_xml/struct.Reader.html

use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::io::{BufReader, Read, Write};
use std::iter::Sum;
use std::str::FromStr;

use failure::bail;
use quick_xml::events::Event;
use quick_xml::Reader;

use crate::base::GeneralTaxonomy;
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy;
use crate::Result;

/// Write a Taxonomy object out to a `writer` in the PhyloXML format.
///
/// Completely unimplemented.
pub fn save_phyloxml<'t, T: 't, D: 't, X, W>(
    _tax: &'t X,
    _writer: &mut W,
    _root_node: Option<T>,
) -> Result<()>
where
    W: Write,
    T: Clone + Debug + Display + PartialEq,
    D: Debug + Display + PartialOrd + Sum,
    X: Taxonomy<'t, T, D>,
{
    unimplemented!("This feature isn't complete yet")
}

/// Read PhyloXML format into a Taxonomy object out of a `reader`.
///
/// Still somewhat experimental and may not support all Newick features.
pub fn load_phyloxml<R>(reader: &mut R) -> Result<GeneralTaxonomy>
where
    R: Read,
{
    let buf_reader = BufReader::new(reader);
    let mut xml_reader = Reader::from_reader(buf_reader);
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
            Event::Eof => bail!("No valid phyloxml taxonomy found"),
            _ => continue,
        }
        buf.clear();
    }

    let mut tax_ids: Vec<String> = vec![];
    let mut names: Vec<Option<String>> = vec![];
    let mut parent_ids: Vec<usize> = vec![];
    let mut ranks: Vec<Option<TaxRank>> = vec![];
    let mut dists: Vec<f32> = vec![];

    let mut cur_lineage: Vec<usize> = Vec::new();

    let mut current_tag: Vec<u8> = b"".to_vec();
    loop {
        match xml_reader.read_event(&mut buf)? {
            Event::Start(ref e) => {
                match e.name() {
                    b"phylogeny" => bail!("Nested phylogeny not permitted"),
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
                                .parse()?,
                        );
                        ranks.push(None);
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
                    b"branch_length" => *dists.last_mut().unwrap() = text.parse()?,
                    b"rank" => *ranks.last_mut().unwrap() = TaxRank::from_str(&text).ok(),
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
        Some(
            names
                .into_iter()
                .map(|n| n.unwrap_or_else(|| "".to_string()))
                .collect(),
        )
    };

    // if not everything has a tax id, should we try to use the names instead?
    // or should we should have a fallback that makes up an ID for this (and
    // for Newick trees) for anything that doesn't have any kind of identifier?

    GeneralTaxonomy::new(tax_ids, parent_ids, cleaned_names, Some(ranks), Some(dists))
}

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use super::*;

    #[test]
    fn test_load_phyloxml() -> Result<()> {
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
        let tax = load_phyloxml(&mut text_cursor)?;
        assert_eq!(Taxonomy::<&str, f32>::len(&tax), 5);

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
        let tax = load_phyloxml(&mut text_cursor)?;
        assert_eq!(Taxonomy::<&str, f32>::len(&tax), 5);
        Ok(())
    }

    #[test]
    fn test_no_valid_phyloxml() -> Result<()> {
        let text_xml = r#"
        <document></document>
        "#;
        let mut text_cursor = Cursor::new(text_xml);
        assert!(load_phyloxml(&mut text_cursor).is_err());
        Ok(())
    }
}
