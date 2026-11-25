use std::io::{Read, Write};

use crate::base::GeneralTaxonomy;
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::Taxonomy;

use memchr::memchr2;
use std::collections::VecDeque;
use std::fmt::{Debug, Display};

/// NewickToken is used as an intermediate during the tokenization of a Newick string.
#[derive(PartialEq)]
enum NewickToken {
    Start,
    End,
    Delim,
    NameDist(String),
}

impl NewickToken {
    fn as_bytes(&self) -> &[u8] {
        match self {
            NewickToken::Start => b"(",
            NewickToken::End => b")",
            NewickToken::Delim => b",",
            NewickToken::NameDist(s) => s.as_bytes(),
        }
    }
}

pub fn save<'t, T: 't, W: Write>(
    writer: &mut W,
    taxonomy: &'t impl Taxonomy<'t, T>,
    root_node: Option<T>,
) -> TaxonomyResult<()>
where
    T: Clone + Debug + Display + PartialEq,
{
    let root_node = root_node.unwrap_or_else(|| taxonomy.root());
    let mut out_buf = VecDeque::new();

    for (node, pre) in taxonomy.traverse(root_node)? {
        if pre {
            out_buf.push_back(NewickToken::Start);
        } else {
            out_buf.push_back(NewickToken::End);
            let mut name = taxonomy
                .name(node.clone())?
                .replace(char::is_whitespace, "_");
            if let Some((_, dist)) = taxonomy.parent(node)? {
                if dist > 0.0 {
                    name.push(':');
                    name.push_str(&dist.to_string());
                }
            }
            out_buf.push_back(NewickToken::NameDist(name));
            out_buf.push_back(NewickToken::Delim);
        }
    }

    let mut skip_next = false;
    while let Some(token) = out_buf.pop_front() {
        if skip_next {
            skip_next = false;
            continue;
        }
        let next_token = out_buf.front();
        // turn `()`s into `` (remove lists with no children)
        if token == NewickToken::Start && next_token == Some(&NewickToken::End) {
            skip_next = true;
            continue;
        }
        // remove terminal commas from lists
        if token == NewickToken::Delim
            && (next_token == Some(&NewickToken::End) || next_token.is_none())
        {
            continue;
        }
        writer.write_all(token.as_bytes())?;
    }
    writer.write_all(b";")?;
    Ok(())
}

/// Read Newick format into a Taxonomy object out of a `reader`.
///
/// Node labels are stored as names. Newick labels are optional and often repeated, so each
/// node's tax ID is its position in the tree (the root is "0").
///
/// Still somewhat experimental and may not support all Newick features.
pub fn load<R: Read>(reader: &mut R) -> TaxonomyResult<GeneralTaxonomy> {
    // TODO: handle quoted labels
    // TODO: handle ending `;`
    let mut buffer: Vec<u8> = Vec::new();
    reader.read_to_end(&mut buffer)?;

    // make the root node explicitly first
    let mut names: Vec<String> = vec!["".to_string()];
    let mut parent_ids: Vec<usize> = vec![0];
    let mut dists: Vec<f32> = vec![1.];
    let mut cur_lineage: Vec<usize> = Vec::new();

    let mut cur_pos = 0;
    while cur_pos < buffer.len() {
        match buffer[cur_pos] {
            b'(' => {
                // make a new node
                names.push("".to_string());
                parent_ids.push(*cur_lineage.last().unwrap_or(&0));
                dists.push(1.);
                cur_lineage.push(names.len() - 1);
                cur_pos += 1;
            }
            b',' => {
                // finish up with the current node
                cur_lineage.pop();
                // make a new node
                names.push("".to_string());
                parent_ids.push(*cur_lineage.last().unwrap_or(&0));
                dists.push(1.);
                cur_lineage.push(names.len() - 1);
                cur_pos += 1;
            }
            b')' => {
                // just go back to the parent
                cur_lineage.pop();
                cur_pos += 1;
            }
            _ => {
                // set the name and dist on the cur_node
                let cur_node = *cur_lineage.last().unwrap_or(&0);
                let pos = memchr2(b',', b')', &buffer[cur_pos..])
                    .map(|x| x + cur_pos)
                    .unwrap_or_else(|| buffer.len());
                let name_dist = &buffer[cur_pos..pos];
                let mut chunk_iter = std::str::from_utf8(name_dist)
                    .map_err(|_| {
                        Error::new(ErrorKind::ImportError {
                            line: 0,
                            msg: format!(
                                "Could not parse name/dist \"{}\" as unicode",
                                String::from_utf8_lossy(name_dist)
                            ),
                        })
                    })?
                    .trim_end_matches(';')
                    .splitn(2, |x| x == ':');
                names[cur_node] = chunk_iter.next().unwrap_or("").to_string();
                dists[cur_node] = chunk_iter.next().unwrap_or("1").parse().map_err(|e| {
                    Error::new(ErrorKind::ImportError {
                        line: 0,
                        msg: format!("Could not parse distance as a number: {}", e),
                    })
                })?;
                cur_pos = pos;
            }
        }
    }

    let tax_ids = (0..names.len()).map(|i| i.to_string()).collect();
    GeneralTaxonomy::from_arrays(tax_ids, parent_ids, Some(names), None, Some(dists), None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::taxonomy::Taxonomy;

    const HOMARIDAE: &[u8] = b"(Nephrops_norvegicus:0.1,Homarinus_capensis:0.2,(Homarus_americanus:0.3,Homarus_gammarus:0.4)Homarus:0.5)Homaridae;";

    fn id_of<'t>(tax: &'t GeneralTaxonomy, name: &str) -> &'t str {
        let ids = tax.find_all_by_name(name);
        assert_eq!(ids.len(), 1, "expected exactly one node named {}", name);
        ids[0]
    }

    #[test]
    fn test_write_newick() {
        let tax = load(&mut HOMARIDAE.as_ref()).unwrap();

        let mut bytes = Vec::new();
        save(&mut bytes, &tax, Some(id_of(&tax, "Homarus"))).unwrap();
        assert_eq!(
            bytes,
            b"(Homarus_gammarus:0.4,Homarus_americanus:0.3)Homarus:0.5;"
        );
    }

    #[test]
    fn test_write_newick_replaces_whitespace_in_names() {
        let tax = GeneralTaxonomy::from_arrays(
            vec!["0".to_string(), "1".to_string()],
            vec![0, 0],
            Some(vec![
                "Homaridae".to_string(),
                "Homarus americanus".to_string(),
            ]),
            None,
            None,
            None,
        )
        .unwrap();

        let mut bytes = Vec::new();
        save(&mut bytes, &tax, None::<&str>).unwrap();
        assert_eq!(bytes, b"(Homarus_americanus:1)Homaridae;");
    }

    #[test]
    fn test_load_newick() {
        let tax = load(&mut b"(())".as_ref()).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 3);

        let tax = load(
            &mut b"(Homarus_americanus,Homarus_gammarus,(Nephrops_norvegicus,Metanephrops_japonicus));"
                .as_ref(),
        )
        .unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 6);

        let tax = load(&mut HOMARIDAE.as_ref()).unwrap();
        assert_eq!(Taxonomy::<&str>::root(&tax), "0");
        assert_eq!(tax.name("0").unwrap(), "Homaridae");
        assert!(!tax.contains("Homarus_gammarus"));
        assert_eq!(
            tax.parent(id_of(&tax, "Homarus_gammarus")).unwrap(),
            Some((id_of(&tax, "Homarus"), 0.4))
        );
        assert_eq!(
            tax.parent(id_of(&tax, "Homarus")).unwrap(),
            Some((id_of(&tax, "Homaridae"), 0.5))
        );
    }

    #[test]
    fn test_load_newick_gives_unlabelled_and_duplicate_nodes_unique_ids() {
        let tax = load(
            &mut b"((Homarus_americanus,Homarus_americanus),(Nephrops_norvegicus,Metanephrops_japonicus))Homaridae;"
                .as_ref(),
        )
        .unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 7);
        assert_eq!(tax.find_all_by_name("").len(), 2);
        assert_eq!(tax.find_all_by_name("Homarus_americanus").len(), 2);
        tax.validate_uniqueness().unwrap();
    }
}
