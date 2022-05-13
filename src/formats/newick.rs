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
            let mut name: String = format!("{}", node);
            if let Some((_, dist)) = taxonomy.parent(node)? {
                if dist > 0.0 {
                    name.push(':');
                    name.push_str(&format!("{}", dist));
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
            && (next_token == Some(&NewickToken::End) || next_token == None)
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
/// Still somewhat experimental and may not support all Newick features.
pub fn load<R: Read>(reader: &mut R) -> TaxonomyResult<GeneralTaxonomy> {
    // TODO: handle empty tax ids?
    // TODO: handle quoted labels
    // TODO: handle ending `;`
    let mut buffer: Vec<u8> = Vec::new();
    reader.read_to_end(&mut buffer)?;

    // make the root node explicitly first
    let mut tax_ids: Vec<String> = vec!["".to_string()];
    let mut parent_ids: Vec<usize> = vec![0];
    let mut dists: Vec<f32> = vec![1.];
    let mut cur_lineage: Vec<usize> = Vec::new();

    let mut cur_pos = 0;
    while cur_pos < buffer.len() {
        match buffer[cur_pos] {
            b'(' => {
                // make a new node
                tax_ids.push("".to_string());
                parent_ids.push(*cur_lineage.last().unwrap_or(&0));
                dists.push(1.);
                cur_lineage.push(tax_ids.len() - 1);
                cur_pos += 1;
            }
            b',' => {
                // finish up with the current node
                cur_lineage.pop();
                // make a new node
                tax_ids.push("".to_string());
                parent_ids.push(*cur_lineage.last().unwrap_or(&0));
                dists.push(1.);
                cur_lineage.push(tax_ids.len() - 1);
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
                tax_ids[cur_node] = chunk_iter.next().unwrap_or("").to_string();
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

    GeneralTaxonomy::from_arrays(tax_ids, parent_ids, None, None, Some(dists), None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::taxonomy::Taxonomy;

    #[test]
    fn test_write_newick() {
        let newick_str = b"((((((((56812:1)62322:1)22:1)135622:1,((((765909:1)61598:1)53452:1)1046:1)135613:1)1236:1)1224:1)2:1)131567:1)1;";
        let tax = load(&mut newick_str.as_ref()).unwrap();

        // The leaves are swapped, the tree is still ok but the output doesn't match the original string
        // let mut bytes = Vec::new();
        // save(&mut bytes, &tax, None).unwrap();
        // println!("{:?}", std::str::from_utf8(&bytes));
        // assert_eq!(bytes, newick_str);

        // try saving a subtree
        let mut bytes = Vec::new();
        save(&mut bytes, &tax, Some("53452")).unwrap();
        assert_eq!(bytes, b"((765909:1)61598:1)53452:1;");
    }

    #[test]
    fn test_load_newick() {
        let newick_str = b"(())";
        let tax = load(&mut newick_str.as_ref()).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 3);

        let newick_str = b"(A,B,(C,D));";
        let tax = load(&mut newick_str.as_ref()).unwrap();
        assert_eq!(Taxonomy::<&str>::len(&tax), 6);

        let newick_str = b"(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;";
        let tax = load(&mut newick_str.as_ref()).unwrap();
        assert_eq!(tax.parent("D").unwrap(), Some(("E", 0.4)));
        assert_eq!(tax.parent("E").unwrap(), Some(("F", 0.5)));
    }
}
