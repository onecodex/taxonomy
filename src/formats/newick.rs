use std::collections::VecDeque;
use std::fmt::{Debug, Display};
use std::io::{Read, Write};
use std::iter::Sum;
use std::str;

use memchr::memchr2;

use crate::base::GeneralTaxonomy;
use crate::taxonomy::Taxonomy;
use crate::Result;

// TODO: add "use name intead of id" options?
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

// TODO: this allocates a buffer for the *entire* taxonomy; it should be possible
// to stream the writing in the traversal loop, but we need to buffer enough to
// determine whether to drop `()`s and extra `,`s
/// Write a Taxonomy object formatted as Newick to a `writer`.
///
/// Takes a `root_node` to allow writing only a subset of the total taxonomy.
///
/// Still somewhat experimental and may not support all Newick features.
pub fn save_newick<'t, T: 't, D: 't, X, W>(
    tax: &'t X,
    writer: &mut W,
    root_node: Option<T>,
) -> Result<()>
where
    W: Write,
    T: Clone + Debug + Display + PartialEq,
    D: Debug + Display + PartialOrd + Sum,
    X: Taxonomy<'t, T, D>,
{
    let root_id = if let Some(tid) = root_node {
        tid
    } else {
        tax.root()
    };

    let mut out_buf: VecDeque<NewickToken> = VecDeque::new();
    let zero: D = Vec::new().into_iter().sum();

    for (node, pre) in tax.traverse(root_id)? {
        if pre {
            out_buf.push_back(NewickToken::Start);
        } else {
            out_buf.push_back(NewickToken::End);
            let mut name: String = format!("{}", node);
            if let Some((_, dist)) = tax.parent(node)? {
                if dist > zero {
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
pub fn load_newick<R>(reader: &mut R) -> Result<GeneralTaxonomy>
where
    R: Read,
{
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
                let mut chunk_iter = str::from_utf8(name_dist)?
                    .trim_end_matches(';')
                    .splitn(2, |x| x == ':');
                tax_ids[cur_node] = chunk_iter.next().unwrap_or("").to_string();
                dists[cur_node] = chunk_iter.next().unwrap_or("1").parse()?;
                cur_pos = pos;
            }
        }
    }

    Ok(GeneralTaxonomy::new(
        tax_ids,
        parent_ids,
        None,
        None,
        Some(dists),
    ))
}

#[test]
fn test_write_newick() {
    use crate::taxonomy::test::MockTax;
    let tax = MockTax;
    let mut s: Vec<u8> = Vec::new();
    save_newick(&tax, &mut s, None).unwrap();
    assert_eq!(s, b"(((((((((765909:1)61598:1)53452:1)1046:1)135613:1,(((56812:1)62322:1)22:1)135622:1)1236:1)1224:1)2:1)131567:1)1;".to_vec());
}

#[test]
fn test_load_newick() {
    use crate::taxonomy::Taxonomy;

    let newick_str = b"(())";
    let tax = load_newick(&mut newick_str.as_ref());

    let newick_str = b"(A,B,(C,D));";
    let tax = load_newick(&mut newick_str.as_ref());

    let newick_str = b"(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;";
    let tax = load_newick(&mut newick_str.as_ref()).unwrap();
    assert_eq!(tax.parent("D").unwrap(), Some(("E", 0.4)));
    assert_eq!(tax.parent("E").unwrap(), Some(("F", 0.5)));
}
