# Taxonomy

[![PyPI version](https://badge.fury.io/py/taxonomy.svg)](https://pypi.org/project/taxonomy/)
[![Crates version](https://img.shields.io/crates/v/taxonomy.svg)](https://crates.io/crates/taxonomy)
![CI](https://github.com/onecodex/taxonomy/workflows/CI/badge.svg)

This is a Rust library for reading, writing, and editing biological taxonomies. There are associated Python bindings for accessing most of the functionality from Python.

This library was developed initially as a component in One Codex's metagenomic classification pipeline before being refactored out, expanded, and open-sourced. It is designed such that it can be used *as is* with a number of taxonomic formats *or* the Taxonomy trait it provides can be used to add last common ancestor, traversal, etc. methods to a downstream package's taxonomy implementation.

The library ships with a number of features:
 - [X] Common support for taxonomy handling across Rust and Python
 - [X] Fast and low(er) memory usage
 - [X] NCBI taxonomy, JSON ("tree" and "node_link_data" formats), Newick, and PhyloXML support
 - [X] Easily extensible (in Rust) to support other formats and operations


## Installation

### Rust
This library can be added to an existing Cargo.toml file and installed straight from crates.io.

### Python
You can install the Python bindings directly from PyPI (binaries are only built for select architectures) with:

```bash
pip install taxonomy
```

## Python Usage

The Python taxonomy API can open and manipulate all of the formats from the Rust library.
Note that Taxonomy IDs in NCBI format are integers, but they're converted to strings on import. We find working with "string taxonomy IDs" greatly simplifies inter-operation between different taxonomy systems.

### Loading a taxonomy

Taxonomy can be loaded from a variety of sources.

1. `Taxonomy.from_newick(value: str)`: loads a Taxonomy from a Newick-encoded string.

2. `Taxonomy.from_ncbi(ncbi_filder: str)`: loads a Taxonomy from a pair of NCBI dump files. The folder needs to contain the individual files in the NCBI taxonomy directory (e.g. nodes.dmp and names.dmp).

3. `Taxonomy.from_json(value: str, /, json_pointer: str)`: loads a Taxonomy from a JSON-encoded string. The format can either be
of the tree or node_link_data types and will be automatically detected. If `path` is specified, the JSON will be traversed to that sub-object before being parsed as a taxonomy.

4. `Taxonomy.from_phyloxml(value: &str)`: loads a Taxonomy from a PhyloXML-encoded string. **Experimental**

### Exporting a taxonomy

Assuming that the taxonomy has been instantiated as a variable named `tax`.

1. `tax.to_newick()`: exports a Taxonomy as a Newick-encoded byte string.
2. `tax.to_json_tree()`: exports a Taxonomy as a JSON-encoded byte string in a tree format
3. `tax.to_json_node_links()`: exports a Taxonomy as a JSON-encoded byte string in a node links format

### Using a taxonomy

Assuming that the taxonomy has been instantiated as a variable named `tax`. Note that `TaxonomyNode` is a class with
the following schema:

```python
class TaxonomyNode:
    id: str
    name: str
    parent: Optional[str]
    rank: str
```

Note that tax_id in parameters passed in functions described below are string but for example in the case of NCBI need
to be essentially quoting integers: `562 -> "562"`. 
If you loaded a taxonomy via JSON and you had additional data in your file, you can access it via indexing, `node["readcount"]` for example.

#### `tax.root -> TaxonomyNode`
Points to the root of the taxonomy

#### `tax.parent(tax_id: str, /, at_rank: str) -> Optional[TaxonomyNode]`
Return the immediate parent TaxonomyNode of the node id.

If `at_rank` is provided, scan all the nodes in the node's lineage and return
the parent id at that rank.

Examples:

```py
parent = tax.parent("612")
parent = tax.parent("612", at_rank="species")
parent = tax.parent("612")
# Both variables will be `None` if we can't find the parent
parent = tax.parent("unknown")
```

#### `tax.parent_with_distance(tax_id: str, /, at_rank: str) -> (Optional[TaxonomyNode], Optional[float])`
Same as `parent` but return the distance in addition, as a `(TaxonomyNode, float)` tuple.

#### `tax.node(tax_id: str) -> Optional[TaxonomyNode]`

Returns the node at that id. Returns `None` if not found.
You can also use indexing to accomplish that: `tax["some_id"]` but this will raise an exception if the node
is not found.

#### `tax.find_by_name(name: str) -> Optional[TaxonomyNode]`

Returns the node with that name. Returns `None` if not found.
In NCBI, it only accounts for *scientific names* and not synonyms.

#### `tax.children(tax_id: str) -> List[TaxonomyNode]`

Returns all nodes below the given tax id.

#### `tax.lineage(tax_id: str) -> List[TaxonomyNode]`

Returns all nodes above the given tax id, including itself.

#### `tax.parents(tax_id: str) -> List[TaxonomyNode]`

Returns all nodes above the given tax id.

#### `tax.lca(id1: str, id2: str) -> Optional[TaxonomyNode]`

Returns the [lowest common ancestor](https://en.wikipedia.org/wiki/Lowest_common_ancestor) for the 2 given nodes.

#### `tax.prune(keep: List[str], remove: List[str])-> Taxonomy`

Return a copy of the taxonomy containing:

- only the nodes in `keep` and their parents if provided
- all of the nodes except those in remove and their children if provided

#### `tax.remove_node(tax_id: str)`

Remove the node from the tree, re-attaching parents as needed: only a single node is removed.

#### `tax.add_node(parent_tax_id: str, new_tax_id: str)`

Add a new node to the tree at the parent provided.

#### `edit_node(tax_id: str, /, name: str, rank: str, parent_id: str, parent_dist: float)`

Edit properties on a taxonomy node.

### Exceptions
Only one exception is raised intentionally by the library: `TaxonomyError`.
If you get a `pyo3_runtime.PanicException` (or anything with `pyo3` in its name), this is a bug in the underlying Rust library, please open an issue.

## Development

### Rust
There is a test suite runable with `cargo test`. To test the Python-bindings you need to use the additional `python_test` feature: `cargo test --features python_test`.

### Python
To work on the Python library on a Mac OS X/Unix system (requires Python 3):
```bash
# you need the nightly version of Rust installed
curl https://sh.rustup.rs -sSf | sh

# finally, install the library in the local virtualenv
maturin develop --cargo-extra-args="--features=python"

# or using pip
pip install .
```

#### Building binary wheels and pushing to PyPI

```
# The Mac build requires switching through a few different python versions
maturin build --cargo-extra-args="--features=python" --release --strip

# The linux build is automated through cross-compiling in a docker image
docker run --rm -v $(pwd):/io ghcr.io/pyo3/maturin build --cargo-extra-args="--features=python" --release --strip
twine upload target/wheels/*
```

## Other Taxonomy Libraries

There are taxonomic toolkits for other programming languages that offer different features and provided some inspiration for this library:

*ETE Toolkit (http://etetoolkit.org/)* A Python taxonomy library

*Taxize (https://ropensci.github.io/taxize-book/)* An R toolkit for working with taxonomic data
