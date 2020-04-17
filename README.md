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


# Python Usage

The Python taxonomy API can open and manipulate all of the formats from the Rust library:

```python
from taxonomy import Taxonomy

tax = Taxonomy.from_newick('(A,(B,C)D)E;')
assert tax.parent('A') == 'E'
assert tax.parent('B') == 'D'
```

If you have the NCBI taxonomy locally ([found on their FTP](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)), you can use that too:

```python
ncbi_tax = Taxonomy.from_ncbi('./nodes.dmp', './names.dmp')
assert tax.name('562') == 'Escherichia coli'
assert tax.rank('562') == 'species'
```

Note that Taxonomy IDs in NCBI format are integers, but they're converted to strings on import. We find working with "string taxonomy IDs" greatly simplifies interoperation between different taxonomy systems.


# Installation

## Rust
This library can be added to an existing Cargo.toml file and installed straight from crates.io.

## Python
You can install the Python bindings directly from PyPI (binaries are only built for select architextures) with:
```bash
pip install taxonomy
```


# Development

## Rust
There is a test suite runable with `cargo test`. To test the Python-bindings you need to use the additional `python_test` feature: `cargo test --features python_test`.

## Python
To work on the Python library on a Mac OS X/Unix system (requires Python 3):
```bash
# you need the nightly version of Rust installed
curl https://sh.rustup.rs -sSf | sh
rustup default nightly

# finally, install the library
maturin install --cargo-extra-args="--features=python"
```

### Building binary wheels and pushing to PyPI

```
# The Mac build requires switching through a few different python versions
maturin build --cargo-extra-args="--features=python" --release --strip

# The linux build is automated through cross-compiling in a docker image
docker run --rm -v $(pwd):/io konstin2/maturin:master build --cargo-extra-args="--features=python" --release --strip
twine upload target/wheels/*
```

# Other Taxonomy Libraries

There are taxonomic toolkits for other programming languages that offer different features and provided some inspiration for this library:

*ETE Toolkit (http://etetoolkit.org/)* A Python taxonomy library

*Taxize (https://ropensci.github.io/taxize-book/)* An R toolkit for working with taxonomic data
