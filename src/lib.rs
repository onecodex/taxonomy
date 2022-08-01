#[cfg(any(feature = "python", feature = "python_test"))]
extern crate pyo3;

mod base;
mod edit;
pub mod errors;
mod formats;
#[cfg(any(feature = "python", feature = "python_test"))]
mod python;
mod rank;
mod taxonomy;

pub use crate::taxonomy::Taxonomy;
pub use base::GeneralTaxonomy;
pub use edit::{prune_away, prune_to};
pub use errors::{Error, ErrorKind};
pub use formats::gtdb;
pub use formats::json;
pub use formats::ncbi;
pub use formats::newick;
pub use formats::phyloxml;
pub use rank::TaxRank;
