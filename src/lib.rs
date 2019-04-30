#![cfg_attr(feature = "python", feature(specialization))]
#[cfg(feature = "python")]
#[macro_use]
extern crate pyo3;

use std::result::Result as StdResult;

use failure::{Error, Fail};

/// A wrapper type for taxonomy results.
/// Alias for `Result<T, failure::Error>`
pub type Result<T> = StdResult<T, Error>;

/// A taxonomically-specific error
#[derive(Debug, Fail)]
pub enum TaxonomyError {
    /// The given key could not be found in the taxonomy.
    #[fail(display = "Key is not valid: {}", key)]
    NoSuchKey { key: String },
    /// A string could not be parsed into a TaxRank.
    #[fail(display = "Rank is not supported: {}", rank)]
    UnrecognizedRank { rank: String },
    /// A string could not be parsed into a TaxRank.
    #[fail(display = "Error in file {} at line {}: {}", file, line, msg)]
    ImportError {
        file: String,
        line: usize,
        msg: String,
    },
}

mod base;
pub mod distances;
pub mod edit;
pub mod formats;
#[cfg(feature = "python")]
pub mod python;
mod rank;
mod taxonomy;

pub use crate::base::GeneralTaxonomy;
pub use crate::rank::TaxRank;
pub use crate::taxonomy::Taxonomy;
