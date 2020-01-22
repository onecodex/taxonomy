#[cfg(any(feature = "python", feature = "python_test"))]
extern crate pyo3;

use std::fmt;
use std::result::Result as StdResult;

/// A wrapper type for taxonomy results.
/// Alias for `Result<T, failure::Error>`
pub type Result<T> = StdResult<T, TaxonomyError>;

/// A taxonomically-specific error
#[derive(Debug)]
pub enum TaxonomyError {
    /// The given key could not be found in the taxonomy.
    // #[fail(display = "Key is not valid: {}", key)]
    NoSuchKey {
        key: String,
    },
    /// A string could not be parsed into a TaxRank.
    // #[fail(display = "Rank is not supported: {}", rank)]
    UnrecognizedRank {
        rank: String,
    },
    /// A string could not be parsed into a TaxRank.
    // #[fail(display = "Error in file {} at line {}: {}", file, line, msg)]
    ImportError {
        line: usize,
        msg: String,
    },
    /// The taxonomy has an unrooted node or some other issue that breaks
    /// the assumption its a tree
    // #[fail(display = "Taxonomy is not a tree; broken at {}", tax_id)]
    MalformedTree {
        tax_id: String,
    },
    // #[fail(display = "Each taxa must have only one {}", field)]
    CreationFailed {
        field: String,
    },
}

impl fmt::Display for TaxonomyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TaxonomyError::NoSuchKey { key } => write!(f, "Key is not valid: {}", key),
            TaxonomyError::UnrecognizedRank { rank } => {
                write!(f, "Rank is not supported: {}", rank)
            }
            TaxonomyError::ImportError { line, msg } => {
                if *line != 0 {
                    write!(f, "Error importing at line {}: {}", line, msg)
                } else {
                    write!(f, "Error importing: {}", msg)
                }
            }
            TaxonomyError::MalformedTree { tax_id } => {
                write!(f, "Taxonomy is not a tree; broken at {}", tax_id)
            }
            TaxonomyError::CreationFailed { field } => {
                write!(f, "Each taxa must have only one {}", field)
            }
        }
    }
}

impl From<serde_json::error::Error> for TaxonomyError {
    fn from(error: serde_json::error::Error) -> Self {
        TaxonomyError::ImportError {
            line: error.line(),
            msg: error.to_string(),
        }
    }
}

impl From<quick_xml::Error> for TaxonomyError {
    fn from(error: quick_xml::Error) -> Self {
        TaxonomyError::ImportError {
            line: 0,
            msg: error.to_string(),
        }
    }
}

impl From<std::io::Error> for TaxonomyError {
    fn from(error: std::io::Error) -> Self {
        TaxonomyError::ImportError {
            line: 0,
            msg: error.to_string(),
        }
    }
}

mod base;
pub mod edit;
pub mod formats;
#[cfg(any(feature = "python", feature = "python_test"))]
pub mod python;
mod rank;
mod taxonomy;
pub mod weights;

pub use crate::base::GeneralTaxonomy;
pub use crate::rank::TaxRank;
pub use crate::taxonomy::Taxonomy;
