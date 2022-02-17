use std::fmt;

#[derive(Debug, PartialEq, Eq)]
pub enum ErrorKind {
    UnknownRank(String),
    ImportError { line: usize, msg: String },
    InvalidTaxonomy(String),
    NoSuchTaxId(String),
    OperationNotAllowed(String),
}

#[derive(Debug)]
pub struct Error {
    kind: ErrorKind,
    source: Option<Box<dyn std::error::Error + Send + Sync>>,
}

impl Error {
    pub fn new(kind: ErrorKind) -> Self {
        Self { kind, source: None }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            ErrorKind::NoSuchTaxId(s) => write!(f, "Tax ID {} not found in taxonomy", s),
            ErrorKind::UnknownRank(r) => write!(f, "Rank {} is unknown", r),
            ErrorKind::ImportError { line, msg } => {
                write!(f, "Failed to import taxonomy {} at line {}", msg, line)
            }
            ErrorKind::InvalidTaxonomy(s) => write!(f, "Invalid taxonomy: {}", s),
            ErrorKind::OperationNotAllowed(s) => {
                write!(f, "Operation on taxonomy not allowed: {}", s)
            }
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source.as_ref().map(|e| e.as_ref() as _)
    }
}

impl From<serde_json::error::Error> for Error {
    fn from(error: serde_json::error::Error) -> Self {
        let mut err = Error::new(ErrorKind::ImportError {
            line: error.line(),
            msg: error.to_string(),
        });
        err.source = Some(Box::new(error));
        err
    }
}

impl From<std::io::Error> for Error {
    fn from(error: std::io::Error) -> Self {
        let mut err = Error::new(ErrorKind::ImportError {
            line: 0,
            msg: "Failed to read data".to_owned(),
        });
        err.source = Some(Box::new(error));
        err
    }
}

impl From<quick_xml::Error> for Error {
    fn from(error: quick_xml::Error) -> Self {
        let mut err = Error::new(ErrorKind::ImportError {
            line: 0,
            msg: "Error parsing XML".to_owned(),
        });
        err.source = Some(Box::new(error));
        err
    }
}

pub type TaxonomyResult<T> = Result<T, Error>;
