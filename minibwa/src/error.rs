use std::path::PathBuf;

/// Result alias for this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors returned by the safe minibwa API.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// The index could not be loaded from disk. The path that was attempted and
    /// a diagnostic message are included.
    #[error("failed to load index at {path}: {msg}")]
    IndexLoad { path: PathBuf, msg: String },

    /// Index construction failed. Includes the FASTA path and a diagnostic
    /// message from the C layer.
    #[error("failed to build index from {fasta}: {msg}")]
    IndexBuild { fasta: PathBuf, msg: String },

    /// An [`Opts`](crate::Opts) operation failed, e.g. an unknown preset name.
    #[error("invalid options: {0}")]
    InvalidOpts(String),

    /// A caller-supplied value was rejected: empty sequence, NUL in name, etc.
    #[error("invalid input: {0}")]
    InvalidInput(String),

    /// An I/O error propagated from the standard library.
    #[error("I/O: {0}")]
    Io(#[from] std::io::Error),
}
