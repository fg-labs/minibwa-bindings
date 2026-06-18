//! Safe, idiomatic Rust bindings for minibwa — a fast short-read, adaptive,
//! and Hi-C DNA aligner.
//!
//! # Workflow
//!
//! 1. **Build an index** once from a FASTA with [`Index::build_from_fasta`].
//! 2. **Load the index** with [`Index::load`]; the loaded index is `Send + Sync`
//!    and can be shared across threads.
//! 3. **Create an [`Aligner`]** from `&Index` and `&Opts` — cheap, no allocation.
//! 4. **Create one [`ThreadBuf`] per worker thread** — holds the per-thread scratch
//!    arena; it is `Send` but not `Sync`.
//! 5. **Call [`Aligner::map`] or [`Aligner::map_pair`]** to align reads and get
//!    owned [`Hit`] vectors.
//!
//! # Methylation support
//!
//! Pass [`Meth::C2T`] or [`Meth::G2A`] to [`Aligner::map`] for single-read
//! bisulfite alignment, or enable [`Opts::set_methylation`] together with
//! [`Aligner::map_pair`] for paired-end bisulfite alignment.
//!
//! # Paired-end support
//!
//! Use [`Aligner::map_pair`] for paired-end alignment. Insert-size parameters
//! can be tuned with [`Opts::set_pe_insert_size`].
//!
//! # Parallelism
//!
//! Thread-level parallelism is caller-owned: create one [`ThreadBuf`] per thread
//! and call [`Aligner::map`] / [`Aligner::map_pair`] concurrently. The shared
//! [`Index`] is read-only during mapping.
//!
//! # Example
//!
//! ```no_run
//! use minibwa::{Index, Opts, Aligner, ThreadBuf, Meth};
//! Index::build_from_fasta("ref.fa", "ref", false, 4)?;
//! let idx = Index::load("ref", false)?;
//! let opts = Opts::new();
//! let aligner = Aligner::new(&idx, &opts);
//! let mut buf = ThreadBuf::new();
//! for hit in aligner.map(&mut buf, "read1", b"ACGT...", Meth::None)? {
//!     println!(
//!         "{} {}..{} {}",
//!         hit.contig.as_deref().unwrap_or("*"),
//!         hit.ref_start,
//!         hit.ref_end,
//!         hit.cigar_string(),
//!     );
//! }
//! # Ok::<(), minibwa::Error>(())
//! ```
#![forbid(unsafe_op_in_unsafe_fn)]

mod align;
mod error;
mod hit;
mod index;
mod meth;
mod opts;
#[cfg(any(test, feature = "test-util"))]
pub mod test_util;

pub use align::{Aligner, ThreadBuf};
pub use error::{Error, Result};
pub use hit::{CigarKind, CigarOp, Hit, Strand};
pub use index::Index;
pub use meth::Meth;
pub use opts::Opts;
