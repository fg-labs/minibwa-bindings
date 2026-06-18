/// Methylation conversion type for a read, mapped to minibwa's `mt` parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Meth {
    /// Unmethylated / normal alignment.
    None,
    /// Read 1: C-to-T converted strand.
    C2T,
    /// Read 2: G-to-A converted strand.
    G2A,
}

impl Meth {
    pub(crate) fn as_mt(self) -> i32 {
        match self {
            Meth::None => 0,
            Meth::C2T => 1,
            Meth::G2A => 2,
        }
    }
}
