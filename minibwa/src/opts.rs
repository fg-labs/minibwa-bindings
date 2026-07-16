use crate::error::{Error, Result};
use minibwa_sys as sys;
use std::ffi::CString;

/// Alignment options, wrapping minibwa's `mb_opt_t`.
#[derive(Clone, Debug)]
pub struct Opts(sys::mb_opt_t);

impl Opts {
    /// Default options (minibwa's `mb_opt_init`, "adap" preset).
    pub fn new() -> Self {
        // SAFETY: mb_opt_init fully initializes a zeroed mb_opt_t.
        let mut raw: sys::mb_opt_t = unsafe { std::mem::zeroed() };
        unsafe { sys::mb_opt_init(&mut raw) };
        Opts(raw)
    }

    /// Options initialized from a preset: `"sr"`, `"adap"`, or `"lr"`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidOpts`] if the preset name is unknown or contains
    /// a NUL byte.
    pub fn with_preset(preset: &str) -> Result<Self> {
        let mut opts = Opts::new();
        let c = CString::new(preset).map_err(|_| Error::InvalidOpts("preset has NUL".into()))?;
        // SAFETY: valid opt ptr + NUL-terminated preset string.
        let rc = unsafe { sys::mb_opt_preset(&mut opts.0, c.as_ptr()) };
        if rc != 0 {
            return Err(Error::InvalidOpts(format!("unknown preset {preset:?}")));
        }
        Ok(opts)
    }

    fn set_flag(&mut self, flag: u64, on: bool) {
        if on {
            self.0.flag |= flag;
        } else {
            self.0.flag &= !flag;
        }
    }

    /// Enable/disable paired-end mode (`MB_F_PE`).
    pub fn set_paired(mut self, on: bool) -> Self {
        self.set_flag(sys::MB_F_PE as u64, on);
        self
    }

    /// Enable/disable methylation mode (`MB_F_METH`).
    pub fn set_methylation(mut self, on: bool) -> Self {
        self.set_flag(sys::MB_F_METH as u64, on);
        self
    }

    /// Minimum seed length.
    pub fn set_min_seed_len(mut self, v: i32) -> Self {
        self.0.min_len = v;
        self
    }

    /// Max number of secondary alignments to output.
    pub fn set_out_n(mut self, v: i32) -> Self {
        self.0.out_n = v;
        self
    }

    /// Max number of secondary hits reported in the `XA` tag.
    ///
    /// A hit is eligible for `XA` only if there are at most this many
    /// suboptimal hits scoring above [`set_xa_ratio`](Self::set_xa_ratio) of
    /// the best hit's score. minibwa's default is `5`.
    ///
    /// Setting `v <= 0` disables `XA`-tag output entirely (equivalent to bwa's
    /// `-h 0`); the value is otherwise passed through to minibwa unchecked.
    pub fn set_xa_max(mut self, v: i32) -> Self {
        self.0.xa_max = v;
        self
    }

    /// Score-ratio threshold for `XA`-tag inclusion.
    ///
    /// A suboptimal hit is emitted to the `XA` tag only if its score is at
    /// least `ratio * best_score`. minibwa's default is `0.8`. `ratio` is a
    /// fraction of the best hit's score, so only values in `[0.0, 1.0]` are
    /// meaningful: `0.0` includes every suboptimal hit and `1.0` only those
    /// tied with the best. Non-finite values (`NaN`/`inf`) make the comparison
    /// never hold and silently suppress all `XA` hits. The value is passed
    /// through to minibwa unchecked.
    pub fn set_xa_ratio(mut self, ratio: f32) -> Self {
        self.0.xa_ratio = ratio;
        self
    }

    /// Match score (Smith-Waterman `a` parameter).
    pub fn set_match_score(mut self, score: i32) -> Self {
        self.0.a = score;
        self
    }

    /// Mismatch penalty (Smith-Waterman `b` parameter).
    pub fn set_mismatch_penalty(mut self, penalty: i32) -> Self {
        self.0.b = penalty;
        self
    }

    /// Gap-open penalty (`q` parameter).
    pub fn set_gap_open(mut self, open: i32) -> Self {
        self.0.q = open;
        self
    }

    /// Gap-extend penalty (`e` parameter).
    pub fn set_gap_extend(mut self, extend: i32) -> Self {
        self.0.e = extend;
        self
    }

    /// Paired-end insert-size parameters.
    ///
    /// * `avg` — expected insert-size mean
    /// * `std` — expected insert-size standard deviation
    /// * `lo`  — lower insert-size bound for proper-pair classification
    /// * `hi`  — upper insert-size bound for proper-pair classification
    pub fn set_pe_insert_size(mut self, avg: i32, std: i32, lo: i32, hi: i32) -> Self {
        self.0.pe_avg = avg;
        self.0.pe_std = std;
        self.0.pe_lo = lo;
        self.0.pe_hi = hi;
        self
    }

    pub fn is_paired(&self) -> bool {
        self.0.flag & sys::MB_F_PE as u64 != 0
    }

    pub fn is_methylation(&self) -> bool {
        self.0.flag & sys::MB_F_METH as u64 != 0
    }

    /// A copy of the raw options with paired-end (`MB_F_PE`) forced on.
    pub(crate) fn paired_copy(&self) -> minibwa_sys::mb_opt_t {
        let mut o = self.0;
        o.flag |= minibwa_sys::MB_F_PE as u64;
        o
    }

    pub(crate) fn as_ptr(&self) -> *const sys::mb_opt_t {
        &self.0
    }
}

impl Default for Opts {
    fn default() -> Self {
        Opts::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preset_and_flags() {
        let o = Opts::with_preset("sr").unwrap();
        let o = o.set_paired(true).set_methylation(true);
        assert!(o.is_paired());
        assert!(o.is_methylation());
        assert!(Opts::with_preset("nope").is_err());
    }

    #[test]
    fn preset_nul_byte_returns_err() {
        assert!(matches!(
            Opts::with_preset("sr\0x"),
            Err(crate::Error::InvalidOpts(_))
        ));
    }

    #[test]
    fn all_presets_succeed() {
        assert!(Opts::with_preset("adap").is_ok());
        assert!(Opts::with_preset("lr").is_ok());
        assert!(Opts::with_preset("sr").is_ok());
        // "lr" is a long-read preset — paired-end is irrelevant and should be off by default.
        assert!(!Opts::with_preset("lr").unwrap().is_paired());
    }

    #[test]
    fn xa_defaults_and_setters() {
        // minibwa's mb_opt_init defaults (options.c): xa_max = 5, xa_ratio = 0.8.
        let o = Opts::new();
        assert_eq!(o.0.xa_max, 5);
        assert!((o.0.xa_ratio - 0.8).abs() < 1e-6);
        let o = o.set_xa_max(12).set_xa_ratio(0.5);
        assert_eq!(o.0.xa_max, 12);
        assert!((o.0.xa_ratio - 0.5).abs() < 1e-6);
    }

    #[test]
    fn flags_toggle_off() {
        let o = Opts::new().set_paired(true).set_methylation(true);
        assert!(o.is_paired());
        assert!(o.is_methylation());
        let o = o.set_paired(false).set_methylation(false);
        assert!(!o.is_paired());
        assert!(!o.is_methylation());
    }
}
