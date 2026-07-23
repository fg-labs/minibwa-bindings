#ifndef MINIBWA_SHIM_H
#define MINIBWA_SHIM_H

#include <stdint.h>
#include "minibwa.h" /* mb_idx_t (opaque) */

#ifdef __cplusplus
extern "C" {
#endif

/* Build a minibwa index from a FASTA. Returns 0 on success, nonzero on error.
 * On is_meth != 0, also writes <prefix>.meth.mbw.
 * NOTE: minibwa's main_index() may abort() on grossly invalid input
 * (e.g. unreadable FASTA); callers must pre-validate paths in Rust. */
int mb_index_build(const char *fasta_path, const char *prefix, int is_meth, int n_threads);

/* Last error message for the current thread (empty string if none). */
const char *mb_shim_last_error(void);

/* --- Parity-oracle surface + FM-index buffer accessors --- */
/* These expose minibwa's internals so an independent implementation can be diffed against the C
 * code path-for-path, and so throughput can be measured against minibwa's own thread pool. */

/* Two-round seeding oracle (raw mb_seed_intv; no dedup). Returns the true seed count. */
int mb_shim_seed_intv(const mb_idx_t *idx, const uint8_t *seq, int len,
                      int min_len, int max_sub_occ,
                      uint32_t *out_qs, uint32_t *out_qe,
                      uint64_t *out_size, uint64_t *out_x0, int max_out);

/* CPU seeding-throughput baseline: run C two-round mb_seed_intv over the whole read batch
 * using minibwa's own thread pool (kt_for, per-thread kalloc pools — no malloc contention).
 * Reads are nt4-encoded, concatenated (seq) with per-read off/len. Returns elapsed
 * nanoseconds; *out_total_seeds = total seeds (a correctness sanity). This is the fair CPU
 * baseline an alternative seeder should be measured against. */
int64_t mb_shim_seed_bench(const mb_idx_t *idx, int n_reads,
                           const uint8_t *seq, const uint64_t *off, const uint32_t *len,
                           int min_len, int max_sub_occ, int n_threads, uint64_t *out_total_seeds);

/* End-to-end mapping baseline: full C minibwa pipeline (seed→chain→extend→pair) over the whole
 * read batch via mb_map + kt_for, per-thread mb_tbuf. Reads are raw ASCII (concatenated, with
 * per-read off/len). Returns elapsed nanoseconds; *out_total_hits = total alignments,
 * *out_total_cigar = total CIGAR ops. The end-to-end speed baseline for an alternative pipeline.
 * `fset` is OR'd into opt.flag AFTER mb_opt_init, so the caller can select a narrower scope:
 * MB_F_NO_ALN (--chain-only) stops before base alignment, which lets a partial implementation be
 * compared like-for-like rather than against the full pipeline.
 * Returns -1 on allocation failure. */
int64_t mb_shim_map_bench(const mb_idx_t *idx, int n_reads,
                          const char *seq, const uint64_t *off, const uint32_t *len,
                          int n_threads, uint64_t fset,
                          uint64_t *out_total_hits, uint64_t *out_total_cigar);

/* Smith-Waterman extension oracle: C ksw_extz2_sse (single-gap affine, banded, z-drop) on an
 * nt4 query/target pair; 5x5 matrix from (a, b, b_ts, b_ambi). out[8] = {max, max_q, max_t, mqe,
 * mqe_t, mte, mte_q, score}; *out_flags = zdropped | reach_end<<1. Returns n_cigar (cigar[] filled
 * up to max_cigar). Bit-exact reference for the Rust extz2 port. */
int mb_shim_extz2(const uint8_t *query, int qlen, const uint8_t *target, int tlen,
                  int8_t a, int8_t b, int8_t b_ts, int8_t b_ambi,
                  int8_t gapo, int8_t gape, int w, int zdrop, int end_bonus, int flag,
                  int32_t out[8], uint32_t *out_cigar, int max_cigar, int32_t *out_flags);

/* Dual-gap SW extension oracle: C ksw_extd2_sse (gap = min(q+len*e, q2+len*e2)). Same out/return
 * convention as mb_shim_extz2. Bit-exact reference for the Rust extd2 port. */
int mb_shim_extd2(const uint8_t *query, int qlen, const uint8_t *target, int tlen,
                  int8_t a, int8_t b, int8_t b_ts, int8_t b_ambi,
                  int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2,
                  int w, int zdrop, int end_bonus, int flag,
                  int32_t out[8], uint32_t *out_cigar, int max_cigar, int32_t *out_flags);

/* Anchor-generation oracle (mb_anchor, non-methylation path). Consumes SA-interval seeds
 * (x0/size/info), returns anchors (sid/len/qpos/tpos). Returns the anchor count, or -1 if
 * n_seed is negative/oversized or the seed buffer cannot be allocated. */
int mb_shim_anchor(const mb_idx_t *idx,
                   int n_seed, const uint64_t *in_x0, const uint64_t *in_size, const uint64_t *in_info,
                   int min_len, int qlen, const uint8_t *qseq, int max_occ,
                   int32_t *out_sid, int32_t *out_len, int32_t *out_qpos, int64_t *out_tpos, int max_a);

/* Chaining DP oracle (mb_lchain_dp). Feeds identical anchors, returns chains: out_u[i] =
 * (score<<32)|count, plus the compacted anchor array. Returns #chains; *out_n_a = #anchors.
 * Returns -1 if n is negative/oversized or the anchor buffer cannot be allocated. */
int mb_shim_lchain_dp(const mb_idx_t *idx,
                      int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter,
                      int min_sc, float chn_pen_gap,
                      int64_t n, const int32_t *in_sid, const int32_t *in_len,
                      const int32_t *in_qpos, const int64_t *in_tpos,
                      uint64_t *out_u, int max_u,
                      int32_t *out_sid, int32_t *out_len, int32_t *out_qpos, int64_t *out_tpos,
                      int max_a, int *out_n_a);

/* End-to-end alignment-record oracle: full C mb_map per read, returning the PRIMARY hit as a
 * record. Per read i: out_tid[i] (-1 if unmapped); out_pos[i*2..]=ts,te; out_q[i*2..]=qs,qe;
 * out_misc[i*4..]=flags(rev|pri<<1),mapq,dp_score,n_cigar; out_cigar[i*cig_stride..]=CIGAR ops.
 * Reads are raw ASCII concat + off/len. Returns 0, or -1 if any read's CIGAR exceeds cig_stride. */
int mb_shim_map(const mb_idx_t *idx, int n_reads,
                const char *seq, const uint64_t *off, const uint32_t *len, int n_threads,
                int64_t *out_tid, int64_t *out_pos, int32_t *out_q, int32_t *out_misc,
                uint32_t *out_cigar, int cig_stride);

/* Per-region oracle (single read): full C mb_map, returning ALL n_hit regions (primary +
 * secondaries), not just the primary. Diagnostic for the long-read align1 bisect. Per region j:
 * out_tid[j]; out_ts[j]/out_te[j]; out_qs[j]/out_qe[j]; out_dp[j]=p->dp_score; out_dpmax[j]=p->dp_max;
 * out_mapq[j]; out_flags[j]=(rev|sam_pri<<1); out_parent[j]; out_ncigar[j];
 * out_cigar[j*cig_stride..]=CIGAR ops. Returns n_hit, or -1 on overflow (region count > max_regions
 * or any CIGAR > cig_stride). */
int mb_shim_map_regions(const mb_idx_t *idx, const char *seq, int len,
                        int64_t *out_tid, int64_t *out_ts, int64_t *out_te,
                        int32_t *out_qs, int32_t *out_qe, int32_t *out_dp, int32_t *out_dpmax,
                        int32_t *out_mapq, int32_t *out_flags, int32_t *out_parent,
                        int32_t *out_ncigar, uint32_t *out_cigar, int cig_stride, int max_regions);

/* Paired-end oracle: C mb_map on 2*n_frag interleaved reads + mb_pestat + mb_pair per fragment
 * (incl. mate rescue). Per mate m (2*n_frag): out_tid[m], out_pos[m*2]=ts,te, out_q[m*2]=qs,qe,
 * out_misc[m*4]={flags(rev|proper<<1|pri<<2),mapq,dp_score,n_cigar}, out_cigar[m*cig_stride].
 * Returns 0, or -1 on failure/overflow. */
int mb_shim_pair(const mb_idx_t *idx, int n_frag,
                 const char *seq, const uint64_t *off, const uint32_t *len,
                 int64_t *out_tid, int64_t *out_pos, int32_t *out_q, int32_t *out_misc,
                 uint32_t *out_cigar, int cig_stride);

/* Insert-size estimation oracle: C mb_map on 2*n_frag interleaved reads + mb_pestat. Per
 * orientation d in 0..4: out_i[d*3+{0,1,2}]={failed,lo,hi}, out_d[d*2+{0,1}]={avg,std}.
 * out_i >= 12 int32, out_d >= 8 double. Returns 0, or -1 on failure. */
int mb_shim_pestat(const mb_idx_t *idx, int n_frag,
                   const char *seq, const uint64_t *off, const uint32_t *len,
                   int32_t *out_i, double *out_d);

/* SAM-header oracle: C mb_fmt_sam_hdr (no read-group / command line), version = MB_VERSION.
 * Writes header text into out_buf; returns its length, or -1 on overflow. */
int mb_shim_sam_hdr(const mb_idx_t *idx, char *out_buf, uint64_t out_cap);

/* SAM-output oracle: run C mb_map per read and format its hits with mb_fmt_sam (single-end),
 * concatenating per-read SAM text into out_buf; out_off[i..i+1) bounds read i's block (out_off
 * has n_reads+1 entries). `names[i]` is the QNAME. Returns 0, or -1 on out_buf overflow. */
int mb_shim_map_sam(const mb_idx_t *idx, int n_reads,
                    const char *const *names, const char *seq, const uint64_t *off, const uint32_t *len,
                    int paf, int64_t extra_flag,
                    const char *preset, const int64_t *ov, int n_ov, int has_pri, float pri,
                    uint64_t fset, uint64_t fclr,
                    char *out_buf, uint64_t out_cap, uint64_t *out_off);

/* ksw_ll oracle: C ksw_ll_qinit + ksw_ll_{u8,i16}_core (size 1 or 2). out[5] = {score, te, qe,
 * score2, te2}. Bit-exact reference for the Rust ksw_ll port. */
void mb_shim_ksw_ll(int size, const uint8_t *query, int qlen, const uint8_t *target, int tlen,
                    int8_t a, int8_t b, int8_t b_ts, int8_t b_ambi,
                    int gapo, int gape, int xtra, int32_t out[5]);

/* Paired-end SAM oracle: mb_map + mb_pestat + mb_pair per fragment, then format both mates with
 * mb_fmt_sam (n_seg=2). Per fragment f, both mates' SAM lands at out_buf[out_off[f]..out_off[f+1]);
 * out_off has n_frag+1 entries. Returns 0, or -1 on failure/overflow. */
int mb_shim_pair_sam(const mb_idx_t *idx, int n_frag,
                     const char *const *names, const char *seq, const uint64_t *off, const uint32_t *len,
                     const char *preset, const int64_t *ov, int n_ov, int has_pri, float pri,
                     uint64_t fset, uint64_t fclr,
                     char *out_buf, uint64_t out_cap, uint64_t *out_off);

/* Methylation paired-end oracle: mb_map_sai(mt) per mate on a meth index + mb_pestat + mb_pair
 * (MB_F_METH; mate rescue uses the meth path). Per mate m: out_tid/pos/q/misc/cigar. 0 or -1. */
int mb_shim_pair_meth(const mb_idx_t *idx, int n_frag,
                      const char *seq, const uint64_t *off, const uint32_t *len,
                      int64_t *out_tid, int64_t *out_pos, int32_t *out_q, int32_t *out_misc,
                      uint32_t *out_cigar, int cig_stride);

/* Methylation mapping oracle: seed the C2T/G2A-converted read on a meth index, then mb_map_sai on
 * the original read + mt. Returns the primary hit (out_tid/pos/q/misc/cigar). Returns 0 or -1. */
int mb_shim_map_meth(const mb_idx_t *idx, const char *read, int len, int mt, int cig_stride,
                     int64_t *out_tid, int64_t *out_pos, int32_t *out_q, int32_t *out_misc,
                     uint32_t *out_cigar);

/* SE methylation SAM oracle: mb_shim_map_sam on a meth index — C2T-convert (mt=1), seed converted,
 * mb_map_sai on the original read + mt, format with mb_fmt_sam. MB_F_METH forced on; same override/
 * flag handling and out_buf/out_off layout as mb_shim_map_sam. Returns 0, or -1 on overflow. */
int mb_shim_map_meth_sam(const mb_idx_t *idx, int n_reads,
                         const char *const *names, const char *seq, const uint64_t *off, const uint32_t *len,
                         const char *preset, const int64_t *ov, int n_ov, int has_pri, float pri,
                         uint64_t fset, uint64_t fclr,
                         char *out_buf, uint64_t out_cap, uint64_t *out_off);

/* PE methylation SAM oracle: mb_shim_pair_sam on a meth index — R1->C2T/R2->G2A convert, seed
 * converted, mb_map_sai + mb_pestat + mb_pair (MB_F_METH), format both mates with mb_fmt_sam
 * (n_seg=2). Honors MB_F_NO_PAIRING. Same signature/layout as mb_shim_pair_sam. Returns 0, or -1. */
int mb_shim_pair_meth_sam(const mb_idx_t *idx, int n_frag,
                          const char *const *names, const char *seq, const uint64_t *off, const uint32_t *len,
                          const char *preset, const int64_t *ov, int n_ov, int has_pri, float pri,
                          uint64_t fset, uint64_t fclr,
                          char *out_buf, uint64_t out_cap, uint64_t *out_off);

/* Debug oracle: dump ALL hits of one read. Per hit k: out[k*8+{0..7}] = tid, ts, te, qs, qe,
 * flags(rev|inv<<1|sam_pri<<2|(parent==id)<<3), mapq, dp_score. `out` is filled for the first
 * max_hits hits; returns the true n_hit (may exceed max_hits), or -1 on failure. */
int mb_shim_map_hits(const mb_idx_t *idx, const char *seq, int len, int max_hits, int64_t *out);

/* Length-adapted option oracle: dump the exact options the C pipeline uses for a read of length
 * `qlen` (mb_opt_init + mb_opt_adap). See the .c for the fixed field order. `oi` >=29 int32,
 * `of` >=3 float, `ol` >=2 int64. */
void mb_shim_align_opt(int qlen, int32_t *oi, float *of, int64_t *ol);

/* Override-aware mb_shim_align_opt: apply CLI overrides (preset -x, the OV_* scalar array with
 * MB_SHIM_OVR_KEEP=INT64_MIN sentinels, the -p float when has_pri, and flag set/clear masks) onto
 * mb_opt_init before mb_opt_adap. `n_ov` must be >= OV_N (15) for the array to be applied. Same
 * oi/of/ol output layout as mb_shim_align_opt. Layout mirrors cli::OptOverride::encode in Rust. */
void mb_shim_align_opt_ovr(int qlen, const char *preset, const int64_t *ov, int n_ov,
                           int has_pri, float pri, uint64_t fset, uint64_t fclr,
                           int32_t *oi, float *of, int64_t *ol);

/* Reference-sequence extraction oracle: C l2b_getseq (nt4, N=4) for contig tid, local [st,en),
 * into out; returns bases written. Bit-exact reference for the Rust L2b::getseq port. */
int64_t mb_shim_getseq(const mb_idx_t *idx, int64_t tid, int64_t st, int64_t en, uint8_t *out);

/* Chaining gap-penalty oracle: exactly comput_sc's `(int)(chn_pen_gap*dd + .5f*mb_log2(dd+1))`
 * (lchain.c:148-150), the single FMA-sensitive site. Bit-exact reference for the Rust
 * `mbwa_core::chain::gap_penalty` and the no-FMA device gate. dd==0 -> 0. */
int32_t mb_shim_gap_penalty(float chn_pen_gap, int64_t dd);

/* FM-index buffer accessors: hand out the C-owned buffers so a caller can copy them elsewhere
 * (e.g. to another device) without reloading or duplicating the index. */
const uint64_t *mb_shim_bwt_data(const mb_idx_t *idx);
uint64_t mb_shim_bwt_data_len(const mb_idx_t *idx);
uint64_t mb_shim_bwt_primary(const mb_idx_t *idx);
uint64_t mb_shim_bwt_seq_len(const mb_idx_t *idx);
const uint64_t *mb_shim_bwt_sa(const mb_idx_t *idx);
uint64_t mb_shim_bwt_n_sa(const mb_idx_t *idx);
uint32_t mb_shim_bwt_sa_bit(const mb_idx_t *idx);
const uint32_t *mb_shim_bwt_cnt_table(const mb_idx_t *idx);
void mb_shim_bwt_l2(const mb_idx_t *idx, uint64_t out[5]);

#ifdef __cplusplus
}
#endif

#endif
