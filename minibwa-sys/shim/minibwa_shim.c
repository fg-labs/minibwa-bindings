/* clock_gettime / struct timespec / CLOCK_MONOTONIC are POSIX; glibc hides them under
 * -std=c99 unless a feature-test macro is set. macOS exposes them by default, and a strict
 * _POSIX_C_SOURCE there instead hides C99 (snprintf), so scope this to non-Apple. Must
 * precede all includes. */
#if !defined(__APPLE__)
#define _POSIX_C_SOURCE 200809L
#endif
#include "minibwa_shim.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#define HAVE_KALLOC  /* make ksw2.h use the real kalloc funcs (kalloc.h is included below) */
#include "mbpriv.h"  /* internal: mb_idx_s layout, mb_seed_intv, mb_bwt_t, mb_sai_v */
#include "kthread.h" /* kt_for */
#include "kalloc.h"  /* km_init / km_destroy / kfree (per-thread pools) */
#include "ksw2.h"    /* ksw_extz2_sse, ksw_gen_nt4_mat, ksw_extz_t (SW extension oracle) */
#include "l2bit.h"   /* l2b_getseq (reference-extraction oracle) */

/* Override layout for the *_ovr shim entry points: mb_opt_t base scalars a user can set from the
 * CLI, applied after mb_opt_init(+preset) and before mb_opt_adap. MB_SHIM_OVR_KEEP = keep the
 * init/preset default. This index order is shared verbatim with the Rust cli::OptOverride encoder
 * (cli.rs: EncodedOpt::ov); keep the two in sync. Declared here (before the oracle fns that use it)
 * so C99 sees OV_N + the shim_apply_ovr prototype; the definition is further down. */
#define MB_SHIM_OVR_KEEP INT64_MIN
enum { OV_A = 0, OV_B, OV_MIN_LEN, OV_MAX_OCC, OV_MIN_CHAIN, OV_BEST_N, OV_PEN_UNPAIR,
       OV_MAX_SR_LEN, OV_MAX_GAP, OV_BW, OV_BW_LONG, OV_MIN_DP_MAX, OV_OUT_N, OV_MAX_RESCUE,
       OV_XA_MAX, OV_N };
static void shim_apply_ovr(mb_opt_t *opt, const char *preset, const int64_t *ov,
                           int has_pri, float pri, uint64_t fset, uint64_t fclr);
#include "bseq.h"    /* mb_bseq1_t (SAM formatting oracle) */

/* Provided by minibwa's index.c (not declared in minibwa.h). */
extern int main_index(int argc, char **argv);

static __thread char g_err[512];

static void shim_set_err(const char *msg) {
    snprintf(g_err, sizeof(g_err), "%s", msg ? msg : "");
}

const char *mb_shim_last_error(void) { return g_err; }

/* Fill tbuf[0..n) with fresh thread buffers for the kt_for worker pool. On failure, destroy the
 * partial set and return -1 — the workers index tbuf[tid] unconditionally, so a NULL slot would
 * fault. mb_tbuf_destroy dereferences its argument, hence the bound at `t` rather than `n`.
 * (The km_init pools in mb_shim_seed_bench need no such guard: km == 0 is minibwa's "no pool"
 * sentinel and kalloc falls back to malloc, so a failed init degrades instead of faulting.) */
static int shim_tbuf_alloc(mb_tbuf_t **tbuf, int n) {
    for (int t = 0; t < n; ++t) {
        tbuf[t] = mb_tbuf_init(0);
        if (!tbuf[t]) {
            for (int u = 0; u < t; ++u) mb_tbuf_destroy(tbuf[u]);
            return -1;
        }
    }
    return 0;
}

int mb_index_build(const char *fasta_path, const char *prefix, int is_meth, int n_threads) {
    g_err[0] = '\0';
    if (!fasta_path || !prefix) { shim_set_err("null fasta_path or prefix"); return 2; }

    char nthr[32];
    snprintf(nthr, sizeof(nthr), "%d", n_threads > 0 ? n_threads : 1);

    /* argv: index -t <n> [--meth] <fasta> <prefix> */
    char *argv[8];
    int argc = 0;
    argv[argc++] = "index";
    argv[argc++] = "-t";
    argv[argc++] = nthr;
    if (is_meth) argv[argc++] = "--meth";
    argv[argc++] = (char *)fasta_path;
    argv[argc++] = (char *)prefix;
    argv[argc] = NULL;

    int rc = main_index(argc, argv);
    if (rc != 0) shim_set_err("minibwa index build failed");
    return rc;
}

/* ---- Parity-oracle surface + FM-index buffer accessors (see header) ---- */

/* CPU seeding-throughput baseline (see header). Per-read fresh mb_sai_v allocated from a
 * per-thread kalloc pool, freed back after each read — mirrors production threading. */
typedef struct {
    const mb_idx_t *idx;
    const uint8_t *seq;
    const uint64_t *off;
    const uint32_t *len;
    int min_len, max_sub_occ;
    void **km;              /* one kalloc pool per thread */
    uint64_t *seed_counts;  /* per read */
} seed_bench_t;

static void seed_bench_worker(void *data, long i, int tid) {
    seed_bench_t *d = (seed_bench_t *)data;
    mb_sai_v v;
    v.n = 0; v.m = 0; v.a = 0;
    mb_seed_intv(d->km[tid], d->idx->bwt, (int32_t)d->len[i], d->seq + d->off[i],
                 d->min_len, d->max_sub_occ, &v);
    d->seed_counts[i] = (uint64_t)v.n;
    kfree(d->km[tid], v.a); /* return to the thread's pool for reuse */
}

int64_t mb_shim_seed_bench(const mb_idx_t *idx, int n_reads,
                           const uint8_t *seq, const uint64_t *off, const uint32_t *len,
                           int min_len, int max_sub_occ, int n_threads, uint64_t *out_total_seeds) {
    if (n_threads < 1) n_threads = 1;
    seed_bench_t d;
    d.idx = idx; d.seq = seq; d.off = off; d.len = len;
    d.min_len = min_len; d.max_sub_occ = max_sub_occ;
    d.seed_counts = (uint64_t *)calloc(n_reads, sizeof(uint64_t));
    d.km = (void **)calloc(n_threads, sizeof(void *));
    if (!d.seed_counts || !d.km) { /* OOM: free what we got and signal failure (-1 ns). */
        free(d.seed_counts); free(d.km);
        if (out_total_seeds) *out_total_seeds = 0;
        shim_set_err("mb_shim_seed_bench: calloc failed");
        return -1;
    }
    for (int t = 0; t < n_threads; ++t) d.km[t] = km_init();

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    kt_for(n_threads, seed_bench_worker, &d, (long)n_reads);
    clock_gettime(CLOCK_MONOTONIC, &t1);

    uint64_t total = 0;
    for (int i = 0; i < n_reads; ++i) total += d.seed_counts[i];
    if (out_total_seeds) *out_total_seeds = total;

    for (int t = 0; t < n_threads; ++t) km_destroy(d.km[t]);
    free(d.km);
    free(d.seed_counts);
    return (int64_t)(t1.tv_sec - t0.tv_sec) * 1000000000LL + (int64_t)(t1.tv_nsec - t0.tv_nsec);
}

/* End-to-end mapping baseline: run the FULL C minibwa pipeline (seed→chain→extend→pair) per
 * read via mb_map, over the whole batch through minibwa's thread pool (kt_for) with a per-thread
 * mb_tbuf. Reads are raw ASCII (mb_map nt4-encodes internally), concatenated with per-read
 * off/len. Returns elapsed nanoseconds; *out_total_hits = total alignment records, *out_total_cigar
 * = total CIGAR ops (a coarse correctness signal). This is the end-to-end speed baseline an
 * alternative pipeline is measured against; the per-stage cost split and the full alignment-record
 * oracle are separate entry points. */
typedef struct {
    const mb_opt_t *opt;
    const mb_idx_t *idx;
    const char *seq;       /* raw ASCII, concatenated */
    const uint64_t *off;
    const uint32_t *len;
    mb_tbuf_t **tbuf;      /* one thread buffer per worker thread */
    uint64_t *n_hits;      /* per read */
    uint64_t *n_cigar;     /* per read: sum of CIGAR ops over that read's hits */
} map_bench_t;

static void map_bench_worker(void *data, long i, int tid) {
    map_bench_t *d = (map_bench_t *)data;
    int32_t n_hit = 0;
    mb_hit_t *hit = mb_map(d->opt, d->idx, (int32_t)d->len[i], d->seq + d->off[i],
                           L2B_METH_NONE, &n_hit, d->tbuf[tid], 0);
    uint64_t tc = 0;
    for (int32_t j = 0; j < n_hit; ++j) {
        if (hit[j].p) { tc += (uint64_t)hit[j].p->n_cigar; free(hit[j].p); }
    }
    d->n_hits[i] = (uint64_t)n_hit;
    d->n_cigar[i] = tc;
    free(hit); /* mb_map returns malloc'd hit[] + malloc'd .p (see map-main.c free convention) */
}

int64_t mb_shim_map_bench(const mb_idx_t *idx, int n_reads,
                          const char *seq, const uint64_t *off, const uint32_t *len,
                          int n_threads, uint64_t fset,
                          uint64_t *out_total_hits, uint64_t *out_total_cigar) {
    if (n_threads < 1) n_threads = 1;
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = n_threads;
    opt.flag |= (int64_t)fset; /* e.g. MB_F_NO_ALN for the --chain-only scope */

    map_bench_t d;
    d.opt = &opt; d.idx = idx; d.seq = seq; d.off = off; d.len = len;
    d.n_hits = (uint64_t *)calloc(n_reads, sizeof(uint64_t));
    d.n_cigar = (uint64_t *)calloc(n_reads, sizeof(uint64_t));
    d.tbuf = (mb_tbuf_t **)calloc(n_threads, sizeof(mb_tbuf_t *));
    if (!d.n_hits || !d.n_cigar || !d.tbuf) {
        free(d.n_hits); free(d.n_cigar); free(d.tbuf);
        if (out_total_hits) *out_total_hits = 0;
        if (out_total_cigar) *out_total_cigar = 0;
        shim_set_err("mb_shim_map_bench: calloc failed");
        return -1;
    }
    if (shim_tbuf_alloc(d.tbuf, n_threads) < 0) {
        free(d.n_hits); free(d.n_cigar); free(d.tbuf);
        if (out_total_hits) *out_total_hits = 0;
        if (out_total_cigar) *out_total_cigar = 0;
        shim_set_err("mb_shim_map_bench: tbuf init failed");
        return -1;
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    kt_for(n_threads, map_bench_worker, &d, (long)n_reads);
    clock_gettime(CLOCK_MONOTONIC, &t1);

    uint64_t th = 0, tc = 0;
    for (int i = 0; i < n_reads; ++i) { th += d.n_hits[i]; tc += d.n_cigar[i]; }
    if (out_total_hits) *out_total_hits = th;
    if (out_total_cigar) *out_total_cigar = tc;

    for (int t = 0; t < n_threads; ++t) mb_tbuf_destroy(d.tbuf[t]);
    free(d.tbuf); free(d.n_hits); free(d.n_cigar);
    return (int64_t)(t1.tv_sec - t0.tv_sec) * 1000000000LL + (int64_t)(t1.tv_nsec - t0.tv_nsec);
}

/* Smith-Waterman extension oracle: run C ksw_extz2_sse (single-gap affine, banded, z-drop) on an
 * nt4-encoded query/target pair, building the 5x5 score matrix via ksw_gen_nt4_mat (non-meth).
 * Returns n_cigar (true count; cigar[] filled up to max_cigar). The scalar fields land in out[8] =
 * {max, max_q, max_t, mqe, mqe_t, mte, mte_q, score}; *out_flags = zdropped | reach_end<<1. This is
 * the bit-exact reference for the Rust extz2 port. */
int mb_shim_extz2(const uint8_t *query, int qlen, const uint8_t *target, int tlen,
                  int8_t a, int8_t b, int8_t b_ts, int8_t b_ambi,
                  int8_t gapo, int8_t gape, int w, int zdrop, int end_bonus, int flag,
                  int32_t out[8], uint32_t *out_cigar, int max_cigar, int32_t *out_flags) {
    int8_t mat[25];
    ksw_gen_nt4_mat(mat, a, b, b_ts, b_ambi, 0);
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ez));
    ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gapo, gape, w, zdrop, end_bonus, flag, &ez);
    out[0] = ez.max; out[1] = ez.max_q; out[2] = ez.max_t; out[3] = ez.mqe;
    out[4] = ez.mqe_t; out[5] = ez.mte; out[6] = ez.mte_q; out[7] = ez.score;
    if (out_flags) *out_flags = (int32_t)ez.zdropped | ((int32_t)ez.reach_end << 1);
    int n = ez.n_cigar;
    for (int i = 0; i < n && i < max_cigar; ++i) out_cigar[i] = ez.cigar[i];
    free(ez.cigar); /* km=0 path uses malloc/realloc, so free() is correct */
    return n;
}

/* Dual-gap SW extension oracle: C ksw_extd2_sse (gap = min(q+len*e, q2+len*e2), banded, z-drop).
 * Same out/return convention as mb_shim_extz2. Bit-exact reference for the Rust extd2 port. */
int mb_shim_extd2(const uint8_t *query, int qlen, const uint8_t *target, int tlen,
                  int8_t a, int8_t b, int8_t b_ts, int8_t b_ambi,
                  int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2,
                  int w, int zdrop, int end_bonus, int flag,
                  int32_t out[8], uint32_t *out_cigar, int max_cigar, int32_t *out_flags) {
    int8_t mat[25];
    ksw_gen_nt4_mat(mat, a, b, b_ts, b_ambi, 0);
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ez));
    ksw_extd2_sse(0, qlen, query, tlen, target, 5, mat, gapo, gape, gapo2, gape2, w, zdrop, end_bonus, flag, &ez);
    out[0] = ez.max; out[1] = ez.max_q; out[2] = ez.max_t; out[3] = ez.mqe;
    out[4] = ez.mqe_t; out[5] = ez.mte; out[6] = ez.mte_q; out[7] = ez.score;
    if (out_flags) *out_flags = (int32_t)ez.zdropped | ((int32_t)ez.reach_end << 1);
    int n = ez.n_cigar;
    for (int i = 0; i < n && i < max_cigar; ++i) out_cigar[i] = ez.cigar[i];
    free(ez.cigar);
    return n;
}

/* Two-round seeding oracle (raw mb_seed_intv output; no dedup — matches C exactly).
 * Fills up to max_out seeds; returns the true total (may exceed max_out). */
int mb_shim_seed_intv(const mb_idx_t *idx, const uint8_t *seq, int len,
                      int min_len, int max_sub_occ,
                      uint32_t *out_qs, uint32_t *out_qe,
                      uint64_t *out_size, uint64_t *out_x0, int max_out) {
    mb_sai_v v;
    v.n = 0; v.m = 0; v.a = 0;
    mb_seed_intv(0, idx->bwt, len, seq, min_len, max_sub_occ, &v);
    int n = (int)v.n;
    for (int i = 0; i < n && i < max_out; ++i) {
        out_qs[i] = (uint32_t)(v.a[i].info >> 32);
        out_qe[i] = (uint32_t)(v.a[i].info & 0xffffffffu);
        out_size[i] = v.a[i].size;
        out_x0[i] = v.a[i].x[0];
    }
    free(v.a); /* km=0 path uses malloc/realloc, so free() is correct */
    return n;
}

/* Anchor-generation oracle: run mb_anchor (non-meth path) on caller-supplied SA-interval
 * seeds and return the resulting anchors. mb_anchor mutates `u` in place (sort/dedup). */
int mb_shim_anchor(const mb_idx_t *idx,
                   int n_seed, const uint64_t *in_x0, const uint64_t *in_size, const uint64_t *in_info,
                   int min_len, int qlen, const uint8_t *qseq, int max_occ,
                   int32_t *out_sid, int32_t *out_len, int32_t *out_qpos, int64_t *out_tpos, int max_a) {
    if (n_seed < 0 || (uint64_t)n_seed > SIZE_MAX / sizeof(mb_sai_t)) {
        shim_set_err("mb_shim_anchor: invalid seed count");
        return -1;
    }
    mb_sai_v u;
    u.n = n_seed; u.m = n_seed; u.a = 0;
    if (n_seed > 0) {
        u.a = (mb_sai_t *)malloc((size_t)n_seed * sizeof(mb_sai_t));
        if (!u.a) { shim_set_err("mb_shim_anchor: malloc failed"); return -1; }
        for (int i = 0; i < n_seed; ++i) {
            u.a[i].x[0] = in_x0[i]; u.a[i].x[1] = 0;
            u.a[i].size = in_size[i]; u.a[i].info = in_info[i];
        }
    }
    mb_anchor_v v; v.n = 0; v.m = 0; v.a = 0;
    (void)mb_anchor(0, idx, &u, min_len, qlen, qseq, L2B_METH_NONE, max_occ, &v);
    int n = (int)v.n;
    for (int i = 0; i < n && i < max_a; ++i) {
        out_sid[i] = v.a[i].sid; out_len[i] = v.a[i].len;
        out_qpos[i] = v.a[i].qpos; out_tpos[i] = v.a[i].tpos;
    }
    free(u.a); free(v.a); /* km=0 path uses malloc, so free() is correct */
    return n;
}

/* Chaining DP oracle: build an mb_anchor_t[] from the input columns, run mb_lchain_dp
 * (km=0 → malloc/free), and return the compacted chains + anchors. */
int mb_shim_lchain_dp(const mb_idx_t *idx,
                      int max_dist_x, int max_dist_y, int bw, int max_skip, int max_iter,
                      int min_sc, float chn_pen_gap,
                      int64_t n, const int32_t *in_sid, const int32_t *in_len,
                      const int32_t *in_qpos, const int64_t *in_tpos,
                      uint64_t *out_u, int max_u,
                      int32_t *out_sid, int32_t *out_len, int32_t *out_qpos, int64_t *out_tpos,
                      int max_a, int *out_n_a) {
    *out_n_a = 0;
    if (n < 0 || (uint64_t)n > SIZE_MAX / sizeof(mb_anchor_t)) {
        shim_set_err("mb_shim_lchain_dp: invalid anchor count");
        return -1;
    }
    /* n == 0 stays a NULL `a`: mb_lchain_dp contracts on `n == 0 || a == 0` and returns an
     * empty chain set, so there is nothing to allocate. */
    mb_anchor_t *a = 0;
    if (n > 0) {
        a = (mb_anchor_t *)malloc((size_t)n * sizeof(mb_anchor_t));
        if (!a) { shim_set_err("mb_shim_lchain_dp: malloc failed"); return -1; }
        for (int64_t i = 0; i < n; ++i) {
            memset(&a[i], 0, sizeof(a[i]));
            a[i].sid = in_sid[i]; a[i].len = in_len[i];
            a[i].qpos = in_qpos[i]; a[i].tpos = in_tpos[i];
        }
    }
    int n_u = 0;
    uint64_t *u = 0;
    mb_anchor_t *b = mb_lchain_dp(0, idx->l2b, max_dist_x, max_dist_y, bw, max_skip, max_iter,
                                  min_sc, chn_pen_gap, n, a, &n_u, &u); /* frees a */
    int64_t n_a = 0;
    for (int i = 0; i < n_u; ++i) n_a += (uint32_t)u[i];
    for (int i = 0; i < n_u && i < max_u; ++i) out_u[i] = u[i];
    for (int64_t i = 0; i < n_a && i < max_a; ++i) {
        out_sid[i] = b[i].sid; out_len[i] = b[i].len;
        out_qpos[i] = b[i].qpos; out_tpos[i] = b[i].tpos;
    }
    *out_n_a = (int)n_a;
    free(u); free(b); /* both malloc'd via km=0; free(NULL) is safe when n_u==0 */
    return n_u;
}

/* End-to-end alignment-record oracle: run the full C mb_map pipeline per read and return, per
 * read, the PRIMARY hit (sam_pri) as an alignment record — the bit-exact target the Rust SW/align
 * stage is validated against. Reads are raw ASCII concat + off/len. Per read i:
 *   out_tid[i]  = target id, or -1 if unmapped
 *   out_pos[i]  = [ts, te]         (2 int64 per read)
 *   out_q[i]    = [qs, qe]         (2 int32 per read)
 *   out_misc[i] = [flags(rev|pri<<1|rev2..), mapq, dp_score, n_cigar]  (4 int32 per read)
 *   out_cigar[i*cig_stride ..]     = CIGAR ops (len<<4|op), up to cig_stride
 * Threaded via kt_for with per-thread mb_tbuf. Returns 0, or -1 if any read's CIGAR exceeds
 * cig_stride (grow it and retry). */
typedef struct {
    const mb_opt_t *opt;
    const mb_idx_t *idx;
    const char *seq;
    const uint64_t *off;
    const uint32_t *len;
    mb_tbuf_t **tbuf;
    int64_t *out_tid, *out_pos;
    int32_t *out_q, *out_misc;
    uint32_t *out_cigar;
    int cig_stride;
    int overflow;
} map_rec_t;

static void map_rec_worker(void *data, long i, int tid) {
    map_rec_t *d = (map_rec_t *)data;
    int32_t n_hit = 0;
    mb_hit_t *hit = mb_map(d->opt, d->idx, (int32_t)d->len[i], d->seq + d->off[i],
                           L2B_METH_NONE, &n_hit, d->tbuf[tid], 0);
    int pri = -1;
    for (int32_t j = 0; j < n_hit; ++j)
        if (hit[j].sam_pri) { pri = j; break; }
    if (pri < 0 && n_hit > 0) pri = 0; // fall back to the first hit
    if (pri < 0) {
        d->out_tid[i] = -1;
        d->out_misc[i * 4 + 3] = 0;
    } else {
        mb_hit_t *h = &hit[pri];
        d->out_tid[i] = h->tid;
        d->out_pos[i * 2] = h->ts; d->out_pos[i * 2 + 1] = h->te;
        d->out_q[i * 2] = h->qs; d->out_q[i * 2 + 1] = h->qe;
        int32_t nc = h->p ? h->p->n_cigar : 0;
        d->out_misc[i * 4] = (int32_t)(h->rev | (h->sam_pri << 1));
        d->out_misc[i * 4 + 1] = h->mapq;
        d->out_misc[i * 4 + 2] = h->p ? h->p->dp_score : 0;
        d->out_misc[i * 4 + 3] = nc;
        if (nc > d->cig_stride) { d->overflow = 1; nc = d->cig_stride; }
        for (int32_t k = 0; k < nc; ++k) d->out_cigar[(size_t)i * d->cig_stride + k] = h->p->cigar[k];
    }
    for (int32_t j = 0; j < n_hit; ++j) free(hit[j].p);
    free(hit);
}

int mb_shim_map(const mb_idx_t *idx, int n_reads,
                const char *seq, const uint64_t *off, const uint32_t *len, int n_threads,
                int64_t *out_tid, int64_t *out_pos, int32_t *out_q, int32_t *out_misc,
                uint32_t *out_cigar, int cig_stride) {
    if (n_threads < 1) n_threads = 1;
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = n_threads;
    map_rec_t d;
    memset(&d, 0, sizeof(d));
    d.opt = &opt; d.idx = idx; d.seq = seq; d.off = off; d.len = len;
    d.out_tid = out_tid; d.out_pos = out_pos; d.out_q = out_q; d.out_misc = out_misc;
    d.out_cigar = out_cigar; d.cig_stride = cig_stride;
    d.tbuf = (mb_tbuf_t **)calloc(n_threads, sizeof(mb_tbuf_t *));
    if (!d.tbuf) { shim_set_err("mb_shim_map: calloc failed"); return -1; }
    if (shim_tbuf_alloc(d.tbuf, n_threads) < 0) {
        free(d.tbuf);
        shim_set_err("mb_shim_map: tbuf init failed");
        return -1;
    }
    kt_for(n_threads, map_rec_worker, &d, (long)n_reads);
    for (int t = 0; t < n_threads; ++t) mb_tbuf_destroy(d.tbuf[t]);
    free(d.tbuf);
    return d.overflow ? -1 : 0;
}

/* Paired-end oracle: run C mb_map on 2*n_frag interleaved reads, mb_pestat over the batch, then
 * mb_pair per fragment (INCLUDING mate rescue). Returns, per mate m (2*n_frag total), the primary
 * (sam_pri) hit: out_tid[m] (-1 if none); out_pos[m*2]=ts,te; out_q[m*2]=qs,qe;
 * out_misc[m*4]={flags(rev|proper_pair<<1|sam_pri<<2), mapq, dp_score, n_cigar};
 * out_cigar[m*cig_stride]. Returns 0, or -1 on failure / cigar overflow. */
/* Per-region oracle (single read): return ALL regions from C mb_map. See header for layout. */
int mb_shim_map_regions(const mb_idx_t *idx, const char *seq, int len,
                        int64_t *out_tid, int64_t *out_ts, int64_t *out_te,
                        int32_t *out_qs, int32_t *out_qe, int32_t *out_dp, int32_t *out_dpmax,
                        int32_t *out_mapq, int32_t *out_flags, int32_t *out_parent,
                        int32_t *out_ncigar, uint32_t *out_cigar, int cig_stride, int max_regions) {
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    /* Optional: enable C's built-in AF (anchor) / AD (gap-fill block) trace to stderr for the
     * long-read align1 bisect. MB_DBG_AN_POS = 0x20. */
    extern int kom_dbg_flag;
    int saved_dbg = kom_dbg_flag;
    if (getenv("MBWA_DBG_ANPOS")) kom_dbg_flag |= 0x20;
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_map_regions: tbuf init failed"); return -1; }
    int32_t n_hit = 0;
    mb_hit_t *hit = mb_map(&opt, idx, len, seq, L2B_METH_NONE, &n_hit, b, 0);
    kom_dbg_flag = saved_dbg;
    int rc = n_hit;
    if (n_hit > max_regions) { shim_set_err("mb_shim_map_regions: region overflow"); rc = -1; }
    else {
        for (int32_t j = 0; j < n_hit; ++j) {
            mb_hit_t *h = &hit[j];
            out_tid[j] = h->tid;
            out_ts[j] = h->ts; out_te[j] = h->te;
            out_qs[j] = h->qs; out_qe[j] = h->qe;
            out_dp[j] = h->p ? h->p->dp_score : 0;
            out_dpmax[j] = h->p ? h->p->dp_max : 0;
            out_mapq[j] = h->mapq;
            out_flags[j] = (int32_t)(h->rev | (h->sam_pri << 1));
            out_parent[j] = h->parent;
            int32_t nc = h->p ? h->p->n_cigar : 0;
            out_ncigar[j] = nc;
            if (nc > cig_stride) { shim_set_err("mb_shim_map_regions: cigar overflow"); rc = -1; break; }
            for (int32_t k = 0; k < nc; ++k) out_cigar[(size_t)j * cig_stride + k] = h->p->cigar[k];
        }
    }
    for (int32_t j = 0; j < n_hit; ++j) free(hit[j].p);
    free(hit);
    mb_tbuf_destroy(b);
    return rc;
}

int mb_shim_pair(const mb_idx_t *idx, int n_frag,
                 const char *seq, const uint64_t *off, const uint32_t *len,
                 int64_t *out_tid, int64_t *out_pos, int32_t *out_q, int32_t *out_misc,
                 uint32_t *out_cigar, int cig_stride) {
    int n_reads = n_frag * 2;
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_pair: tbuf init failed"); return -1; }
    mb_hit_t **hit = (mb_hit_t **)calloc(n_reads, sizeof(mb_hit_t *));
    int32_t *n_hit = (int32_t *)calloc(n_reads, sizeof(int32_t));
    int32_t *seg_off = (int32_t *)calloc(n_frag, sizeof(int32_t));
    int32_t *seg_cnt = (int32_t *)calloc(n_frag, sizeof(int32_t));
    if (!hit || !n_hit || !seg_off || !seg_cnt) {
        free(hit); free(n_hit); free(seg_off); free(seg_cnt); mb_tbuf_destroy(b);
        shim_set_err("mb_shim_pair: calloc failed"); return -1;
    }
    for (int j = 0; j < n_reads; ++j)
        hit[j] = mb_map(&opt, idx, (int32_t)len[j], seq + off[j], L2B_METH_NONE, &n_hit[j], b, 0);
    for (int i = 0; i < n_frag; ++i) { seg_off[i] = i * 2; seg_cnt[i] = 2; }
    mb_pestat_t pes[4];
    mb_pestat(mb_tbuf_km(b), &opt, n_frag, seg_off, seg_cnt, n_hit, hit, pes);

    int ret = 0;
    for (int i = 0; i < n_frag; ++i) {
        int32_t nh[2] = { n_hit[i * 2], n_hit[i * 2 + 1] };
        mb_hit_t *hh[2] = { hit[i * 2], hit[i * 2 + 1] };
        int32_t ql[2] = { (int32_t)len[i * 2], (int32_t)len[i * 2 + 1] };
        char *qs[2] = { (char *)(seq + off[i * 2]), (char *)(seq + off[i * 2 + 1]) };
        mb_pair(mb_tbuf_km(b), &opt, idx->l2b, nh, hh, pes, ql, qs);
        hit[i * 2] = hh[0]; n_hit[i * 2] = nh[0];
        hit[i * 2 + 1] = hh[1]; n_hit[i * 2 + 1] = nh[1];
        for (int m = 0; m < 2; ++m) {
            int mm = i * 2 + m;
            int pri = -1;
            for (int32_t j = 0; j < nh[m]; ++j)
                if (hh[m][j].sam_pri) { pri = j; break; }
            if (pri < 0 && nh[m] > 0) pri = 0;
            if (pri < 0) { out_tid[mm] = -1; out_misc[mm * 4 + 3] = 0; continue; }
            mb_hit_t *h = &hh[m][pri];
            out_tid[mm] = h->tid;
            out_pos[mm * 2] = h->ts; out_pos[mm * 2 + 1] = h->te;
            out_q[mm * 2] = h->qs; out_q[mm * 2 + 1] = h->qe;
            int32_t nc = h->p ? h->p->n_cigar : 0;
            out_misc[mm * 4] = (int32_t)(h->rev | (h->proper_pair << 1) | (h->sam_pri << 2));
            out_misc[mm * 4 + 1] = h->mapq;
            out_misc[mm * 4 + 2] = h->p ? h->p->dp_score : 0;
            out_misc[mm * 4 + 3] = nc;
            if (nc > cig_stride) { ret = -1; nc = cig_stride; }
            for (int32_t k = 0; k < nc; ++k) out_cigar[(size_t)mm * cig_stride + k] = h->p->cigar[k];
        }
    }
    for (int j = 0; j < n_reads; ++j) {
        for (int32_t k = 0; k < n_hit[j]; ++k) free(hit[j][k].p);
        free(hit[j]);
    }
    free(hit); free(n_hit); free(seg_off); free(seg_cnt);
    mb_tbuf_destroy(b);
    return ret;
}

/* Insert-size estimation oracle: run C mb_map on 2*n_frag interleaved reads (R1,R2 per fragment),
 * then C mb_pestat over the batch. Reads are raw ASCII (concat + off/len), read 2i / 2i+1 are the
 * mates of fragment i. Per orientation d in 0..4: out_i[d*3 + {0,1,2}] = {failed, lo, hi};
 * out_d[d*2 + {0,1}] = {avg, std}. Returns 0, or -1 on allocation failure. */
int mb_shim_pestat(const mb_idx_t *idx, int n_frag,
                   const char *seq, const uint64_t *off, const uint32_t *len,
                   int32_t *out_i, double *out_d) {
    int n_reads = n_frag * 2;
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_pestat: tbuf init failed"); return -1; }
    mb_hit_t **hit = (mb_hit_t **)calloc(n_reads, sizeof(mb_hit_t *));
    int32_t *n_hit = (int32_t *)calloc(n_reads, sizeof(int32_t));
    int32_t *seg_off = (int32_t *)calloc(n_frag, sizeof(int32_t));
    int32_t *seg_cnt = (int32_t *)calloc(n_frag, sizeof(int32_t));
    if (!hit || !n_hit || !seg_off || !seg_cnt) {
        free(hit); free(n_hit); free(seg_off); free(seg_cnt); mb_tbuf_destroy(b);
        shim_set_err("mb_shim_pestat: calloc failed");
        return -1;
    }
    for (int j = 0; j < n_reads; ++j)
        hit[j] = mb_map(&opt, idx, (int32_t)len[j], seq + off[j], L2B_METH_NONE, &n_hit[j], b, 0);
    for (int i = 0; i < n_frag; ++i) { seg_off[i] = i * 2; seg_cnt[i] = 2; }

    mb_pestat_t pes[4];
    mb_pestat(mb_tbuf_km(b), &opt, n_frag, seg_off, seg_cnt, n_hit, hit, pes);
    for (int d = 0; d < 4; ++d) {
        out_i[d * 3 + 0] = pes[d].failed;
        out_i[d * 3 + 1] = pes[d].lo;
        out_i[d * 3 + 2] = pes[d].hi;
        out_d[d * 2 + 0] = pes[d].avg;
        out_d[d * 2 + 1] = pes[d].std;
    }
    for (int j = 0; j < n_reads; ++j) {
        for (int32_t k = 0; k < n_hit[j]; ++k) free(hit[j][k].p);
        free(hit[j]);
    }
    free(hit); free(n_hit); free(seg_off); free(seg_cnt);
    mb_tbuf_destroy(b);
    return 0;
}

/* SAM-header oracle: C mb_fmt_sam_hdr with no read-group and no command line, version = MB_VERSION.
 * Writes the header text into out_buf; returns its length, or -1 on overflow. */
int mb_shim_sam_hdr(const mb_idx_t *idx, char *out_buf, uint64_t out_cap) {
    kstring_t s = {0, 0, 0};
    mb_fmt_sam_hdr(&s, idx->l2b, 0, MB_VERSION, 0, 0);
    int ret = -1;
    if (s.l <= out_cap) { memcpy(out_buf, s.s, s.l); ret = (int)s.l; }
    free(s.s);
    return ret;
}

/* SAM-output oracle: run the full C mb_map pipeline per read and format its hits with the real
 * mb_fmt_sam (single-end, n_seg=1), replicating minibwa's main output loop (a line per primary
 * region, plus up to opt.out_n secondaries; an unmapped line when there are no hits). The
 * concatenated SAM text for read i lands at out_buf[out_off[i] .. out_off[i+1]); out_off has
 * n_reads+1 entries. `names[i]` is the QNAME; reads have no qual (QUAL='*'). Single-threaded,
 * in read order. Returns 0, or -1 on out_buf overflow (grow out_cap). */
int mb_shim_map_sam(const mb_idx_t *idx, int n_reads,
                    const char *const *names, const char *seq, const uint64_t *off, const uint32_t *len,
                    int paf, int64_t extra_flag,
                    const char *preset, const int64_t *ov, int n_ov, int has_pri, float pri,
                    uint64_t fset, uint64_t fclr,
                    char *out_buf, uint64_t out_cap, uint64_t *out_off) {
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    if (paf) opt.flag |= MB_F_PAF; /* 0x1: PAF instead of SAM */
    opt.flag |= (uint64_t)extra_flag; /* e.g. MB_F_EQX / MB_F_WRITE_CS/DS/MD */
    /* CLI overrides: must match the Rust mapping path (derive_read_opt) exactly for oracle parity. */
    shim_apply_ovr(&opt, preset, (ov && n_ov >= OV_N) ? ov : NULL, has_pri, pri, fset, fclr);
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_map_sam: tbuf init failed"); return -1; }
    kstring_t s = {0, 0, 0};
    uint64_t pos = 0;
    int ret = 0;
    for (int i = 0; i < n_reads; ++i) {
        out_off[i] = pos;
        int32_t n_hit = 0;
        /* The real binary seeds the tie-break hash with the read name (map-main.c:88 ->
           map-algo.c:573). Passing 0 here made the oracle disagree with `minibwa map` on ~5-7% of
           named reads -- equal-scoring hits, different winner. Pass the name so the oracle matches
           the binary; the Rust mirror passes align::hash_str(name) to match. */
        mb_hit_t *hit = mb_map(&opt, idx, (int32_t)len[i], seq + off[i], L2B_METH_NONE, &n_hit, b, names[i]);
        mb_bseq1_t t;
        memset(&t, 0, sizeof(t));
        t.l_seq = len[i];
        t.name = (char *)names[i];
        t.seq = (char *)(seq + off[i]);
        s.l = 0;
        mb_hit_t *hit_arr[1] = { hit };
        if (n_hit > 0) {
            int32_t n_sec = 0;
            for (int32_t j = 0; j < n_hit; ++j) {
                mb_hit_t *h = &hit[j];
                if (h->parent == h->id || n_sec < opt.out_n)
                    mb_format(mb_tbuf_km(b), &s, idx->l2b, &t, 1, &n_hit, hit_arr, j, &opt, 0, 0);
                n_sec += (h->parent != h->id);
            }
        } else {
            mb_format(mb_tbuf_km(b), &s, idx->l2b, &t, 1, &n_hit, hit_arr, -1, &opt, 0, 0);
        }
        if (pos + s.l > out_cap) { ret = -1; }
        else { memcpy(out_buf + pos, s.s, s.l); pos += s.l; }
        for (int32_t j = 0; j < n_hit; ++j) free(hit[j].p);
        free(hit);
        if (ret) break;
    }
    out_off[n_reads] = pos;
    free(s.s);
    mb_tbuf_destroy(b);
    return ret;
}

/* ksw_ll oracle: C ksw_ll_qinit + ksw_ll_u8_core (size=1) or ksw_ll_i16_core (size=2) on an nt4
 * query/target pair with a 5x5 matrix from (a,b,b_ts,b_ambi). out[5] = {score, te, qe, score2,
 * te2}. Bit-exact reference for the Rust ksw_ll port. */
void mb_shim_ksw_ll(int size, const uint8_t *query, int qlen, const uint8_t *target, int tlen,
                    int8_t a, int8_t b, int8_t b_ts, int8_t b_ambi,
                    int gapo, int gape, int xtra, int32_t out[5]) {
    int8_t mat[25];
    ksw_gen_nt4_mat(mat, a, b, b_ts, b_ambi, 0);
    void *qp = ksw_ll_qinit(0, size, qlen, query, 5, mat); /* km=0: libc alloc (b[] grows via realloc) */
    ksw_llrst_t r = size > 1 ? ksw_ll_i16_core(qp, tlen, target, gapo, gape, xtra)
                             : ksw_ll_u8_core(qp, tlen, target, gapo, gape, xtra);
    out[0] = r.score; out[1] = r.te; out[2] = r.qe; out[3] = r.score2; out[4] = r.te2;
    free(qp);
}

/* Paired-end SAM oracle: run mb_map on 2*n_frag interleaved reads, mb_pestat, mb_pair per
 * fragment (incl. rescue), then format BOTH mates' records with mb_fmt_sam (n_seg=2), matching
 * minibwa's PE output loop. Per fragment f, both mates' SAM text lands at
 * out_buf[out_off[f] .. out_off[f+1]); out_off has n_frag+1 entries. qname=0 for the seed hash
 * (matches Rust). Returns 0, or -1 on failure/overflow. */
int mb_shim_pair_sam(const mb_idx_t *idx, int n_frag,
                     const char *const *names, const char *seq, const uint64_t *off, const uint32_t *len,
                     const char *preset, const int64_t *ov, int n_ov, int has_pri, float pri,
                     uint64_t fset, uint64_t fclr,
                     char *out_buf, uint64_t out_cap, uint64_t *out_off) {
    int n_reads = n_frag * 2;
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    /* CLI overrides: must match the Rust mapping path (derive_read_opt) exactly for oracle parity. */
    shim_apply_ovr(&opt, preset, (ov && n_ov >= OV_N) ? ov : NULL, has_pri, pri, fset, fclr);
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_pair_sam: tbuf init failed"); return -1; }
    mb_hit_t **hit = (mb_hit_t **)calloc(n_reads, sizeof(mb_hit_t *));
    int32_t *n_hit = (int32_t *)calloc(n_reads, sizeof(int32_t));
    int32_t *seg_off = (int32_t *)calloc(n_frag, sizeof(int32_t));
    int32_t *seg_cnt = (int32_t *)calloc(n_frag, sizeof(int32_t));
    if (!hit || !n_hit || !seg_off || !seg_cnt) {
        free(hit); free(n_hit); free(seg_off); free(seg_cnt); mb_tbuf_destroy(b);
        shim_set_err("mb_shim_pair_sam: calloc failed"); return -1;
    }
    for (int j = 0; j < n_reads; ++j)
        hit[j] = mb_map(&opt, idx, (int32_t)len[j], seq + off[j], L2B_METH_NONE, &n_hit[j], b, names[j]);
    for (int i = 0; i < n_frag; ++i) { seg_off[i] = i * 2; seg_cnt[i] = 2; }
    /* -P / --hic (MB_F_NO_PAIRING): map-main.c gates the whole PE block (pestat + mb_pair) on
     * !(flag & MB_F_NO_PAIRING); with it set, each mate keeps its SE hits and only the n_seg=2
     * formatting differs. Match that here so the oracle stays faithful to native `map -P`. */
    int do_pairing = !(opt.flag & MB_F_NO_PAIRING);
    mb_pestat_t pes[4];
    if (do_pairing)
        mb_pestat(mb_tbuf_km(b), &opt, n_frag, seg_off, seg_cnt, n_hit, hit, pes);

    kstring_t s = {0, 0, 0};
    uint64_t pos = 0;
    int ret = 0;
    for (int f = 0; f < n_frag; ++f) {
        out_off[f] = pos;
        int32_t nh[2] = { n_hit[f * 2], n_hit[f * 2 + 1] };
        mb_hit_t *hh[2] = { hit[f * 2], hit[f * 2 + 1] };
        int32_t ql[2] = { (int32_t)len[f * 2], (int32_t)len[f * 2 + 1] };
        char *qs[2] = { (char *)(seq + off[f * 2]), (char *)(seq + off[f * 2 + 1]) };
        if (do_pairing)
            mb_pair(mb_tbuf_km(b), &opt, idx->l2b, nh, hh, pes, ql, qs);
        hit[f * 2] = hh[0]; n_hit[f * 2] = nh[0];
        hit[f * 2 + 1] = hh[1]; n_hit[f * 2 + 1] = nh[1];
        s.l = 0;
        for (int seg = 0; seg < 2; ++seg) {
            int32_t mate_qlen = ql[1 - seg];
            mb_bseq1_t t;
            memset(&t, 0, sizeof(t));
            t.l_seq = ql[seg];
            t.name = (char *)names[f * 2 + seg];
            t.seq = qs[seg];
            if (nh[seg] > 0) {
                int32_t n_sec = 0;
                for (int32_t j = 0; j < nh[seg]; ++j) {
                    if (hh[seg][j].parent == hh[seg][j].id || n_sec < opt.out_n)
                        mb_format(mb_tbuf_km(b), &s, idx->l2b, &t, 2, nh, hh, j, &opt, seg, mate_qlen);
                    n_sec += (hh[seg][j].parent != hh[seg][j].id);
                }
            } else { /* unmapped mate: still emitted unless MB_F_NO_UNMAP */
                mb_format(mb_tbuf_km(b), &s, idx->l2b, &t, 2, nh, hh, -1, &opt, seg, mate_qlen);
            }
        }
        if (pos + s.l > out_cap) { ret = -1; }
        else { memcpy(out_buf + pos, s.s, s.l); pos += s.l; }
        if (ret) break;
    }
    out_off[n_frag] = pos;
    free(s.s);
    for (int j = 0; j < n_reads; ++j) {
        for (int32_t k = 0; k < n_hit[j]; ++k) free(hit[j][k].p);
        free(hit[j]);
    }
    free(hit); free(n_hit); free(seg_off); free(seg_cnt);
    mb_tbuf_destroy(b);
    return ret;
}

/* Methylation paired-end oracle: like mb_shim_pair but on a meth index. Each read j is
 * nt4-encoded, converted (even=C2T/R1, odd=G2A/R2), seeded converted, then mb_map_sai(mt) on the
 * original read; then mb_pestat + mb_pair (MB_F_METH set, so mate rescue uses the meth path) per
 * fragment. Returns per mate m: out_tid/pos/q/misc/cigar (as mb_shim_pair). Returns 0 or -1. */
int mb_shim_pair_meth(const mb_idx_t *idx, int n_frag,
                      const char *seq, const uint64_t *off, const uint32_t *len,
                      int64_t *out_tid, int64_t *out_pos, int32_t *out_q, int32_t *out_misc,
                      uint32_t *out_cigar, int cig_stride) {
    int n_reads = n_frag * 2;
    mb_opt_t opt, opt_adap;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    opt.flag |= MB_F_METH;
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_pair_meth: tbuf init failed"); return -1; }
    void *km = mb_tbuf_km(b);
    mb_hit_t **hit = (mb_hit_t **)calloc(n_reads, sizeof(mb_hit_t *));
    int32_t *n_hit = (int32_t *)calloc(n_reads, sizeof(int32_t));
    int32_t *seg_off = (int32_t *)calloc(n_frag, sizeof(int32_t));
    int32_t *seg_cnt = (int32_t *)calloc(n_frag, sizeof(int32_t));
    if (!hit || !n_hit || !seg_off || !seg_cnt) {
        free(hit); free(n_hit); free(seg_off); free(seg_cnt); mb_tbuf_destroy(b);
        shim_set_err("mb_shim_pair_meth: calloc failed"); return -1;
    }
    for (int j = 0; j < n_reads; ++j) {
        l2b_meth_t m = (j & 1) == 0 ? L2B_METH_C2T : L2B_METH_G2A;
        int l = (int)len[j];
        uint8_t *seq4 = (uint8_t *)kmalloc(km, l);
        for (int i = 0; i < l; ++i) seq4[i] = kom_nt4_table[(uint8_t)(seq + off[j])[i]];
        l2b_meth_convert(m, l, seq4);
        mb_sai_v u = {0, 0, 0};
        mb_seed_intv(km, idx->bwt, l, seq4, opt.min_len, opt.max_sub_occ, &u);
        kfree(km, seq4);
        mb_opt_adap(&opt, l, &opt_adap);
        hit[j] = mb_map_sai(&opt_adap, idx, l, seq + off[j], m, &u, &n_hit[j], b, 0);
    }
    for (int i = 0; i < n_frag; ++i) { seg_off[i] = i * 2; seg_cnt[i] = 2; }
    mb_pestat_t pes[4];
    mb_pestat(km, &opt, n_frag, seg_off, seg_cnt, n_hit, hit, pes);

    int ret = 0;
    for (int i = 0; i < n_frag; ++i) {
        int32_t nh[2] = { n_hit[i * 2], n_hit[i * 2 + 1] };
        mb_hit_t *hh[2] = { hit[i * 2], hit[i * 2 + 1] };
        int32_t ql[2] = { (int32_t)len[i * 2], (int32_t)len[i * 2 + 1] };
        char *qs[2] = { (char *)(seq + off[i * 2]), (char *)(seq + off[i * 2 + 1]) };
        mb_pair(km, &opt, idx->l2b, nh, hh, pes, ql, qs);
        hit[i * 2] = hh[0]; n_hit[i * 2] = nh[0];
        hit[i * 2 + 1] = hh[1]; n_hit[i * 2 + 1] = nh[1];
        for (int mm2 = 0; mm2 < 2; ++mm2) {
            int mm = i * 2 + mm2, pri = -1;
            for (int32_t j = 0; j < nh[mm2]; ++j)
                if (hh[mm2][j].sam_pri) { pri = j; break; }
            if (pri < 0 && nh[mm2] > 0) pri = 0;
            if (pri < 0) { out_tid[mm] = -1; out_misc[mm * 4 + 3] = 0; continue; }
            mb_hit_t *h = &hh[mm2][pri];
            out_tid[mm] = h->tid;
            out_pos[mm * 2] = h->ts; out_pos[mm * 2 + 1] = h->te;
            out_q[mm * 2] = h->qs; out_q[mm * 2 + 1] = h->qe;
            int32_t nc = h->p ? h->p->n_cigar : 0;
            out_misc[mm * 4] = (int32_t)(h->rev | (h->proper_pair << 1) | (h->sam_pri << 2));
            out_misc[mm * 4 + 1] = h->mapq;
            out_misc[mm * 4 + 2] = h->p ? h->p->dp_score : 0;
            out_misc[mm * 4 + 3] = nc;
            if (nc > cig_stride) { ret = -1; nc = cig_stride; }
            for (int32_t k = 0; k < nc; ++k) out_cigar[(size_t)mm * cig_stride + k] = h->p->cigar[k];
        }
    }
    for (int j = 0; j < n_reads; ++j) {
        for (int32_t k = 0; k < n_hit[j]; ++k) free(hit[j][k].p);
        free(hit[j]);
    }
    free(hit); free(n_hit); free(seg_off); free(seg_cnt);
    mb_tbuf_destroy(b);
    return ret;
}

/* Methylation mapping oracle: replicate the meth path for ONE read on a meth index. nt4-encode the
 * read, l2b_meth_convert it (mt=1 C2T / mt=2 G2A), seed the CONVERTED read, then mb_map_sai on the
 * ORIGINAL read + those seeds + mt. Returns the PRIMARY hit: out_tid (-1 if none); out_pos[0..1]=
 * ts,te; out_q[0..1]=qs,qe; out_misc[0..3]={flags(rev|pri<<1),mapq,dp_score,n_cigar};
 * out_cigar[0..n_cigar). Returns 0, or -1 on cigar overflow (>cig_stride). */
int mb_shim_map_meth(const mb_idx_t *idx, const char *read, int len, int mt, int cig_stride,
                     int64_t *out_tid, int64_t *out_pos, int32_t *out_q, int32_t *out_misc,
                     uint32_t *out_cigar) {
    mb_opt_t opt, opt_adap;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    if (mt) opt.flag |= MB_F_METH;
    mb_opt_adap(&opt, len, &opt_adap);
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_map_meth: tbuf init failed"); return -1; }
    void *km = mb_tbuf_km(b);
    uint8_t *seq4 = (uint8_t *)kmalloc(km, len);
    for (int i = 0; i < len; ++i) seq4[i] = kom_nt4_table[(uint8_t)read[i]];
    l2b_meth_t m = mt == 1 ? L2B_METH_C2T : mt == 2 ? L2B_METH_G2A : L2B_METH_NONE;
    if (m != L2B_METH_NONE) l2b_meth_convert(m, len, seq4);
    mb_sai_v u = {0, 0, 0};
    mb_seed_intv(km, idx->bwt, len, seq4, opt.min_len, opt.max_sub_occ, &u);
    kfree(km, seq4);
    int32_t n_hit = 0;
    mb_hit_t *hit = mb_map_sai(&opt_adap, idx, len, read, m, &u, &n_hit, b, 0);
    int ret = 0, pri = -1;
    for (int32_t j = 0; j < n_hit; ++j)
        if (hit[j].sam_pri) { pri = j; break; }
    if (pri < 0 && n_hit > 0) pri = 0;
    if (pri < 0) {
        out_tid[0] = -1; out_misc[3] = 0;
    } else {
        mb_hit_t *h = &hit[pri];
        out_tid[0] = h->tid;
        out_pos[0] = h->ts; out_pos[1] = h->te;
        out_q[0] = h->qs; out_q[1] = h->qe;
        int32_t nc = h->p ? h->p->n_cigar : 0;
        out_misc[0] = (int32_t)(h->rev | (h->sam_pri << 1));
        out_misc[1] = h->mapq;
        out_misc[2] = h->p ? h->p->dp_score : 0;
        out_misc[3] = nc;
        if (nc > cig_stride) { ret = -1; nc = cig_stride; }
        for (int32_t k = 0; k < nc; ++k) out_cigar[k] = h->p->cigar[k];
    }
    for (int32_t j = 0; j < n_hit; ++j) free(hit[j].p);
    free(hit);
    mb_tbuf_destroy(b);
    return ret;
}

/* SE methylation SAM oracle: like mb_shim_map_sam but on a meth index. Each read is nt4-encoded,
 * C2T-converted (single-end meth uses mt=1, matching map-main.c's (j&1)==0 -> C2T), seeded on the
 * converted read, mb_map_sai'd on the ORIGINAL read + mt, then formatted with mb_fmt_sam (n_seg=1).
 * MB_F_METH is forced on. Same override/flag handling and out_buf/out_off layout as mb_shim_map_sam.
 * qname=0 for the seed hash (matches Rust read_hash(0,...)). Returns 0, or -1 on overflow. */
int mb_shim_map_meth_sam(const mb_idx_t *idx, int n_reads,
                         const char *const *names, const char *seq, const uint64_t *off, const uint32_t *len,
                         const char *preset, const int64_t *ov, int n_ov, int has_pri, float pri,
                         uint64_t fset, uint64_t fclr,
                         char *out_buf, uint64_t out_cap, uint64_t *out_off) {
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    opt.flag |= MB_F_METH;
    shim_apply_ovr(&opt, preset, (ov && n_ov >= OV_N) ? ov : NULL, has_pri, pri, fset, fclr);
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_map_meth_sam: tbuf init failed"); return -1; }
    void *km = mb_tbuf_km(b);
    kstring_t s = {0, 0, 0};
    uint64_t pos = 0;
    int ret = 0;
    for (int i = 0; i < n_reads; ++i) {
        out_off[i] = pos;
        int l = (int)len[i];
        uint8_t *seq4 = (uint8_t *)kmalloc(km, l);
        for (int c = 0; c < l; ++c) seq4[c] = kom_nt4_table[(uint8_t)(seq + off[i])[c]];
        l2b_meth_convert(L2B_METH_C2T, l, seq4);
        mb_sai_v u = {0, 0, 0};
        mb_seed_intv(km, idx->bwt, l, seq4, opt.min_len, opt.max_sub_occ, &u);
        kfree(km, seq4);
        mb_opt_t opt_adap;
        mb_opt_adap(&opt, l, &opt_adap);
        int32_t n_hit = 0;
        mb_hit_t *hit = mb_map_sai(&opt_adap, idx, l, seq + off[i], L2B_METH_C2T, &u, &n_hit, b, names[i]);
        mb_bseq1_t t;
        memset(&t, 0, sizeof(t));
        t.l_seq = l;
        t.name = (char *)names[i];
        t.seq = (char *)(seq + off[i]);
        s.l = 0;
        mb_hit_t *hit_arr[1] = { hit };
        if (n_hit > 0) {
            int32_t n_sec = 0;
            for (int32_t j = 0; j < n_hit; ++j) {
                mb_hit_t *h = &hit[j];
                if (h->parent == h->id || n_sec < opt.out_n)
                    mb_format(km, &s, idx->l2b, &t, 1, &n_hit, hit_arr, j, &opt, 0, 0);
                n_sec += (h->parent != h->id);
            }
        } else {
            mb_format(km, &s, idx->l2b, &t, 1, &n_hit, hit_arr, -1, &opt, 0, 0);
        }
        if (pos + s.l > out_cap) { ret = -1; }
        else { memcpy(out_buf + pos, s.s, s.l); pos += s.l; }
        for (int32_t j = 0; j < n_hit; ++j) free(hit[j].p);
        free(hit);
        if (ret) break;
    }
    out_off[n_reads] = pos;
    free(s.s);
    mb_tbuf_destroy(b);
    return ret;
}

/* PE methylation SAM oracle: like mb_shim_pair_sam but on a meth index. Each mate is nt4-encoded and
 * converted by segment (R1/seg0 -> C2T mt=1, R2/seg1 -> G2A mt=2, matching map-main.c's (j&1)),
 * seeded converted, mb_map_sai'd on the original read + mt; then mb_pestat + mb_pair per fragment
 * with MB_F_METH set (mate rescue uses the meth path), and both mates formatted with mb_fmt_sam
 * (n_seg=2). Honors MB_F_NO_PAIRING (skip pairing) like mb_shim_pair_sam. Returns 0, or -1. */
int mb_shim_pair_meth_sam(const mb_idx_t *idx, int n_frag,
                          const char *const *names, const char *seq, const uint64_t *off, const uint32_t *len,
                          const char *preset, const int64_t *ov, int n_ov, int has_pri, float pri,
                          uint64_t fset, uint64_t fclr,
                          char *out_buf, uint64_t out_cap, uint64_t *out_off) {
    int n_reads = n_frag * 2;
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    opt.flag |= MB_F_METH;
    shim_apply_ovr(&opt, preset, (ov && n_ov >= OV_N) ? ov : NULL, has_pri, pri, fset, fclr);
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_pair_meth_sam: tbuf init failed"); return -1; }
    void *km = mb_tbuf_km(b);
    mb_hit_t **hit = (mb_hit_t **)calloc(n_reads, sizeof(mb_hit_t *));
    int32_t *n_hit = (int32_t *)calloc(n_reads, sizeof(int32_t));
    int32_t *seg_off = (int32_t *)calloc(n_frag, sizeof(int32_t));
    int32_t *seg_cnt = (int32_t *)calloc(n_frag, sizeof(int32_t));
    if (!hit || !n_hit || !seg_off || !seg_cnt) {
        free(hit); free(n_hit); free(seg_off); free(seg_cnt); mb_tbuf_destroy(b);
        shim_set_err("mb_shim_pair_meth_sam: calloc failed"); return -1;
    }
    for (int j = 0; j < n_reads; ++j) {
        l2b_meth_t m = (j & 1) == 0 ? L2B_METH_C2T : L2B_METH_G2A;
        int l = (int)len[j];
        uint8_t *seq4 = (uint8_t *)kmalloc(km, l);
        for (int i = 0; i < l; ++i) seq4[i] = kom_nt4_table[(uint8_t)(seq + off[j])[i]];
        l2b_meth_convert(m, l, seq4);
        mb_sai_v u = {0, 0, 0};
        mb_seed_intv(km, idx->bwt, l, seq4, opt.min_len, opt.max_sub_occ, &u);
        kfree(km, seq4);
        mb_opt_t opt_adap;
        mb_opt_adap(&opt, l, &opt_adap);
        hit[j] = mb_map_sai(&opt_adap, idx, l, seq + off[j], m, &u, &n_hit[j], b, names[j]);
    }
    for (int i = 0; i < n_frag; ++i) { seg_off[i] = i * 2; seg_cnt[i] = 2; }
    int do_pairing = !(opt.flag & MB_F_NO_PAIRING);
    mb_pestat_t pes[4];
    if (do_pairing)
        mb_pestat(km, &opt, n_frag, seg_off, seg_cnt, n_hit, hit, pes);

    kstring_t s = {0, 0, 0};
    uint64_t pos = 0;
    int ret = 0;
    for (int f = 0; f < n_frag; ++f) {
        out_off[f] = pos;
        int32_t nh[2] = { n_hit[f * 2], n_hit[f * 2 + 1] };
        mb_hit_t *hh[2] = { hit[f * 2], hit[f * 2 + 1] };
        int32_t ql[2] = { (int32_t)len[f * 2], (int32_t)len[f * 2 + 1] };
        char *qs[2] = { (char *)(seq + off[f * 2]), (char *)(seq + off[f * 2 + 1]) };
        if (do_pairing)
            mb_pair(km, &opt, idx->l2b, nh, hh, pes, ql, qs);
        hit[f * 2] = hh[0]; n_hit[f * 2] = nh[0];
        hit[f * 2 + 1] = hh[1]; n_hit[f * 2 + 1] = nh[1];
        s.l = 0;
        for (int seg = 0; seg < 2; ++seg) {
            int32_t mate_qlen = ql[1 - seg];
            mb_bseq1_t t;
            memset(&t, 0, sizeof(t));
            t.l_seq = ql[seg];
            t.name = (char *)names[f * 2 + seg];
            t.seq = qs[seg];
            if (nh[seg] > 0) {
                int32_t n_sec = 0;
                for (int32_t j = 0; j < nh[seg]; ++j) {
                    if (hh[seg][j].parent == hh[seg][j].id || n_sec < opt.out_n)
                        mb_format(km, &s, idx->l2b, &t, 2, nh, hh, j, &opt, seg, mate_qlen);
                    n_sec += (hh[seg][j].parent != hh[seg][j].id);
                }
            } else {
                mb_format(km, &s, idx->l2b, &t, 2, nh, hh, -1, &opt, seg, mate_qlen);
            }
        }
        if (pos + s.l > out_cap) { ret = -1; }
        else { memcpy(out_buf + pos, s.s, s.l); pos += s.l; }
        if (ret) break;
    }
    out_off[n_frag] = pos;
    free(s.s);
    for (int j = 0; j < n_reads; ++j) {
        for (int32_t k = 0; k < n_hit[j]; ++k) free(hit[j][k].p);
        free(hit[j]);
    }
    free(hit); free(n_hit); free(seg_off); free(seg_cnt);
    mb_tbuf_destroy(b);
    return ret;
}

/* Debug oracle: run C mb_map on ONE read and dump ALL hits (not just the primary). Per hit k:
 * out[k*8 + {0..7}] = tid, ts, te, qs, qe, flags(rev|inv<<1|sam_pri<<2|(parent==id)<<3), mapq,
 * dp_score. `out` is filled for the first max_hits hits; the return is the TRUE hit count, which
 * may exceed max_hits (grow `out` and retry). Returns -1 on failure. */
int mb_shim_map_hits(const mb_idx_t *idx, const char *seq, int len, int max_hits, int64_t *out) {
    mb_opt_t opt;
    mb_opt_init(&opt);
    opt.n_thread = 1;
    mb_tbuf_t *b = mb_tbuf_init(0);
    if (!b) { shim_set_err("mb_shim_map_hits: tbuf init failed"); return -1; }
    int32_t n_hit = 0;
    mb_hit_t *hit = mb_map(&opt, idx, len, seq, L2B_METH_NONE, &n_hit, b, 0);
    int n = n_hit < max_hits ? n_hit : max_hits;
    for (int k = 0; k < n; ++k) {
        mb_hit_t *h = &hit[k];
        out[k * 8 + 0] = h->tid;
        out[k * 8 + 1] = h->ts;
        out[k * 8 + 2] = h->te;
        out[k * 8 + 3] = h->qs;
        out[k * 8 + 4] = h->qe;
        out[k * 8 + 5] = h->rev | (h->inv << 1) | (h->sam_pri << 2) | ((h->id == h->parent) << 3);
        out[k * 8 + 6] = h->mapq;
        out[k * 8 + 7] = h->p ? h->p->dp_score : 0;
    }
    for (int32_t k = 0; k < n_hit; ++k) free(hit[k].p);
    free(hit);
    mb_tbuf_destroy(b);
    return n_hit;
}

/* Length-adapted option oracle: run mb_opt_init + mb_opt_adap(qlen) and dump the exact options
 * the C pipeline uses for a read of length `qlen`. The `mb_map` pipeline adapts bw/zdrop/max_gap/
 * min_chain_score/min_dp_max/best_n to the read length (see mb_opt_adap), so the Rust seed→chain→
 * align port must use these exact values to stay bit-exact. Fixed field order:
 *   oi[0..8)  = a, b, b_ts, b_ambi, q, e, q2, e2
 *   oi[8..18) = bw, bw_long, max_gap, min_len, min_ksw_len, min_chain_score, min_dp_max, end_bonus,
 *               zdrop, zdrop_inv
 *   oi[18..25)= max_occ, max_sub_occ, max_chain_skip, max_chain_iter, mask_len, best_n, is_sr
 *   oi[25..29)= pen_unpair, max_rescue, max_sr_len, out_n     (PE-pairing + output knobs)
 *   of[0..3)  = chain_gap_scale, pri_ratio, mask_level         (floats)
 *   ol[0..2)  = max_sw_mat, (int64)flag
 * `oi` must hold >=29 int32, `of` >=3 float, `ol` >=2 int64. */
/* Apply CLI overrides onto a freshly-mb_opt_init'd base opt: -x preset first, then scalar overrides
 * (sentinel-guarded), the -p float, then the flag set/clear masks — mirroring main_map's order.
 * (Declared near the top of the file; layout enum/macro live there too.) */
static void shim_apply_ovr(mb_opt_t *opt, const char *preset, const int64_t *ov,
                           int has_pri, float pri, uint64_t fset, uint64_t fclr) {
    if (preset && *preset) mb_opt_preset(opt, preset);
    if (ov) {
        if (ov[OV_A]          != MB_SHIM_OVR_KEEP) opt->a               = (int32_t)ov[OV_A];
        if (ov[OV_B]          != MB_SHIM_OVR_KEEP) opt->b               = (int32_t)ov[OV_B];
        if (ov[OV_MIN_LEN]    != MB_SHIM_OVR_KEEP) opt->min_len         = (int32_t)ov[OV_MIN_LEN];
        if (ov[OV_MAX_OCC]    != MB_SHIM_OVR_KEEP) opt->max_occ         = (int32_t)ov[OV_MAX_OCC];
        if (ov[OV_MIN_CHAIN]  != MB_SHIM_OVR_KEEP) opt->min_chain_score = (int32_t)ov[OV_MIN_CHAIN];
        if (ov[OV_BEST_N]     != MB_SHIM_OVR_KEEP) opt->best_n          = (int32_t)ov[OV_BEST_N];
        if (ov[OV_PEN_UNPAIR] != MB_SHIM_OVR_KEEP) opt->pen_unpair      = (int32_t)ov[OV_PEN_UNPAIR];
        if (ov[OV_MAX_SR_LEN] != MB_SHIM_OVR_KEEP) opt->max_sr_len      = (int32_t)ov[OV_MAX_SR_LEN];
        if (ov[OV_MAX_GAP]    != MB_SHIM_OVR_KEEP) opt->max_gap         = (int32_t)ov[OV_MAX_GAP];
        if (ov[OV_BW]         != MB_SHIM_OVR_KEEP) opt->bw              = (int32_t)ov[OV_BW];
        if (ov[OV_BW_LONG]    != MB_SHIM_OVR_KEEP) opt->bw_long         = (int32_t)ov[OV_BW_LONG];
        if (ov[OV_MIN_DP_MAX] != MB_SHIM_OVR_KEEP) opt->min_dp_max      = (int32_t)ov[OV_MIN_DP_MAX];
        if (ov[OV_OUT_N]      != MB_SHIM_OVR_KEEP) opt->out_n           = (int32_t)ov[OV_OUT_N];
        if (ov[OV_MAX_RESCUE] != MB_SHIM_OVR_KEEP) opt->max_rescue      = (int32_t)ov[OV_MAX_RESCUE];
        if (ov[OV_XA_MAX]     != MB_SHIM_OVR_KEEP) opt->xa_max          = (int32_t)ov[OV_XA_MAX];
    }
    if (has_pri) opt->pri_ratio = pri;
    opt->flag |= fset;
    opt->flag &= ~fclr;
}

/* Override-aware variant of mb_shim_align_opt: build the base opt via mb_opt_init + CLI overrides,
 * then mb_opt_adap for qlen, and export the same oi/of/ol arrays. `n_ov` guards a short/absent
 * override array (must be >= OV_N to be applied). */
void mb_shim_align_opt_ovr(int qlen, const char *preset, const int64_t *ov, int n_ov,
                           int has_pri, float pri, uint64_t fset, uint64_t fclr,
                           int32_t *oi, float *of, int64_t *ol) {
    mb_opt_t opt0, opt;
    mb_opt_init(&opt0);
    shim_apply_ovr(&opt0, preset, (ov && n_ov >= OV_N) ? ov : NULL, has_pri, pri, fset, fclr);
    mb_opt_adap(&opt0, qlen, &opt);
    oi[0] = opt.a;  oi[1] = opt.b;  oi[2] = opt.b_ts;  oi[3] = opt.b_ambi;
    oi[4] = opt.q;  oi[5] = opt.e;  oi[6] = opt.q2;    oi[7] = opt.e2;
    oi[8] = opt.bw; oi[9] = opt.bw_long; oi[10] = opt.max_gap; oi[11] = opt.min_len;
    oi[12] = opt.min_ksw_len; oi[13] = opt.min_chain_score; oi[14] = opt.min_dp_max;
    oi[15] = opt.end_bonus; oi[16] = opt.zdrop; oi[17] = opt.zdrop_inv;
    oi[18] = opt.max_occ; oi[19] = opt.max_sub_occ; oi[20] = opt.max_chain_skip;
    oi[21] = opt.max_chain_iter; oi[22] = opt.mask_len; oi[23] = opt.best_n;
    oi[24] = mb_is_sr_mode(&opt, qlen);
    oi[25] = opt.pen_unpair; oi[26] = opt.max_rescue; oi[27] = opt.max_sr_len; oi[28] = opt.out_n;
    of[0] = opt.chain_gap_scale; of[1] = opt.pri_ratio; of[2] = opt.mask_level;
    ol[0] = opt.max_sw_mat; ol[1] = (int64_t)opt.flag;
}

void mb_shim_align_opt(int qlen, int32_t *oi, float *of, int64_t *ol) {
    mb_shim_align_opt_ovr(qlen, NULL, NULL, 0, 0, 0.0f, 0, 0, oi, of, ol);
}

/* Reference-sequence extraction oracle: C l2b_getseq over the loaded index's l2bit store. Writes
 * nt4 (with N=4) for contig tid, contig-local [st,en), into out; returns bases written. Bit-exact
 * reference for the Rust L2b::getseq port. */
int64_t mb_shim_getseq(const mb_idx_t *idx, int64_t tid, int64_t st, int64_t en, uint8_t *out) {
    return l2b_getseq(idx->l2b, tid, st, en, out);
}

/* Chaining gap-penalty oracle: exactly comput_sc's penalty (lchain.c:148-150). The single
 * FMA-sensitive site; bit-exact reference for mbwa_core::chain::gap_penalty and the device gate. */
int32_t mb_shim_gap_penalty(float chn_pen_gap, int64_t dd) {
    if (dd == 0) return 0;
    float lin_pen = chn_pen_gap * (float)dd;
    float log_pen = dd >= 1 ? mb_log2((float)(dd + 1)) : 0.0f;
    return (int32_t)(lin_pen + 0.5f * log_pen);
}

/* Index-buffer accessors: expose the C-owned FM-index buffers so a caller can copy them
 * elsewhere without reloading or duplicating the index. */
const uint64_t *mb_shim_bwt_data(const mb_idx_t *idx) { return idx->bwt->data; }
uint64_t mb_shim_bwt_data_len(const mb_idx_t *idx) { return idx->bwt->data_len; }
uint64_t mb_shim_bwt_primary(const mb_idx_t *idx) { return idx->bwt->primary; }
uint64_t mb_shim_bwt_seq_len(const mb_idx_t *idx) { return idx->bwt->seq_len; }
const uint64_t *mb_shim_bwt_sa(const mb_idx_t *idx) { return idx->bwt->sa; }
uint64_t mb_shim_bwt_n_sa(const mb_idx_t *idx) { return idx->bwt->n_sa; }
uint32_t mb_shim_bwt_sa_bit(const mb_idx_t *idx) { return idx->bwt->sa_bit; }
const uint32_t *mb_shim_bwt_cnt_table(const mb_idx_t *idx) { return idx->bwt->cnt_table; }
void mb_shim_bwt_l2(const mb_idx_t *idx, uint64_t out[5]) {
    for (int i = 0; i < 5; ++i) out[i] = idx->bwt->L2[i];
}
