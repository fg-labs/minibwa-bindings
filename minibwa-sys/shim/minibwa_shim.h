#ifndef MINIBWA_SHIM_H
#define MINIBWA_SHIM_H

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

#ifdef __cplusplus
}
#endif

#endif
