#include "minibwa_shim.h"
#include <stdio.h>
#include <string.h>

/* Provided by minibwa's index.c (not declared in minibwa.h). */
extern int main_index(int argc, char **argv);

static __thread char g_err[512];

static void shim_set_err(const char *msg) {
    snprintf(g_err, sizeof(g_err), "%s", msg ? msg : "");
}

const char *mb_shim_last_error(void) { return g_err; }

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
