/* Rebin a *_ylm_2pt_full_jackblocks_*.dat file (jack_block_size=1, fixed
 * per-block pair order matching full_pairs_for_l() in symmetry_test.py)
 * into n_blocks/group coarser jackknife blocks, summing sum_re/sum_im.
 * Relies on the C++ writer emitting exactly n_pairs lines per block in a
 * fixed, identical order every block (verified against the file header).
 *
 * Usage: rebin_jackblocks <in.dat> <group> <out.dat>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s <in.dat> <group> <out.dat>\n", argv[0]);
        return 1;
    }
    FILE *fin = fopen(argv[1], "r");
    if (!fin) { perror("fopen in"); return 1; }
    long group = atol(argv[2]);
    FILE *fout = fopen(argv[3], "w");
    if (!fout) { perror("fopen out"); return 1; }

    char line[512];
    if (!fgets(line, sizeof(line), fin)) { fprintf(stderr, "empty file\n"); return 1; }
    long n_pairs = 0, jbs = 0, n_blocks = 0;
    if (sscanf(line, "# n_pairs=%ld jack_block_size=%ld n_blocks=%ld", &n_pairs, &jbs, &n_blocks) != 3) {
        fprintf(stderr, "bad header: %s\n", line);
        return 1;
    }

    long n_target = (n_blocks + group - 1) / group;
    double *sumre = calloc((size_t)n_target * n_pairs, sizeof(double));
    double *sumim = calloc((size_t)n_target * n_pairs, sizeof(double));
    long *bucketcount = calloc((size_t)n_target, sizeof(long));
    char **labels = calloc((size_t)n_pairs, sizeof(char *));

    long count = 0;
    long k;
    char l_s[8], mm_s[8], mp_s[8];
    double sre, sim;

    while (fgets(line, sizeof(line), fin)) {
        if (sscanf(line, "%ld %7s %7s %7s %lf %lf", &k, l_s, mm_s, mp_s, &sre, &sim) != 6) {
            fprintf(stderr, "bad line at count=%ld: %s\n", count, line);
            return 1;
        }
        long pair_idx = count % n_pairs;
        long block_idx = count / n_pairs;
        if (block_idx != k) {
            fprintf(stderr, "order mismatch at count=%ld: expected block %ld, got k=%ld\n",
                    count, block_idx, k);
            return 1;
        }
        if (labels[pair_idx] == NULL) {
            char buf[32];
            snprintf(buf, sizeof(buf), "%s %s %s", l_s, mm_s, mp_s);
            labels[pair_idx] = strdup(buf);
        }
        long bucket = block_idx / group;
        sumre[bucket * n_pairs + pair_idx] += sre;
        sumim[bucket * n_pairs + pair_idx] += sim;
        if (pair_idx == 0) bucketcount[bucket]++;
        count++;
    }
    fclose(fin);

    if (count != n_blocks * n_pairs) {
        fprintf(stderr, "warning: read %ld lines, expected %ld\n", count, n_blocks * n_pairs);
    }

    fprintf(fout, "# n_pairs=%ld jack_block_size=%ld n_blocks=%ld\n", n_pairs, jbs * group, n_target);
    for (long b = 0; b < n_target; b++) {
        for (long p = 0; p < n_pairs; p++) {
            fprintf(fout, "%ld %s %.16e %.16e %ld\n",
                    b, labels[p], sumre[b * n_pairs + p], sumim[b * n_pairs + p], bucketcount[b]);
        }
    }
    fclose(fout);
    fprintf(stderr, "wrote %ld blocks x %ld pairs to %s\n", n_target, n_pairs, argv[3]);
    return 0;
}
