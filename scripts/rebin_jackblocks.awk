#!/usr/bin/awk -f
# Rebin a *_ylm_2pt_full_jackblocks_*.dat file from many small jackknife
# blocks into `n_target` coarser blocks, by summing sum_re/sum_im across
# every `group` consecutive raw blocks. Assumes n_k=1 per raw block (true
# for jack_block_size=1 runs), so the rebinned block's sample count is
# just the count of raw block indices folded into it.
#
# Usage: awk -v group=500 -f rebin_jackblocks.awk in.dat > out_body.dat
NR == 1 { next }  # header handled separately by the caller
{
    k = $1; l = $2; mm = $3; mp = $4; sre = $5; sim = $6
    b = int(k / group)
    key = b SUBSEP l SUBSEP mm SUBSEP mp
    sumre[key] += sre
    sumim[key] += sim
    bk = b SUBSEP k
    if (!(bk in seenk)) { seenk[bk] = 1; bucketcount[b]++ }
}
END {
    for (key in sumre) {
        split(key, a, SUBSEP)
        b = a[1]
        printf "%s %s %s %s %.16e %.16e %d\n", b, a[2], a[3], a[4], sumre[key], sumim[key], bucketcount[b]
    }
}
