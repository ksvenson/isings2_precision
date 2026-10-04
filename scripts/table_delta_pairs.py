#!/usr/bin/env python3
"""Writes a plain-text/CSV table of Delta_{l,pair}(a) predictions (the
single-recursion-step scaling-dimension estimate, exact value 1/8) per
n_refine, per l pair, for naive and equal_area -- the tabulated companion
to plot_delta_pairs_mesh_compare.py's plots. Reuses that script's
point_recursion() (same compute_F_l/delta_l_pair_formula reduction) so the
numbers are guaranteed identical to what the plots show.

Usage:
  scripts/table_delta_pairs.py <naive_dir> <eqarea_dir> [--l_max 8] [-o out.csv]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_delta_pairs_mesh_compare import collect, point_recursion  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("naive_dir")
    ap.add_argument("eqarea_dir")
    ap.add_argument("--l_max", type=int, default=8)
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    naive_by_n = collect(args.naive_dir)
    eqarea_by_n = collect(args.eqarea_dir)
    ns = sorted(set(naive_by_n) & set(eqarea_by_n))

    out = args.out or os.path.join(
        os.path.dirname(args.naive_dir.rstrip("/")), "analysis_full", "delta_l_pair_table.csv"
    )
    os.makedirs(os.path.dirname(out), exist_ok=True)

    rows = []
    for n in ns:
        naive_res = point_recursion(sorted(naive_by_n[n]), args.l_max)
        eqarea_res = point_recursion(sorted(eqarea_by_n[n]), args.l_max)
        for l in range(1, args.l_max + 1):
            rows.append(dict(
                n_refine=n, l_minus_1=l - 1, l=l,
                naive_delta=naive_res["delta_l_pair"][l],
                naive_err=naive_res["delta_l_pair_err"][l],
                eqarea_delta=eqarea_res["delta_l_pair"][l],
                eqarea_err=eqarea_res["delta_l_pair_err"][l],
            ))

    header = ("n_refine,l_minus_1,l,naive_delta_l_pair,naive_err,"
              "eqarea_delta_l_pair,eqarea_err,naive_sigma_from_1_8,eqarea_sigma_from_1_8")
    with open(out, "w") as f:
        f.write(header + "\n")
        for r in rows:
            naive_sigma = (r["naive_delta"] - 0.125) / r["naive_err"] if r["naive_err"] > 0 else float("nan")
            eqarea_sigma = (r["eqarea_delta"] - 0.125) / r["eqarea_err"] if r["eqarea_err"] > 0 else float("nan")
            f.write(
                f"{r['n_refine']},{r['l_minus_1']},{r['l']},"
                f"{r['naive_delta']:.6f},{r['naive_err']:.6f},"
                f"{r['eqarea_delta']:.6f},{r['eqarea_err']:.6f},"
                f"{naive_sigma:.2f},{eqarea_sigma:.2f}\n"
            )
    print(f"wrote {out} ({len(rows)} rows)")

    # also print a human-readable summary to stdout: largest n_refine only
    n_max = max(ns)
    print(f"\nDelta_l_pair(a) at n_refine={n_max} (exact value: 0.125000):")
    print(f"{'l-1->l':>8} {'naive':>10} {'+/-':>8} {'sigma':>6}   {'eqarea':>10} {'+/-':>8} {'sigma':>6}")
    for r in rows:
        if r["n_refine"] != n_max:
            continue
        naive_sigma = (r["naive_delta"] - 0.125) / r["naive_err"]
        eqarea_sigma = (r["eqarea_delta"] - 0.125) / r["eqarea_err"]
        print(
            f"{r['l_minus_1']:>3}->{r['l']:<3} "
            f"{r['naive_delta']:>10.6f} {r['naive_err']:>8.6f} {naive_sigma:>6.2f}   "
            f"{r['eqarea_delta']:>10.6f} {r['eqarea_err']:>8.6f} {eqarea_sigma:>6.2f}"
        )


if __name__ == "__main__":
    main()
