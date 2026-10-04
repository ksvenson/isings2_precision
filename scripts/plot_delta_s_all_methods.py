#!/usr/bin/env python3
"""Combined Delta_sigma(a) vs 1/n_refine comparison plot across every
method this campaign has tried: push7's unwindowed l-space Delta_s(a)
(plateaus ~0.33-0.41, the original anomaly), the real-space z-binned
method's jackknifed production ladder (job 7311120, extrapolates to
~0.125), and the real-space orbit-averaged ladder (job 7311202, plotted
with whatever shards are available -- safe to rerun once it completes for
a full comparison). Hardcodes this session's specific `campaign_runs/`
directory names; not a general-purpose reusable tool, more a snapshot of
the state of the Delta-extraction investigation as of 2026-08-25 -- see
journal.md for the narrative these numbers support.

Usage: scripts/plot_delta_s_all_methods.py
"""
import sys, os, glob, re
sys.path.insert(0, os.path.dirname(__file__))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares
from fit_real_space_ladder import analyze_point
from cft_symmetry_test import analyze as cft_analyze

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# --- 1. push7 l-space unwindowed Delta_s(a) (naive mesh) ---
push7_dir = os.path.join(ROOT, "campaign_runs/push7_2026-08-25_dual_naive_eqarea")
ladder7 = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]
k7, ds7, ds7_err = [], [], []
for k in ladder7:
    paths = sorted(glob.glob(os.path.join(push7_dir, f"q5k{k}", "shard_*", f"q5k{k}", "*_jackblocks_*.dat")))
    if not paths:
        continue
    res = cft_analyze(paths)
    k7.append(k)
    ds7.append(res["delta_s"])
    ds7_err.append(res["delta_s_err"])
k7 = np.array(k7, dtype=float)
ds7 = np.array(ds7)
ds7_err = np.array(ds7_err)

# --- 2. real-space z-binned method, full jackknifed ladder ---
rs_dir = os.path.join(ROOT, "campaign_runs/real_space_ladder_2026-08-25")
theta_min, theta_max = 0.3, 2.8
point_dirs = sorted(glob.glob(os.path.join(rs_dir, "q5k*")), key=lambda p: int(re.search(r"q5k(\d+)", p).group(1)))
k_rs, d_rs, e_rs = [], [], []
for pd in point_dirs:
    k = int(re.search(r"q5k(\d+)", pd).group(1))
    n_shards, d_full, d_bar, d_err = analyze_point(pd, theta_min, theta_max)
    k_rs.append(k); d_rs.append(d_bar); e_rs.append(d_err)
k_rs = np.array(k_rs, dtype=float)
d_rs = np.array(d_rs)
e_rs = np.array(e_rs)

# fit real-space to Delta_inf + c/n_refine (validated best power=1 earlier)
def resid(params):
    d_inf, c = params
    return (d_inf + c / k_rs - d_rs) / e_rs
sol = least_squares(resid, (0.125, 0.0))
chi2 = np.sum(sol.fun**2)
dof = len(k_rs) - 2
cov = np.linalg.inv(sol.jac.T @ sol.jac) * (chi2 / dof)
d_inf, d_inf_err = sol.x[0], np.sqrt(cov[0, 0])

# --- 3. orbit-averaged ladder, if available (partial ladder OK) ---
orbit_dir = os.path.join(ROOT, "campaign_runs/real_space_ladder_orbit_2026-08-25")
k_orb, d_orb, e_orb = [], [], []
for pd in sorted(glob.glob(os.path.join(orbit_dir, "q5k*")), key=lambda p: int(re.search(r"q5k(\d+)", p).group(1))):
    k = int(re.search(r"q5k(\d+)", pd).group(1))
    n_shards_avail = len(glob.glob(os.path.join(pd, "shard_*", "*.dat")))
    if n_shards_avail < 4:
        continue
    try:
        n_shards, d_full, d_bar, d_err = analyze_point(pd, theta_min, theta_max)
    except Exception:
        continue
    k_orb.append(k); d_orb.append(d_bar); e_orb.append(d_err)
k_orb = np.array(k_orb, dtype=float)
d_orb = np.array(d_orb)
e_orb = np.array(e_orb)

# --- plot ---
fig, ax = plt.subplots(figsize=(9, 6.5))

ax.errorbar(1.0 / k7, ds7, yerr=ds7_err, fmt="^", color="#C44E52", ms=6, capsize=3,
            label="l-space Delta_s(a)=2F1/(F1+F0), unwindowed, full sphere (push7)")

ax.errorbar(1.0 / k_rs, d_rs, yerr=e_rs, fmt="o", color="#4C72B0", ms=6, capsize=3,
            label="real-space fit Delta(a), z-binned (job 7311120, jackknifed)")
xx = np.linspace(0, (1.0 / k_rs).max() * 1.1, 100)
ax.plot(xx, sol.x[0] + sol.x[1] * xx, "-", color="#4C72B0", lw=1.5, alpha=0.7,
        label=f"real-space extrap: Delta_inf={d_inf:.4f}({int(round(d_inf_err*1e4))})  chi2/dof={chi2/dof:.2f}")

if len(k_orb) > 0:
    ax.errorbar(1.0 / k_orb, d_orb, yerr=e_orb, fmt="s", color="#55A868", ms=6, capsize=3,
                label="real-space fit Delta(a), orbit-averaged (job 7311202, partial ladder, jackknifed)")

ax.axhline(0.125, color="gray", ls=":", lw=1.5, label="exact CFT value: 1/8")
ax.set_xlabel("1/n_refine")
ax.set_ylabel("Delta_sigma(a)")
ax.set_title("Delta_sigma(a) vs 1/n_refine: l-space (unwindowed) vs. real-space methods")
ax.legend(fontsize=8, loc="center right")
ax.grid(alpha=0.3)
fig.tight_layout()
out = os.path.join(ROOT, "campaign_runs", "delta_s_all_methods_2026-08-25.png")
fig.savefig(out, dpi=150)
print("wrote", out)
print(f"push7 l-space points: {len(k7)}, real-space z-binned points: {len(k_rs)}, orbit-avg points ready: {len(k_orb)}")
