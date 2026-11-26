"""Figure fig:adv: the approach flow ahead of the Phillips inlet against the ADV
measurements of Goharrokhi et al. (2019).

(a) CFD vertical profile 5 mm ahead of the tip with the digitised ADV profiles (their Fig. 3);
(b) CFD axial velocity along the inlet axis at the four sweep velocities, with each bore mean RU.
Also prints the axis-cell comparison against U used in the text and Table tab:labR.
Run from cfd/: python geometry/plot_adv.py
Needs the v4 sweep slices in results/mesh_slices/ and results/goharrokhi2019_fig3_digitised.csv.

The paper figure is drawn in MATLAB (paper/figures/matlab/make_fig_adv.m) from the CSV files
this script writes to paper/figures/data/:
- fig_adv_profile.csv (z_mm, u_cfd, u_scaled);
- fig_adv_meas.csv (flume_cm_s, u_over_U, z_mm);
- fig_adv_axis.csv (U, s_mm, u_over_U).
 (fig_adv_bore.csv, used by no figure)
The matplotlib version is only a preview, written to paper/figures/preview/ with --preview.
"""
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "geometry")
from goharrokhi_adv import TIP, cell_mean, interpolator, late_R, plane  # noqa: E402

INK, MUTED = "#2c3e50", "#7f8c8d"
BLUES = {0.154: "#74a9cf", 0.3: "#3690c0", 0.585: "#0570b0", 1.0: "#034e7b"}
C_CFD = "#2a78d6"
SWEEP = [0.154, 0.3, 0.585, 1.0]
MATCH = {15: 0.154, 30: 0.3, 60: 0.585}       # ADV flume velocity (cm/s) -> CFD case
DIG = Path("results/goharrokhi2019_fig3_digitised.csv")
OUT = Path("../paper/figures/preview/fig_adv")
DATA = Path("../paper/figures/data")


def slice_path(U):
    return f"results/mesh_slices/v4_s1.00_U{U:.3f}_field.npz"


def ts_path(U):
    return f"results/urans_v4/phillips2000_s1.00_sst-gamma/U{U:.3f}_timeseries.csv"


def main():
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6,
                         "font.family": "DejaVu Sans"})
    dig = [dict(s=int(r["series_cm_s"]), u=float(r["u_cm_s"]), z=float(r["z_mm_from_axis"]))
           for r in csv.DictReader(open(DIG))]
    axis_z = {s: [d["z"] for d in dig if d["s"] == s and -5 <= d["z"] <= 2.5] for s in MATCH}
    zs_axis = sorted(set(round(z, 1) for zz in axis_z.values() for z in zz))  # their axis cells

    fields = {U: interpolator(slice_path(U), U)[0] for U in SWEEP}
    # the v4 0.3 m/s time series is not in the results; its late mean is used instead
    R = {U: late_R(ts_path(U)) if Path(ts_path(U)).exists() else {0.3: 0.1857}[U] for U in SWEEP}

    fig, (prof, axl) = plt.subplots(1, 2, figsize=(6.2, 2.8),
                                    gridspec_kw=dict(width_ratios=[1, 1.35]))

    # (a) vertical profile 5 mm ahead; the CFD shape is nearly the same at every U
    ratios = []
    for cm, U in MATCH.items():
        pts = [d for d in dig if d["s"] == cm]
        cz = plane(fields[U], TIP - 0.005, np.array([d["z"] for d in pts]) * 1e-3)[:, 0]
        ratios += list(np.array([d["u"] / cm for d in pts]) / cz)
    k = float(np.mean(ratios))
    z = np.linspace(-12, 12, 241)
    u3 = plane(fields[0.3], TIP - 0.005, z * 1e-3)[:, 0]
    prof.plot(u3, z, color=C_CFD, lw=1.5, label="CFD")
    prof.plot(k * u3, z, color=C_CFD, lw=1.2, ls=(0, (4, 2)), label=f"CFD × {k:.2f}")
    for cm, mk in ((15, "D"), (30, "^"), (60, "o")):
        pts = [d for d in dig if d["s"] == cm]
        prof.plot([d["u"] / cm for d in pts], [d["z"] for d in pts], mk, ms=3.4, mfc="white",
                  mec=INK, mew=0.9, ls="none", label=f"ADV {cm}")
    prof.axhline(0, color=MUTED, lw=0.5, ls=(0, (2, 2)))
    prof.set_xlim(0.3, 1.08); prof.set_ylim(-12, 12)
    prof.set_xlabel("streamwise velocity / U")
    prof.set_ylabel("height above the inlet axis (mm)")
    prof.legend(loc="lower left", fontsize=6.0, frameon=False, handlelength=1.8,
                borderaxespad=0.2)
    prof.set_title("(a) profile 5 mm ahead of the tip", loc="left", fontsize=8.5, color=INK)

    # (b) centreline approach
    s = np.linspace(-60, 20, 801)
    axl.axvspan(0, 20, color="#eef1f3", lw=0)
    axl.text(10, 1.0, "inside the\ninlet tube", ha="center", va="top", fontsize=7, color=MUTED)
    axl.axvline(-5, color=MUTED, lw=0.6, ls=(0, (2, 2)))
    axl.text(-6.5, 0.5, "ADV cells\n5 mm ahead", ha="right", va="center", fontsize=7,
             color=MUTED)
    for U in SWEEP:
        u = plane(fields[U], TIP + s * 1e-3, np.zeros_like(s))[:, 0]
        axl.plot(s[s <= 0], u[s <= 0], color=BLUES[U], lw=1.5)
        axl.plot([0.5, 19.5], [R[U]] * 2, color=BLUES[U], lw=1.5, ls=(0, (3, 1.5)))
        axl.text(20.8, R[U], f"{U:.3g}" if U < 1 else "1.0", va="center", fontsize=6.5,
                 color=BLUES[U])
    axl.text(20.8, max(R.values()) + 0.07, "U (m/s);\nbore RU", va="bottom", fontsize=6.5,
             color=MUTED)
    axl.set_xlim(-60, 30); axl.set_ylim(0, 1.02)
    axl.set_xlabel("distance from the inlet tip along the axis (mm)")
    axl.set_ylabel("streamwise velocity / U")
    axl.set_title("(b) approach along the inlet axis", loc="left", fontsize=8.5, color=INK)

    for ax in (prof, axl):
        ax.grid(color="#e5e8ea", lw=0.5)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout(w_pad=1.5)
    if "--preview" in sys.argv:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(OUT.with_suffix(".pdf"))
        fig.savefig(OUT.with_suffix(".png"), dpi=300)

    DATA.mkdir(parents=True, exist_ok=True)
    with open(DATA / "fig_adv_profile.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(("z_mm", "u_cfd", "u_scaled"))
        w.writerows((f"{a:.4g}", f"{b:.5f}", f"{k * b:.5f}") for a, b in zip(z, u3))
    with open(DATA / "fig_adv_meas.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(("flume_cm_s", "u_over_U", "z_mm"))
        w.writerows((d["s"], f"{d['u'] / d['s']:.5f}", f"{d['z']:.3f}") for d in dig
                    if d["s"] in MATCH)
    with open(DATA / "fig_adv_axis.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(("U", "s_mm", "u_over_U"))
        for U in SWEEP:
            u = plane(fields[U], TIP + s * 1e-3, np.zeros_like(s))[:, 0]
            w.writerows((f"{U:g}", f"{a:.4g}", f"{b:.5f}") for a, b in zip(s, u) if a <= 0)

    # axis-cell comparison against U (text and Table tab:labR)
    def stack(U, dx):
        return np.mean([cell_mean(fields[U], TIP + dx, zz * 1e-3, 0) for zz in zs_axis])

    meas = {s: np.mean([d["u"] for d in dig if d["s"] == s and -5 <= d["z"] <= 2.5]) / s
            for s in (15, 21, 30, 35, 60)}
    print("their axis cells (mm):", zs_axis)
    print("U, CFD axis-cell mean (5 mm), band 4-8 mm, bore R")
    for U in SWEEP:
        v4, v5, v8 = stack(U, -0.004), stack(U, -0.005), stack(U, -0.008)
        print(f"  {U:.3f}  {v5:.3f}  {min(v4, v8):.3f}-{max(v4, v8):.3f}  {R[U]:.4f}")
    print(f"profile scale ADV/CFD k = {k:.3f} (sd {np.std(ratios):.3f}, n {len(ratios)})")
    print("ADV measured axis-stack mean / U:", {c: round(v, 3) for c, v in meas.items()})


if __name__ == "__main__":
    main()
