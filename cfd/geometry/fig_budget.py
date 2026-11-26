"""Figure fig:budget: where the free-stream dynamic pressure goes on its way through
each sampler at 0.3 m/s.

From the symmetry-plane field slices (results/mesh_slices/*_U0.300_field.npz, written on the HPC by
fluent/field_slice.py from each production run's final case), the total pressure is averaged
across the bore at stations along the flow path: the inlet mouth (A, 2 mm inside the tip), the
inner end of the inlet bore (B, 2 mm short of it), the chamber (C, mid-length), 2 mm into the
outlet bore (D) and the outlet mouth (E, 2 mm short of the tip). On the plane y = 0 an
axisymmetric cross-section average weights each face by its distance from the bore axis, and
total pressure is weighted by mass flux as well. The free-stream dynamic pressure is split into:

  approach      1 - Cpt(A)          lost before the mouth
  inlet bore    Cpt(A) - Cpt(B)     friction (and the insert contraction for Perks)
  expansion     Cpt(B) - Cp(C)      jet into the chamber
  outlet entry  Cp(C) - Cpt(D)      entry into the outlet bore
  outlet bore   Cpt(D) - Cpt(E)     friction, fittings and the elbow (Elliott)
  exit energy   Cpt(E) - Cp(E)      kinetic energy leaving the outlet
  base          Cp(E)               outlet static pressure; negative (suction) adds to the drive

so the segments sum to 1. Cp = (p - p_inf) / (rho U^2 / 2), with p_inf the inlet-boundary static
pressure of the run. The meter variant has no station B or D (the meter fills the path), so its
inlet bore and expansion are one segment, as are its outlet entry and bore.

Run from cfd/: python geometry/fig_budget.py   (writes results/budget_U0.300.csv and its copy
../paper/figures/data/fig_budget.csv)

The paper figure is drawn in MATLAB (paper/figures/matlab/make_fig_budget.m) from
fig_budget.csv. The matplotlib version is only a preview, written to paper/figures/preview/
with --preview.
"""
import csv
import json
import shutil
import statistics as st
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from tims_geometry import SPECS

RHO, U = 998.2, 0.3
# Elliott's outlet mouth sits in the vertical 4 mm tube below the elbow, where the y = 0 slice
# does not resolve the bore cleanly (static pressure jumps from ~9 to ~2 Pa across the station).
# There the static pressure is the Fluent surface report and the kinetic energy leaving is
# ALPHA_EXIT * rho * ubar^2 / 2, ubar from the reported mass flow. ALPHA_EXIT = 2.0 is what the
# slices give at the outlet mouths of the axial designs (1.96-2.17); 1.3-2 would move about 0.03
# of the budget between "outlet bore" and "kinetic energy leaving".
ALPHA_EXIT = 2.0
INK, MUTED, GRID = "#2c3e50", "#7f8c8d", "#e5e8ea"
SLICES = Path("results/mesh_slices")
DESIGNS = [  # label, run, slice file, results folder, geometry version
    ("Phillips", "phillips2000", "phillips2000_v5_s1.00", "urans_v5", "v5"),
    ("Foets", "foets2020", "foets2020_v5bore_s1.00", "foets", "v5bore"),
    ("Perks", "perks2014", "perks2014_v5_s1.00", "urans_v5", "v5"),
    ("Perks + meter", "perks2014_flowmeter", "perks2014_flowmeter_v5_s1.00", "urans_v5", "v5"),
    ("Elliott", "elliott2017", "elliott2017_v5_s1.00", "urans_v5", "v5"),
    ("McDonald", "mcdonald2010", "mcdonald2010_v5_s1.00", "urans_v5bore", "v5bore"),
]
SEGMENTS = [  # key, label, colour (one hue per stage of the path, light to dark along it)
    ("approach", "approach to the mouth", "#c9d6e3"),
    ("inlet", "inlet bore", "#2a78d6"),
    ("expansion", "expansion into the chamber", "#8fb8ea"),
    ("out_entry", "outlet entry", "#eb6834"),
    ("out_bore", "outlet bore and fittings", "#a8441d"),
    ("exit", "kinetic energy leaving", "#f3b38f"),
    ("base", "outlet static pressure", "#9aa5ad"),
]


def load_slice(stem):
    z = np.load(SLICES / f"{stem}_U0.300_field.npz")
    off, xz = z["off"], z["xz"]
    n = np.diff(off)
    return dict(x=np.add.reduceat(xz[:, 0], off[:-1]) / n, z=np.add.reduceat(xz[:, 1], off[:-1]) / n,
                a=z["area"], u=z["u"], w=z["w"], p=z["p"],
                q=0.5 * RHO * (z["u"] ** 2 + z["v"] ** 2 + z["w"] ** 2))


def section(m, x0, r, zc=0.0, h=6e-4, vertical=False):
    """Mass-weighted total and area-weighted static pressure across a bore of radius r, on the
    plane y = 0: faces within h of the station, weighted by distance from the bore axis."""
    if vertical:                              # bore along z at x = x0, station at height zc
        s = (np.abs(m["z"] - zc) < h) & (np.abs(m["x"] - x0) < r)
        rad, flux = np.abs(m["x"][s] - x0), np.abs(m["w"][s])
    else:
        s = (np.abs(m["x"] - x0) < h) & (np.abs(m["z"] - zc) < r)
        rad, flux = np.abs(m["z"][s] - zc), np.abs(m["u"][s])
    if s.sum() < 10:
        raise ValueError(f"only {s.sum()} faces at x={x0}, z={zc}")
    wa = rad * m["a"][s]
    pt = m["p"][s] + m["q"][s]
    return float((pt * flux * wa).sum() / (flux * wa).sum()), float((m["p"][s] * wa).sum() / wa.sum())


def stations(run, ver):
    zones = json.load(open(f"geometry/out/{ver}/scale_1.00/{run}_zones.json"))
    x = zones["stations_m"]
    s = SPECS[run]
    r_in = (s.nozzle_id or s.tube_id) / 2               # bore at the inner end of the inlet
    r_tip = s.tube_id / 2
    out_d = getattr(s, "outlet_id", 0) or s.tube_id
    st_ = dict(A=dict(x0=x["mouth_up"], r=r_tip), B=dict(x0=x["in_inner"] - 0.002, r=r_in),
               C=dict(x0=0.5 * (x["cav0"] + x["cav1"]), r=s.body_id / 2),
               D=dict(x0=x["out_inner"] + 0.002, r=(s.nozzle_id or out_d) / 2),
               E=dict(x0=x["mouth_down"], r=out_d / 2))
    if run == "elliott2017":                            # through the 6.35 mm fitting, elbow, 4 mm tube down
        st_["D"]["r"] = s.pass_d / 2
        st_["E"] = dict(report=True, r=s.tube_id / 2)     # see ALPHA_EXIT
    if run == "mcdonald2010":                           # 3 mm hole through the plug, off axis
        st_["D"] = dict(x0=x["out_inner"] + 0.002, r=s.hole_d / 2, zc=s.hole_z)
        st_["E"] = dict(x0=x["mouth_down"], r=s.hole_d / 2, zc=s.hole_z)
    if run == "perks2014_flowmeter":
        st_.pop("B"), st_.pop("D")
    return st_


def budget(label, run, stem, folder, ver):
    m = load_slice(stem)
    base = Path("results") / folder / f"{run}_s1.00_sst-gamma" / "U0.300"
    summ = json.load(open(f"{base}_summary.json"))
    rows = [r for r in csv.DictReader(open(f"{base}_timeseries.csv"))
            if float(r["t"]) >= summ["average_from"] - 1e-9]
    p_inf = st.mean(float(r["p_inlet"]) for r in rows)
    p_up_report = st.mean(float(r["p_mouth_up"]) for r in rows)
    qinf = 0.5 * RHO * U**2
    sts = stations(run, ver)
    v = {k: section(m, **kw) for k, kw in sts.items() if not kw.get("report")}
    if sts["E"].get("report"):
        pE = st.mean(float(r["p_mouth_down"]) for r in rows)
        ubar = abs(summ["mean"]["mdot_mouth_down"]) / RHO / (np.pi * sts["E"]["r"] ** 2 / 2)
        v["E"] = (pE + ALPHA_EXIT * 0.5 * RHO * ubar**2, pE)
    cp = lambda p: (p - p_inf) / qinf                   # noqa: E731
    ptA, pA = v["A"]
    ptE, pE = v["E"]
    pC = v["C"][1]
    b = dict(approach=1 - cp(ptA), exit=cp(ptE) - cp(pE), base=cp(pE))
    if "B" in v:
        b.update(inlet=cp(ptA) - cp(v["B"][0]), expansion=cp(v["B"][0]) - cp(pC),
                 out_entry=cp(pC) - cp(v["D"][0]), out_bore=cp(v["D"][0]) - cp(ptE))
    else:
        b.update(inlet=cp(ptA) - cp(pC), expansion=0.0, out_entry=cp(pC) - cp(ptE), out_bore=0.0)
    b["R"] = summ["mean"]["R"]
    b["check_p_mouth_up"] = (pA - p_up_report) / qinf   # slice against the Fluent report
    return b


def main():
    rows = [(label, budget(label, run, stem, folder, ver)) for label, run, stem, folder, ver in DESIGNS]
    keys = [k for k, _, _ in SEGMENTS]
    with open("results/budget_U0.300.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design", "R"] + keys + ["sum", "check_p_mouth_up"])
        for label, b in rows:
            w.writerow([label, f"{b['R']:.4f}"] + [f"{b[k]:.4f}" for k in keys]
                       + [f"{sum(b[k] for k in keys):.4f}", f"{b['check_p_mouth_up']:+.4f}"])
            print(f"{label:14s} R {b['R']:.4f} " + " ".join(f"{k} {b[k]:+.3f}" for k in keys)
                  + f"  sum {sum(b[k] for k in keys):.3f}  slice-report {b['check_p_mouth_up']:+.4f}")

    data = Path("../paper/figures/data")
    data.mkdir(parents=True, exist_ok=True)
    shutil.copyfile("results/budget_U0.300.csv", data / "fig_budget.csv")
    if "--preview" not in sys.argv:
        return

    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6,
                         "font.family": "DejaVu Sans"})
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    order = rows[::-1]                                   # Figure 3 order, top to bottom
    for i, (label, b) in enumerate(order):
        left = 0.0
        merged = b["expansion"] == 0 and b["out_bore"] == 0    # meter variant
        for k, _, colour in SEGMENTS:
            val = b[k]
            if val <= 0:
                continue
            ax.barh(i, val, left=left, height=0.62, color=colour, edgecolor="white", lw=0.8)
            if merged and k in ("inlet", "out_entry"):
                ax.text(left + val / 2, i, "inlet bore, meter, expansion" if k == "inlet"
                        else "outlet entry and bore", ha="center", va="center", fontsize=6,
                        color="white")
            left += val
        if b["base"] < 0:                               # outlet suction: drive beyond rho U^2/2
            ax.barh(i, -b["base"], left=1.0, height=0.62, color="none", edgecolor=INK,
                    hatch="////", lw=0)
        ax.text(-0.02, i, f"{label}\n$R$ = {b['R']:.3f}", ha="right", va="center", fontsize=7,
                color=INK)
    ax.axvline(1.0, color=INK, lw=0.7, ls=(0, (2, 2)))
    ax.set_yticks([])
    ax.set_xlim(0, 1.13)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel(r"share of the free-stream dynamic pressure $\rho U^2/2$ ($U$ = 0.3 m s$^{-1}$)")
    ax.grid(axis="x", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, ec="white") for _, _, c in SEGMENTS]
    labels = [lab for _, lab, _ in SEGMENTS]
    handles.append(plt.Rectangle((0, 0), 1, 1, fc="white", ec=INK, hatch="////", lw=0))
    labels.append("drive from outlet suction")
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=6.8,
               bbox_to_anchor=(0.5, 0.0), handlelength=1.6, columnspacing=1.2)
    fig.tight_layout(rect=(0.1, 0.17, 1, 1))
    out = Path("../paper/figures/preview/fig_budget")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=300)


if __name__ == "__main__":
    main()
