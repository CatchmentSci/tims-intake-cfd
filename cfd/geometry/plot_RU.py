"""Figure fig:RU: intake ratio against free-stream velocity for the six designs.

(a) R(U); (b) the local power-law slope d ln R / d ln U between adjacent velocities, plotted at
their geometric mean. Linear axes; --log gives log-log axes, on which (b) is
the gradient of the curves in (a).
Slope 1 is the friction-only limit (R proportional to U); slope 0 is the limit where only
inertial (u_in^2) losses act and R no longer depends on U.
Run from cfd/: python geometry/plot_RU.py   (writes ../paper/figures/fig_RU.{pdf,png})
       python geometry/plot_RU.py --log      (log-log axes, writes fig_RU_log.{pdf,png})
Reads every U*_summary.json of the medium-mesh (s1.00) sweep in results/; velocities still
running are skipped, so rerun it as the sweep completes. Also prints R and the slopes.
Panel (a) carries the 1-D pipe-network screening band of each design (screening_band.py; drive
dCp 0.89-1.08, the range of the CFD's own drive; continuous friction law), shaded in the design's colour; --no-band leaves it out.

The paper figure is drawn in MATLAB (paper/figures/matlab/make_fig_RU.m, in the house style of
the group's other papers) from the CSV files this script writes to paper/figures/data/:
fig_RU.csv (design, U, R), fig_RU_slope.csv (design, U at the geometric mean, slope) and
fig_q.csv (intake flow, below). The 1-D screening band (R_lo, R_hi) is no longer written to the paper data (the model lines come from screening_band.py's fig_RU_model_own.csv). The matplotlib version is only a preview, written to
paper/figures/preview/ with --preview.

fig_q.csv: intake flow q = R U A_in of the whole sampler (each design has one inlet tube; A_in is the
bore of the inlet tube at the R plane, tube_id in tims_geometry.py, as the driver uses; R is already for
the full sampler, the half model's flow doubled), in mL/min. Uncertainty from R: the medium mesh reads
low by the design's mesh error at 0.3 m/s (tab:verification; one-sided, taken as velocity-independent,
since the fine mesh at 1.5 m/s agrees with the medium within 0.6 %), so q_hi = q / (1 + e) is the
mesh-corrected value; time_pct is the run's |drift| between the halves of its averaging window (% of q).
The exponent of q in U is 1 + d ln R / d ln U (fig_RU_slope.csv).
"""
import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from screening_band import DCP, PATHS, band  # noqa: E402
from table_R import metrics  # noqa: E402
from tims_geometry import SPECS  # noqa: E402

INK, MUTED, GRID = "#2c3e50", "#7f8c8d", "#e5e8ea"
SWEEP = [0.154, 0.3, 0.585, 1.0, 1.25, 1.5]
# Colours: three hues plus near-black, which separate under protan/deutan simulation for every
# pair (OKLab dE >= 9.2). The two variants share their parent's colour with a dashed line and
# open marker, so no design is told apart by colour alone.
DESIGNS = [  # label, results folder, run, colour, line style, marker, filled
    ("Phillips", "urans_v5", "phillips2000", "#2c2c2c", "-", "o", True),
    ("Foets", "foets", "foets2020", "#2c2c2c", (0, (4, 2)), "o", False),
    ("Perks", "urans_v5", "perks2014", "#eb6834", "-", "s", True),
    ("Perks + meter", "urans_v5", "perks2014_flowmeter", "#eb6834", (0, (4, 2)), "s", False),
    ("Elliott", "urans_v5", "elliott2017", "#2a78d6", "-", "^", True),
    ("McDonald", "urans_v5bore", "mcdonald2010", "#1baf7a", "-", "D", True),
]
# medium-mesh error at 0.3 m/s (tab:verification), fraction; negative = medium reads low
MESH_ERR = {"Phillips": -0.025, "Foets": -0.027, "Perks": -0.0055, "Perks + meter": -0.0036,
            "Elliott": -0.014, "McDonald": -0.016}
OUT = Path("../paper/figures/preview/fig_RU")
DATA = Path("../paper/figures/data")


def joined_R(cont, U):
    """Mean R of a run continued from a checkpoint: part 1 to t_ckpt, then part 2 shifted by
    t_ckpt, averaged over the second half of the full run as for every other run."""
    info = json.load(open(cont / f"U{U:.3f}_continuation.json"))
    t0, total = info["t_ckpt"], info["sim_time_total"]
    rows = [(float(r["t"]), float(r["R"]))
            for r in csv.DictReader(open(cont / f"U{U:.3f}_part1_timeseries.csv"))
            if float(r["t"]) <= t0 + 1e-6]
    rows += [(float(r["t"]) + t0, float(r["R"]))
             for r in csv.DictReader(open(cont / f"U{U:.3f}_timeseries.csv"))]
    return float(np.mean([R for t, R in rows if total / 2 - 1e-6 <= t <= total + 1e-6]))


def load(folder, run):
    base = Path("results") / folder / f"{run}_s1.00_sst-gamma"
    cont = Path("results") / f"{folder}_cont" / f"{run}_s1.00_sst-gamma"
    R = {}
    for U in SWEEP:
        summary = base / f"U{U:.3f}_summary.json"
        if summary.exists():
            R[U] = json.load(open(summary))["mean"]["R"]
        elif (cont / f"U{U:.3f}_continuation.json").exists():
            R[U] = joined_R(cont, U)
    return R


def write_csv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows([r[0], *(f"{v:.6g}" for v in r[1:])] for r in rows)


def slopes(R):
    U = sorted(R)
    return [(np.sqrt(a * b), np.log(R[b] / R[a]) / np.log(b / a)) for a, b in zip(U, U[1:])]


def main():
    linear = "--log" not in sys.argv
    out = OUT if linear else OUT.with_name(OUT.name + "_log")
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6,
                         "font.family": "DejaVu Sans"})
    fig, (axR, axs) = plt.subplots(1, 2, figsize=(6.2, 3.3))
    handles = []
    rows_R, rows_s, rows_b, rows_q = [], [], [], []
    for label, folder, run, colour, ls, mk, filled in DESIGNS:
        if "--no-band" not in sys.argv and run in PATHS:
            lo, hi = zip(*band(run, SWEEP))
            axR.fill_between(SWEEP, lo, hi, color=colour, alpha=0.11, lw=0, zorder=0)
            rows_b += [(label, u, a, b) for u, a, b in zip(SWEEP, lo, hi)]
        R = load(folder, run)
        style = dict(color=colour, ls=ls, lw=1.2, marker=mk, ms=4.2, mew=1.0, mec=colour,
                     mfc=colour if filled else "white")
        U = sorted(R)
        h, = axR.plot(U, [R[u] for u in U], label=label, **style)
        handles.append(h)
        s = slopes(R)
        rows_R += [(label, u, R[u]) for u in U]
        d = SPECS[run].tube_id
        base = Path("results") / folder / f"{run}_s1.00_sst-gamma"
        for u in U:
            q = R[u] * u * np.pi * d**2 / 4 * 6e7            # m3/s -> mL/min
            m = metrics(base, u)
            rows_q.append((label, u, R[u], 1e3 * d, q, q / (1 + MESH_ERR[label]),
                           abs(m["halves"]) if m else float("nan")))
        rows_s += [(label, v, m) for v, m in s]
        if s:
            axs.plot(*zip(*s), **style)
        print(f"{label:14s} R: " + "  ".join(f"{u:g}={R[u]:.4f}" for u in U))
        print(f"{'':14s} slope: " + "  ".join(f"{m:.3f}@{v:.3f}" for m, v in
                                               ((m, v) for v, m in s)))

    ticks = [0.154, 0.3, 0.585, 1.0, 1.5]
    for ax in (axR, axs):
        if linear:
            lt = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
            ax.set_xticks(lt, [f"{t:g}" for t in lt])
            ax.set_xlim(0, 1.6)
        else:
            ax.set_xscale("log")
            ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
            ax.minorticks_off()
            ax.set_xlim(0.135, 1.7)
        ax.set_xlabel("free-stream velocity $U$ (m s$^{-1}$)")
        ax.grid(color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

    if linear:
        axR.set_ylim(0, 0.5)
    else:
        axR.set_yscale("log")
        yt = [0.02, 0.05, 0.1, 0.2, 0.5]
        axR.set_yticks(yt, [f"{t:g}" for t in yt])
        axR.yaxis.set_minor_locator(plt.NullLocator())
        axR.set_ylim(0.018, 0.6)
    axR.set_ylabel("intake ratio $R$")
    axR.set_title("(a) intake ratio", loc="left", fontsize=8.5, color=INK)

    for y, text in ((1, "friction only: $R \\propto U$"),
                    (0, "inertial losses only: $R$ independent of $U$")):
        axs.axhline(y, color=MUTED, lw=0.7, ls=(0, (2, 2)))
        axs.text(axs.get_xlim()[0] + 0.03, y + 0.025, text, ha="left", va="bottom",
                 fontsize=6.5, color=MUTED)
    axs.set_ylim(-0.05, 1.12)
    axs.set_ylabel("local slope d ln $R$ / d ln $U$")
    axs.set_title("(b) sensitivity to velocity", loc="left", fontsize=8.5, color=INK)

    if "--no-band" not in sys.argv:
        handles.append(Patch(color=MUTED, alpha=0.3, lw=0,
                             label=f"1-D model, $\Delta C_p$ {DCP[0]:g}–{DCP[1]:g}"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=7,
               handlelength=2.6, columnspacing=1.4, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(w_pad=2.0, rect=(0, 0.13, 1, 1))
    if linear:
        write_csv(DATA / "fig_RU.csv", ("design", "U", "R"), rows_R)
        write_csv(DATA / "fig_RU_slope.csv", ("design", "U", "slope"), rows_s)
        write_csv(DATA / "fig_q.csv", ("design", "U", "R", "d_in_mm", "q_mL_min", "q_hi_mL_min", "time_pct"),
                  rows_q)
    if "--preview" in sys.argv:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out.with_suffix(".pdf"))
        fig.savefig(out.with_suffix(".png"), dpi=300)


if __name__ == "__main__":
    main()
