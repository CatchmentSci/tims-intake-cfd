"""Views of the Fluent volume mesh on the symmetry plane (y = 0), from slices written on the HPC by
fluent/sym_slice.py. Shows the cell-size distribution, the prism layers at the tube lip and the
rear-cap corner, and the inlet tip at the three GCI scales.
Usage: python plot_mesh.py [--slices ../results/mesh_slices] [--out phillips2000_mesh_v5.png]
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon

from tims_geometry import SPECS, x_stations
from plot_design import sampler_polygons, SOLID, EDGE, LABEL_BOX

CELLS = {"1.65": "1.73M", "1.00": "3.97M", "0.64": "9.46M"}   # from the Fluent Meshing logs
EDGE_INK = "#2c3e50"


def load(path):
    d = np.load(path)
    off, xz = d["off"], d["xz"].astype(float) * 1000            # mm
    n = np.diff(off)
    cen = np.add.reduceat(xz, off[:-1]) / n[:, None]
    return dict(off=off, xz=xz, cen=cen, h=np.sqrt(d["area"]) * 1000)


def polys_in(m, win, pad=2.0):
    x0, x1, z0, z1 = win
    c = m["cen"]
    sel = np.nonzero((c[:, 0] > x0 - pad) & (c[:, 0] < x1 + pad) & (c[:, 1] > z0 - pad) & (c[:, 1] < z1 + pad))[0]
    return [m["xz"][m["off"][i]:m["off"][i + 1]] for i in sel], sel


def solids(ax, s, xs, fc=SOLID):
    for p in sampler_polygons(s, xs):
        ax.add_patch(Polygon(p, closed=True, fc=fc, ec=EDGE, lw=0.6, hatch="////", zorder=3))


def edges_panel(ax, m, s, xs, win, title, lw=0.25):
    P, _ = polys_in(m, win)
    ax.add_collection(PolyCollection(P, facecolors="white", edgecolors=EDGE_INK, linewidths=lw))
    solids(ax, s, xs)
    ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3]); ax.set_aspect("equal")
    ax.tick_params(labelsize=8)
    ax.set_title(title, loc="left", fontsize=10)


def main(slices, out):
    s = SPECS["phillips2000"]
    xs = x_stations(s)
    X = {k: v * 1000 for k, v in xs.items()}
    rb_o, rt_o, rt_i = s.body_od * 500, s.tube_od * 500, s.tube_id * 500
    X["mouth_up"] = X["in_tip"] + rt_i          # flow-rate plane, as in tims_geometry.build (v5)
    M = {k: load(Path(slices) / f"v5_s{k}_sym.npz") for k in CELLS}
    m = M["1.00"]

    fig = plt.figure(figsize=(16, 10.8))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.62, 1.05, 0.78], hspace=0.38, wspace=0.16)

    # (a) cell size over the sampler, medium mesh ------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    win = (X["in_tip"] - 200, X["out_tip"] + 330, -150, 150)
    P, sel = polys_in(m, win, pad=40)
    pc = PolyCollection(P, array=m["h"][sel], cmap="Blues_r", norm=LogNorm(0.4, 40),
                        edgecolors="none", linewidths=0)
    ax.add_collection(pc)
    solids(ax, s, xs, fc="#e5e7e9")
    cb = fig.colorbar(pc, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("cell size √(face area), mm", fontsize=8.5); cb.ax.tick_params(labelsize=8)
    ax.annotate("", (X["in_tip"] - 25, 120), (X["in_tip"] - 150, 120),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2))
    ax.text(X["in_tip"] - 88, 127, "flow", fontsize=9, ha="center")
    ax.text(510, 0, "settling chamber: 1.2 mm jet core along the axis,\ncoarsening towards the walls",
            fontsize=8, ha="center", va="center", color="k", bbox=LABEL_BOX, zorder=5)
    ax.text(X["out_tip"] + 170, -128, "wake region from 150 mm\nbehind the outlet tip (v5)",
            fontsize=8, ha="center", bbox=LABEL_BOX)
    ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3]); ax.set_aspect("equal")
    ax.tick_params(labelsize=8)
    ax.set_title("(a) Medium mesh (v5, s1.00, 3.97M cells): cell size on the symmetry plane "
                 "(dark = fine; mm)", loc="left", fontsize=10.5)

    # (b) inlet: tube in the nose cone ------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    edges_panel(ax, m, s, xs, (X["in_tip"] - 14, X["in_tip"] + 26, -12, 12),
                "(b) Inlet tip, medium mesh (mm)", lw=0.2)
    ax.plot([X["mouth_up"]] * 2, [-rt_i, rt_i], color="#7d3c98", lw=1.6, zorder=4)
    ax.annotate("flow-rate plane\n(2 mm inside the bore)", (X["mouth_up"], -rt_i), (X["in_tip"] + 4, -10.5),
                color="#7d3c98", fontsize=7.5, arrowprops=dict(arrowstyle="-", color="#7d3c98", lw=0.6),
                bbox=LABEL_BOX, zorder=6)
    ax.annotate("nose cone", (X["in_tip"] + 23, 7), (X["in_tip"] + 13, 10.3), fontsize=7.5,
                bbox=LABEL_BOX, zorder=6)

    # (c) prism layers at the lip ------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    edges_panel(ax, m, s, xs, (X["in_tip"] - 1.6, X["in_tip"] + 2.4, 1.2, 3.8),
                "(c) Tube lip close-up: prism layers (mm)", lw=0.35)
    ax.text(X["in_tip"] + 1.3, 2.5, "tube wall\n(6 × 4 mm)", fontsize=8, ha="center", va="center",
            bbox=LABEL_BOX, zorder=6)
    ax.text(X["in_tip"] - 1.5, 1.27, "10 layers, first 25 µm, growth 1.2 (0.65 mm total)\n"
            "on the bore, outer tube and end face", fontsize=7.5, va="bottom", bbox=LABEL_BOX, zorder=6)

    # (d) rear-cap corner -------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 2])
    edges_panel(ax, m, s, xs, (X["rear"] - 8, X["rear"] + 8, rb_o - 7, rb_o + 5),
                "(d) Rear-cap outer corner (v4 fault site; mm)", lw=0.3)
    ax.text(X["rear"] - 7.5, rb_o + 4.4, "prisms from the side wall and rear face\nnow fan round the 90° edge",
            fontsize=7.5, va="top", bbox=LABEL_BOX, zorder=6)

    # (e-g) inlet tip at the three GCI scales ------------------------------------------------
    for j, (k, title) in enumerate((("1.65", "(e) Coarse"), ("1.00", "(f) Medium"), ("0.64", "(g) Fine"))):
        ax = fig.add_subplot(gs[2, j])
        edges_panel(ax, M[k], s, xs, (X["in_tip"] - 4, X["in_tip"] + 8, -1, 6),
                    f"{title}, s{k}: {CELLS[k]} cells (mm)", lw=0.3)
        ax.text(X["in_tip"] + 7.7, 5.7, f"first prism {25 * float(k):.0f} µm", fontsize=8,
                ha="right", va="top", bbox=LABEL_BOX, zorder=6)

    fig.suptitle("Phillips et al. (2000) TIMS: Fluent poly-hexcore mesh, geometry v5 "
                 "(symmetry plane through the axis)", fontsize=13, x=0.01, ha="left", y=0.995,
                 weight="bold")
    fig.text(0.01, 0.005, "Hatched grey = solid (not meshed). Face sizes scale with s: bore 0.25 s mm, tube 0.4 s, "
             "tip 1.5 s, jet core 1.2 s, cone 2 s, body 3 s, near region 6 s, wake 10 s, max 40 s mm. "
             "Refinement ratios 1.33 (coarse→medium) and 1.34 (medium→fine).", fontsize=8.5, color="0.3")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    here = Path(__file__).parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--slices", default=str(here.parent / "results/mesh_slices"))
    ap.add_argument("--out", default=str(here / "phillips2000_mesh_v5.png"))
    a = ap.parse_args()
    main(a.slices, a.out)
