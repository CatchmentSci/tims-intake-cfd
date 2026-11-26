"""Flow field on the symmetry plane for one solved case, from a slice written on the HPC by
fluent/field_slice.py (one value per face, from the cell behind it).
Usage: python plot_flow.py [--field ../results/mesh_slices/v4_s1.00_U0.300_field.npz] [--U 0.3]
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize, TwoSlopeNorm
from matplotlib.path import Path as MPath
from scipy.interpolate import griddata

from tims_geometry import SPECS, x_stations
from plot_design import sampler_polygons, LABEL_BOX
from plot_mesh import polys_in, solids

SEQ = LinearSegmentedColormap.from_list("seq", ["#f4f8fc", "#c6dbef", "#6baed6", "#2171b5", "#08306b"])
DIV = LinearSegmentedColormap.from_list("div", ["#b03a2e", "#e8a598", "#e6e6e6", "#9ab8dc", "#1f4e99"])
INK = "#2c3e50"


def load(path, U):
    d = np.load(path)
    off, xz = d["off"], d["xz"].astype(float) * 1000            # mm
    cen = np.add.reduceat(xz, off[:-1]) / np.diff(off)[:, None]
    u, v, w = (d[k].astype(float) / U for k in "uvw")
    return dict(off=off, xz=xz, cen=cen, ux=u, uz=w, mag=np.sqrt(u**2 + v**2 + w**2),
                h=np.sqrt(d["area"]) * 1000)


def field_panel(ax, m, s, xs, win, val, cmap, norm, title, pad=5):
    P, sel = polys_in(m, win, pad=pad)
    pc = PolyCollection(P, array=val[sel], cmap=cmap, norm=norm, edgecolors="face", linewidths=0.15)
    ax.add_collection(pc)
    solids(ax, s, xs, fc="#d5d8dc")
    ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3]); ax.set_aspect("equal")
    ax.tick_params(labelsize=8)
    ax.set_title(title, loc="left", fontsize=10)
    return pc


def streamlines(ax, m, s, xs, win, n=(360, 220), density=1.4, **kw):
    """Streamlines of the in-plane velocity, interpolated from face centroids onto a grid."""
    x = np.linspace(win[0], win[1], n[0]); z = np.linspace(win[2], win[3], n[1])
    X, Z = np.meshgrid(x, z)
    c = m["cen"]
    pad = 0.1 * max(win[1] - win[0], win[3] - win[2])
    k = (c[:, 0] > win[0] - pad) & (c[:, 0] < win[1] + pad) & (c[:, 1] > win[2] - pad) & (c[:, 1] < win[3] + pad)
    UX = griddata(c[k], m["ux"][k], (X, Z), method="linear")
    UZ = griddata(c[k], m["uz"][k], (X, Z), method="linear")
    inside = np.zeros(X.shape, bool)
    for p in sampler_polygons(s, xs):
        inside |= MPath(p).contains_points(np.c_[X.ravel(), Z.ravel()]).reshape(X.shape)
    UX[inside] = np.nan; UZ[inside] = np.nan
    ax.streamplot(X, Z, UX, UZ, color=INK, linewidth=0.5, arrowsize=0.6, density=density, zorder=2, **kw)


def main(field, U, out):
    s = SPECS["phillips2000"]
    xs = x_stations(s)
    X = {k: v * 1000 for k, v in xs.items()}
    rb_o, rt_i = s.body_od * 500, s.tube_id * 500
    m = load(field, U)

    fig = plt.figure(figsize=(16, 13))
    gs = fig.add_gridspec(4, 2, height_ratios=[0.8, 1.15, 0.42, 0.75], hspace=0.42, wspace=0.12)

    # (a) overview ------------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    win = (-450, 2300, -260, 260)
    pc = field_panel(ax, m, s, xs, win, m["mag"], SEQ, Normalize(0, 1.4),
                     "(a) Speed |U|/U₀ around the sampler (mm)", pad=60)
    cb = fig.colorbar(pc, ax=ax, fraction=0.02, pad=0.01); cb.set_label("|U| / U₀", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    ax.text(-430, 225, "flow →", fontsize=9, va="top")
    ax.text(-200, 95, "slow flow round the tube tip\nand nose cone (see b)", fontsize=8, ha="center",
            bbox=LABEL_BOX, zorder=6)
    ax.annotate("fastest at the cone–body shoulder (1.38 U₀)", (0, 57), (330, 120), fontsize=8,
                arrowprops=dict(arrowstyle="-", color="0.3", lw=0.6), bbox=LABEL_BOX, zorder=6)
    ax.text(1750, 95, "wake of the blunt rear cap: still 0.8 U₀ on the axis at 2.3 m", fontsize=8,
            ha="center", bbox=LABEL_BOX, zorder=6)

    # (b) inlet tip with streamlines ------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    win = (X["in_tip"] - 30, X["in_tip"] + 30, -18, 18)
    pc = field_panel(ax, m, s, xs, win, m["mag"], SEQ, Normalize(0, 1.2),
                     "(b) Inlet tip: speed and streamlines (mm)")
    streamlines(ax, m, s, xs, win)
    cb = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02); cb.ax.tick_params(labelsize=8)
    ax.text(X["in_tip"] + 29, -16.5, "flow slows ahead of the mouth; only the central\nstreamtube enters the bore. The flow separates\nwhere the tube meets the cone",
            fontsize=8, ha="right", va="bottom", bbox=LABEL_BOX, zorder=6)

    # (c) outlet in the wake: axial velocity ------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    win = (X["rear"] - 30, X["rear"] + 150, -75, 75)
    pc = field_panel(ax, m, s, xs, win, m["ux"], DIV, TwoSlopeNorm(0, -0.4, 1.4),
                     "(c) Rear cap and outlet: axial velocity uₓ/U₀ (mm)")
    streamlines(ax, m, s, xs, win, density=1.1)
    cb = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02); cb.ax.tick_params(labelsize=8)
    ax.text(X["rear"] + 148, -73, "red = reversed flow (down to −0.37 U₀):\nrecirculation bubble 137 mm long behind\n"
            "the cap; the outlet tip (1140) sits inside it", fontsize=8, ha="right", va="bottom", bbox=LABEL_BOX, zorder=6)

    # (d) chamber jet (log scale) -----------------------------------------------------------
    ax = fig.add_subplot(gs[2, :])
    win = (X["cav0"] - 5, X["cav1"] + 5, -rb_o + 3, rb_o - 3)
    pc = field_panel(ax, m, s, xs, win, np.maximum(m["mag"], 1e-5), SEQ, LogNorm(1e-4, 0.3),
                     "(d) Settling chamber: speed |U|/U₀, log scale (mm)", pad=3)
    cb = fig.colorbar(pc, ax=ax, fraction=0.02, pad=0.01); cb.ax.tick_params(labelsize=8)
    ax.text(500, 30, "the inlet jet spreads across the chamber within ~250 mm; beyond 300 mm, 99 % of the plane is below 10⁻³ U₀",
            fontsize=8, ha="center", bbox=LABEL_BOX, zorder=6)

    # (e) speed along the axis ------------------------------------------------------------
    ax = fig.add_subplot(gs[3, :])
    c = m["cen"]
    k = np.abs(c[:, 1]) < np.maximum(0.15 * rt_i, 0.5 * m["h"])     # faces nearest the axis
    o = np.argsort(c[k, 0])
    ax.plot(c[k, 0][o], m["ux"][k][o], color="#2171b5", lw=1.2)
    ax.axhline(1, color="0.6", lw=0.6, ls=":")
    R = 0.1857
    for a, b in ((X["in_tip"], X["in_inner"]), (X["out_inner"], X["out_tip"])):
        ax.axvspan(a, b, color="0.92", zorder=0)
    ax.axvspan(X["cav0"], X["cav1"], color="#f4f8fc", zorder=0)
    ax.text(X["in_tip"] + 75, 0.92, "inlet tube", fontsize=8, ha="center", color="0.35")
    ax.text(X["out_inner"] + 75, 0.92, "outlet tube", fontsize=8, ha="center", color="0.35")
    ax.text(500, 0.92, "settling chamber", fontsize=8, ha="center", color="0.35")
    ax.annotate(f"centreline in the bore 0.345 U₀ = 1.86 × the mean (R = {R:.3f});\n"
                "laminar, Re ≈ 220 (fully developed Poiseuille flow would give 2×)",
                (X["in_tip"] + 70, 0.345), (X["in_tip"] + 260, 0.55),
                fontsize=8, arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6), bbox=LABEL_BOX)
    ax.set_xlim(-450, 2300); ax.set_ylim(-0.25, 1.1)
    ax.set_xlabel("x along the axis (mm)", fontsize=9); ax.set_ylabel("uₓ / U₀ on the axis", fontsize=9)
    ax.tick_params(labelsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", color="0.9", lw=0.6)
    ax.set_title("(e) Axial velocity along the sampler axis", loc="left", fontsize=10)

    fig.suptitle(f"Phillips et al. (2000) TIMS: flow at U₀ = {U} m/s (URANS SST-γ, medium mesh, "
                 "t = 20 s, symmetry plane)", fontsize=13, x=0.01, ha="left", y=0.995, weight="bold")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    here = Path(__file__).parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--field", default=str(here.parent / "results/mesh_slices/v4_s1.00_U0.300_field.npz"))
    ap.add_argument("--U", type=float, default=0.3)
    ap.add_argument("--out", default=str(here / "phillips2000_flow_U0.300.png"))
    a = ap.parse_args()
    main(a.field, a.U, a.out)
