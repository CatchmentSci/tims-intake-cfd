"""Annotated design drawing of a TIMS case, drawn from the same spec used to build the CFD geometry.

Dimensions taken from the literature are labelled in dark blue; assumed dimensions in orange.
The bottom-right panel plots the actual generated triangulation on the symmetry plane.
Usage: python plot_design.py phillips2000 [--mesh out/scale_1.00/phillips2000]
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon, Rectangle

from tims_geometry import SPECS, DomainSpec, x_stations, R_JET, R_TIP, R_NEAR, R_WAKE
from check_geometry import read_fluent_boundary_mesh

SRC, ASM, SOLID, EDGE = "#1f3a93", "#d35400", "#c8ccd2", "#4a4f57"
LABEL_BOX = dict(fc="white", ec="none", alpha=0.85, pad=0.6)


def sampler_polygons(s, xs):
    """Solid parts in the x-z section (mm). Each item is a list of (x, z) vertices."""
    mm = 1000
    rb_o, rb_i, rt_o, rt_i = (v * mm / 2 for v in (s.body_od, s.body_id, s.tube_od, s.tube_id))
    x = {k: v * mm for k, v in xs.items()}
    polys = []
    for sgn in (1, -1):
        polys.append([(0, sgn * rb_i), (x["rear"], sgn * rb_i), (x["rear"], sgn * rb_o), (0, sgn * rb_o)])
        polys.append([(0, sgn * rt_o), (x["cav0"], sgn * rt_o), (x["cav0"], sgn * rb_i), (0, sgn * rb_i)])
        polys.append([(x["cav1"], sgn * rt_o), (x["rear"], sgn * rt_o), (x["rear"], sgn * rb_i),
                      (x["cav1"], sgn * rb_i)])
        polys.append([(0, sgn * rt_o), (0, sgn * rb_o), (x["cone_tip"], sgn * rt_o)])
        if s.tail_cone_len > 0:
            polys.append([(x["rear"], sgn * rt_o), (x["rear"], sgn * rb_o), (x["tail_tip"], sgn * rt_o)])
        polys.append([(x["in_tip"], sgn * rt_i), (x["in_inner"], sgn * rt_i),
                      (x["in_inner"], sgn * rt_o), (x["in_tip"], sgn * rt_o)])
        polys.append([(x["out_inner"], sgn * rt_i), (x["out_tip"], sgn * rt_i),
                      (x["out_tip"], sgn * rt_o), (x["out_inner"], sgn * rt_o)])
    return polys


def draw_sampler(ax, s, xs):
    for p in sampler_polygons(s, xs):
        ax.add_patch(Polygon(p, closed=True, fc=SOLID, ec=EDGE, lw=0.6, hatch="////"))
    ax.axhline(0, color="0.55", lw=0.5, ls=(0, (8, 3, 2, 3)))


def dim_h(ax, x0, x1, z, text, color, above=True, fs=8):
    ax.annotate("", (x0, z), (x1, z), arrowprops=dict(arrowstyle="<->", color=color, lw=0.8,
                                                      shrinkA=0, shrinkB=0))
    ax.text(0.5 * (x0 + x1), z + (1.2 if above else -1.2), text, color=color, fontsize=fs,
            ha="center", va="bottom" if above else "top", bbox=LABEL_BOX)


def dim_v(ax, x, z0, z1, text, color, side="right", fs=8):
    ax.annotate("", (x, z0), (x, z1), arrowprops=dict(arrowstyle="<->", color=color, lw=0.8,
                                                      shrinkA=0, shrinkB=0))
    ax.text(x + (1.5 if side == "right" else -1.5), 0.5 * (z0 + z1), text, color=color,
            fontsize=fs, ha="left" if side == "right" else "right", va="center", rotation=0)


def main(case, mesh_stem, out):
    s, dom = SPECS[case], DomainSpec()
    xs = x_stations(s)
    X = {k: v * 1000 for k, v in xs.items()}
    rb_o, rb_i, rt_o, rt_i = (v * 500 for v in (s.body_od, s.body_id, s.tube_od, s.tube_id))

    fig = plt.figure(figsize=(16, 10.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.55, 1.1, 1.2], hspace=0.28, wspace=0.18)

    # (a) full sampler, true scale --------------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    draw_sampler(ax, s, xs)
    ax.annotate("", (X["in_tip"] - 25, 0), (X["in_tip"] - 140, 0),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.2))
    ax.text(X["in_tip"] - 82, 7, "flow", fontsize=9, ha="center")
    dim_h(ax, X["cav0"], X["cav1"], rb_o + 22, f"settling chamber {s.body_len*1000:.0f} mm (internal)", SRC)
    dim_h(ax, X["cone_tip"], 0, rb_o + 22, f"nose cone {s.cone_len*1000:.0f} mm", ASM)
    dim_h(ax, X["in_tip"], X["in_inner"], -rb_o - 22, f"inlet tube {s.inlet_len*1000:.0f} mm",
          SRC, above=False)
    dim_h(ax, X["out_inner"], X["out_tip"], -rb_o - 22, f"outlet tube {s.outlet_len*1000:.0f} mm",
          SRC, above=False)
    dim_v(ax, 800, -rb_i, rb_i, f"i.d. {s.body_id*1000:.0f} mm", SRC)
    dim_v(ax, 300, -rb_o, rb_o, f"o.d. {s.body_od*1000:.0f} mm", ASM, side="left")
    ax.text(X["out_tip"] + 15, 14, "outlet in the\nwake of a blunt\nrear cap", fontsize=8, va="bottom")
    ax.text(510, 14, "chamber volume 7.54 L   ·   area ratio (98/4)² = 600", fontsize=8.5,
            ha="center", va="bottom", color=SRC)
    ax.set_xlim(X["in_tip"] - 170, X["out_tip"] + 120)
    ax.set_ylim(-rb_o - 45, rb_o + 45)
    ax.set_aspect("equal")
    ax.set_title(f"(a) Phillips et al. (2000) TIMS: longitudinal section through the axis "
                 f"(true scale, mm)", loc="left", fontsize=10.5)
    ax.set_yticks([-55, 0, 55]); ax.tick_params(labelsize=8)

    # (b) inlet detail ------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    draw_sampler(ax, s, xs)
    ax.annotate("", (X["in_tip"] + 2, 0), (X["in_tip"] - 22, 0),
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1))
    dim_h(ax, X["in_tip"], X["cone_tip"], 9, "20 mm\nbeyond tip", ASM)
    dim_h(ax, X["cav0"], X["in_inner"], 9, "20 mm\ninto chamber", SRC)
    lead = dict(arrowstyle="-", lw=0.6)
    ax.annotate("bore 4 mm i.d.", (X["in_tip"] + 14, 0.8), (X["in_tip"] + 26, 27), color=SRC,
                fontsize=8, arrowprops=dict(lead, color=SRC))
    ax.annotate("tube 6 mm o.d.", (X["in_tip"] + 8, -rt_o), (X["in_tip"] - 22, -20), color=ASM,
                fontsize=8, arrowprops=dict(lead, color=ASM))
    ax.annotate("interior_mouth_upstream (flow-rate plane, 2 mm inside)", (X["in_tip"] + rt_i, -rt_i), (X["in_tip"] - 22, -29),
                color="#7d3c98", fontsize=7.5, arrowprops=dict(lead, color="#7d3c98"))
    dim_h(ax, 0, X["cav0"], 30, "cap 10 mm", ASM)
    ax.text(-40, -18, "sealed funnel\n(solid in model)", fontsize=8, color=ASM, ha="center", bbox=LABEL_BOX)
    ax.plot([X["in_tip"] + rt_i] * 2, [-rt_i, rt_i], color="#7d3c98", lw=2)
    ax.set_xlim(X["in_tip"] - 25, X["in_inner"] + 22); ax.set_ylim(-33, 40)
    ax.set_aspect("equal"); ax.tick_params(labelsize=8)
    ax.set_title("(b) Inlet: tube, nose cone and upstream cap (mm)", loc="left", fontsize=10)

    # (c) outlet detail -----------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    draw_sampler(ax, s, xs)
    dim_h(ax, X["out_inner"], X["cav1"], 9, "20 mm\ninto chamber", SRC)
    dim_h(ax, X["rear"], X["out_tip"], 9, f"{X['out_tip'] - X['rear']:.0f} mm beyond rear cap", SRC)
    dim_h(ax, X["cav1"], X["rear"], 30, "cap 10 mm", ASM)
    ax.plot([X["out_tip"] - rt_i] * 2, [-rt_i, rt_i], color="#7d3c98", lw=2)
    ax.text(X["out_tip"] - rt_i, -6, "interior_mouth_downstream", fontsize=7.5, ha="center", va="top", color="#7d3c98")
    ax.set_xlim(X["out_inner"] - 20, X["out_tip"] + 22); ax.set_ylim(-33, 40)
    ax.set_aspect("equal"); ax.tick_params(labelsize=8)
    ax.set_title("(c) Outlet: tube and blunt rear cap (mm)", loc="left", fontsize=10)

    # (d) legend / key -------------------------------------------------------------------
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    ax.set_title("Key and sources", loc="left", fontsize=10)
    ax.add_patch(Rectangle((0.02, 0.86), 0.08, 0.07, fc=SOLID, ec=EDGE, hatch="////", transform=ax.transAxes))
    ax.text(0.13, 0.895, "solid (not meshed); white = water", transform=ax.transAxes, va="center", fontsize=9)
    ax.text(0.02, 0.76, "■ dark blue", color=SRC, transform=ax.transAxes, fontsize=9, weight="bold")
    ax.text(0.25, 0.76, "published: Phillips et al. 2000; Russell et al. 2000", transform=ax.transAxes, fontsize=8.5)
    ax.text(0.02, 0.68, "■ orange", color=ASM, transform=ax.transAxes, fontsize=9, weight="bold")
    ax.text(0.25, 0.68, "assumed: confirm or replace", transform=ax.transAxes, fontsize=8.5)
    notes = [
        "Assumptions",
        "• body o.d. 110 mm: standard UK PVC pipe",
        "• end caps 10 mm, flush with body (Smith & Owens report",
        "   8.68 L vs 7.54 L; caps with recesses may explain it)",
        "• funnel: 100 mm cone, body o.d. → tube o.d., sealed",
        "• tube o.d. 6 mm (6 × 4 mm nylon pneumatic tube)",
        "• mounting uprights omitted",
        "",
        "Model",
        "• half domain, symmetry plane through the axis",
        "• water 20 °C; uniform inflow U = 0.10–2.5 m/s",
        "• intake driven by the pressure field (not imposed)",
    ]
    for i, t in enumerate(notes):
        ax.text(0.02, 0.57 - i * 0.052, t, transform=ax.transAxes, fontsize=8.3,
                weight="bold" if t in ("Assumptions", "Model") else "normal",
                color=ASM if t.startswith("•") and i < 7 else "k")

    # (e) domain and sizing regions -------------------------------------------------------
    ax = fig.add_subplot(gs[2, :2])
    x0 = (xs["in_tip"] - dom.upstream) * 1000
    x1 = (xs["out_tip"] + dom.downstream) * 1000
    H = dom.half_width * 1000
    ax.add_patch(Rectangle((x0, -H), x1 - x0, 2 * H, fc="#f4f8fb", ec="k", lw=1))
    regions = [
        ("near-body region (6 mm faces)", X["in_tip"] - 150, X["out_tip"] + 150, R_NEAR, "#2e86c1"),
        ("wake region (10 mm)", X["out_tip"] + 150, X["out_tip"] + 1650, R_WAKE, "#17a589"),
    ]
    for label, a, b, r, c in regions:
        r *= 1000
        ax.add_patch(Rectangle((a, -r), b - a, 2 * r, fill=False, ec=c, lw=1, ls="--"))
        ax.text(a + 10, r + 12, label, color=c, fontsize=8)
    for p in sampler_polygons(s, xs):
        ax.add_patch(Polygon(p, closed=True, fc="#555", ec="#555", lw=0.3))
    ax.annotate("velocity inlet U\nTI 5 %, μt/μ 10", (x0, 0), (x0 + 260, 330), fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6))
    ax.text(x1 - 20, 0, "pressure\noutlet\n0 Pa", fontsize=8.5, ha="right", va="center")
    ax.text(0.5 * (x0 + x1), H + 25, "slip (symmetry) far-field boundaries", fontsize=8.5, ha="center")
    dim_h(ax, x0, X["in_tip"], -H - 60, "1.0 m", "k", above=False, fs=8.5)
    dim_h(ax, X["out_tip"], x1, -H - 60, "3.0 m", "k", above=False, fs=8.5)
    dim_v(ax, x1 + 40, -H, H, "1.1 m", "k")
    ax.text(0.5 * (x0 + x1), -H - 330, "half model: 0.55 m deep behind this symmetry plane  ·  "
            "blockage 0.8 %\ninner regions (jet core r = 12 mm, tip spheres r = 30 mm) are too small to show at this scale",
            fontsize=8, ha="center", color="0.3")
    ax.set_xlim(x0 - 60, x1 + 260); ax.set_ylim(-H - 440, H + 90)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("(e) Computational domain and mesh-sizing regions (true scale, mm)", loc="left", fontsize=10)

    # (f) actual triangulation on the symmetry plane near the inlet tip --------------------
    ax = fig.add_subplot(gs[2, 2])
    if mesh_stem:
        nodes, zones, names = read_fluent_boundary_mesh(Path(mesh_stem).with_suffix(".msh"))
        zid = next(z for z, (n, _) in names.items() if n == "symmetry")
        tri = zones[zid]
        P = nodes[:, [0, 2]] * 1000
        win = (X["in_tip"] - 20, X["in_tip"] + 70, -28, 28)
        keep = np.all((P[tri, 0] > win[0] - 5) & (P[tri, 0] < win[1] + 5) &
                      (P[tri, 1] > win[2] - 5) & (P[tri, 1] < win[3] + 5), axis=1)
        segs = np.concatenate([P[tri[keep]][:, [0, 1]], P[tri[keep]][:, [1, 2]], P[tri[keep]][:, [2, 0]]])
        ax.add_collection(LineCollection(segs, colors="#34495e", linewidths=0.2))
        for p in sampler_polygons(s, xs):
            ax.add_patch(Polygon(p, closed=True, fc="#e5e7e9", ec="#c0392b", lw=0.8))
        ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3]); ax.set_aspect("equal")
        ax.tick_params(labelsize=8)
        ax.set_title("(f) Generated surface mesh, symmetry plane at the inlet\n"
                     "(medium resolution; Fluent Meshing remeshes from this)", loc="left", fontsize=10)

    fig.suptitle("TIMS baseline CFD geometry (v5): Phillips et al. (2000) sampler", fontsize=13,
                 x=0.01, ha="left", y=0.995, weight="bold")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=sorted(SPECS))
    ap.add_argument("--mesh", default=str(Path(__file__).parent / "out/v5/scale_1.00/phillips2000"))
    ap.add_argument("--out", default=str(Path(__file__).parent / "phillips2000_design.png"))
    a = ap.parse_args()
    main(a.case, a.mesh, a.out)
