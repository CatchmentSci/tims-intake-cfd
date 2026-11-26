"""Data for Fig. 1 (the modelled designs), drawn in MATLAB (make_fig_designs.m).

The outlines are taken from the generated boundary mesh, not redrawn: the symmetry-plane triangles
(zone `symmetry`) cover the water on y = 0, so their free edges are the water/solid boundaries.
Edges on the domain boundary are dropped; the remaining free edges are chained into closed loops,
one per solid part as cut by the plane (each loop encloses solid; everything else is water).
Collinear points are removed (tolerance 0.005 mm).

Writes to paper/figures/data/:
- fig_designs_outline.csv: design, loop, x_mm, z_mm. x from the inlet tip (positive downstream),
  z up. Each loop is closed (last point = first). Plot each loop as a filled polygon.
- fig_designs_assumed.csv: design, part, shape (poly / box), loop of vertices x_mm, z_mm: the parts
  marked assumed in Table 1 (Phillips nose cone; Foets outlet tube o.d.; Elliott outlet elbow and
  down-tube; McDonald funnel; meter position), for hatching or a label.
- fig_designs_stations.csv: design, station, x_mm (inlet tip, nose start, body front/rear faces,
  outlet tip), for alignment and labels.

Meshes: production medium meshes (v5; McDonald and Foets v5bore).
Run from cfd/:  python geometry/fig_designs_data.py
"""
import csv
import json
from collections import defaultdict
from pathlib import Path

import sys

import numpy as np

sys.path.insert(0, "geometry")
from check_geometry import read_fluent_boundary_mesh  # noqa: E402

OUT = Path("../paper/figures/data")
CASES = [("Phillips", "v5", "phillips2000"), ("Perks", "v5", "perks2014"),
         ("Perks + meter", "v5", "perks2014_flowmeter"), ("Elliott", "v5", "elliott2017"),
         ("McDonald", "v5bore", "mcdonald2010"), ("Foets", "v5bore", "foets2020")]


def loops_from_triangles(P, tri, keep_edge):
    cnt = defaultdict(int)
    for t in tri:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            cnt[(min(a, b), max(a, b))] += 1
    edges = [e for e, c in cnt.items() if c == 1 and keep_edge(P[e[0]], P[e[1]])]
    nbr = defaultdict(list)
    for a, b in edges:
        nbr[a].append(b)
        nbr[b].append(a)
    used, loops = set(), []
    for a, b in edges:
        if (a, b) in used or (b, a) in used:
            continue
        loop, prev, cur = [a], a, b
        used.add((a, b))
        while cur != a:
            loop.append(cur)
            nxt = [n for n in nbr[cur] if n != prev and (cur, n) not in used and (n, cur) not in used]
            if not nxt:
                break
            used.add((cur, nxt[0]))
            prev, cur = cur, nxt[0]
        loops.append(loop + [a])
    return loops


def simplify(xy, tol=0.005):
    out = [xy[0]]
    for i in range(1, len(xy) - 1):
        p, q, r = np.array(out[-1]), np.array(xy[i]), np.array(xy[i + 1])
        d = r - p
        n = np.hypot(*d)
        if n == 0 or abs(d[0] * (q - p)[1] - d[1] * (q - p)[0]) / n > tol:
            out.append(xy[i])
    out.append(xy[-1])
    return out


def main():
    rows, asm, sta = [], [], []
    for name, ver, case in CASES:
        stem = Path(f"geometry/out/{ver}/scale_1.00/{case}")
        nodes, zones, names = read_fluent_boundary_mesh(stem.with_suffix(".msh"))
        meta = json.loads(stem.with_name(stem.name + "_zones.json").read_text())
        sp, X = meta["spec"], {k: v * 1000 for k, v in meta["stations_m"].items()}
        x0 = X["in_tip"]
        zid = next(z for z, (n, _) in names.items() if n == "symmetry")
        P = nodes[:, [0, 2]] * 1000
        P[:, 0] -= x0
        tri = zones[zid]
        lo, hi = P.min(0), P.max(0)
        on_domain = lambda p, q: any(abs(c - b) < 1e-6 for c in (p[0], q[0]) for b in (lo[0], hi[0])) or \
            any(abs(c - b) < 1e-6 for c in (p[1], q[1]) for b in (lo[1], hi[1]))
        loops = loops_from_triangles(P, tri, lambda p, q: not on_domain(p, q))
        loops = [simplify([tuple(P[i]) for i in L]) for L in loops if len(L) > 3]
        npts = 0
        for k, L in enumerate(sorted(loops, key=lambda L: min(p[0] for p in L))):
            for x, z in L:
                rows.append([name, k, f"{x:.4f}", f"{z:.4f}"])
            npts += len(L)
        allx = [p[0] for L in loops for p in L]
        rt_o = sp["tube_od"] * 500
        # nose start: first solid point behind the tip lying outside the tube's outer radius
        nose = np.inf
        for L in loops:
            for (xa, za), (xb, zb) in zip(L[:-1], L[1:]):
                za, zb = abs(za), abs(zb)
                th = rt_o + 0.2
                if max(za, zb) > th and max(xa, xb) > 0.5:
                    xc = xa if za > th else (xb if za == zb else xa + (th - za) / (zb - za) * (xb - xa))
                    if min(za, zb) > th:
                        xc = min(xa, xb)
                    nose = min(nose, xc)
        st = {"inlet tip": 0.0, "nose start": nose, "body front face": -x0,
              "body rear face": X["rear"] - x0, "outlet tip": X["out_tip"] - x0,
              "aft end": max(allx)}
        for k, v in st.items():
            sta.append([name, k, f"{v:.2f}"])
        rb_o = sp["body_od"] * 500
        if name == "Phillips":
            for sg in (1, -1):
                asm.append([name, "nose cone", "poly", [(X["cone_tip"] - x0, sg * rt_o), (-x0, sg * rt_o), (-x0, sg * rb_o)]])
        if name == "Foets":
            ri, ro = sp["outlet_id"] * 500, sp["outlet_od"] * 500
            for sg in (1, -1):
                asm.append([name, "outlet tube (o.d. assumed)", "poly",
                            [(X["out_inner"] - x0, sg * ri), (X["out_tip"] - x0, sg * ri), (X["out_tip"] - x0, sg * ro),
                             (X["out_inner"] - x0, sg * ro)]])
        if name == "Elliott":
            asm.append([name, "outlet elbow and down-tube", "box",
                        [(X["rear"] - x0, sp["out_tip_z"] * 1000), (max(allx), sp["out_tip_z"] * 1000),
                         (max(allx), 15.0), (X["rear"] - x0, 15.0)]])
        if name == "McDonald":
            asm.append([name, "funnel (hollow, dimensions assumed)", "box",
                        [(0.0, -rb_o), (-x0 + sp["collar_len"] * 1000 + 10, -rb_o),
                         (-x0 + sp["collar_len"] * 1000 + 10, rb_o), (0.0, rb_o)]])
        if name == "Perks + meter":
            x1 = X["in_end"] - x0 - sp["dev_port_len"] * 1000        # meter body ends before its outlet port
            asm.append([name, "meter position (on the axis)", "box",
                        [(x1 - sp["dev_len"] * 1000, -sp["dev_h"] * 1000 + sp["dev_top"] * 1000),
                         (x1, -sp["dev_h"] * 1000 + sp["dev_top"] * 1000), (x1, sp["dev_top"] * 1000),
                         (x1 - sp["dev_len"] * 1000, sp["dev_top"] * 1000)]])
        print(f"{name:14s} {len(loops)} solid loops, {npts} points; nose starts {nose:.1f} mm behind the tip; "
              f"rear face {X['rear'] - x0:.1f}, aft end {max(allx):.1f} mm")
    OUT.mkdir(exist_ok=True)
    with (OUT / "fig_designs_outline.csv").open("w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["design", "loop", "x_mm", "z_mm"]); w.writerows(rows)
    with (OUT / "fig_designs_assumed.csv").open("w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["design", "part", "shape", "vertex", "x_mm", "z_mm"])
        for d, part, shape, V in asm:
            for i, (x, z) in enumerate(V + [V[0]]):
                w.writerow([d, part, shape, i, f"{x:.3f}", f"{z:.3f}"])
    with (OUT / "fig_designs_stations.csv").open("w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["design", "station", "x_mm"]); w.writerows(sta)
    print(f"wrote fig_designs_outline.csv ({len(rows)} rows), fig_designs_assumed.csv, fig_designs_stations.csv")


if __name__ == "__main__":
    main()
