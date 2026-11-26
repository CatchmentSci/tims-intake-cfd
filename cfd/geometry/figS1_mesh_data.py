"""Data for Supplemental Figure S1: computational domain, mesh-sizing regions and mesh
close-ups (Phillips, medium mesh s = 1, half model, symmetry plane y = 0).

Writes to paper/figures/data/ (x from the inlet tip, z up, mm):
- figS1_domain.csv: region, shape (rect / circle), x0_mm, x1_mm, z0_mm, z1_mm, xc_mm, zc_mm, r_mm,
  face_mm (target face size at s = 1). Domain box, near-body and wake cylinders and the jet core
  (rectangles on the plane), tip spheres (circles). From tims_geometry.py (DomainSpec, R_*; the region
  construction in build()); the sampler outline itself is in fig_designs_outline.csv.
- figS1_mesh_<inlet|outlet>.csv: x1_mm, z1_mm, x2_mm, z2_mm. Edges of the volume mesh's faces on the
  symmetry plane (the poly-hexcore cells as cut by y = 0, prism layers included), from the field slice
  of the production run, inside the Fig. 2 windows (inlet: -60..120 x +-30 mm from the tip; outlet:
  -90..360 x +-75 mm from the rear outer face, given here from the tip).

Run from cfd/:  python geometry/figS1_mesh_data.py
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "geometry")
from tims_geometry import DomainSpec, R_JET, R_TIP, R_NEAR, R_WAKE  # noqa: E402

OUT = Path("../paper/figures/data")
GEO = Path("geometry/out/v5/scale_1.00/phillips2000_zones.json")
SLICE = Path("results/mesh_slices/phillips2000_v5_s1.00_U0.300_field.npz")
FACE = {"far field": 40, "wake": 10, "near body": 6, "tip sphere": 1.5, "jet core": 1.2}   # mm at s = 1


def main():
    X = {k: v * 1000 for k, v in json.loads(GEO.read_text())["stations_m"].items()}
    dom = DomainSpec()
    t = X["in_tip"]
    x = {k: v - t for k, v in X.items()}
    H = dom.half_width * 1000
    x_in, x_out = -dom.upstream * 1000, x["out_tip"] + dom.downstream * 1000
    near0, near1 = x.get("front", 0.0) - 150, x["out_tip"] + 150
    wake1 = near1 + min(1500, x_out - near1 - 200)
    x_ie = x.get("in_end", x["in_inner"])
    jet0 = x_ie - min(10, 0.5 * (x_ie - x["cav0"]))
    jet1 = x["out_inner"] + min(10, 0.5 * (x["cav1"] - x["out_inner"]))
    rows = [
        ("domain (far field)", "rect", x_in, x_out, -H, H, "", "", "", FACE["far field"]),
        ("near body", "rect", near0, near1, -R_NEAR * 1000, R_NEAR * 1000, "", "", "", FACE["near body"]),
        ("wake", "rect", near1, wake1, -R_WAKE * 1000, R_WAKE * 1000, "", "", "", FACE["wake"]),
        ("jet core", "rect", jet0, jet1, -R_JET * 1000, R_JET * 1000, "", "", "", FACE["jet core"]),
        ("inlet tip sphere", "circle", "", "", "", "", 0.0, 0.0, R_TIP * 1000, FACE["tip sphere"]),
        ("outlet tip sphere", "circle", "", "", "", "", x["out_tip"], 0.0, R_TIP * 1000, FACE["tip sphere"]),
    ]
    OUT.mkdir(exist_ok=True)
    with (OUT / "figS1_domain.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["region", "shape", "x0_mm", "x1_mm", "z0_mm", "z1_mm", "xc_mm", "zc_mm", "r_mm", "face_mm"])
        for r in rows:
            w.writerow([f"{v:.2f}" if isinstance(v, float) else v for v in r])
    print(f"domain x {x_in:.0f}..{x_out:.0f} mm, z +-{H:.0f}; near {near0:.0f}..{near1:.0f}; wake ..{wake1:.0f}; "
          f"jet {jet0:.0f}..{jet1:.0f} mm")

    d = np.load(SLICE)
    off, xz = d["off"], d["xz"].astype(float) * 1000
    xz[:, 0] -= t
    names = [str(n) for n in d["names"]]
    plane = np.array([n == "symmetry" or n.startswith("symmetry:") for n in names])[d["zone"]]
    rear = x["rear"]
    for tag, (wx, wz) in (("inlet", ((-60, 120), (-30, 30))), ("outlet", ((rear - 90, rear + 360), (-75, 75)))):
        edges = set()
        for i in np.flatnonzero(plane):
            P = xz[off[i]:off[i + 1]]
            if P[:, 0].max() < wx[0] or P[:, 0].min() > wx[1] or P[:, 1].max() < wz[0] or P[:, 1].min() > wz[1]:
                continue
            Q = np.round(P, 4)
            for a, b in zip(Q, np.roll(Q, -1, axis=0)):
                e = (tuple(a), tuple(b)) if tuple(a) < tuple(b) else (tuple(b), tuple(a))
                edges.add(e)
        with (OUT / f"figS1_mesh_{tag}.csv").open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["x1_mm", "z1_mm", "x2_mm", "z2_mm"])
            for (a, b) in sorted(edges):
                w.writerow([f"{a[0]:.4f}", f"{a[1]:.4f}", f"{b[0]:.4f}", f"{b[1]:.4f}"])
        print(f"{tag}: {len(edges)} edges in x {wx[0]:.0f}..{wx[1]:.0f}, z {wz[0]}..{wz[1]} mm (from the tip)")


if __name__ == "__main__":
    main()
