"""Data for Fig. 2 (flow at 0.3 m/s), drawn in MATLAB (make_fig_flow.m).

From the symmetry-plane slices (fluent/field_slice.py) of the final fields of the production
medium runs at 0.3 m/s. Face-centre values are interpolated linearly onto regular grids; grid
points farther than 0.75 face sizes from any fluid face centre are inside a solid and are NaN.

Writes to paper/figures/data/:
- fig_flow_inlet_<key>.csv  (Phillips, Perks, Elliott, McDonald): common inlet window, x from the
  inlet tip; columns x_mm, z_mm, umag, u, w (all over U0). Grid step INLET["step"].
- fig_flow_wake_<key>.csv   (same four): common outlet/wake window, x from the body's rear outer
  face (x_tip_mm = x from the inlet tip as well); columns x_mm, x_tip_mm, z_mm, umag, u, w.
- fig_flow_axis.csv (all six designs): design, x_mm (from the inlet tip), u (u/U0 on z = 0; NaN
  inside solids), from 1 m ahead of the tip to 2 m behind the rear face.
Window sizes: inlet -60..120 x +-30 mm, wake
-90..360 x +-75 mm, both 3:1); change them here.

Run from cfd/:  python geometry/fig_flow_data.py
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import griddata
from matplotlib.path import Path as MPath

sys.path.insert(0, "geometry")
from plot_flow import load  # noqa: E402

U0 = 0.3
SL = Path("results/mesh_slices")
OUT = Path("../paper/figures/data")
INLET = dict(x=(-60.0, 120.0), z=(-30.0, 30.0), step=0.5)
WAKE = dict(x=(-90.0, 360.0), z=(-75.0, 75.0), step=1.5)
AXIS = dict(ahead=1000.0, behind=2000.0, step=1.0)
# key, label, zones json (geometry), slice
CASES = [
    ("phillips", "Phillips", "v5/scale_1.00/phillips2000", "phillips2000_v5_s1.00"),
    ("perks", "Perks", "v5/scale_1.00/perks2014", "perks2014_v5_s1.00"),
    ("perks_meter", "Perks + meter", "v5/scale_1.00/perks2014_flowmeter", "perks2014_flowmeter_v5_s1.00"),
    ("elliott", "Elliott", "v5/scale_1.00/elliott2017", "elliott2017_v5_s1.00"),
    ("mcdonald", "McDonald", "v5bore/scale_1.00/mcdonald2010", "mcdonald2010_v5bore_s1.00"),
    ("foets", "Foets", "v5bore/scale_1.00/foets2020", "foets2020_v5bore_s1.00"),
]
WINDOWED = ("phillips", "perks", "elliott", "mcdonald")


def stations(geo):
    meta = json.loads(Path(f"geometry/out/{geo}_zones.json").read_text())
    return {k: v * 1000 for k, v in meta["stations_m"].items()}


def solids(label):
    rows = [r for r in csv.DictReader(open(OUT / "fig_designs_outline.csv")) if r["design"] == label]
    loops = {}
    for r in rows:
        loops.setdefault(r["loop"], []).append((float(r["x_mm"]), float(r["z_mm"])))
    return [MPath(np.array(v)) for v in loops.values()]


def sample(m, xs, zs, polys, tip, keys=("mag", "ux", "uz")):
    """Values at points (xs, zs) in model mm; NaN inside the solid outlines (tip-referenced)."""
    c = m["cen"]
    pts = np.c_[xs, zs]
    pad = 10.0
    k = (c[:, 0] > xs.min() - pad) & (c[:, 0] < xs.max() + pad) & (c[:, 1] > zs.min() - pad) & (c[:, 1] < zs.max() + pad)
    q = np.c_[xs - tip, zs]
    solid = np.zeros(len(xs), bool)
    for P in polys:
        solid |= P.contains_points(q)
    out = []
    for key in keys:
        v = griddata(c[k], m[key][k], pts, method="linear")
        v[solid] = np.nan
        out.append(v)
    return out


def write(name, header, cols, fmt):
    with (OUT / name).open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for row in zip(*cols):
            w.writerow([("" if not np.isfinite(v) else f) % v if f else v for v, f in zip(row, fmt)])


def main():
    OUT.mkdir(exist_ok=True)
    axis_rows = []
    for key, label, geo, sl in CASES:
        f = SL / f"{sl}_U{U0:.3f}_field.npz"
        if not f.exists():
            print(f"{label}: slice {f} not here yet, skipped")
            continue
        X = stations(geo)
        m = load(f, U0)
        tip, rear = X["in_tip"], X["rear"]
        polys = solids(label)
        if key in WINDOWED:
            for tag, W, x_ref in (("inlet", INLET, tip), ("wake", WAKE, rear)):
                gx = np.arange(W["x"][0], W["x"][1] + 1e-9, W["step"])
                gz = np.arange(W["z"][0], W["z"][1] + 1e-9, W["step"])
                GX, GZ = np.meshgrid(gx, gz)
                xs, zs = GX.ravel(), GZ.ravel()
                mag, u, w = sample(m, xs + x_ref, zs, polys, tip)
                if tag == "inlet":
                    write(f"fig_flow_{tag}_{key}.csv", ["x_mm", "z_mm", "umag", "u", "w"],
                          [xs, zs, mag, u, w], ["%.2f", "%.2f", "%.5f", "%.5f", "%.5f"])
                else:
                    write(f"fig_flow_{tag}_{key}.csv", ["x_mm", "x_tip_mm", "z_mm", "umag", "u", "w"],
                          [xs, xs + rear - tip, zs, mag, u, w], ["%.2f", "%.2f", "%.2f", "%.5f", "%.5f", "%.5f"])
                print(f"{label}: {tag} window {len(gx)} x {len(gz)}, solid {np.isnan(u).mean():.0%}")
        ax = np.arange(tip - AXIS["ahead"], rear + AXIS["behind"] + 1e-9, AXIS["step"])
        (u,) = sample(m, ax, np.zeros_like(ax), polys, tip, keys=("ux",))
        axis_rows += [[label, f"{x - tip:.1f}", "" if not np.isfinite(v) else f"{v:.5f}"] for x, v in zip(ax, u)]
        far = u[(ax - rear > 1500)]
        print(f"{label}: axis {len(ax)} points; u at -60 mm {u[np.argmin(abs(ax - tip + 60))]:.3f}, "
              f"-5 mm {u[np.argmin(abs(ax - tip + 5))]:.3f}, min behind rear {np.nanmin(u[ax > rear]):.3f}, "
              f"1.5-2 m behind {np.nanmean(far):.3f}")
    with (OUT / "fig_flow_axis.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["design", "x_mm", "u"])
        w.writerows(axis_rows)
    print(f"wrote fig_flow_axis.csv ({len(axis_rows)} rows)")


if __name__ == "__main__":
    main()
