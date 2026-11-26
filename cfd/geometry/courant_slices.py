"""Convective Courant number of the production runs, from the symmetry-plane slices (
external 4: 'give maximum and representative Courant numbers').

For each face of the y = 0 slice (one value per cell behind it), CFL = |u| dt / Lx + |w| dt / Lz,
with Lx and Lz the face's extent along x and z (so thin prism cells along a wall are not penalised
for their small wall-normal size when the flow runs along them); |v| is zero on the plane. Reported: the maximum, the
99th percentile, and the area-weighted mean over the body region (|x| within the body plus 0.2 m,
|z| < 0.2 m) and over the bores (cells with h < 0.4 mm, which only the bore sizing produces), and
the cell that gives the maximum. dt is the production time step at each velocity (0.006 m / U).

Reads results/mesh_slices/regime/regime_<case>_U<U>_field.npz (0.3-1.5 m/s; Phillips, Elliott,
Perks, Perks + meter) and results/mesh_slices/<case>_v5*_s1.00_U0.300_field.npz (all six at 0.3 m/s).
Writes results/checks/courant.csv. Run from cfd/: python geometry/courant_slices.py
"""
import csv
import glob
import re
from pathlib import Path

import numpy as np

DT = {0.154: 0.04, 0.3: 0.02, 0.585: 0.01, 1.0: 0.006, 1.25: 0.0048, 1.5: 0.004}
S = Path("results/mesh_slices")


def stats(path, case, U):
    d = np.load(path, allow_pickle=True)
    off, xz = d["off"], d["xz"]
    c = np.add.reduceat(xz, off[:-1]) / np.diff(off)[:, None]
    h = np.sqrt(d["area"])
    lx = np.maximum.reduceat(xz[:, 0], off[:-1]) - np.minimum.reduceat(xz[:, 0], off[:-1])
    lz = np.maximum.reduceat(xz[:, 1], off[:-1]) - np.minimum.reduceat(xz[:, 1], off[:-1])
    cfl = np.abs(d["u"]) * DT[U] / lx + np.abs(d["w"]) * DT[U] / lz
    body = (np.abs(c[:, 1]) < 0.2) & (c[:, 0] > -0.25) & (c[:, 0] < 1.45)
    bore = h < 4e-4
    i = int(np.argmax(cfl))
    w = lambda m: float(np.average(cfl[m], weights=d["area"][m])) if m.any() else float("nan")
    return dict(case=case, U=U, dt=DT[U], cfl_max=round(float(cfl.max()), 2),
                cfl_p99=round(float(np.percentile(cfl, 99)), 2),
                cfl_mean_body=round(w(body), 3), cfl_max_bore=round(float(cfl[bore].max()), 2) if bore.any() else "",
                cfl_mean_bore=round(w(bore), 2), at_x_m=round(float(c[i, 0]), 4), at_z_m=round(float(c[i, 1]), 4),
                at_h_mm=round(float(h[i] * 1e3), 3))


rows = []
for p in sorted(glob.glob(str(S / "regime" / "regime_*_field.npz"))):
    m = re.search(r"regime_(.+)_U([\d.]+)_field", p)
    rows.append(stats(p, m.group(1), float(m.group(2))))
for case in ("phillips2000", "perks2014", "perks2014_flowmeter", "elliott2017", "mcdonald2010", "foets2020"):
    if any(r["case"] == case and r["U"] == 0.3 for r in rows):
        continue
    f = sorted(glob.glob(str(S / f"{case}_v5*_s1.00_U0.300_field.npz")))
    if f:
        rows.append(stats(f[-1], case, 0.3))
out = Path("results/checks/courant.csv")
with out.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(f"{r['case']:20s} {r['U']:5} max {r['cfl_max']:6} p99 {r['cfl_p99']:6} body {r['cfl_mean_body']:6} "
          f"bore max {r['cfl_max_bore']} mean {r['cfl_mean_bore']}  at x {r['at_x_m']} z {r['at_z_m']} h {r['at_h_mm']} mm")
print("wrote", out)
