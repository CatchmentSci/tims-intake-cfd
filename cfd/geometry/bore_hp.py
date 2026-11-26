"""Bore pressure gradient, friction excess over Hagen-Poiseuille and centreline/mean ratio
in the Phillips inlet tube, from symmetry-plane slices.
Run from cfd/: python geometry/bore_hp.py [skip_from_tip=0.09] [skip_from_inner_end=0.02]   (metres)
The default window starts 90 mm in, past the ~45 mm entrance length at Re ~ 220."""
import sys, json, numpy as np
sys.path.insert(0, "geometry")
import tims_geometry as tg
s = tg.SPECS["phillips2000"]; X = tg.x_stations(s)
r = s.tube_id / 2; mu = 1.002e-3; U = 0.3
runs = {"coarse": ("phillips2000_v5_s1.65", "urans_v5/phillips2000_s1.65_sst-gamma"),
        "medium": ("phillips2000_v5_s1.00", "urans_v5/phillips2000_s1.00_sst-gamma"),
        "fine":   ("phillips2000_v5_s0.64", "urans_v5/phillips2000_s0.64_sst-gamma"),
        "bore0.5": ("phillips2000_v5bore_s1.00", "urans_v5bore/phillips2000_s1.00_sst-gamma"),
        # the two finest bore meshes
        "bore0.25": ("phillips2000_v5bore4_s1.00", "rerun_v5bore4/phillips2000_s1.00_sst-gamma"),
        "fine bore0.5": ("phillips2000_v5bore_s0.64", "rerun_v5bore/phillips2000_s0.64_sst-gamma")}
import csv
from pathlib import Path
for name, (sl, res) in runs.items():
    if not Path(f"results/mesh_slices/{sl}_U0.300_field.npz").exists():
        print(f"{name:8s} no slice yet"); continue
    d = np.load(f"results/mesh_slices/{sl}_U0.300_field.npz")
    off, xz = d["off"], d["xz"].astype(float)
    cen = np.add.reduceat(xz, off[:-1]) / np.diff(off)[:, None]
    x, z = cen[:, 0], cen[:, 1]
    rows = [q for q in csv.DictReader(open(f"results/{res}/U0.300_timeseries.csv")) if float(q["t"]) >= 6]
    R = np.mean([float(q["R"]) for q in rows[-10:]])
    a, b = X["in_tip"] + (float(sys.argv[1]) if len(sys.argv) > 1 else 0.09), X["in_inner"] - (float(sys.argv[2]) if len(sys.argv) > 2 else 0.02)           # developed part of the tube
    inb = (x > a) & (x < b) & (np.abs(z) < r)
    ax = inb & (np.abs(z) < 0.25 * r)
    g = np.polyfit(x[ax], d["p"][ax], 1)[0]
    g_hp = -8 * mu * R * U / r**2
    cl = d["u"][inb & (np.abs(z) < 0.08 * r)].mean() / (R * U)
    print(f"{name:8s} R {R:.5f}  dp/dx {g:7.1f}  HP {g_hp:7.1f}  excess {100*(g/g_hp-1):4.1f}%  "
          f"HP-implied R {-g*r**2/(8*mu*U):.4f}  centreline/mean {cl:.3f}")
