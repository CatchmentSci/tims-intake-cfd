"""Turbulence intensity that reaches the sampler nose, from a symmetry-plane slice written on the HPC
by fluent/field_slice.py. The inlet value decays over the 1 m fetch, so the sensitivity runs are
reported by the inlet value and by the value arriving on the tube axis ahead of the inlet tip.
TI = sqrt(2k/3) / U, with U the free-stream velocity.
Usage: python arriving_ti.py --case phillips2000 --U 0.3 SLICE.npz [SLICE.npz ...]
"""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
AHEAD = (0.9, 0.5, 0.2, 0.1, 0.05)        # m upstream of the inlet tip (the inlet is 1.0 m ahead)
BAND = 0.015                             # m either side of the axis (a narrower band is sparse far upstream)


def axis_profile(path, x_tip, z_axis=0.0):
    d = np.load(path)
    off, xz = d["off"], d["xz"].astype(float)
    cen = np.add.reduceat(xz, off[:-1]) / np.diff(off)[:, None]
    sel = (np.abs(cen[:, 1] - z_axis) < BAND) & (cen[:, 0] < x_tip)
    return cen[sel, 0], d["k"].astype(float)[sel]


def ti_at(x, k, x0, U, half=0.01):
    m = np.abs(x - x0) < half
    return np.sqrt(2 * k[m].mean() / 3) / U if m.any() else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slices", nargs="+")
    ap.add_argument("--case", default="phillips2000")
    ap.add_argument("--U", type=float, default=0.3)
    a = ap.parse_args()
    st = json.loads((HERE / f"out/v5/scale_1.00/{a.case}_zones.json").read_text())["stations_m"]
    x_tip = st["in_tip"]
    print("slice, " + ", ".join(f"TI {s:g} m ahead" for s in AHEAD))
    for p in a.slices:
        x, k = axis_profile(p, x_tip)
        row = [ti_at(x, k, x_tip - s, a.U) for s in AHEAD]
        print(Path(p).name + ", " + ", ".join(f"{100 * t:.2f} %" for t in row))


if __name__ == "__main__":
    main()
