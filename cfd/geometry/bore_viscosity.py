"""Turbulent viscosity ratio nu_t/nu along the inlet path, from symmetry-plane slices written on
the HPC by fluent/field_slice.py (needs its mu_t and mu_lam fields). The check for the inlet-
turbulence runs: ambient nu_t carried into the laminar inlet bore would add friction that is a
model artefact, not a turbulence effect.
Usage: python bore_viscosity.py [--case phillips2000] SLICE.npz [SLICE.npz ...]
"""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent


def regions(st, rt):
    """(label, x0, x1, half-width in z) along the inlet axis z = 0."""
    tip, inner = st["in_tip"], st["in_inner"]
    return [("free stream 0.1 m ahead", tip - 0.11, tip - 0.09, 0.015),
            ("2 mm ahead of the mouth", tip - 0.003, tip - 0.001, rt),
            ("bore, first 10 mm", tip, tip + 0.010, rt),
            ("bore, 10-50 mm", tip + 0.010, tip + 0.050, rt),
            ("bore, rest to the cavity", tip + 0.050, inner, rt)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slices", nargs="+")
    ap.add_argument("--case", default="phillips2000")
    a = ap.parse_args()
    rep = json.loads((HERE / f"out/v5/scale_1.00/{a.case}_zones.json").read_text())
    st, rt = rep["stations_m"], rep["spec"]["tube_id"] / 2
    for p in a.slices:
        d = np.load(p)
        off, xz = d["off"], d["xz"].astype(float)
        cen = np.add.reduceat(xz, off[:-1]) / np.diff(off)[:, None]
        r = d["mu_t"].astype(float) / d["mu_lam"].astype(float)
        print(Path(p).name)
        for lab, x0, x1, hw in regions(st, rt):
            m = (cen[:, 0] > x0) & (cen[:, 0] < x1) & (np.abs(cen[:, 1]) < hw)
            if m.any():
                print(f"  {lab:26s} n={m.sum():5d}  nu_t/nu mean {r[m].mean():8.2f}  max {r[m].max():8.2f}")


if __name__ == "__main__":
    main()
