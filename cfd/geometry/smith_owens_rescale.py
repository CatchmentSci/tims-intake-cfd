"""Rescale the Smith & Owens (2014) flume mass efficiencies from Eq. 1 (Phillips et al. 2000) to the
open-domain CFD intake ratio of the Phillips model (Discussion).

Their case (ESPL 39: 1197-1207, pp. 1200-1203, Table I): Phillips build, 98 mm x 1 m body, 8.68 L
internal volume, 4 mm x 150 mm tubes. Blockage 4-5 % (110 mm body o.d. in the 0.70 m wide flume with
35 cm of water: 3.9 % rectangular, 4.9 % semicircular base); this script uses the open-domain R, and the
Discussion's 0.18-0.19 adds the blockage by linear scaling of the 0.3 m flume (+16-20 % at 10.6 %). U = 0.28 m/s,
6 h trials, sampler filled with flume water at the start. Efficiency = collected / expected mass,
expected mass from SSC_i (0 h grab) and SSC_s (mean of the 2, 4, 6 h grabs) and the volumes.

Their method text (pp. 1201-1202) allows two readings of the expected mass:
  (a) initial fill + inflow:  M = SSC_i V + SSC_s Q T
  (b) inflow only, first flush at SSC_i:  M = SSC_i V + SSC_s (Q T - V)
Collected mass is unchanged by the intake ratio, so E_cfd = E_pub * M_eq1 / M_cfd.

Run from the repo root:  python cfd/geometry/smith_owens_rescale.py
"""
import csv

import numpy as np
from scipy.interpolate import PchipInterpolator

U, T, V = 0.28, 6 * 3600.0, 8.68e-3            # m/s, s, m3
A_IN = np.pi * 0.002 ** 2                      # 4 mm bore
# Table I: run, SSC at 0 h, mean 2-6 h (mg/L); published grab-sample efficiencies, trials 2-3 mean
TRIALS = {"kaolinite": [(2, 132.8, 136.8), (3, 137.0, 136.3)],
          "river sediment": [(2, 84.2, 17.8), (3, 72.6, 18.1)]}
E_PUB = {"kaolinite": 43.0, "river sediment": 87.0}


def r_eq1(u):
    """Eq. 1, cm/s: log u_in = -2.182 + 2.074 log U."""
    return 10 ** (-2.182 + 2.074 * np.log10(100 * u)) / (100 * u)


def r_cfd():
    rows = [r for r in csv.DictReader(open("paper/figures/data/fig_bands.csv")) if r["design"] == "Phillips"]
    u = np.log([float(r["U"]) for r in rows])
    out = {}
    for k in ("R", "R_lo_total", "R_hi_total"):
        out[k] = float(np.exp(PchipInterpolator(u, np.log([float(r[k]) for r in rows]))(np.log(U))))
    return out


def expected(ssc_i, ssc_s, q, reading):
    qt = q * T
    return ssc_i * V + ssc_s * (qt if reading == "a" else qt - V)


def main():
    re = r_eq1(U)
    rc = r_cfd()
    print(f"R at {U} m/s: Eq. 1 {re:.4f}; CFD open domain {rc['R']:.4f} "
          f"(total band {rc['R_lo_total']:.4f}-{rc['R_hi_total']:.4f}); ratio {re / rc['R']:.3f}")
    for lab, r in (("Eq. 1", re), ("CFD", rc["R"])):
        q = r * U * A_IN
        print(f"  {lab}: Q {q * 1e6:.3f} mL/s, {q * T * 1e3:.2f} L in 6 h = {q * T / V:.2f} sampler volumes")
    for reading in ("a", "b"):
        print(f"reading ({reading}):")
        for sed, trials in TRIALS.items():
            f = {}
            for key in ("R", "R_lo_total", "R_hi_total"):
                f[key] = np.mean([expected(si, ss, re * U * A_IN, reading) / expected(si, ss, rc[key] * U * A_IN, reading)
                                  for _, si, ss in trials])
            e = E_PUB[sed]
            print(f"  {sed:14s} factor {f['R']:.3f} ({f['R_hi_total']:.3f}-{f['R_lo_total']:.3f}); "
                  f"E {e:.0f} % -> {e * f['R']:.0f} % ({e * f['R_hi_total']:.0f}-{e * f['R_lo_total']:.0f})")
    print(f"pure flow ratio (no initial fill, one SSC): {re / rc['R']:.3f}")


if __name__ == "__main__":
    main()
