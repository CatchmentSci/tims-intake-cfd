"""Expected effect of ambient (streamwise) turbulence on the intake ratio R, from the mean-flow CFD.

The river eddies (streamwise integral scale 0.25-0.75 m at the sampler, below) are much larger than the intake, so each mouth
responds to them quasi-steadily (as a pitot tube does). Two limits bracket the effect on mean R:
  A  the bore flow follows the fluctuation (slower than the bore response time):
     <Q> = Q(U) [1 + C sigma^2 / 2],  C = n(n-1) + dn/dlnU from ln Q against ln U (sweep)
  B  the fluctuation is faster than the bore can follow: each mouth rectifies it,
     <dp> = dp(U) [1 + k sigma^2], k = 1 (both mouths scale with <U^2>) to Cp_up/dCp (upstream only),
     and Q responds to the mean with s = d ln Q / d ln dp = n/2,  so <R>/R = 1 + s k sigma^2.
The bore is a first-order low-pass: tau = inertance / (d dp / dQ), fc = 1 / (2 pi tau).
A von Karman spectrum (Taylor's hypothesis) gives the fraction phi of the variance below fc, and
the blended estimate is phi*A + (1 - phi)*B.
Changes to the body boundary layer and wake base pressure are not included.
Run from cfd/: python geometry/turbulence_bias.py
"""
import glob
import json

import numpy as np
import pandas as pd
from scipy.integrate import quad

RHO, U = 998.2, 0.3
# inlet + outlet bore length (m) and bore area for the inertance; meter: 7.94 mm tubing only
GEOM = {"phillips2000": (0.30, np.pi * 0.002**2), "perks2014_flowmeter": (0.323, np.pi * 0.00397**2)}
# Streamwise integral scale at the sampler axis, Nezu (1977, thesis, sec. 4.6.4, Eqs. 4.48-4.49):
# L_x/h = B (y/h)^0.5 for y/h <= 0.6, B3 for y/h >= 0.6, with B, B3 = 1.0-1.1 (independent of Fr;
# roughness matters only below y/h ~ 0.1). Axis height y 0.155 m, depths 0.25-3 m.
Z_AXIS, DEPTHS, B_RANGE = 0.155, (0.25, 0.5, 1.0, 3.0), (1.0, 1.1)


def nezu_lx(h, B):
    return B * h * (np.sqrt(Z_AXIS / h) if Z_AXIS / h <= 0.6 else 1.0)


LX = (min(nezu_lx(h, B_RANGE[0]) for h in DEPTHS), max(nezu_lx(h, B_RANGE[1]) for h in DEPTHS))
SIGMA = (0.05, 0.10, 0.20)               # streamwise intensity u'/U


def vk(f, L):
    """von Karman longitudinal spectrum, normalised to unit variance."""
    return 4 * (L / U) / (1 + 70.8 * (f * L / U) ** 2) ** (5 / 6)


def curvature(case):
    pts = sorted((j["U"], j["mean"]["R"]) for j in
                 (json.load(open(f)) for f in glob.glob(f"results/urans_v5/{case}_s1.00_sst-gamma/U*_summary.json")))
    u, r = map(np.array, zip(*pts))
    sel = u <= 0.6
    p = np.polyfit(np.log(u[sel]), np.log(r[sel] * u[sel]), 2)
    x = np.log(U)
    n = 2 * p[0] * x + p[1]
    return n, n * (n - 1) + 2 * p[0]


def main():
    q = 0.5 * RHO * U**2
    print("L_x (Nezu 1977) at y = %.3f m, B = 1.0: " % Z_AXIS + ", ".join(f"h {h:g} m: {nezu_lx(h, 1.0):.3f}" for h in DEPTHS)
          + f"; range used {LX[0]:.2f}-{LX[1]:.2f} m (B 1.0-1.1)")
    for case, (L, A) in GEOM.items():
        h = pd.read_csv(f"results/urans_v5/{case}_s1.00_sst-gamma/U0.300_timeseries.csv")
        h = h[h.t > 10].mean()
        cp_up, dcp = (h.p_mouth_up - h.p_inlet) / q, (h.p_mouth_up - h.p_mouth_down) / q
        n, C = curvature(case)
        s = n / 2
        Q = h.R * U * A
        tau = (RHO * L / A) / ((h.p_mouth_up - h.p_mouth_down) / (s * Q))
        fc = 1 / (2 * np.pi * tau)
        a, b_lo, b_hi = C / 2, s, s * cp_up / dcp
        print(f"{case}: n {n:.2f}, C {C:.2f}, s {s:.2f}, Cp_up {cp_up:.3f}, dCp {dcp:.3f}, "
              f"tau {tau:.2f} s, fc {fc:.2f} Hz; A {a:.2f} sig^2, B {b_lo:.2f}-{b_hi:.2f} sig^2")
        for Lx in LX:
            phi = quad(lambda f: vk(f, Lx), 0, fc)[0]
            for sg in SIGMA:
                lo, hi = (100 * sg**2 * (phi * a + (1 - phi) * b) for b in (b_lo, b_hi))
                print(f"  Lx {Lx:.2f} m, phi {phi:.2f}, u'/U {sg:.2f}: A {100 * a * sg**2:.2f} %, "
                      f"B {100 * b_lo * sg**2:.2f}-{100 * b_hi * sg**2:.2f} %, blended {lo:.2f}-{hi:.2f} %")


if __name__ == "__main__":
    main()
