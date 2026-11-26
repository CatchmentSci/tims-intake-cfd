"""Geometry sensitivity of the intake flow, by the 1-D pipe-network model
For each design, the bore-path dimensions that are assumed, derived or
subject to manufacturing tolerance are set to the ends of a plausible range (Table S1), and the
change in intake flow q is computed at the design's own CFD drive (screening_band.cfd_drive).

The 1-D model reads every design 1-19 % low in absolute terms (Section S3), but a change in one
dimension moves both runs alike, so the relative change is the useful output. Dimensions that act
only on the external flow (body and tube outside diameters, nose and cap shapes, the McDonald
funnel) are outside the model and are listed as such.

Run from cfd/:  python geometry/geometry_sensitivity_1d.py  -> results/checks/geometry_sensitivity_1d.txt
"""
import copy
import math
import sys
from pathlib import Path

from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import screening_band as sb  # noqa: E402

US = (0.154, 0.3, 1.0, 1.5)


def head_ratio(R, U, path):
    """sb.head_ratio with the inlet diameter taken from the path itself."""
    d_in = path[0][1]
    u_in = R * U
    total, alpha = 0.0, 1.0
    for kind, D, x in path:
        s = (d_in / D) ** 2
        u = u_in * s
        if kind == "entry":
            total += x * s**2
        elif kind == "tube":
            f, alpha = sb.tube_loss(u, D, x)
            total += f * s**2
        elif kind == "contract":
            s2 = (d_in / x) ** 2
            total += 0.5 * (1 - (x / D) ** 2) * s2**2
        elif kind == "expand":
            total += alpha * s**2 if x is None else (1 - (D / x) ** 2) ** 2 * s**2
            alpha = 1.0
    return total


def q_of(path, U, dcp):
    R = brentq(lambda R: R**2 * head_ratio(R, U, path) - dcp, 1e-7, 3.0)
    return R * U * math.pi * path[0][1] ** 2 / 4


def edit(path, changes):
    """changes: {element index: (D, x)} with None keeping the value."""
    p = copy.deepcopy(path)
    for i, (D, x) in changes.items():
        k, D0, x0 = p[i]
        p[i] = (k, D0 if D is None else D, x0 if x is None else x)
    return p


def scale_d(path, d_old, d_new):
    """Every element at diameter d_old (and contractions/expansions to it) moved to d_new."""
    out = []
    for k, D, x in path:
        D2 = d_new if abs(D - d_old) < 1e-9 else D
        x2 = d_new if (isinstance(x, float) and k in ("contract", "expand") and abs(x - d_old) < 1e-9) else x
        out.append((k, D2, x2))
    return out


P = sb.PATHS
# (design, item, range text, low path, high path). Element indices refer to sb.PATHS.
CASES = [
    ("phillips2000", "bore diameter, tolerance", "4.0 +- 0.1 mm",
     scale_d(P["phillips2000"], 0.004, 0.0039), scale_d(P["phillips2000"], 0.004, 0.0041)),
    ("foets2020", "outlet bore diameter, tolerance", "2.0 +- 0.1 mm",
     scale_d(P["foets2020"], 0.002, 0.0019), scale_d(P["foets2020"], 0.002, 0.0021)),
    ("foets2020", "outlet tube length (assumed as Phillips)", "150 +- 20 mm",
     edit(P["foets2020"], {4: (None, 0.130)}), edit(P["foets2020"], {4: (None, 0.170)})),
    ("perks2014", "tubing length beyond the funnels (author, ~20 mm)", "+-10 mm each end",
     edit(P["perks2014"], {1: (None, 0.1130 - 0.010), 8: (None, 0.1130 - 0.010)}),
     edit(P["perks2014"], {1: (None, 0.1130 + 0.010), 8: (None, 0.1130 + 0.010)})),
    ("perks2014", "insert bore diameter, tolerance", "6.35 +- 0.15 mm",
     scale_d(P["perks2014"], 0.00635, 0.0062), scale_d(P["perks2014"], 0.00635, 0.0065)),
    ("elliott2017", "bore diameter (3/8 in nylon bored to 4 mm), tolerance", "4.0 +- 0.1 mm",
     scale_d(P["elliott2017"], 0.004, 0.0039), scale_d(P["elliott2017"], 0.004, 0.0041)),
    ("elliott2017", "exposed outlet leg (assumed)", "31 (modelled) to 46 mm",
     P["elliott2017"], edit(P["elliott2017"], {7: (None, 0.046)})),
    ("elliott2017", "exposed outlet leg, text's 150 mm tube", "31 to 150 mm",
     P["elliott2017"], edit(P["elliott2017"], {7: (None, 0.150)})),
    ("elliott2017", "elbow loss coefficient (model input)", "K 0.75-1.5",
     edit(P["elliott2017"], {5: (None, 0.75)}), edit(P["elliott2017"], {5: (None, 1.5)})),
    ("mcdonald2010", "inlet bore diameter, tolerance", "2.0 +- 0.1 mm",
     scale_d(P["mcdonald2010"], 0.002, 0.0019), scale_d(P["mcdonald2010"], 0.002, 0.0021)),
    ("mcdonald2010", "outlet hole length (plug thickness, assumed)", "5-15 mm",
     edit(P["mcdonald2010"], {4: (None, 0.005)}), edit(P["mcdonald2010"], {4: (None, 0.015)})),
    ("mcdonald2010", "outlet hole diameter", "3.0 +- 0.2 mm",
     scale_d(P["mcdonald2010"], 0.003, 0.0028), scale_d(P["mcdonald2010"], 0.003, 0.0032)),
]
OUTSIDE = [
    "Phillips body o.d. (100-114 mm), end caps (5-15 mm, flush or recessed), tube o.d. (6-8 mm), nose cone (90-110 mm): external flow only",
    "Perks end caps: chamber and external flow only",
    "McDonald body wall (2.5-3.5 mm) and funnel (45-75 deg, spout 15-25 mm): external flow only; outlet position: CFD run E, +2.6 %",
    "Foets outlet tube o.d. (3-6 mm): external flow at the outlet tip only",
]


def main():
    lines = ["Change in intake flow q (%) at each design's own CFD drive, 1-D model; low end / high end of each range",
             "design         item                                                    range                     " + "  ".join(f"U={U:<5g}" for U in US)]
    for case, item, rng, lo, hi in CASES:
        dcp = sb.cfd_drive(sb.LABEL[case])
        cells = []
        for U in US:
            q0 = q_of(P[case], U, dcp)
            a, b = (100 * (q_of(p, U, dcp) / q0 - 1) for p in (lo, hi))
            cells.append(f"{a:+5.1f}/{b:+5.1f}")
        lines.append(f"{sb.LABEL[case]:14s} {item[:55]:55s} {rng[:25]:25s} " + "  ".join(cells))
    lines.append("Outside the 1-D model:")
    lines += ["  " + s for s in OUTSIDE]
    out = "\n".join(lines)
    print(out)
    Path("results/checks/geometry_sensitivity_1d.txt").write_text(out + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
