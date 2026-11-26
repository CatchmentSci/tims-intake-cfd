"""Pipe-network screening band for Fig. 3a: R(U) from a 1-D series circuit per design.

Each design is a chain of bores between the nose (pressure coefficient Cp_in) and the outlet
(Cp_out), driven by dCp = Cp_in - Cp_out. Losses are in velocity heads of the local bore:
laminar developing-flow friction (Shah 1978) or Blasius, entry and contraction losses, the
Borda-Carnot loss of each sudden expansion, and the exit kinetic energy (alpha, from 1 to 2 as
the laminar profile develops). Bore dimensions are those of the CFD geometry
(tims_geometry.SPECS), so the band and the CFD share every assumption except the 1-D treatment.
The meter variant has no band: its internal meter is not a 1-D circuit.

Drive and friction law:
- Drive. The band spans dCp = 0.89-1.08, the range of the CFD's own drive over the six designs at
  0.3 m/s: 1 - approach - base from the pressure budget (fig_budget.csv). The earlier 0.98-1.2
  (mouth stagnation to an assumed base pressure of -0.2) sat above every CFD drive. Each design's
  own drive is cfd_drive(name); the per-design model error at that drive is printed and written to
  results/screening_model_errors.csv.
- Friction. Above Re 2300 the bore loss is the larger of developing-laminar (Shah) and Blasius,
  compared on friction + exit energy, so the loss is continuous in Re. The former hard switch to
  Blasius at Re 2300 lowered the loss of these short bores (L/d 25-45) and made the band jump
  upwards between 1.25 and 1.5 m/s.

Run from cfd/:  python geometry/screening_band.py   (prints R at the sweep velocities and the
per-design error against the CFD)
"""
import csv
import math
import sys
from pathlib import Path

from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pipe_network_model as _pnm  # noqa: E402

DATA = Path(__file__).resolve().parents[2] / "paper" / "figures" / "data"
DCP = (0.89, 1.08)       # range of the CFD drive, 1 - approach - base (fig_budget.csv), 0.3 m/s


def _laminar(u, D, L):
    """Developing-laminar loss (Shah 1978) and exit energy, in velocity heads.
    Shah's apparent friction gives the static-pressure drop from a uniform entry profile, which
    includes the rise in kinetic energy as the profile develops (alpha from 1 to 2; with the excess
    wall shear it makes up K(inf) = 1.25). That rise leaves the bore as the exit energy alpha,
    counted at the next expansion, so it is taken off here and the loss is the dissipation only:
    4 f_app L/D - (alpha - 1)."""
    Re = max(u * D / _pnm.NU, 1e-9)
    xp = L / (D * Re)
    fapp_re = (3.44 / math.sqrt(xp)
               + (1.25 / (4 * xp) + 16 - 3.44 / math.sqrt(xp)) / (1 + 2.1e-4 / xp**2))
    alpha = 1 + (1 - math.exp(-xp / 0.03))
    return 4 * fapp_re / Re * L / D - (alpha - 1), alpha


def tube_loss(u, D, L):
    """Bore loss continuous in Re: developing laminar below Re 2300; above it, the larger of that
    and Blasius (alpha 1.05), compared on dissipation + exit energy."""
    f_l, a_l = _laminar(u, D, L)
    Re = u * D / _pnm.NU
    if Re < 2300:
        return f_l, a_l
    f_b, a_b = 0.316 * Re**-0.25 * L / D, 1.05
    return (f_l, a_l) if f_l + a_l >= f_b + a_b else (f_b, a_b)


def cfd_drive(label):
    """The CFD drive 1 - approach - base of design `label` (fig_budget.csv, 0.3 m/s)."""
    for r in csv.DictReader(open(DATA / "fig_budget.csv")):
        if r["design"] == label:
            return 1 - float(r["approach"]) - float(r["base"])
    raise KeyError(label)


LABEL = {"phillips2000": "Phillips", "foets2020": "Foets", "perks2014": "Perks",
         "elliott2017": "Elliott", "mcdonald2010": "McDonald"}
K_ENTRY, K_REENTRANT, K_ELBOW = 0.5, 0.8, 1.0

# Elements: ("entry", D, K) minor loss at D; ("tube", D, L) friction, remembered for the next
# exit; ("expand", D1, D2) Borda-Carnot loss from D1 into D2 (D2 = None: into the chamber or the
# free stream, i.e. the exit kinetic energy alpha of the last tube); ("contract", D1, D2).
PATHS = {
    "phillips2000": [("entry", 0.004, K_ENTRY), ("tube", 0.004, 0.150), ("expand", 0.004, None),
                     ("entry", 0.004, K_REENTRANT), ("tube", 0.004, 0.150), ("expand", 0.004, None)],
    "foets2020":    [("entry", 0.004, K_ENTRY), ("tube", 0.004, 0.150), ("expand", 0.004, None),
                     ("entry", 0.002, K_REENTRANT), ("tube", 0.002, 0.150), ("expand", 0.002, None)],
    # tubing 7.94 mm, then the 6.35 mm panel insert (48.5 mm) through each cap
    "perks2014":    [("entry", 0.00794, K_ENTRY), ("tube", 0.00794, 0.1615 - 0.0485),
                     ("contract", 0.00794, 0.00635), ("tube", 0.00635, 0.0485), ("expand", 0.00635, None),
                     ("entry", 0.00635, K_REENTRANT), ("tube", 0.00635, 0.0485),
                     ("expand", 0.00635, 0.00794), ("tube", 0.00794, 0.1615 - 0.0485),
                     ("expand", 0.00794, None)],
    # inlet 4 mm x 178 mm; outlet through the 6.35 mm fitting (42 mm), an elbow, ~31 mm of 4 mm tube
    "elliott2017":  [("entry", 0.004, K_ENTRY), ("tube", 0.004, 0.178418), ("expand", 0.004, None),
                     ("entry", 0.00635, K_REENTRANT), ("tube", 0.00635, 0.042418),
                     ("entry", 0.00635, K_ELBOW), ("contract", 0.00635, 0.004),
                     ("tube", 0.004, 0.031), ("expand", 0.004, None)],
    # inlet 2 mm x 96 mm behind the open funnel; outlet a 3 mm hole through the 10 mm plug
    "mcdonald2010": [("entry", 0.002, K_ENTRY), ("tube", 0.002, 0.096162), ("expand", 0.002, None),
                     ("entry", 0.003, K_ENTRY), ("tube", 0.003, 0.010), ("expand", 0.003, None)],
}
D_IN = {k: v[0][1] for k, v in PATHS.items()}


def head_ratio(R, U, path):
    """Total loss in inlet velocity heads, sum K_i (u_i/u_in)^2, at inlet velocity R*U."""
    u_in, d_in = R * U, D_IN[next(iter(k for k, v in PATHS.items() if v is path))]
    total, alpha = 0.0, 1.0
    for kind, D, x in path:
        s = (d_in / D) ** 2                       # u_i / u_in
        u = u_in * s
        if kind == "entry":
            total += x * s**2
        elif kind == "tube":
            f, alpha = tube_loss(u, D, x)
            total += f * s**2
        elif kind == "contract":                  # loss on the downstream (smaller) velocity
            s2 = (d_in / x) ** 2
            total += 0.5 * (1 - (x / D) ** 2) * s2**2
        elif kind == "expand":
            if x is None:
                total += alpha * s**2
            else:
                total += (1 - (D / x) ** 2) ** 2 * s**2
            alpha = 1.0
    return total


def solve_R(U, case, dcp):
    path = PATHS[case]
    return brentq(lambda R: R**2 * head_ratio(R, U, path) - dcp, 1e-7, 3.0)


def band(case, velocities):
    """(R at dCp low, R at dCp high) for each velocity."""
    return [(solve_R(U, case, DCP[0]), solve_R(U, case, DCP[1])) for U in velocities]


if __name__ == "__main__":
    US = (0.154, 0.3, 0.585, 1.0, 1.25, 1.5)
    for case in PATHS:
        print(case, "  ".join(f"{U}: {lo:.4f}-{hi:.4f}" for U, (lo, hi) in zip(US, band(case, US))))
    # Model at each design's own CFD drive against the CFD sweep (fig_RU.csv)
    cfd = {(r["design"], float(r["U"])): float(r["R"]) for r in csv.DictReader(open(DATA / "fig_RU.csv"))}
    rows = []
    for case, lab in LABEL.items():
        d = cfd_drive(lab)
        errs = []
        for U in US:
            rm, rc = solve_R(U, case, d), cfd.get((lab, U))
            if rc:
                errs.append(100 * (rm / rc - 1))
                rows.append((lab, U, f"{d:.3f}", f"{rm:.5f}", f"{rc:.5f}", f"{errs[-1]:+.1f}"))
        print(f"{lab:9s} drive {d:.3f}: model/CFD - 1 = " + " ".join(f"{e:+.0f}" for e in errs) + " %")
    out = Path("results") / "screening_model_errors.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(("design", "U", "dCp_cfd", "R_model", "R_cfd", "error_pct"))
        w.writerows(rows)
    print("wrote", out)
    # Fig. 3a: each design's model at its own CFD drive, as a line over the sweep (40 velocities)
    grid = [0.154 * (1.5 / 0.154) ** (i / 39) for i in range(40)]
    own = DATA / "fig_RU_model_own.csv"
    with own.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(("design", "U", "dCp", "R_model"))
        for case, lab in LABEL.items():
            d = cfd_drive(lab)
            for U in grid:
                w.writerow((lab, f"{U:.4f}", f"{d:.3f}", f"{solve_R(U, case, d):.5f}"))
    print("wrote", own)
