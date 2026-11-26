"""Grid convergence at 0.3 m/s for the paper (Table 2 / tab:verification, Table S3): one source for the
medium-mesh error e, the fine-mesh GCI and the extrapolated R, used by sensitivity_bands.py and
supplement_tables.py.
-1.4 % in the other, against -1.43 % from the data).

Method: Celik et al. (2008). h ~ N^(-1/3) with N the cell count; r21 = (N_fine / N_medium)^(1/3),
r32 = (N_medium / N_coarse)^(1/3); apparent order p by fixed-point iteration, capped at 2 (the
schemes' formal order); R_ext = (r21^p R_fine - R_medium) / (r21^p - 1); e = R_medium / R_ext - 1
(negative: the medium mesh reads low); GCI_fine = 1.25 |R_fine - R_medium| / R_fine / (r21^p - 1).
Phillips is not in the asymptotic range; its e is taken against the upper end of its physical bracket
for R_inf (0.1896-0.1915, best 0.1905, Hagen-Poiseuille bore check), not from the GCI.

Cell counts: results/checks/mesh_min_orthogonal_quality.csv (from the meshing logs). R: the
U0.300 summary.json of each mesh. Run from cfd/: python geometry/mesh_gci.py
"""
import csv
import json
import math
from pathlib import Path

RES = Path("results")
QUAL = RES / "checks" / "mesh_min_orthogonal_quality.csv"
# design -> (mesh name in QUAL, results folder, run, (coarse, medium, fine) scales)
SETS = {
    "Phillips": ("v5/phillips2000", "urans_v5", "phillips2000", ("1.65", "1.00", "0.64")),
    "Perks": ("v5/perks2014", "urans_v5", "perks2014", ("1.65", "1.00", "0.64")),
    "Perks + meter": ("v5/perks2014_flowmeter", "urans_v5", "perks2014_flowmeter", ("1.65", "1.00", "0.64")),
    "Elliott": ("v5/elliott2017", "urans_v5", "elliott2017", ("1.65", "1.00", "0.64")),
    "McDonald": ("v5bore/mcdonald2010", "urans_v5bore", "mcdonald2010", ("2.20", "1.00", "0.45")),
    "Foets": ("v5bore/foets2020", "foets", "foets2020", ("1.65", "1.00", "0.64")),
}
PHILLIPS_RINF = (0.1896, 0.1905, 0.1915)   # lower, best, upper (Table S3 note)


def cells():
    """{(mesh, scale): cells} from the meshing-log table."""
    return {(r["mesh"], f"{float(r['scale']):.2f}"): int(r["cells"]) for r in csv.DictReader(open(QUAL))}


def R_at(folder, run, scale, U=0.3):
    f = RES / folder / f"{run}_s{scale}_sst-gamma" / f"U{U:.3f}_summary.json"
    return json.load(open(f))["mean"]["R"]


def gci(name):
    """Grid-convergence numbers for one design (dict)."""
    mesh, folder, run, scales = SETS[name]
    N = cells()
    n3, n2, n1 = (N[(mesh, s)] for s in scales)
    f3, f2, f1 = (R_at(folder, run, s) for s in scales)
    r21, r32 = (n1 / n2) ** (1 / 3), (n2 / n3) ** (1 / 3)
    e21, e32 = f2 - f1, f3 - f2
    sgn = math.copysign(1.0, e32 / e21)
    p = 2.0
    for _ in range(500):
        p = abs(math.log(abs(e32 / e21)) + math.log((r21**p - sgn) / (r32**p - sgn))) / math.log(r21)
    pc = min(p, 2.0)
    r_ext = (r21**pc * f1 - f2) / (r21**pc - 1)
    out = dict(cells=(n3, n2, n1), R=(f3, f2, f1), r21=r21, r32=r32, p=p, p_used=pc,
               monotone=sgn > 0, R_ext=r_ext, e=f2 / r_ext - 1,
               gci_fine=1.25 * abs(f1 - f2) / f1 / (r21**pc - 1),
               ratio=(f2 / r_ext - 1) / (f2 / f1 - 1))       # e over the fine-medium difference
    if name == "Phillips":                                     # not asymptotic: physical bracket
        out.update(R_ext=PHILLIPS_RINF[1], e=f2 / PHILLIPS_RINF[2] - 1, gci_fine=None,
                   ratio=(f2 / PHILLIPS_RINF[2] - 1) / (f2 / f1 - 1))
    return out


if __name__ == "__main__":
    for name in SETS:
        g = gci(name)
        gf = "-" if g["gci_fine"] is None else f"{100 * g['gci_fine']:.2f} %"
        print(f"{name:14s} cells {'/'.join(f'{n / 1e6:.2f}' for n in g['cells'])} M  "
              f"R {'/'.join(f'{r:.5f}' for r in g['R'])}  r21 {g['r21']:.3f} r32 {g['r32']:.3f}  "
              f"p {g['p']:.2f} (used {g['p_used']:.0f})  R_ext {g['R_ext']:.4f}  e {100 * g['e']:+.2f} %  "
              f"GCI_fine {gf}  e / (medium/fine - 1) {g['ratio']:.2f}")
