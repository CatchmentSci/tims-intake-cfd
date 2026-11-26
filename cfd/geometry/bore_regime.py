"""Flow regime in the inlet bores across the velocity sweep: does SST-gamma make the bore laminar
or turbulent where it should?

From symmetry-plane slices with mu_t, mu_lam and gamma (fluent/field_slice.py), for each design
and velocity:
  - intermittency gamma and nu_t/nu in the bore (axisymmetric average, faces weighted by |z|);
  - the static-pressure drop between two stations along the straight bore, in velocity heads,
    against the laminar developing-flow loss (Shah) and the turbulent Blasius loss over the
    same length (analysis/pipe_network_model.tube_loss, without its Re switch).
Run from cfd/: python geometry/bore_regime.py   (slices in results/mesh_slices/regime/)
       python geometry/bore_regime.py sst   (plain-SST runs at high velocity: slices
       results/mesh_slices/sst_<case>_U*_field.npz, R from results/sst_highU)
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from fig_budget import RHO, section  # noqa: E402

NU = 1.002e-3 / RHO
SST = len(sys.argv) > 1 and sys.argv[1] == "sst"
SL = Path("results/mesh_slices") if SST else Path("results/mesh_slices/regime")
PATTERN, RUNS = ("sst_{}_U*_field.npz", "sst_highU/{}_s1.00_sst") if SST else     ("regime_{}_U*_field.npz", "urans_v5/{}_s1.00_sst-gamma")
# case: (bore diameter, station 1 and 2 as distances from the inlet tip, m)
CASES = {"phillips2000": (0.004, 0.03, 0.14), "elliott2017": (0.004, 0.03, 0.13),
         "perks2014": (0.00794, 0.03, 0.10), "perks2014_flowmeter": (0.00794, 0.03, 0.10)}  # ends 13 mm before the 6.35 mm insert (noz_in)
# 2 mm bores: case -> (label, bore diameter,
# station key in <case>_zones.json, station 1 and 2 downstream of it, m). McDonald's 2 mm inlet runs
# 96 mm from the tip; Foets's 2 mm bore is the outlet tube, from out_inner (0.99 m) to its tip (1.14 m).
# Both production sweeps are on the v5bore meshes (results urans_v5bore and foets).
CASES2 = {"mcdonald2010": [("inlet 2 mm", 0.002, "in_tip", 0.02, 0.08)],
          "foets2020": [("inlet 4 mm", 0.004, "in_tip", 0.03, 0.14), ("outlet 2 mm", 0.002, "out_inner", 0.03, 0.14)]}
RUNS2 = {"mcdonald2010": "urans_v5bore/mcdonald2010_s1.00_sst-gamma", "foets2020": "foets/foets2020_s1.00_sst-gamma"}


def k_lam(u, D, L):
    """Laminar developing-flow loss coefficient from the entrance to L (Shah 1978)."""
    Re = u * D / NU
    xp = L / (D * Re)
    fre = 3.44 / math.sqrt(xp) + (1.25 / (4 * xp) + 16 - 3.44 / math.sqrt(xp)) / (1 + 2.1e-4 / xp**2)
    return 4 * fre / Re * L / D


def load(f):
    d = np.load(f)
    off, xz = d["off"], d["xz"]
    n = np.diff(off)
    m = dict(x=np.add.reduceat(xz[:, 0], off[:-1]) / n, z=np.add.reduceat(xz[:, 1], off[:-1]) / n,
             a=d["area"], u=d["u"], w=d["w"], p=d["p"],
             q=0.5 * RHO * (d["u"] ** 2 + d["v"] ** 2 + d["w"] ** 2))
    m["nut"] = d["mu_t"] / d["mu_lam"]
    m["gamma"] = d["gamma"] if "gamma" in d.files else np.full_like(m["nut"], np.nan)
    return m


def main():
    print("case, U, Re_b, gamma bore mean (core max), nu_t/nu bore mean (max), "
          "K CFD, K laminar, K turbulent, K CFD / K laminar")
    for case, (D, s1, s2) in CASES.items():
        st = json.loads(Path(f"geometry/out/v5/scale_1.00/{case}_zones.json").read_text())["stations_m"]
        tip = st["in_tip"]
        for f in sorted(SL.glob(PATTERN.format(case))):
            U = float(f.name.split("_U")[1][:5])
            R = json.loads(Path(f"results/{RUNS.format(case)}/U{U:.3f}_summary.json")
                           .read_text())["mean"]["R"]
            ub = R * U
            Re = ub * D / NU
            m = load(f)
            sel = (m["x"] > tip + s1) & (m["x"] < tip + s2) & (np.abs(m["z"]) < D / 2)
            w = np.abs(m["z"][sel]) * m["a"][sel]
            core = sel & (np.abs(m["z"]) < D / 4)
            g, nt = m["gamma"][sel], m["nut"][sel]
            # Shah's apparent friction is a static-pressure drop (wall friction plus the momentum
            # gained as the profile develops), so compare the area-averaged static pressure
            p1, p2 = section(m, tip + s1, D / 2)[1], section(m, tip + s2, D / 2)[1]
            K = (p1 - p2) / (0.5 * RHO * ub**2)
            Kl = k_lam(ub, D, s2) - k_lam(ub, D, s1)
            Kt = 0.316 * Re**-0.25 * (s2 - s1) / D
            print(f"{case}, {U:.3f}, {Re:5.0f}, {np.average(g, weights=w):.3f} "
                  f"({m['gamma'][core].max():.3f}), {np.average(nt, weights=w):7.2f} ({nt.max():7.2f}), "
                  f"{K:6.3f}, {Kl:6.3f}, {Kt:6.3f}, {K / Kl:5.2f}")


def main2():
    """The 2 mm bores (McDonald inlet, Foets outlet; Foets inlet for comparison). The bore velocity is
    the mean over the bore cross-section from R and the inlet area (both bores carry the same flow)."""
    print("\n2 mm bores: case, bore, U, Re_b, gamma bore mean (core max), nu_t/nu bore mean (max), "
          "K CFD, K laminar, K turbulent, K CFD / K laminar")
    for case, bores in CASES2.items():
        st = json.loads(Path(f"geometry/out/v5bore/scale_1.00/{case}_zones.json").read_text())["stations_m"]
        d_in = bores[0][1]
        for f in sorted(SL.glob(f"regime_{case}_U*_field.npz")):
            U = float(f.name.split("_U")[1][:5])
            R = json.loads(Path(f"results/{RUNS2[case]}/U{U:.3f}_summary.json").read_text())["mean"]["R"]
            m = load(f)
            for label, D, key, s1, s2 in bores:
                ub = R * U * (d_in / D) ** 2
                Re = ub * D / NU
                x0 = st[key]
                sel = (m["x"] > x0 + s1) & (m["x"] < x0 + s2) & (np.abs(m["z"]) < D / 2)
                w = np.abs(m["z"][sel]) * m["a"][sel]
                core = sel & (np.abs(m["z"]) < D / 4)
                g, nt = m["gamma"][sel], m["nut"][sel]
                p1, p2 = section(m, x0 + s1, D / 2)[1], section(m, x0 + s2, D / 2)[1]
                K = (p1 - p2) / (0.5 * RHO * ub**2)
                Kl = k_lam(ub, D, s2) - k_lam(ub, D, s1)
                Kt = 0.316 * Re**-0.25 * (s2 - s1) / D
                print(f"{case}, {label}, {U:.3f}, {Re:5.0f}, {np.average(g, weights=w):.3f} "
                      f"({m['gamma'][core].max():.3f}), {np.average(nt, weights=w):7.2f} ({nt.max():7.2f}), "
                      f"{K:6.3f}, {Kl:6.3f}, {Kt:6.3f}, {K / Kl:5.2f}")


if __name__ == "__main__":
    main()
    if not SST:
        main2()
