"""Pipe-loss upper bounds of the Phillips sampler (Table 4(a), Fig. 5) and the flume drive rise
(Results, "Blockage in the 0.3 m flume raises the pressure difference between the mouths by 19-21 %").

hard-coded RIG_BOUND / PIPE_BOUND in rig_reading.py).

Bound: the largest R that a drive of DCP dynamic pressures can push through the two 4 mm x 150 mm
tubes, with the least loss laminar flow can have: fully developed (Hagen-Poiseuille) friction along
both tubes, the 0.5 velocity head of a sharp entry and the one velocity head the jet loses on leaving
each tube into a much larger space (the inlet jet into the body; the exit jet into the wake counted
once, so 1.5 heads in all is a lower bound on the minor losses). Real (developing) friction and the
expansion losses only add to this, so R cannot exceed it.
- open water: DCP = 1.2 (the mouth total pressure, about 1.0, less the wake pressure, about -0.2);
- rig in the flume: DCP = 1.2 x 1.21 (the flume raises the mouth pressure difference by up to
  21 %, below), plus Hagen-Poiseuille friction in the 579 mm x 7 mm glass bore at its own velocity.
Water at 20 C as in the CFD (mu 1.002e-3 Pa s, rho 998.2 kg/m3).

Flume drive rise: the static pressure difference between the inlet and outlet mouths, averaged over
each run's window, flume (rerun_flume_v5) over open domain (urans_v5), Phillips at 0.154-0.585 m/s.
Run from cfd/: python geometry/pipe_bounds.py   (writes results/checks/pipe_bounds.txt)
"""
import csv
import json
from pathlib import Path

from scipy.optimize import brentq

RHO, MU = 998.2, 1.002e-3
NU = MU / RHO
D, L = 0.004, 0.150                 # Phillips tubes
DG, LG = 0.007, 0.579               # calibration glass bore (SPECS phillips2000_rig, x_stations)
MINOR = 1.5                         # entry 0.5 + one jet velocity head
DCP_OPEN, FLUME = 1.2, 1.21


def k_hp(u, d, length):
    """Hagen-Poiseuille friction in velocity heads of u."""
    return 64 * NU / (u * d) * length / d


def bound(U, rig=False):
    dcp = DCP_OPEN * (FLUME if rig else 1.0)

    def res(R):
        u = R * U
        k = 2 * k_hp(u, D, L) + MINOR
        if rig:
            ug = u * (D / DG) ** 2
            k += k_hp(ug, DG, LG) * (ug / u) ** 2
        return k * R**2 - dcp
    return brentq(res, 1e-6, 5.0)


def window_mean(folder, run, U, keys):
    base = Path("results") / folder / f"{run}_s1.00_sst-gamma" / f"U{U:.3f}"
    a = json.load(open(f"{base}_summary.json"))["average_from"]
    rows = [r for r in csv.DictReader(open(f"{base}_timeseries.csv")) if float(r["t"]) >= a - 1e-9]
    return {k: sum(float(r[k]) for r in rows) / len(rows) for k in keys}


def drive_rise(U):
    keys = ("p_mouth_up", "p_mouth_down")
    o = window_mean("urans_v5", "phillips2000", U, keys)
    f = window_mean("rerun_flume_v5", "phillips2000_flume", U, keys)
    return (f["p_mouth_up"] - f["p_mouth_down"]) / (o["p_mouth_up"] - o["p_mouth_down"]) - 1


OPEN_BOUND = [(U, bound(U)) for U in (0.154, 0.3, 0.585)]
RIG_BOUND = [(U, bound(U, rig=True)) for U in (0.154, 0.3, 0.585)]

if __name__ == "__main__":
    lines = ["Pipe-loss upper bounds of the Phillips sampler (geometry/pipe_bounds.py; Table 4(a), Fig. 5)",
             "U_m_s,open_water,rig_in_flume"]
    lines += [f"{U:g},{o:.3f},{r:.3f}" for (U, o), (_, r) in zip(OPEN_BOUND, RIG_BOUND)]
    lines += ["", "published points (Phillips et al. 2000 Fig. 2, digitised) against the rig bound:"]
    for U, Rp in [(0.1544, 0.1103), (0.2132, 0.1973), (0.2917, 0.2220), (0.3557, 0.3714), (0.5631, 0.4435),
                  (0.5880, 0.5142)]:
        lines.append(f"U {U:.4f}: published {Rp:.4f}, rig bound {bound(U, rig=True):.3f}, ratio {Rp / bound(U, rig=True):.3f}")
    lines += ["", "Flume drive rise (static mouth-to-mouth pressure difference, flume / open - 1):"]
    lines += [f"U {U:g}: {100 * drive_rise(U):+.1f} %" for U in (0.154, 0.3, 0.585)]
    out = "\n".join(lines) + "\n"
    Path("results/checks/pipe_bounds.txt").write_text(out)
    print(out)
