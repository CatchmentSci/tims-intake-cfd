"""Wake and external-flow numbers of the paper that are not produced by another script:

1. the peak speed at the shoulder (Results: "accelerates to 1.25-1.38 U at the shoulder");
2. the rear-face (base) pressure coefficient of the Phillips sampler at 0.3 m/s (Discussion, against
   the -0.19 of Pavia et al. 2019);
3. the one-sided wake estimate: the rise in R if shedding lowered the outlet-mouth pressure by the gap
   between the measured and computed base pressures, from the one-dimensional model at each design's
   own CFD drive (geometry/screening_band.py);
4. the inflow-kick test (results/full_kick): side-force decay and R;
5. the first-cell y+ on the bed of the Goharrokhi et al. flume (Supplement S2, "10-50"), from the
   log law at the first cell of the symmetry-plane slices (results/mesh_slices/gflume_U*_final_field.npz);
6. the rise speed of a small air bubble in the glass tube (Supplement S4, "a 0.5-1 mm bubble rises at
   0.06-0.3 m/s"): terminal velocity of a bubble of that DIAMETER with a rigid (contaminated) interface,
   Schiller-Naumann drag, and a clean interface, Mei, Klausner and Lawrence (1994) drag.
Run from cfd/: python geometry/number_checks.py   (writes results/checks/number_checks.txt)
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, "geometry")
import fig_budget as fb          # noqa: E402  (slices, SPECS)
import screening_band as sb      # noqa: E402  (one-dimensional model)

RHO = 998.2
CP_BASE_MEASURED = -0.19         # Pavia et al. (2019), L/D 5, tripped
STEMS = {"Phillips": ("phillips2000", "phillips2000_v5_s1.00", "v5"), "Perks": ("perks2014", "perks2014_v5_s1.00", "v5"),
         "Elliott": ("elliott2017", "elliott2017_v5_s1.00", "v5"), "McDonald": ("mcdonald2010", "mcdonald2010_v5bore_s1.00", "v5bore")}
FOLDER = {"Phillips": "urans_v5", "Perks": "urans_v5", "Elliott": "urans_v5", "McDonald": "urans_v5bore"}


def p_inf(folder, run, U=0.3):
    base = Path("results") / folder / f"{run}_s1.00_sst-gamma" / f"U{U:.3f}"
    a = json.load(open(f"{base}_summary.json"))["average_from"]
    rows = [r for r in csv.DictReader(open(f"{base}_timeseries.csv")) if float(r["t"]) >= a - 1e-9]
    return float(np.mean([float(r["p_inlet"]) for r in rows]))


def shoulder_speed(label, U=0.3):
    """Largest |u|/U outside the body (radius beyond the body's outer radius), between the inlet tip
    and the rear face: the acceleration over the shoulder."""
    run, stem, ver = STEMS[label]
    m = fb.load_slice(stem)
    st = json.load(open(f"geometry/out/{ver}/scale_1.00/{run}_zones.json"))["stations_m"]
    ro = fb.SPECS[run].body_od / 2
    s = (np.abs(m["z"]) > ro + 5e-4) & (m["x"] > st["in_tip"]) & (m["x"] < st["rear"])
    speed = np.sqrt(2 * m["q"][s] / RHO) / U
    j = int(np.argmax(speed))
    return float(speed[j]), float(m["x"][s][j] - st["in_tip"])


def base_cp(label="Phillips", U=0.3):
    """Area-weighted pressure coefficient on the rear face, from the first cell layer behind it on the
    symmetry plane (each radial bin's face nearest the face), between the outlet tube and the body edge."""
    run, stem, ver = STEMS[label]
    m = fb.load_slice(stem)
    st = json.load(open(f"geometry/out/{ver}/scale_1.00/{run}_zones.json"))["stations_m"]
    sp = fb.SPECS[run]
    r0, r1 = sp.tube_od / 2, sp.body_od / 2
    s = (m["x"] > st["rear"]) & (m["x"] < st["rear"] + 0.003) & (np.abs(m["z"]) > r0) & (np.abs(m["z"]) < r1)
    z, x, p = np.abs(m["z"][s]), m["x"][s], m["p"][s]
    vals, rads = [], []
    for lo, hi in zip(np.linspace(r0, r1, 41)[:-1], np.linspace(r0, r1, 41)[1:]):
        k = (z >= lo) & (z < hi)
        if k.any():
            j = int(np.argmin(np.where(k, x, np.inf)))
            vals.append(p[j]), rads.append(0.5 * (lo + hi))
    vals, rads = np.array(vals), np.array(rads)
    q = 0.5 * RHO * U**2
    return float(((vals - p_inf(FOLDER[label], run, U)) * rads).sum() / rads.sum() / q)


def wake_rise(gap, velocities=(0.154, 0.3, 0.585)):
    """{design: [% rise in R per velocity]} when the outlet-mouth pressure falls by `gap` dynamic pressures."""
    out = {}
    for case, lab in sb.LABEL.items():
        d = sb.cfd_drive(lab)
        out[lab] = [100 * (sb.solve_R(U, case, d + gap) / sb.solve_R(U, case, d) - 1) for U in velocities]
    return out


def kick():
    base = Path("results/full_kick/phillips2000_full_s1.00_sst-gamma/U0.300")
    rows = list(csv.DictReader(open(f"{base}_timeseries.csv")))
    t = np.array([float(r["t"]) for r in rows])
    fy = np.abs(np.array([float(r["side_y"]) for r in rows]))
    R = np.array([float(r["R"]) for r in rows])
    jp = int(np.argmax(fy))
    after = t > t[jp]
    below = t[after & (fy < 1e-4)]
    r_ref = json.load(open("results/rerun_full_v5/phillips2000_full_s1.00_sst-gamma/U0.300_summary.json"))["mean"]["R"]
    period = 0.110 / (0.2 * 0.3)        # St_D 0.2, D 0.110 m, U 0.3 m/s
    t_dec = below[0] if below.size else float("nan")
    return dict(peak=fy[jp], t_peak=t[jp], t_1e4=t_dec, periods=(t_dec - 2.0) / period, R_end=R[-1], R_ref=r_ref,
                t_end=t[-1])


def flume_bed_yplus(U):
    """First-cell y+ on the flume bed upstream of the sampler (x from -1.0 to -0.2 m), log law
    u_p / u_tau = ln(E y+) / kappa (kappa 0.41, E 9.793) at each first-cell centroid."""
    d = np.load(f"results/mesh_slices/gflume_U{U:.3f}_final_field.npz")
    off, xz = d["off"], d["xz"]
    n = np.diff(off)
    x, z = np.add.reduceat(xz[:, 0], off[:-1]) / n, np.add.reduceat(xz[:, 1], off[:-1]) / n
    nu, zb = 1.002e-3 / RHO, -0.145                    # bed 0.145 m below the inlet axis
    s = (x > -1.0) & (x < -0.2)
    first = s & (z < z[s].min() + 2e-4)
    yp = []
    for y, u in zip(z[first] - zb, np.abs(d["u"][first])):
        ut = brentq(lambda t: u / t - np.log(9.793 * y * t / nu) / 0.41, 1e-6, 1.0)
        yp.append(y * ut / nu)
    return float(np.min(yp)), float(np.median(yp)), float(np.max(yp))


def bubble_rise(d, clean):
    """Terminal rise speed (m/s) of an air bubble of diameter d (m) in water at 20 C."""
    mu = 1.002e-3

    def f(v):
        Re = RHO * v * d / mu
        if clean:
            cd = 16 / Re * (1 + 1 / (8 / Re + 0.5 * (1 + 3.315 / np.sqrt(Re))))
        else:
            cd = 24 / Re * (1 + 0.15 * Re**0.687)
        return 0.75 * cd * v * v / d - 9.81
    return brentq(f, 1e-6, 2.0)


if __name__ == "__main__":
    out = ["Paper A number checks (geometry/number_checks.py; final numbers check, 26 Nov)", ""]
    out.append("1. Peak speed over the shoulder, |u|/U (distance from the inlet tip, m):")
    for lab in STEMS:
        v, x = shoulder_speed(lab)
        out.append(f"   {lab}: {v:.3f} at {x:.3f}")
    cpb = base_cp()
    gap = CP_BASE_MEASURED - cpb
    out += ["", f"2. Phillips rear-face Cp at 0.3 m/s (area-weighted, first cell layer): {cpb:.3f}",
            f"   gap to the measured {CP_BASE_MEASURED}: {-gap:.3f}", "",
            f"3. Rise in R if the outlet-mouth pressure fell by {-gap:.3f} (1-D model, own drive), % at 0.154 / 0.3 / 0.585 m/s:"]
    for lab, v in wake_rise(-gap).items():
        out.append(f"   {lab}: " + " / ".join(f"{x:.1f}" for x in v))
    out += ["   (for comparison, a gap of 0.07:)"]
    for lab, v in wake_rise(0.07).items():
        out.append(f"   {lab}: " + " / ".join(f"{x:.1f}" for x in v))
    k = kick()
    out += ["", "4. Inflow kick (5 deg for 2 s, full domain, Phillips 0.3 m/s):",
            f"   side force peak {k['peak']:.4f} N at t = {k['t_peak']:.1f} s; below 1e-4 N from t = {k['t_1e4']:.1f} s,"
            f" {k['t_1e4'] - 2.0:.1f} s after the kick ended ({k['periods']:.1f} shedding periods at St_D 0.2)",
            f"   R at t = {k['t_end']:.1f} s {k['R_end']:.6f} against the unperturbed full-domain {k['R_ref']:.6f}"
            f" ({100 * (k['R_end'] / k['R_ref'] - 1):+.4f} %)"]
    out += ["", "5. Flume bed first-cell y+ (log law; min / median / max), Goharrokhi flume:"]
    for U in (0.15, 0.3, 0.6):
        lo, med, hi = flume_bed_yplus(U)
        out.append(f"   U {U:g}: {lo:.0f} / {med:.0f} / {hi:.0f}")
    out += ["", "6. Air-bubble rise speed (diameter; rigid interface / clean interface), m/s:"]
    for dmm in (0.5, 1.0):
        slow = bubble_rise(dmm / 1000, False)
        out.append(f"   d {dmm} mm: {slow:.3f} / {bubble_rise(dmm / 1000, True):.3f};"
                   f" time to rise 3.5 mm at the slower: {0.0035 / slow:.3f} s")
    txt = "\n".join(out) + "\n"
    Path("results/checks/number_checks.txt").write_text(txt)
    print(txt)
