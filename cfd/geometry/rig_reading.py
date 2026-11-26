"""What Phillips et al.'s (2000) bubble method would have read: the calibration-rig CFD.

The rig runs (out/rig_v5: the Phillips sampler with the 7 mm x 600 mm glass tube sleeved over its
inlet, in the 0.3 m flume) give, for each ambient velocity U:

  R_rig      intake ratio with the glass tube fitted (bore mean velocity of the 4 mm inlet / U)
  issue B    R_rig / R_flume - 1, the cut in intake caused by the glass tube itself
  u_g        bulk velocity in the glass tube, u_in / (7/4)^2  (continuity)
  bubble     speed of a bubble timed over 500 mm of the glass tube, for three assumptions about
             where the bubble rides: with the bulk flow (k = 1), on the centreline (developing
             towards 2 u_g; from the symmetry-plane slice), or near the top of the bore
  reading    Phillips' conversion, u_in,Ph = 3.0625 x bubble speed, as R_Ph = u_in,Ph / U, against
             his Equation 1. The paper prints 1.75 (the diameter ratio), but J. M. Phillips told
             the author that the printed value was probably an error and that the area
             ratio (7/4)^2 was applied. So R_Ph = k R_rig, with k the bubble speed over the glass bulk
             speed. The printed 1.75 would need k = 2.35 and 3.12 to return Eq. 1 at 0.3 and 0.585 m/s,
             above the laminar centreline (1.97, 1.91); at 0.154 m/s it needs 1.96, just below the
             centreline (1.99). That check is printed at the end as supporting evidence.

The bubble speed over a timed window is the window length over the transit time, so the local
velocity is averaged harmonically along x. The three windows are the upstream, centred and
downstream 500 mm of the 579 mm passage from the glass mouth to the inlet-tube tip; downstream is
the headline.

Run from cfd/:  python geometry/rig_reading.py
Reads  results/rig_v5/.../U*_summary.json (or, with --partial, the last 25 % of the time series
       of a run still in progress) and results/mesh_slices/rig_U*_field[_ckpt].npz
Writes results/rig_reading.csv and ../paper/figures/data/fig_phillips_*.csv (drawn by
       paper/figures/matlab/make_fig_phillips.m).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

from tims_geometry import SPECS, x_stations

RIG = Path("results/rig_v5/phillips2000_rig_flume_s1.00_sst-gamma")
FLUME = Path("results/rerun_flume_v5/phillips2000_flume_s1.00_sst-gamma")
OPEN = Path("results/urans_v5/phillips2000_s1.00_sst-gamma")
SLICES = Path("results/mesh_slices")
DATA = Path("../paper/figures/data")
VELS = ["0.081", "0.154", "0.300", "0.585"]
SPEC = SPECS["phillips2000_rig"]
XS = x_stations(SPEC)
RG = SPEC.glass_id / 2
AREA_RATIO = (SPEC.glass_id / 0.004) ** 2          # 3.0625
PH_FACTOR = AREA_RATIO                              # Phillips' conversion: the area ratio (J. M. Phillips)
PRINTED = 1.75                                      # the factor printed in Phillips et al. (2000), ruled out
PASSAGE = XS["in_tip"] - XS["glass_mouth"]          # 0.579 m of 7 mm bore
WINDOWS = {"upstream": (0.0, 0.5), "centred": ((PASSAGE - 0.5) / 2, (PASSAGE + 0.5) / 2),
           "downstream": (PASSAGE - 0.5, PASSAGE)}  # distance from the glass mouth, m
# pipe-loss upper bounds (Table tab:labR), computed in geometry/pipe_bounds.py (hard-coded here before)
import pipe_bounds  # noqa: E402
PIPE_BOUND = pipe_bounds.OPEN_BOUND
RIG_BOUND = pipe_bounds.RIG_BOUND
# Phillips et al. (2000) Fig. 2: U, u_in in cm/s, digitised from the embedded 300 ppi scan by
# geometry/digitise_phillips_fig2.py (marker centres; reading error +-0.0004 log U, +-0.0010 log
# u_in; the six points refit to y = -2.179 + 2.070 x against the printed -2.182 + 2.074 x).
# Replaces earlier by-eye values (within 1.6 %).
FIG2_CSV = Path("results/phillips2000_fig2_digitised.csv")
FIG2 = [(float(r["U_cm_s"]), float(r["u_in_cm_s"])) for r in
        csv.DictReader(l for l in FIG2_CSV.read_text().splitlines() if not l.startswith("#"))]


def eq1(U):
    """Phillips' Equation 1 as an intake ratio: log u_in = -2.182 + 2.074 log U (cm/s)."""
    return 10 ** (-2.182 + 2.074 * np.log10(100 * U)) / (100 * U)


def run_R(folder, vel, partial=False):
    f = folder / f"U{vel}_summary.json"
    if f.exists():
        return json.loads(f.read_text())["mean"]["R"], "final"
    ts = folder / f"U{vel}_timeseries.csv"
    if partial and ts.exists():
        a = np.genfromtxt(ts, delimiter=",", names=True)
        n = max(1, len(a) // 4)
        return float(np.mean(a[a.dtype.names[1]][-n:])), f"partial (t = {a[a.dtype.names[0]][-1]:.1f} s)"
    return None, "missing"


def glass_profile(vel):
    """Bulk, centreline and near-top velocity along the glass bore, in 5 mm bins from the mouth."""
    f = SLICES / f"rig_U{vel}_field.npz"
    tag = "final"
    if not f.exists():
        f, tag = SLICES / f"rig_U{vel}_field_ckpt.npz", "checkpoint"
    if not f.exists():
        return None, "missing"
    z = np.load(f)
    off = z["off"]
    n = np.diff(off)
    x = np.add.reduceat(z["xz"][:, 0], off[:-1]) / n
    zz = np.add.reduceat(z["xz"][:, 1], off[:-1]) / n
    a, u = z["area"], z["u"]
    # the slice also holds the flume's free-slip floor and side wall (zones symmetry_flume_*),
    # which project onto the same (x, z); keep the y = 0 plane only
    names = [str(nm) for nm in z["names"]]
    plane = np.array(["flume" not in nm for nm in names])[z["zone"]]
    edges = np.arange(0.0, PASSAGE + 1e-9, 0.005)
    s = x - XS["glass_mouth"]
    rows = []
    for s0, s1 in zip(edges[:-1], edges[1:]):
        m = plane & (s >= s0) & (s < s1) & (np.abs(zz) < RG)
        if not m.any():
            continue
        w = np.abs(zz[m]) * a[m]                     # axisymmetric weight on the y = 0 plane
        c = m & (np.abs(zz) < 0.0004)
        if not c.any():
            c = m & (np.abs(zz) < 0.0008)
        t = m & (zz > 0.0025) & (zz < 0.0034)
        rows.append((0.5 * (s0 + s1), np.sum(u[m] * w) / np.sum(w),
                     np.sum(u[c] * a[c]) / np.sum(a[c]), np.sum(u[t] * a[t]) / np.sum(a[t])))
    return np.array(rows), tag


def window_speed(s, v, lo, hi):
    """Speed of something carried at local velocity v(s) over [lo, hi]: length / transit time."""
    m = (s >= lo) & (s <= hi)
    return (hi - lo) / np.sum(0.005 / v[m]) * (np.count_nonzero(m) * 0.005 / (hi - lo))


def main():
    partial = "--partial" in sys.argv
    out = []
    print(f"glass passage {1e3 * PASSAGE:.0f} mm; area ratio {AREA_RATIO:.4f}; Phillips' factor {PH_FACTOR}")
    for vel in VELS:
        U = float(vel)
        R, state = run_R(RIG, vel, partial)
        if R is None:
            print(f"U = {vel}: no rig result yet")
            continue
        Rf, _ = run_R(FLUME, vel)
        Ro, _ = run_R(OPEN, vel)
        ug = R * U / AREA_RATIO
        row = dict(U=U, state=state, R_rig=R, R_flume=Rf, R_open=Ro, R_eq1=eq1(U),
                   issueB=(R / Rf - 1) if Rf else None, u_glass=ug,
                   Re_glass=998.2 * ug * SPEC.glass_id / 1.003e-3,
                   R_Ph_bulk=PH_FACTOR * ug / U)
        prof, ptag = glass_profile(vel)
        row["slice"] = ptag
        if prof is not None:
            s, ub, uc, ut = prof.T
            row["bulk_check"] = float(np.mean(ub[(s > 0.05) & (s < PASSAGE - 0.02)]) / ug)
            for name, (lo, hi) in WINDOWS.items():
                kc = window_speed(s, uc, lo, hi) / ug
                kt = window_speed(s, ut, lo, hi) / ug
                row[f"k_centre_{name}"] = kc
                row[f"k_top_{name}"] = kt
                row[f"R_Ph_centre_{name}"] = PH_FACTOR * kc * ug / U
                row[f"R_Ph_top_{name}"] = PH_FACTOR * kt * ug / U
            row["t_transit_bulk_s"] = 0.5 / ug
            row["t_transit_centre_s"] = 0.5 / (row["k_centre_downstream"] * ug)
            np.savetxt(Path("results") / f"rig_profile_U{vel}.csv", np.column_stack([s, ub / ug, uc / ug, ut / ug]),
                       delimiter=",", header="s_from_mouth_m,bulk_over_ug,centre_over_ug,top_over_ug", comments="")
        out.append(row)
        print(f"\nU = {vel} m/s ({state}; slice {ptag})")
        print(f"  R_rig {R:.4f}   flume bare {Rf if Rf is None else round(Rf, 4)}   open {Ro if Ro is None else round(Ro, 4)}"
              f"   Eq. 1 {eq1(U):.4f}")
        if Rf:
            print(f"  issue B (glass tube): {100 * row['issueB']:+.1f} %")
        print(f"  glass bulk {100 * ug:.3f} cm/s, Re {row['Re_glass']:.0f}; 500 mm at bulk speed takes {0.5 / ug:.0f} s")
        print(f"  reading with a bulk-speed bubble (k = 1): R_Ph {row['R_Ph_bulk']:.4f} = {row['R_Ph_bulk'] / eq1(U):.2f} x Eq. 1")
        if prof is not None:
            print(f"  slice bulk / continuity bulk: {row['bulk_check']:.3f}")
            for name in WINDOWS:
                print(f"  {name:10s}: centreline k {row[f'k_centre_{name}']:.3f} -> R_Ph {row[f'R_Ph_centre_{name}']:.4f}"
                      f" ({row[f'R_Ph_centre_{name}'] / eq1(U):.2f} x Eq. 1);  near-top k {row[f'k_top_{name}']:.3f}"
                      f" -> R_Ph {row[f'R_Ph_top_{name}']:.4f}")
    if not out:
        return
    keys = sorted({k for r in out for k in r}, key=lambda k: (k != "U", k))
    with open("results/rig_reading.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    # slopes between adjacent velocities
    for key in ("R_rig", "R_Ph_centre_downstream"):
        v = [(r["U"], r[key]) for r in out if r.get(key)]
        if len(v) > 1:
            sl = [np.log(b[1] / a[1]) / np.log(b[0] / a[0]) for a, b in zip(v[:-1], v[1:])]
            print(f"\nlocal slopes d ln / d ln U of {key}: " + ", ".join(f"{x:.2f}" for x in sl))
    # k a bubble would need for the method to return Eq. 1 and each published point: R_Ph = k R_rig
    # with the area ratio, (1.75 / 3.0625) k R_rig with the printed factor. R / bound is the k a point
    # needs to be drawable at all.
    cal = [(r["U"], r["R_rig"], r.get("k_centre_downstream")) for r in out if r["U"] >= 0.15]
    if len(cal) > 1:
        def loglin(u, pts):
            lu, lr = np.log([q[0] for q in pts]), np.log([q[1] for q in pts])
            if u > pts[-1][0]:                           # 0.588 m/s: extend the last segment
                return float(np.exp(np.polyval(np.polyfit(lu[-2:], lr[-2:], 1), np.log(u))))
            return float(np.exp(np.interp(np.log(u), lu, lr)))
        print("\nk needed to return Eq. 1 (area ratio / printed 1.75), and the laminar centreline k:")
        for u, rr, kc in cal:
            print(f"  U {u:.3f}: {eq1(u) / rr:.2f} / {eq1(u) / rr * AREA_RATIO / PRINTED:.2f}   centreline {kc:.2f}")
        print("published points (Fig. 2): R, k needed (area ratio / printed 1.75), R over the rig bound")
        implied = []
        for Ucm, ucm in FIG2:
            u, R = Ucm / 100, ucm / Ucm
            k = R / loglin(u, [(c[0], c[1]) for c in cal])
            rb = R / pipe_bounds.bound(u, rig=True)   # exact at each point (was log-interpolated)
            implied.append((u, R, k, k * AREA_RATIO / PRINTED, rb))
            print(f"  U {u:.4f}  R {R:.4f}  k {k:.2f} / {k * AREA_RATIO / PRINTED:.2f}  R / bound {rb:.3f}")
        # beta implied by each published point: shown apart from
        # the literature prior band (1.25-1.4), never used to set it
        if not partial:
            with open(DATA / "fig_phillips_beta.csv", "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["U", "R_published", "beta_area_ratio", "beta_printed_1.75", "R_over_rig_bound"])
                w.writerows([f"{u:.4f}", f"{R:.4f}", f"{k:.3f}", f"{kp:.3f}", f"{rb:.3f}"] for u, R, k, kp, rb in implied)
    if partial:
        print("\n--partial: figure data not written")
        return
    DATA.mkdir(parents=True, exist_ok=True)
    with open(DATA / "fig_phillips_cfd.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["series", "U", "R"])
        for r in out:
            for series, key in (("open", "R_open"), ("flume", "R_flume"), ("rig", "R_rig"),
                                ("reading_bulk", "R_Ph_bulk"), ("reading_centre", "R_Ph_centre_downstream")):
                if r.get(key):
                    w.writerow([series, r["U"], f"{r[key]:.5f}"])
    with open(DATA / "fig_phillips_lab.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["series", "U", "R"])
        for Ucm, ucm in FIG2:
            w.writerow(["phillips_data", f"{Ucm / 100:.4f}", f"{ucm / Ucm:.4f}"])
        for U in np.linspace(0.154, 0.585, 40):
            w.writerow(["phillips_fit", f"{U:.4f}", f"{eq1(U):.5f}"])
        for U in np.linspace(0.15, 0.60, 40):             # Goharrokhi et al. (2019): u = 0.89 U - 1.3 cm/s
            w.writerow(["adv_fit", f"{U:.4f}", f"{0.89 - 1.3 / (100 * U):.5f}"])
        for U, R in PIPE_BOUND:
            w.writerow(["pipe_bound", U, f"{R:.3f}"])
    print("\nwrote results/rig_reading.csv and", DATA / "fig_phillips_{cfd,lab}.csv")


if __name__ == "__main__":
    main()
