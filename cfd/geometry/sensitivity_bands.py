"""Data for fig:bands: R(U) of each design with sensitivity bands.

Writes two CSVs to paper/figures/data/ (the MATLAB figure reads them):
- fig_bands_sources.csv (long): design, U, source, lo_pct, hi_pct, basis. One row per design,
  velocity and source; lo/hi are the % change in R the source allows (lo <= 0 <= hi).
- fig_bands.csv (wide): design, U, R, then R_lo / R_hi for the numerical band (mesh, time step,
  inner iterations, domain), the model band (closure) and the total, and the flume R for Phillips.

Sources and how each is carried to the designs and velocities where it was not run:
- mesh: one-sided at 0.3 m/s. The medium mesh reads low by e (tab:verification), so the true R
  lies between R and R / (1 + e). Velocity-specific for the
  two designs with a fine-mesh run at 1.5 m/s. There e(1.5) = (medium / fine - 1 at 1.5 m/s) x the
  design's own Richardson factor at 0.3 m/s (e / (medium / fine - 1)), keeping its sign:
  Elliott reads HIGH at 1.5 m/s (fine 0.58 % below medium), so its term flips to the low side.
  Where the 1.5 m/s fine-medium difference is no larger than that fine run's time-mean mouth
  imbalance (Perks: 0.14 % against 0.14 %), the sign is unresolved and the term is +-|e(1.5)|.
  The low and high ends are linear in U between 0.3 and 1.5 m/s, and flat below 0.3. The other
  designs have no high-U mesh test: their 0.3 m/s term is carried to every velocity and labelled so.
- closure: velocity-dependent. At 0.3 m/s, Phillips: plain SST +0.5 %,
  realizable k-epsilon -1.8 % against SST-gamma. At high U, where the bore Reynolds number
  passes 2000 and SST-gamma still keeps the bores nearly laminar, the alternative closures make the
  bores turbulent and lower R: for Perks at 1.0 m/s plain SST -3.55 % (sst_highU; extra
  loss in the outlet insert and tubing) and realizable k-epsilon -6.53 % (rke_highU; path
  loss dCp/R^2 +21 %, drive +5.7 %); for Elliott at 1.5 m/s plain SST -5.24 % (turbulent
  inlet bore, less recovery at the jet, less outlet suction). The low side at each anchor velocity
  is the most negative closure change measured there or at any lower anchor (a running minimum:
  the anchors are different designs, and the envelope does not narrow as the bore Reynolds number
  rises), linear in U between anchors, -1.8 % at and below 0.3 m/s, the last anchor's value above.
  Applied to every design as an envelope: conservative for the designs whose bores stay laminar
  (the 2 mm McDonald and Foets bores).
  The high side (forced-laminar bore runs, reading rule fixed before any result): the
  Perks and Perks + meter runs at 0.585, 1.0 and 1.5 m/s with every bore cell forced laminar
  (results/laminar_bore), dR_lam = R_laminar / R_SST-gamma - 1. At each anchor (0.3, 0.585, 1.0,
  1.5 m/s) the high side is the largest of +0.5 % (Phillips, plain SST, 0.3 m/s) and every dR_lam
  at that velocity or below (a running maximum), linear in U between anchors, flat outside them.
  A negative dR_lam enters the low side by the running minimum; it never narrows it.
  All six runs are final; the script stops if any summary.json is missing.
- time step: Elliott at 1.25 m/s, the least steady production run, at three time steps (ratio 2):
  2 dt -0.79 % (rerun_dt2x), dt/2 -0.22 % (rerun_dt_highU) against the production dt. The
  changes alternate in sign (convergence ratio -0.28), so there is no observed order; the term is
  half the range of the triplet, +-0.39 % (oscillatory convergence, Stern et al. 2001; Celik et al.
  2008), falling back to the dt/2 change if the 2 dt run is absent. Carried to every design and
  velocity as a conservative allocation (Phillips at 0.3 m/s, steady, gives -0.001 %). Inner iterations (Elliott fine mesh, 20 against 10,
  +0.007 % over 10-20 s, rr-inner20 with its continuation) and full domain (Phillips -0.19 %, from rerun_full_v5;
  Perks + meter pending, full_meter_v5): read from results/ where present, applied as +- to all.
- domain extent: one-sided. Phillips at 0.3 m/s with every boundary moved out (open_big: side and
  top/bottom at +-1.10 m, cross-section x4; upstream and downstream x1.5)
  reads 0.65 % below the production domain, so the production R is high.
  Taking the blockage error as proportional to frontal area over cross-section (as in the classical
  corrections), the change to an unbounded domain is 4/3 of that, -0.86 %. Carried to every design
  and velocity: only Elliott has a larger frontal area than Phillips (body o.d. 114.3 against 110 mm).
  Like the mesh error, a bias of known sign, so it is added linearly to the quadrature sum.
- flume: Phillips in the 0.3 m flume at 0.154, 0.3, 0.585 m/s. A configuration, not an
  uncertainty of the river case, so it is reported alongside and not added to the total.

Totals: the positive and negative parts of the independent sources are combined in quadrature
(numerical: mesh, dt, inner, domain; model: closure; total: both).

Run from cfd/:  python geometry/sensitivity_bands.py
"""
import csv
import json
import math
from pathlib import Path

RES = Path("results")
DATA = Path("../paper/figures/data")
SWEEP = [0.154, 0.3, 0.585, 1.0, 1.25, 1.5]
# label, results folder, run, medium-mesh error e at 0.3 m/s (tab:verification; medium reads low).
# e comes from geometry/mesh_gci.py (one source with Table S3).
# Elliott -0.014 against -0.0143 from the data). Phillips: against the upper end of R_inf (0.1896-0.1915).
import mesh_gci  # noqa: E402
DESIGNS = [(name, folder, run, mesh_gci.gci(name)["e"]) for name, folder, run in [
    ("Phillips", "urans_v5", "phillips2000"),
    ("Perks", "urans_v5", "perks2014"),
    ("Perks + meter", "urans_v5", "perks2014_flowmeter"),
    ("Elliott", "urans_v5", "elliott2017"),
    ("McDonald", "urans_v5bore", "mcdonald2010"),
    ("Foets", "foets", "foets2020"),
]]
CLOSURE = (-1.8, +0.5)          # % (k-epsilon, plain SST), Phillips 0.3 m/s
# high-U closure anchors: (U, design label, run) with plain SST (sst_highU) and k-epsilon (rke_highU)
CLOSURE_HIGHU = [(1.0, "Perks", "perks2014"), (1.5, "Elliott", "elliott2017")]
# forced-laminar bore runs (U, run)
# holding values were replaced by these finals (git history, 2405e66).
LAMINAR = [(U, r) for U in (0.585, 1.0, 1.5) for r in ("perks2014", "perks2014_flowmeter")]
DT_STEADY = -0.001              # %, Phillips 0.3 m/s dt/2 (rerun_dt0.01), steady: reported, not used
DT = None                       # %, set in main() from the Elliott 1.25 m/s time-step triplet
DT_NOTE = ""
INNER = +0.007                  # %, Elliott fine mesh, 20 vs 10 inner, R over 10-20 s (rr-inner20 + continuation)


def R_of(folder, run, U, scale="s1.00", turb="sst-gamma"):
    f = RES / folder / f"{run}_{scale}_{turb}" / f"U{U:.3f}_summary.json"
    return json.load(open(f))["mean"]["R"] if f.exists() else None


# fine-mesh runs at 1.5 m/s for the velocity-specific mesh term: design -> (run, fine folder at 1.5)
MESH_HIGHU = {"Perks": ("perks2014", "urans_v5_cont2"), "Elliott": ("elliott2017", "urans_v5_cont2")}


def mesh_term(name, run, e, U):
    """(lo %, hi %, note) of the mesh term at U."""
    h = lambda ee: 100 * (1 / (1 + ee) - 1)                   # noqa: E731  % change to the true R
    lo0, hi0 = (0.0, h(e)) if e < 0 else (h(e), 0.0)
    note0 = f"medium-mesh error {100 * e:+.2f} % at 0.3 m/s, one-sided"
    if name not in MESH_HIGHU:
        return lo0, hi0, note0 + "; carried to all U (no high-U mesh test)"
    r, fine_folder = MESH_HIGHU[name]
    d03 = R_of("urans_v5", r, 0.3) / R_of("urans_v5", r, 0.3, scale="s0.64") - 1
    fs = RES / fine_folder / f"{r}_s0.64_sst-gamma" / "U1.500_summary.json"
    mf = json.load(open(fs))["mean"]
    d15 = R_of("urans_v5", r, 1.5) / mf["R"] - 1
    imb = abs((mf["mdot_mouth_up"] + mf["mdot_mouth_down"]) / mf["mdot_mouth_up"])
    e15 = d15 * e / d03
    if abs(d15) <= imb:                                      # sign not resolved by the fine run
        lo1, hi1 = -abs(h(e15)), abs(h(e15))
        how = f"unresolved (|medium/fine - 1| {100 * abs(d15):.2f} % <= imbalance {100 * imb:.2f} %), +-"
    else:
        lo1, hi1 = (0.0, h(e15)) if e15 < 0 else (h(e15), 0.0)
        how = "reads " + ("low" if e15 < 0 else "high")
    if U <= 0.3:
        return lo0, hi0, note0 + "; at and below 0.3 m/s"
    f = min(1.0, (U - 0.3) / 1.2)
    note = (f"{note0}; at 1.5 m/s medium/fine - 1 = {100 * d15:+.2f} %, x Richardson factor "
            f"{e / d03:.2f} -> e {100 * e15:+.2f} % ({how}); lo and hi linear in U from 0.3 to 1.5 m/s")
    return lo0 + f * (lo1 - lo0), hi0 + f * (hi1 - hi0), note


def closure_anchors():
    """(U, low-side %, note) at 0.3 m/s and at each high-U anchor with a result."""
    pts = [(0.3, CLOSURE[0], "Phillips 0.3 m/s: k-epsilon -1.80 %, plain SST +0.50 %")]
    for U, label, run in CLOSURE_HIGHU:
        base = R_of("urans_v5", run, U)
        alts = {m: R_of(folder, run, U, turb=m) for m, folder in (("sst", "sst_highU"), ("rke", "rke_highU"))}
        pct = {m: 100 * (r / base - 1) for m, r in alts.items() if r}
        if pct:
            lo = min(min(pct.values()), pts[-1][1])       # running minimum over the anchors
            pts.append((U, lo, f"{label} {U:g} m/s: " + ", ".join(
                f"{'plain SST' if m == 'sst' else 'k-epsilon'} {v:+.2f} %" for m, v in pct.items())))
    return pts


def laminar_changes():
    """{(U, run): (dR_lam %, final?, note)} for the forced-laminar bore runs."""
    out = {}
    for U, run in LAMINAR:
        lam, base = R_of("laminar_bore", run, U), R_of("urans_v5", run, U)
        if lam and base:
            out[(U, run)] = (100 * (lam / base - 1), True, "final")
        else:
            raise SystemExit(f"forced-laminar result missing: {run} {U} m/s (results/laminar_bore or urans_v5)")
    return out


def closure_hi_anchors(lam):
    """(U, high-side %, note) at 0.3 m/s and each forced-laminar velocity: running maximum."""
    pts = [(0.3, CLOSURE[1], "Phillips 0.3 m/s: plain SST +0.50 %")]
    for U in sorted({u for u, _ in lam}):
        vals = {r: v for (u, r), v in lam.items() if u == U}
        hi = max([pts[-1][1]] + [v[0] for v in vals.values()])
        pts.append((U, hi, f"laminar bores {U:g} m/s: " + ", ".join(
            f"{'Perks + meter' if r.endswith('flowmeter') else 'Perks'} {v[0]:+.2f} %"
            + ("" if v[1] else " (provisional)") for r, v in sorted(vals.items()))))
    return pts


def closure_lo(U, pts):
    """Low side of the closure term at U: linear in U between anchors, flat outside them."""
    if U <= pts[0][0]:
        return pts[0][1]
    for (u0, l0, _), (u1, l1, _) in zip(pts[:-1], pts[1:]):
        if U <= u1:
            return l0 + (l1 - l0) * (U - u0) / (u1 - u0)
    return pts[-1][1]


def quad(vals):
    return math.sqrt(sum(v * v for v in vals))


def main():
    global DT, DT_NOTE
    dt1 = R_of("urans_v5", "elliott2017", 1.25)
    dt2 = R_of("rerun_dt_highU", "elliott2017", 1.25)
    dtx = R_of("rerun_dt2x", "elliott2017", 1.25)
    DT = 100 * (dt2 / dt1 - 1) if dt2 else -0.22
    DT_NOTE = f"Elliott 1.25 m/s (least steady run), dt/2 {DT:+.2f} %; carried conservatively"
    if dt2 and dtx:                  # triplet: half its range (oscillatory, no observed order)
        trip = [dtx, dt1, dt2]
        DT = 50 * (max(trip) - min(trip)) / dt1
        DT_NOTE = (f"Elliott 1.25 m/s (least steady run), 2dt {100 * (dtx / dt1 - 1):+.2f} %, dt/2 "
                   f"{100 * (dt2 / dt1 - 1):+.2f} %: oscillatory, half range {DT:.2f} %; carried conservatively")
    print("time step:", DT_NOTE)
    big = R_of("domain_big", "phillips2000_big", 0.3)
    c_ratio = (1.10 / 0.55) ** 2                     # cross-section, open_big over open
    dom_ext = 100 * (big / R_of("urans_v5", "phillips2000", 0.3) - 1) * c_ratio / (c_ratio - 1) if big else None
    print("domain extent:", f"{dom_ext:+.3f} % (to an unbounded domain)" if big else "not in yet")
    full_ph = R_of("rerun_full_v5", "phillips2000_full", 0.3)
    full_pct = 100 * (full_ph / R_of("urans_v5", "phillips2000", 0.3) - 1) if full_ph else None
    full_m = R_of("full_meter_v5", "perks2014_flowmeter_full", 0.3)

    cl_pts = closure_anchors()
    cl_note = "; ".join(n for _, _, n in cl_pts)
    lam = laminar_changes()
    hi_pts = closure_hi_anchors(lam)
    hi_note = "; ".join(n for _, _, n in hi_pts)
    prov = any(not v[1] for v in lam.values())
    lam_neg = sorted((u, v[0]) for (u, _), v in lam.items() if v[0] < 0)
    cl_low = lambda U: min([closure_lo(U, cl_pts)] + [d for u, d in lam_neg if u <= U + 1e-9])  # noqa: E731
    print("laminar bores (forced):", "; ".join(f"{r} {u:g} m/s {v[0]:+.3f} % ({v[2]})" for (u, r), v in sorted(lam.items())))
    print("closure high side anchors:", ", ".join(f"{u:g} m/s {h:+.2f} %" for u, h, _ in hi_pts),
          "PROVISIONAL" if prov else "")
    print("closure anchors:", ", ".join(f"{u:g} m/s {lo:+.2f} %" for u, lo, _ in cl_pts))
    src_rows, wide = [], []
    for name, folder, run, e in DESIGNS:
        for U in SWEEP:
            R = R_of(folder, run, U)
            if R is None:
                continue
            mesh_lo, mesh_hi, mesh_note = mesh_term(name, run, e, U)
            dom = full_pct if full_pct is not None else None
            if name == "Perks + meter" and full_m:
                dom = 100 * (full_m / R_of("urans_v5", "perks2014_flowmeter", 0.3) - 1)
            rows = [
                ("mesh", mesh_lo, mesh_hi, mesh_note),
                ("closure", cl_low(U), closure_lo(U, hi_pts),
                 ("PROVISIONAL (forced-laminar runs in progress): " if prov else "")
                 + f"low side linear in U between anchors ({cl_note}); high side running maximum, linear in U "
                 f"between anchors ({hi_note}); all designs"),
                ("time step", -abs(DT), abs(DT), DT_NOTE),
                ("inner iterations", -abs(INNER), abs(INNER), "Elliott fine 0.3 m/s, 20 vs 10 inner (10-20 s); carried"),
            ]
            if dom_ext is not None:
                rows.append(("domain extent", dom_ext, 0.0,
                             f"Phillips 0.3 m/s, boundaries moved out {100 * (big / R_of('urans_v5', 'phillips2000', 0.3) - 1):+.2f} %, "
                             "extrapolated to an unbounded domain; one-sided; carried"))
            if dom is not None:
                rows.append(("no symmetry plane", -abs(dom), abs(dom),
                             ("own design" if name in ("Phillips", "Perks + meter") and
                              (name == "Phillips" or full_m) else "Phillips") + " 0.3 m/s; carried"))
            for s, lo, hi, basis in rows:
                src_rows.append(dict(design=name, U=U, source=s, lo_pct=f"{lo:.3f}", hi_pct=f"{hi:.3f}", basis=basis))
            num = [r for r in rows if r[0] in ("mesh", "time step", "inner iterations", "no symmetry plane", "domain extent")]
            mod = [r for r in rows if r[0] == "closure"]
            # The mesh error is a one-sided bias of known sign, so it is added linearly to the
            # quadrature sum of the other sources (int M2)
            # The domain-extent error is one-sided of known sign too, and is added the same way
            ONE = ("mesh", "domain extent")
            band = lambda rs: (-quad([r[1] for r in rs if r[0] not in ONE]) + sum(r[1] for r in rs if r[0] in ONE),
                               quad([r[2] for r in rs if r[0] not in ONE]) + sum(r[2] for r in rs if r[0] in ONE))
            w = dict(design=name, U=U, R=f"{R:.6f}")
            tags = (("num", num), ("model", mod), ("total", num + mod))
            for tag, rs in tags:
                if rs:
                    lo, hi = band(rs)
                    w[f"R_lo_{tag}"], w[f"R_hi_{tag}"] = f"{R * (1 + lo / 100):.6f}", f"{R * (1 + hi / 100):.6f}"
                else:
                    w[f"R_lo_{tag}"] = w[f"R_hi_{tag}"] = ""
            Rf = R_of("rerun_flume_v5", "phillips2000_flume", U) if name == "Phillips" else None
            w["R_flume"] = f"{Rf:.6f}" if Rf else ""
            wide.append(w)
            if Rf:
                src_rows.append(dict(design=name, U=U, source="flume (configuration)", lo_pct="0.000",
                                     hi_pct=f"{100 * (Rf / R - 1):.3f}",
                                     basis="0.3 m flume, own run; not an uncertainty, not in the total"))
    with (DATA / "fig_bands_sources.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(src_rows[0]))
        wr.writeheader(); wr.writerows(src_rows)
    with (DATA / "fig_bands.csv").open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(wide[0]))
        wr.writeheader(); wr.writerows(wide)
    print(f"wrote {DATA / 'fig_bands.csv'} ({len(wide)} rows) and fig_bands_sources.csv ({len(src_rows)} rows)")
    print(f"full domain: Phillips {full_pct:+.3f} %" if full_pct is not None else "full domain: none",
          "| meter full domain:", "in" if full_m else "pending")
    for Up in (1.0, 1.25, 1.5):
        ws = [w for w in wide if w["U"] == Up]
        lo = [100 * (float(w["R_lo_total"]) / float(w["R"]) - 1) for w in ws]
        hi = [100 * (float(w["R_hi_total"]) / float(w["R"]) - 1) for w in ws]
        print(f"  {Up} m/s: total low {min(lo):+.2f} to {max(lo):+.2f} %, high {min(hi):+.2f} to {max(hi):+.2f} %  "
              + ", ".join(f"{w['design']} {l:+.2f}/{h:+.2f}" for w, l, h in zip(ws, lo, hi)))
    for w in wide:
        if w["U"] == 0.3:
            R = float(w["R"])
            print(f"  {w['design']:14s} R {R:.4f}  total {100 * (float(w['R_lo_total']) / R - 1):+.2f} / "
                  f"{100 * (float(w['R_hi_total']) / R - 1):+.2f} %   flume {w['R_flume'] or '-'}")


if __name__ == "__main__":
    main()
