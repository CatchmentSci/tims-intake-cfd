"""Goharrokhi et al. (2019) flume runs (0.15, 0.30 and 0.60 m/s): the CFD at their ADV cells against
the digitised Fig. 3 points, scored by the acceptance criterion set before the results:
F = measured / CFD, within 1.00 +- 0.05.

Both sides are normalised by the flume section-mean velocity U: the measured u / (series cm/s), the
CFD u / U with U the bulk velocity of the precursor inflow.

For each run it prints:
- the bore R (late mean of the time series) and the cell table of goharrokhi_adv.report (axis-cell
  and all-21-cell means over U, point and 7 mm disc);
- F on the axis cells (digitised points with -5 <= z <= 2.5 mm, as plot_adv.py), per point and as the
  ratio of the means;
- F over every digitised point of that series (mean +- s.d.), and by group (above, on, below the
  axis) for the shape test;
- the CFD axial dip (1 - axis / mean of +-7 mm), against the measured dip of the same series
  (measured_dip(): 7.1-14.7 % over the five series).
  not reproducible, final numbers check).
The CFD value at a digitised point is taken at x = tip - 5 mm, as the plane value at the point
(the definition of the open-domain k = 1.165) and as the 1 mm cell / 7 mm disc mean.

Run from cfd/:
  python geometry/goharrokhi_flume_adv.py --U 0.3 \
      --slice results/mesh_slices/gflume_U0.300_field.npz \
      --timeseries results/goharrokhi_flume_v5/goharrokhi2019_gflume_s1.00_sst-gamma/U0.300_timeseries.csv
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "geometry")
from goharrokhi_adv import TIP, cell_mean, interpolator, late_R, plane, report  # noqa: E402

DIG = Path("results/goharrokhi2019_fig3_digitised.csv")
X0 = TIP - 0.005
OPEN_K = 1.165          # open-domain measured/CFD, all 56 points (plot_adv.py, v4 slices)


def measured_dip(cm):
    """Measured axial dip of one digitised series (Goharrokhi et al. Fig. 3), 1 - axis / mean of the
    stacks above and below. The stacks are the digitised points of each ADV position: axis
    -5 <= z <= 2.5 mm (as plot_adv.py), above 2.5 < z <= 9.5 mm, below -12 <= z < -5 mm, z from the
    inlet axis as digitised (the axis stack plots at -4.4 to +1.8 mm, so the windows follow it)."""
    pts = [(float(r["z_mm_from_axis"]), float(r["u_cm_s"])) for r in csv.DictReader(open(DIG))
           if int(r["series_cm_s"]) == cm]
    z, u = np.array([a for a, _ in pts]), np.array([b for _, b in pts])
    ax, ab, be = (z >= -5) & (z <= 2.5), (z > 2.5) & (z <= 9.5), (z < -5) & (z >= -12)
    return 1 - u[ax].mean() / np.mean([u[ab].mean(), u[be].mean()])


def score(U, slice_path, ts, disc):
    res = report(slice_path, ts, U, disc)
    f = interpolator(slice_path, U)[0]
    cm = int(round(U * 100))
    pts = [dict(u=float(r["u_cm_s"]) / cm, z=float(r["z_mm_from_axis"]))
           for r in csv.DictReader(open(DIG)) if int(r["series_cm_s"]) == cm]
    z = np.array([p["z"] for p in pts]) * 1e-3
    meas = np.array([p["u"] for p in pts])
    cfd_pt = plane(f, X0, z)[:, 0]
    cfd_cell = np.array([cell_mean(f, X0, zz, 0) for zz in z])
    cfd_disc = np.array([cell_mean(f, X0, zz, disc) for zz in z])
    ax = (z >= -5e-3) & (z <= 2.5e-3)
    groups = {"above (z > 2.5 mm)": z > 2.5e-3, "axis (-5 to 2.5 mm)": ax, "below (z < -5 mm)": z < -5e-3}

    print(f"\n   ---- F = measured / CFD, series {cm} cm/s, {len(pts)} digitised points "
          f"({ax.sum()} on the axis cells)")
    out = dict(U=U, R=res["R"], n=len(pts), n_axis=int(ax.sum()), axis_point=res["axis_point"],
               axis_disc=res["axis_disc"], all_point=res["all_point"], all_disc=res["all_disc"])
    for tag, c in (("point", cfd_pt), ("cell", cfd_cell), ("disc", cfd_disc)):
        r = meas / c
        fa = meas[ax].mean() / c[ax].mean()
        out[f"F_axis_{tag}"], out[f"F_all_{tag}"], out[f"sd_all_{tag}"] = fa, r.mean(), r.std()
        print(f"   CFD as {tag:5s}: axis cells meas {meas[ax].mean():.3f} / CFD {c[ax].mean():.3f}"
              f" = F {fa:.3f};  all points F {r.mean():.3f} +- {r.std():.3f} (min {r.min():.3f},"
              f" max {r.max():.3f})")
        print("      by group: " + "; ".join(f"{g} {r[k].mean():.3f} (n {k.sum()})"
                                             for g, k in groups.items() if k.any()))
        out[f"ratios_{tag}"] = r
    out["points"] = (z * 1e3, meas, cfd_pt)        # per-point residual (fig_adv_flume_points.csv)
    stack = {zp: np.mean([cell_mean(f, X0, zp + c, 0) for c in np.arange(-3, 4) * 1e-3])
             for zp in (-0.007, 0.0, 0.007)}       # the seven 1 mm cells at each ADV position
    out["dip"] = 1 - stack[0.0] / np.mean([stack[0.007], stack[-0.007]])
    out["axis_excess"] = float(plane(f, TIP - 0.075, np.array([0.0]))[0, 0])
    print(f"   CFD axial dip at x = tip - 5 mm (7-cell stacks, 1 - axis / mean of +-7 mm): "
          f"{100 * out['dip']:.1f} % (measured {100 * measured_dip(cm):.1f} %)")
    print(f"   CFD u/U on the axis 75 mm ahead of the tip: {out['axis_excess']:.3f}"
          f" (precursor axis-height 1.153 / 1.140 / 1.132)")
    fa = out["F_axis_point"]
    verdict = ("meets 1.00 +- 0.05" if abs(fa - 1) <= 0.05 else
               "partly closed (1.05-1.10)" if 1.05 < fa <= 1.10 else
               "gap persists (> 1.10)" if fa > 1.10 else "below 0.95")
    print(f"   criterion (axis cells, point): F = {fa:.3f}: {verdict}; open-domain k {OPEN_K}")
    return out


def axis_interp(path, U):
    """u/U on the symmetry plane near the inlet axis (|z| < 30 mm), over the whole length."""
    d = np.load(path)
    off, xz = d["off"], d["xz"].astype(float)
    cen = np.add.reduceat(xz, off[:-1]) / np.diff(off)[:, None]
    k = np.abs(cen[:, 1]) < 0.03
    from scipy.interpolate import LinearNDInterpolator
    return LinearNDInterpolator(cen[k], d["u"][k].astype(float) / U)


def closure(U, flume, U_open, open_):
    """Split the open-domain gap (k = 1.165) into (a) the centreline excess of the developed flume
    flow, (b) the sampler's blockage of the confined section, (c) the residual F.
    Ratio flume / open of u/U on the inlet axis: far upstream it is the excess (a); its change
    towards the tip, relative to far upstream, is what confinement adds near the sampler (b, an
    upper bound on blockage: it also holds the 114.3 vs 110 mm body and the sheared approach).
    At the ADV cells the same ratio, over (a), is (b) at the cells."""
    ff, fo = axis_interp(flume, U), axis_interp(open_, U_open)
    s = np.array([-700, -500, -300, -200, -100, -60, -40, -20, -10, -5]) * 1e-3
    uf = ff(np.column_stack([TIP + s, np.zeros_like(s)]))
    uo = fo(np.column_stack([TIP + s, np.zeros_like(s)]))
    a = uf[0] / uo[0]
    print(f"\n   ---- closure, flume U {U} against open domain U {U_open} ({open_})")
    print("   axis u/U, mm ahead of the tip: flume / open / ratio")
    for si, x, y in zip(s, uf, uo):
        print(f"     {si * 1e3:5.0f}  {x:.4f}  {y:.4f}  {x / y:.4f}")
    f1, f2 = interpolator(flume, U)[0], interpolator(open_, U_open)[0]
    cells = np.arange(-3, 4) * 1e-3
    ax_f = np.mean([cell_mean(f1, X0, c, 0) for c in cells])
    ax_o = np.mean([cell_mean(f2, X0, c, 0) for c in cells])
    al_f = np.mean([cell_mean(f1, X0, zp + c, 0) for zp in (-0.007, 0, 0.007) for c in cells])
    al_o = np.mean([cell_mean(f2, X0, zp + c, 0) for zp in (-0.007, 0, 0.007) for c in cells])
    b60 = (uf[5] / uo[5]) / a
    b_ax, b_all = ax_f / ax_o / a, al_f / al_o / a
    print(f"   (a) centreline excess, flume / open far upstream (700 mm ahead): {a:.4f}")
    print(f"   (b) confinement near the sampler, over (a): {100 * (b60 - 1):+.2f} % at 60 mm ahead; "
          f"{100 * (b_ax - 1):+.2f} % on the 7 axis cells; {100 * (b_all - 1):+.2f} % on all 21")
    print(f"   ADV cells flume / open: axis {ax_f:.3f} / {ax_o:.3f} = {ax_f / ax_o:.4f}; "
          f"all 21 {al_f:.3f} / {al_o:.3f} = {al_f / al_o:.4f}")
    return dict(a=a, b60=b60, b_axis=b_ax, b_all=b_all, open_axis=ax_o, open_all=al_o)


def figdata(rows, Us, slices):
    """fig:adv data for the flume runs, beside the open-domain files from plot_adv.py:
    fig_adv_flume_profile.csv (U, z_mm, u_cfd) at x = tip - 5 mm; fig_adv_flume_axis.csv
    (U, s_mm, u_over_U) along the inlet axis; fig_adv_flume_F.csv (per-run R and F);
    fig_adv_flume_points.csv (U, z_mm, u_meas, u_cfd, F) at each digitised point;
    fig_adv_flume_bore.csv (U, R, a, R_local), the bore means for panel (b)."""
    data = Path("../paper/figures/data")
    data.mkdir(parents=True, exist_ok=True)
    z = np.linspace(-12, 12, 241)
    s_mm = np.linspace(-60, 20, 321)
    with open(data / "fig_adv_flume_profile.csv", "w", newline="") as fp,             open(data / "fig_adv_flume_axis.csv", "w", newline="") as fa:
        wp, wa = csv.writer(fp), csv.writer(fa)
        wp.writerow(("U", "z_mm", "u_cfd"))
        wa.writerow(("U", "s_mm", "u_over_U"))
        for U, sl in zip(Us, slices):
            f = interpolator(sl, U)[0]
            wp.writerows((f"{U:g}", f"{a:.4g}", f"{b:.5f}") for a, b in zip(z, plane(f, X0, z * 1e-3)[:, 0]))
            wa.writerows((f"{U:g}", f"{a:.4g}", f"{b:.5f}")
                         for a, b in zip(s_mm, plane(f, TIP + s_mm * 1e-3, np.zeros_like(s_mm))[:, 0]))
    keys = ("U", "R", "axis_point", "axis_disc", "all_point", "all_disc", "F_axis_point",
            "F_axis_disc", "F_all_point", "sd_all_point", "n", "dip")
    with open(data / "fig_adv_flume_F.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        w.writerows([f"{r[k]:.5g}" if isinstance(r[k], float) else r[k] for k in keys] for r in rows)
    with open(data / "fig_adv_flume_points.csv", "w", newline="") as fh:
        w = csv.writer(fh)                              # u/U measured (digitised Fig. 3) and CFD
        w.writerow(("U", "z_mm", "u_meas", "u_cfd", "F"))
        for r in rows:
            for zz, m, c in zip(*r["points"]):
                w.writerow((f"{r['U']:g}", f"{zz:.3g}", f"{m:.5f}", f"{c:.5f}", f"{m / c:.5f}"))
    with open(data / "fig_adv_flume_bore.csv", "w", newline="") as fh:
        w = csv.writer(fh)                              # bore mean RU over the section-mean U is R;
        w.writerow(("U", "R", "a", "R_local"))          # a: centreline excess at the inlet height
        for r in rows:                                  # (closure, needs --open); R_local = R / a
            a = r.get("a", float("nan"))
            w.writerow((f"{r['U']:g}", f"{r['R']:.5f}", f"{a:.4f}", f"{r['R'] / a:.5f}"))
    print(f"figure data written to {data}/fig_adv_flume_*.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--U", type=float, action="append", required=True)
    ap.add_argument("--slice", action="append", required=True)
    ap.add_argument("--timeseries", action="append", required=True)
    ap.add_argument("--disc", type=float, default=0.007)
    ap.add_argument("--open", action="append", default=[],
                    help="open-domain slice per run (same order), for the closure split")
    ap.add_argument("--open-U", type=float, action="append", default=[])
    ap.add_argument("--figdata", action="store_true",
                    help="write the flume series for fig:adv to paper/figures/data/")
    a = ap.parse_args()
    rows = [score(U, s, t, a.disc) for U, s, t in zip(a.U, a.slice, a.timeseries)]
    for r, sl, op, uo in zip(rows, a.slice, a.open, a.open_U):
        r.update(closure(r["U"], sl, uo, op))
    if a.open:
        print("\nclosure: U, (a) excess, (b) at 60 mm / axis cells / all 21, (c) F axis / F all,"
              " product (a)(b_axis)(c_axis), open-domain k 1.165")
        for r in rows:
            if "a" in r:
                print(f"  {r['U']:.3f}  (a) {r['a']:.4f}  (b) {r['b60']:.4f} {r['b_axis']:.4f}"
                      f" {r['b_all']:.4f}  (c) {r['F_axis_point']:.3f} {r['F_all_point']:.3f}"
                      f"  product {r['a'] * r['b_axis'] * r['F_axis_point']:.3f}")
    if a.figdata:
        figdata(rows, a.U, a.slice)
    print("\nsummary: U, R, axis u/U (point, disc), all-21 u/U (point, disc), F axis (point, cell,"
          " disc), F all points mean +- sd (point), dip")
    for r in rows:
        print(f"  {r['U']:.3f}  R {r['R']:.4f}  axis {r['axis_point']:.3f} {r['axis_disc']:.3f}"
              f"  all-21 {r['all_point']:.3f} {r['all_disc']:.3f}  F_axis {r['F_axis_point']:.3f}"
              f" {r['F_axis_cell']:.3f} {r['F_axis_disc']:.3f}  F_all {r['F_all_point']:.3f} +-"
              f" {r['sd_all_point']:.3f} (n {r['n']})  dip {100 * r['dip']:.1f} %")
    if len(rows) > 1:
        allr = np.concatenate([r["ratios_point"] for r in rows])
        print(f"  pooled over {len(allr)} points: F {allr.mean():.3f} +- {allr.std():.3f}")


if __name__ == "__main__":
    main()
