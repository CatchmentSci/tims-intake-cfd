"""CFD velocity at the ADV positions of Goharrokhi et al. (2019), Phillips sampler.

Goharrokhi et al. measured with a down-looking Vectrino II profiler 5 mm upstream of the inlet
tube tip, at three positions (on the inlet axis, 7 mm above and 7 mm below), seven 1 mm cells
each, and took that velocity as the inlet-tube velocity (u_in = 0.89 U - 1.3, cm/s). This samples
the CFD at the same 21 cells from symmetry-plane slices (fluent/field_slice.py), alongside the
bore velocity R*U, the axial profile ahead of the mouth and the vertical profile at x = tip - 5 mm
(their Fig. 3).

Cell values are averaged over the 1 mm cell height and, with --disc, over a horizontal disc of
that diameter (the ADV sampling-volume footprint). The disc average assumes the flow near the tip
is axisymmetric about the tube axis, so a point (dx, dy, z) takes the plane value at
(x + dx, sign(z) * sqrt(dy^2 + z^2)).

Run from cfd/:
  python geometry/goharrokhi_adv.py --U 0.3 --slice results/mesh_slices/v4_s1.00_U0.300_field.npz \
      --timeseries results/urans_v4/phillips2000_s1.00_sst-gamma/U0.300_timeseries.csv
Several --slice/--timeseries pairs (same U) may be given; each is reported in turn. A number in
place of a time-series path is taken as the bore R directly.
"""
import argparse
import csv
import sys

import numpy as np
from scipy.interpolate import LinearNDInterpolator

sys.path.insert(0, "geometry")
import tims_geometry as tg  # noqa: E402

X = tg.x_stations(tg.SPECS["phillips2000"])
TIP = X["in_tip"]                       # inlet tube tip, m (flow towards +x, axis at z = 0)
ADV_DX = -0.005                         # sampling volume centre 5 mm upstream of the tip
POSITIONS = {"below": -0.007, "axis": 0.0, "above": 0.007}
CELLS = np.arange(-3, 4) * 1e-3         # seven 1 mm cells, centred on each position


def goharrokhi_R(U):
    """Their fit, u_in = 0.89 U - 1.3 (cm/s), as a ratio."""
    return (0.89 * U * 100 - 1.3) / (U * 100)


def late_R(path):
    rows = [r for r in csv.DictReader(open(path))]
    t_end = float(rows[-1]["t"])
    late = [float(r["R"]) for r in rows if float(r["t"]) > t_end / 2]
    return float(np.mean(late))


def interpolator(path, U, win=(-0.20, -0.08, -0.04, 0.04)):
    d = np.load(path)
    off, xz = d["off"], d["xz"].astype(float)
    cen = np.add.reduceat(xz, off[:-1]) / np.diff(off)[:, None]
    k = (cen[:, 0] > win[0]) & (cen[:, 0] < win[1]) & (cen[:, 1] > win[2]) & (cen[:, 1] < win[3])
    vals = np.column_stack([d["u"][k], d["w"][k]]).astype(float) / U
    return LinearNDInterpolator(cen[k], vals), int(k.sum()), float(np.sqrt(d["area"][k]).min())


def plane(f, x, z):
    return f(np.column_stack([np.broadcast_to(x, np.shape(z)).ravel(), np.ravel(z)]))


def cell_mean(f, x0, zc, disc, n=11):
    """u/U over one 1 mm cell centred at (x0, zc), optionally over a disc of diameter `disc`."""
    dz = np.linspace(-0.5e-3, 0.5e-3, n)
    if disc <= 0:
        return float(np.nanmean(plane(f, x0, zc + dz)[:, 0]))
    r = disc / 2
    g = np.linspace(-r, r, 2 * n + 1)
    DX, DY = np.meshgrid(g, g)
    inside = DX**2 + DY**2 <= r * r
    DX, DY = DX[inside], DY[inside]
    out = []
    for z in zc + dz:
        zz = np.sqrt(DY**2 + z * z)
        sgn = np.sign(z) if z != 0 else 1.0
        u = plane(f, x0 + DX, sgn * zz)[:, 0]
        if z == 0:                                     # on the axis: average both sides
            u = 0.5 * (u + plane(f, x0 + DX, -zz)[:, 0])
        out.append(np.nanmean(u))
    return float(np.mean(out))


def report(path, ts, U, disc):
    f, nfaces, hmin = interpolator(path, U)
    R = late_R(ts) if ts and not ts.replace(".", "", 1).isdigit() else float(ts or "nan")
    x0 = TIP + ADV_DX
    print(f"\n== {path}\n   U = {U} m/s; {nfaces} slice faces in the window, smallest "
          f"{hmin * 1e3:.3f} mm; bore R (late mean) = {R:.4f}; "
          f"Goharrokhi fit at this U = {goharrokhi_R(U):.3f}")
    print(f"   ADV cells at x = tip - 5 mm (u/U; point = cell height only; "
          f"disc = {disc * 1e3:.0f} mm footprint, axisymmetric):")
    pos_means = {}
    for name, zp in POSITIONS.items():
        pts = [cell_mean(f, x0, zp + c, 0) for c in CELLS]
        dsk = [cell_mean(f, x0, zp + c, disc) for c in CELLS] if disc > 0 else pts
        pos_means[name] = (np.mean(pts), np.mean(dsk))
        print(f"   {name:5s} z = {zp * 1e3:+4.0f} mm  cells " +
              " ".join(f"{v:.3f}" for v in pts) +
              f"  | mean point {np.mean(pts):.3f}, disc {np.mean(dsk):.3f}")
    allp = np.mean([v[0] for v in pos_means.values()])
    alld = np.mean([v[1] for v in pos_means.values()])
    ax_p, ax_d = pos_means["axis"]
    print(f"   all 21 cells: point {allp:.3f}, disc {alld:.3f}")
    print(f"   axis position over bore R: point {ax_p / R:.2f}x, disc {ax_d / R:.2f}x; "
          f"axial dip vs mean of +-7 mm: "
          f"{100 * (1 - ax_p / np.mean([pos_means['above'][0], pos_means['below'][0]])):.1f} %")

    print("   centreline u/U ahead of the tip (mm from tip: value):")
    s = np.array([-60, -40, -30, -20, -15, -10, -8, -6, -5, -4, -3, -2, -1, -0.5, 0, 2, 5, 10, 20])
    u = plane(f, TIP + s * 1e-3, np.zeros_like(s, dtype=float))[:, 0]
    print("   " + "  ".join(f"{a:g}: {b:.3f}" for a, b in zip(s, u)))
    fine = np.linspace(-0.08, 0, 1601)
    uc = plane(f, TIP + fine, np.zeros_like(fine))[:, 0]
    for lev in (0.95, 0.9, 0.8, 0.5):
        j = np.nonzero(uc < lev)[0]
        print(f"   centreline first falls below {lev:.2f} U at {fine[j[0]] * 1e3:+.1f} mm"
              if len(j) else f"   centreline stays above {lev:.2f} U")

    print("   vertical profile at x = tip - 5 mm (z mm: u/U), cf. their Fig. 3:")
    zs = np.arange(-15, 16, 1.0)
    uz = plane(f, x0, zs * 1e-3)[:, 0]
    print("   " + "  ".join(f"{a:+.0f}: {b:.3f}" for a, b in zip(zs, uz)))
    return dict(slice=path, U=U, R=R, axis_point=ax_p, axis_disc=ax_d, all_point=allp,
                all_disc=alld, goharrokhi=goharrokhi_R(U))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--U", type=float, required=True)
    ap.add_argument("--slice", action="append", required=True)
    ap.add_argument("--timeseries", action="append", default=[])
    ap.add_argument("--disc", type=float, default=0.007,
                    help="ADV footprint diameter (m); 0 for cell-height averaging only. 7 mm: "
                         "Nortek Vectrino Profiler (= Vectrino II) datasheet, 2020-03-17, sampling "
                         "volume diameter")
    a = ap.parse_args()
    ts = a.timeseries + [None] * (len(a.slice) - len(a.timeseries))
    rows = [report(p, t, a.U, a.disc) for p, t in zip(a.slice, ts)]
    print("\nsummary (u/U): slice, bore R, ADV axis point, axis disc, all-21 point, all-21 disc, "
          "Goharrokhi fit")
    for r in rows:
        print(f"  {r['slice']}: {r['R']:.4f}  {r['axis_point']:.3f}  {r['axis_disc']:.3f}  "
              f"{r['all_point']:.3f}  {r['all_disc']:.3f}  {r['goharrokhi']:.3f}")


if __name__ == "__main__":
    main()
