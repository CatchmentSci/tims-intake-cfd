"""Digitise Phillips et al. (2000) Fig. 2 (log inlet velocity against log ambient velocity).

Source: Phillips, Russell, Walling (2000), Hydrological Processes 14, 2589-2602, p. 2593 (PDF page 5),
Fig. 2. The figure is embedded in the PDF as a single 1-bit CCITT scan, 1212 x 1141 px at 300 ppi.
Extract it unresampled (pdfimages cannot open the original file name, which holds a Unicode dash;
copy the PDF to an ASCII name first):

    pdfimages -png -f 5 -l 5 phillips2000.pdf fig2      ->  fig2-000.png

Method (all automatic, no clicking):
  1. Axis calibration. The tick marks (x: 1.1 ... 1.8 below the x axis; y: 0 ... 1.6 left of the y
     axis) are found as ink runs on pixel rows / columns just outside each axis line; each tick's
     centre is averaged over 7 rows / columns. A straight line pixel -> value is fitted to each axis
     by least squares; its residuals give the calibration error.
  2. Markers. The six filled squares are the only solid blocks in the plot area: a pixel is kept if
     the whole K x K window round it is ink (K = 15 px; the plotted line is ~4 px wide and text
     strokes are < 12 px). Each surviving component is a square shrunk by (K-1)/2 on every side; the
     marker centre is the centre of its bounding box. The original marker's bounding box (the
     component's box grown back by (K-1)/2) is reported as its size.
  3. Check. The plotted regression line is located on columns clear of the markers and the text,
     and compared with the printed equation y = -2.182 + 2.074 x. Agreement tests the calibration
     and the drawing together.

Error estimate: the centre of a square edge-detected on a 1-bit scan is good to +-0.5 px, plus the
calibration residual. Reported in log units and as a percentage of the velocity.

Usage:  python geometry/digitise_phillips_fig2.py <fig2-000.png>   (from cfd/)
Writes results/phillips2000_fig2_digitised.csv.
"""
import csv
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

XT = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]
YT = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]
K = 15
OUT = Path("results/phillips2000_fig2_digitised.csv")


def runs(v, offset=0):
    idx = np.flatnonzero(v)
    if idx.size == 0:
        return []
    groups = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    return [(g[0] + offset + g[-1] + offset) / 2 for g in groups]


def axis_lines(ink):
    """Row of the x axis and column of the y axis (the longest straight ink runs)."""
    row = int(np.argmax(ink.sum(1)))
    col = int(np.argmax(ink.sum(0)))
    rows = np.flatnonzero(ink.sum(1) > 0.8 * ink.sum(1)[row])
    cols = np.flatnonzero(ink.sum(0) > 0.8 * ink.sum(0)[col])
    return rows.min(), rows.max(), cols.min(), cols.max()


def calibrate(ink):
    r0, r1, c0, c1 = axis_lines(ink)
    # x ticks: rows just below the axis; y ticks: columns just left of it
    xs = np.array([runs(ink[r, c0:], c0) for r in range(r1 + 3, r1 + 10)])
    ys = np.array([runs(ink[: r1 + 1, c], 0) for c in range(c0 - 12, c0 - 5)])
    if xs.shape[1] != len(XT) or ys.shape[1] != len(YT):
        raise RuntimeError(f"found {xs.shape[1]} x ticks and {ys.shape[1]} y ticks")
    px, py = xs.mean(0), ys.mean(0)[::-1]                      # y ticks top-down -> 0 first
    fx, fy = np.polyfit(px, XT, 1), np.polyfit(py, YT, 1)
    resx = np.polyval(fx, px) - XT
    resy = np.polyval(fy, py) - YT
    return dict(px=px, py=py, fx=fx, fy=fy, resx=resx, resy=resy, axes=(r0, r1, c0, c1))


def markers(ink, cal):
    r0, r1, c0, c1 = cal["axes"]
    plot = np.zeros_like(ink)
    plot[: r0 - 2, c1 + 2:] = ink[: r0 - 2, c1 + 2:]
    solid = ndimage.binary_erosion(plot, structure=np.ones((K, K), bool))
    lab, n = ndimage.label(solid)
    h = (K - 1) / 2
    out = []
    for sl in ndimage.find_objects(lab):
        rr, cc = sl
        top, bot = rr.start - h, rr.stop - 1 + h
        lef, rig = cc.start - h, cc.stop - 1 + h
        out.append(dict(xc=(lef + rig) / 2, yc=(top + bot) / 2, w=rig - lef + 1, hgt=bot - top + 1))
    return sorted(out, key=lambda m: m["xc"])


def line_check(ink, cal, mk):
    """Locate the plotted regression line on clear columns; compare with the printed equation."""
    r0, r1, c0, c1 = cal["axes"]
    rows = []
    for c in range(c1 + 10, ink.shape[1] - 60, 25):
        if any(abs(c - m["xc"]) < m["w"] / 2 + 6 for m in mk):
            continue
        x = np.polyval(cal["fx"], c)
        yexp = -2.182 + 2.074 * x
        rexp = np.interp(yexp, np.polyval(cal["fy"], np.arange(r0))[::-1], np.arange(r0)[::-1])
        win = ink[max(0, int(rexp) - 25): int(rexp) + 26, c]
        rs = runs(win, max(0, int(rexp) - 25))
        if len(rs) != 1:                       # text or nothing in the window: skip the column
            continue
        rows.append((x, np.polyval(cal["fy"], rs[0]) - yexp))
    return np.array(rows)


def main():
    img = Path(sys.argv[1])
    ink = np.array(Image.open(img)).astype(bool)
    if ink.mean() > 0.5:                       # ink must be the minority value
        ink = ~ink
    cal = calibrate(ink)
    dx, dy = abs(cal["fx"][0]), abs(cal["fy"][0])
    print(f"image {ink.shape[1]} x {ink.shape[0]} px")
    print(f"x axis: {1 / dx:.2f} px per log unit, tick residual rms {np.sqrt(np.mean(cal['resx'] ** 2)):.5f}"
          f" (max {np.abs(cal['resx']).max():.5f})")
    print(f"y axis: {1 / dy:.2f} px per log unit, tick residual rms {np.sqrt(np.mean(cal['resy'] ** 2)):.5f}"
          f" (max {np.abs(cal['resy']).max():.5f})")
    mk = markers(ink, cal)
    print(f"{len(mk)} markers (K = {K} px)")
    lc = line_check(ink, cal, mk)
    print(f"plotted line minus printed equation, {len(lc)} columns: mean {lc[:, 1].mean():+.4f},"
          f" s.d. {lc[:, 1].std():.4f}, range {lc[:, 1].min():+.4f} to {lc[:, 1].max():+.4f} log units")
    # error: +-0.5 px centre plus the calibration residual (rms), in log units
    ex = np.hypot(0.5 * dx, np.sqrt(np.mean(cal["resx"] ** 2)))
    ey = np.hypot(0.5 * dy, np.sqrt(np.mean(cal["resy"] ** 2)))
    pct = lambda e: 100 * (10 ** e - 1)
    print(f"pixel reading error: log U +-{ex:.4f} ({pct(ex):.2f} %), log u_in +-{ey:.4f} ({pct(ey):.2f} %)")
    # Drafting check against what the paper states: the end velocities are given in the text
    # (p. 2592: six ambient velocities "in the range 15.4-58.5 cm/s"), and the line has a printed equation.
    lx_end = [np.polyval(cal["fx"], mk[0]["xc"]), np.polyval(cal["fx"], mk[-1]["xc"])]
    d_end = np.array(lx_end) - np.log10([15.4, 58.5])
    print(f"end markers minus stated 15.4 / 58.5 cm/s: {d_end[0]:+.4f} / {d_end[1]:+.4f} log units")
    # Overall accuracy: the drafting errors seen in the two checks dominate the pixel error
    acc = float(np.sqrt(np.mean(np.r_[d_end, lc[:, 1]] ** 2)))
    acc = max(acc, ex, ey)
    print(f"overall accuracy (rms of the two checks): +-{acc:.4f} log units ({pct(acc):.2f} %) on each"
          f" coordinate; R +-{pct(np.sqrt(2) * acc):.2f} %")
    OUT.parent.mkdir(exist_ok=True)
    with OUT.open("w", newline="") as fh:
        fh.write("# Phillips et al. (2000) Hydrol. Process. 14:2589-2602, Fig. 2 (p. 2593), digitised by\n"
                 "# geometry/digitise_phillips_fig2.py from the embedded 300 ppi scan (marker centres, ticks\n"
                 f"# fitted per axis). Pixel reading error log U +-{ex:.4f}, log u_in +-{ey:.4f}. Drafting\n"
                 f"# checks: end markers vs the stated 15.4/58.5 cm/s {d_end[0]:+.4f}/{d_end[1]:+.4f}; plotted\n"
                 f"# line vs printed equation mean {lc[:, 1].mean():+.4f}, max {np.abs(lc[:, 1]).max():.4f}."
                 f" Overall +-{acc:.4f} log units per coordinate (R +-{pct(np.sqrt(2) * acc):.1f} %).\n")
        w = csv.writer(fh)
        w.writerow(["log_U", "log_u_in", "U_cm_s", "u_in_cm_s", "R", "R_eq1", "marker_px_w", "marker_px_h"])
        print("\n log U    log u_in   U (cm/s)  u_in (cm/s)   R       Eq.1 R   size (px)")
        for m in mk:
            lx, ly = np.polyval(cal["fx"], m["xc"]), np.polyval(cal["fy"], m["yc"])
            U, u = 10 ** lx, 10 ** ly
            req = 10 ** (-2.182 + 2.074 * lx) / U
            print(f" {lx:.4f}  {ly:.4f}   {U:7.2f}   {u:8.3f}   {u / U:.4f}   {req:.4f}   {m['w']:.0f} x {m['hgt']:.0f}")
            w.writerow([f"{lx:.4f}", f"{ly:.4f}", f"{U:.3f}", f"{u:.4f}", f"{u / U:.5f}", f"{req:.5f}",
                        f"{m['w']:.0f}", f"{m['hgt']:.0f}"])
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
