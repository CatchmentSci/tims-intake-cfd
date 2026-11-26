"""Digitise Figure 3 of Goharrokhi et al. (2019): ADV velocity profiles 5 mm upstream of the
inlet of an original-design sampler, at flume velocities of 15, 21, 30, 35 and 60 cm/s.

The figure is embedded in the PDF as a raster (1308 x 550 px). Axes are calibrated on the
gridlines (x: 0 and 60 cm/s at px 166.5 and 949.5; y: 16 and 13 cm at px 14.5 and 436.5).
Markers are found as dark connected components and assigned to series by velocity range and,
where the 30 and 35 cm/s profiles overlap, by fill fraction (open triangles about 0.45, plus
signs about 0.3). Fused or clipped markers are dropped, so some series have fewer than 21
points. Check the output with --overlay.

Usage (from cfd/): python geometry/digitise_goharrokhi_fig3.py ARTICLE.pdf [--overlay out.png]
ARTICLE.pdf is the published article (Hydrological Processes 33: 2048-2056), which users obtain
themselves; the figure is read from its sixth page.
Writes results/goharrokhi2019_fig3_digitised.csv (series_cm_s, u_cm_s, z_cm, z_mm_from_axis).
"""
import argparse
import csv
from pathlib import Path

import numpy as np
import pypdf
from PIL import Image, ImageDraw
from scipy import ndimage

OUT = Path("results/goharrokhi2019_fig3_digitised.csv")
X0, X60, Y16, Y13 = 166.5, 949.5, 14.5, 436.5
AXIS_CM = 14.5                                  # inlet axis height above the flume bed
GRID_COLS = (166, 167, 297, 298, 427, 428, 558, 559, 688, 689, 819, 820, 949, 950)


def u_of(px):
    return (px - X0) * 60 / (X60 - X0)


def z_of(py):
    return 16 - (py - Y16) * 3 / (Y13 - Y16)


def series_of(u, fill, w):
    if 12 < u < 16:
        return 15
    if 17 < u < 22:
        return 21
    if 22 < u < 36.5:
        if u < 28.5:
            return 30
        return 30 if fill >= 0.38 and w >= 10 else 35   # narrow = plus sign cut by a gridline
    if u > 45:
        return 60
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf", help="the published article (PDF)")
    ap.add_argument("--overlay")
    a = ap.parse_args()
    img = pypdf.PdfReader(a.pdf).pages[5].images[0].image.convert("RGB")
    g = np.asarray(img).astype(float).mean(axis=2)
    dark = g < 110
    dark[:, :168] = False; dark[:, 958:] = False; dark[:17, :] = False; dark[434:, :] = False
    for c in GRID_COLS[2:]:
        dark[:, c] = False
    dark[224:228, :] = False                    # 14.5 cm gridline
    lab, _ = ndimage.label(dark, structure=np.ones((3, 3)))
    rows = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        m = lab[sl] == i
        h, w = m.shape
        if m.sum() < 30 or w > 14 or h > 14:     # fragments and fused pairs
            continue
        cx = (sl[1].start + sl[1].stop - 1) / 2
        cy = (sl[0].start + sl[0].stop - 1) / 2
        if h < 9:                               # cut by the 14.5 cm gridline: centre on it
            cy = 225.5 if abs(cy - 225.5) < 8 else cy
        u, z = u_of(cx), z_of(cy)
        s = series_of(u, m.sum() / (w * h), w)
        if s is not None:
            rows.append((s, round(u, 2), round(z, 3), round((z - AXIS_CM) * 10, 2), cx, cy))
    rows.sort()
    with open(OUT, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["series_cm_s", "u_cm_s", "z_cm", "z_mm_from_axis"])
        wr.writerows(r[:4] for r in rows)
    for s in (15, 21, 30, 35, 60):
        pts = [r for r in rows if r[0] == s]
        ax = [r[1] for r in pts if -5 <= r[3] <= 2.5]
        print(f"{s:2d} cm/s: {len(pts):2d} points; axis stack (-5 to +2.5 mm) n={len(ax)} "
              f"mean {np.mean(ax):.2f} cm/s = {np.mean(ax) / s:.3f} U")
    if a.overlay:
        big = img.resize((img.width * 2, img.height * 2))
        d = ImageDraw.Draw(big)
        col = {15: "red", 21: "orange", 30: "green", 35: "blue", 60: "magenta"}
        for s, u, z, _, cx, cy in rows:
            d.ellipse([2 * cx - 9, 2 * cy - 9, 2 * cx + 9, 2 * cy + 9], outline=col[s], width=2)
        big.save(a.overlay)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
