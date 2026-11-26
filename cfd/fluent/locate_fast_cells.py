"""Mesh screen on a saved case-data file: fast cells and the smallest cells, in every zone.

usage: python locate_fast_cells.py <case.cas.h5> [--fast 1.0] [zone ...]
Counts cells with |U| above --fast (default 1 m/s) in each fluid zone, lists the fastest 15
overall with their positions, and the 15 smallest cells overall with their speed. Ends with
one line "SCREEN fast=<n> umax=<m/s> zone=<zone> vmin=<m3>" for the job summary.
"""
import argparse
import logging
import os
import sys

import numpy as np

sys.path.insert(0, os.getcwd())
import run_tims_case_v3 as r  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                    stream=sys.stdout)
LOG = logging.getLogger("screen")
ap = argparse.ArgumentParser()
ap.add_argument("case")
ap.add_argument("zones", nargs="*")
ap.add_argument("--fast", type=float, default=1.0)
a = ap.parse_args()

solver = r.launch("solver", int(os.environ.get("SLURM_NTASKS", 4)))
s = solver.settings
s.file.read(file_type="case-data", file_name=os.path.abspath(a.case))
zones = a.zones or list(s.setup.cell_zone_conditions.fluid.keys())
svd = solver.fields.solution_variable_data


def get(var, z):
    return np.asarray(svd.get_data(variable_name=var, zone_names=[z], domain_name="mixture")[z])


rows, total_fast = [], 0
for z in zones:
    sp = np.sqrt(get("SV_U", z)**2 + get("SV_V", z)**2 + get("SV_W", z)**2)
    xyz = get("SV_CENTROID", z).reshape(-1, 3)
    vol = get("SV_VOLUME", z)
    nf = int((sp > a.fast).sum())
    total_fast += nf
    LOG.info("%-28s %8d cells  |U| > %.2g: %6d  > 10: %5d  > 100: %4d  max %.4g  vmin %.3g", z,
             len(sp), a.fast, nf, int((sp > 10).sum()), int((sp > 100).sum()),
             sp.max() if len(sp) else 0, vol.min() if len(vol) else 0)
    for i in set(np.argsort(-sp)[:15]) | set(np.argsort(vol)[:15]):
        rows.append((z, float(sp[i]), float(vol[i]), *map(float, xyz[i])))


def show(title, items):
    LOG.info(title)
    for z, u, v, x, y, zc in items:
        LOG.info("  %-26s |U|=%9.4g  vol=%.3g  x=%.5f y=%.5f z=%.5f r=%.5f", z, u, v, x, y, zc,
                 (y * y + zc * zc) ** 0.5)


show("fastest cells:", sorted(rows, key=lambda t: -t[1])[:15])
show("smallest cells:", sorted(rows, key=lambda t: t[2])[:15])
top = max(rows, key=lambda t: t[1])
LOG.info("SCREEN fast=%d umax=%.4g zone=%s vmin=%.3g", total_fast, top[1], top[0],
         min(t[2] for t in rows))
logging.shutdown()
sys.stdout.flush()
# non-zero when the screen fails (fast cells), so an afterok-chained run does not start
os._exit(3 if total_fast else 0)
