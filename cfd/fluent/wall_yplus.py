"""First-cell y+ on every wall zone of a Fluent case/data pair (
y+ on the body, bore and outlet walls). Reads the .cas.h5/.dat.h5 with h5py, as field_slice.py
does; no Fluent session and no licence.

Usage: python wall_yplus.py case.cas.h5 data.dat.h5 out.csv [--first-height 2.5e-5]

Every wall carries 10 uniform prism layers with first height 25 um x the mesh size scale
(run_tims_case_v3.py, "Add Boundary Layers"), so the wall-adjacent cell centre sits y_p = h1/2 from
the wall. y+ is Fluent's own face field SV_WALL_YPLUS when the data file stores it (it does for the
production runs). As a check, and where it is missing, y+ is also computed from the wall-adjacent
cell: the tangential speed u_p at the cell's wall distance (SV_WALL_DIST, else h1/2), with tau_w = mu u_p / y_p in the viscous sublayer and the log law
u_p / u_tau = ln(9.793 y+) / 0.4187 above y+ = 11.2 (Fluent's standard constants). y+ is reported
for the cell centre (Fluent's definition) as area-weighted mean, 95th percentile and maximum per wall
zone, with the location of the maximum.
"""
import csv
import re
import sys

import h5py
import numpy as np

RHO, MU = 998.2, 1.002e-3
args = [a for a in sys.argv[1:] if not a.startswith("--")]
cas, dat, out = args[:3]
h1 = float(sys.argv[sys.argv.index("--first-height") + 1]) if "--first-height" in sys.argv else None
if h1 is None:                                    # size scale from the case name, e.g. _s0.64
    m = re.search(r"_s(\d\.\d+)", cas)
    h1 = 2.5e-5 * (float(m.group(1)) if m else 1.0)
y_p = h1 / 2
print(f"first-layer height {h1 * 1e6:.2f} um, cell centre y_p {y_p * 1e6:.2f} um")

m = h5py.File(cas, "r")["meshes/1"]
dres = h5py.File(dat, "r")["results/1/phase-1"]
d = dres["cells"]
print("cell fields:", sorted(d.keys()))
dfaces = dres["faces"] if "faces" in dres else None
if dfaces is not None:
    print("face fields:", sorted(dfaces.keys()))

nt = m["nodes/zoneTopology"]
nmin, nmax = nt["minId"][:].astype(np.int64), nt["maxId"][:].astype(np.int64)
coords = np.zeros((int(nmax.max()) + 1, 3))
for i, key in enumerate(sorted(m["nodes/coords"].keys(), key=int)):
    c = m["nodes/coords"][key][:]
    coords[nmin[i]:nmin[i] + len(c)] = c

ft = m["faces/zoneTopology"]
names = ft["name"][0]
names = (names.decode() if isinstance(names, bytes) else names).split(";")
fmin, fmax = ft["minId"][:].astype(np.int64), ft["maxId"][:].astype(np.int64)
secs = []
for key in m["faces/nodes"].keys():
    g = m["faces/nodes"][key]
    lo = int(g.attrs["minId"][0]) if "minId" in g.attrs else None
    nn = g["nnodes"][:].astype(np.int64)
    secs.append(dict(key=key, lo=lo, nn=nn, off=np.concatenate([[0], np.cumsum(nn)])))
if any(s["lo"] is None for s in secs):
    for s, lo in zip(sorted(secs, key=lambda s: int(s["key"])), fmin):
        s["lo"] = int(lo)


def zone_faces(zi):
    a, b = int(fmin[zi]), int(fmax[zi])
    s = next(s for s in secs if s["lo"] <= a and b < s["lo"] + len(s["nn"]))
    i0, i1 = a - s["lo"], b - s["lo"] + 1
    nodes = m["faces/nodes"][s["key"]]["nodes"][s["off"][i0]:s["off"][i1]].astype(np.int64)
    c0 = m["faces/c0"][s["key"]][i0:i1].astype(np.int64)
    return s["nn"][i0:i1], nodes, c0, a


def face_geometry(nn, nodes):
    """Area vector (Newell) and centroid of each polygon face."""
    o = np.concatenate([[0], np.cumsum(nn)])
    P = coords[nodes]
    nxt = np.arange(len(P)) + 1
    nxt[o[1:] - 1] = o[:-1]                       # wrap each polygon to its first node
    cr = np.cross(P, P[nxt])
    A = 0.5 * np.add.reduceat(cr, o[:-1])
    cen = np.add.reduceat(P, o[:-1]) / nn[:, None]
    return A, cen


ZONE_IDS = ft["id"][:].astype(np.int64) if "id" in ft else None
print("face-zone ids:", None if ZONE_IDS is None else ZONE_IDS.tolist()[:60])
for k in ("SV_WALL_YPLUS", "SV_WALL_SHEAR"):
    if dfaces is not None and k in dfaces:
        g = dfaces[k]
        print(k, "members:", [(n, getattr(g[n], "shape", None), dict(getattr(g[n], "attrs", {})))
                               for n in list(g.keys())[:8]])


def face_field(k, zi, a, n):
    """Values of face field k for the n faces of zone index zi (first global face id a).
    Fluent writes face data in one dataset per face zone, keyed by the zone id; a dataset with
    minId/maxId attributes is matched by face id instead. Returns None if nothing matches n."""
    if dfaces is None or k not in dfaces:
        return None
    g = dfaces[k]
    if ZONE_IDS is not None and str(int(ZONE_IDS[zi])) in g:
        v = g[str(int(ZONE_IDS[zi]))][:]
        if len(v) == n:
            return v
        print(f"  {k} zone {ZONE_IDS[zi]}: {len(v)} values for {n} faces, not used")
    for sub in g.keys():
        ds = g[sub]
        if "minId" in ds.attrs and "maxId" in ds.attrs:
            lo, hi = int(ds.attrs["minId"][0]), int(ds.attrs["maxId"][0])
            if lo <= a and a + n - 1 <= hi:
                return ds[a - lo:a - lo + n]
    return None


U = {k: d[sv]["1"][:] for k, sv in dict(u="SV_U", v="SV_V", w="SV_W").items()}
WD = d["SV_WALL_DIST"]["1"][:] if "SV_WALL_DIST" in d else None
rows = []
for zi, name in enumerate(names):
    if not name.startswith("wall"):
        continue
    nn, nodes, c0, a = zone_faces(zi)
    A, cen = face_geometry(nn, nodes)
    area = np.linalg.norm(A, axis=1)
    nrm = A / np.maximum(area, 1e-30)[:, None]
    cells = c0 - 1
    vel = np.stack([U["u"][cells], U["v"][cells], U["w"][cells]], axis=1)
    ut = np.linalg.norm(vel - np.sum(vel * nrm, axis=1)[:, None] * nrm, axis=1)
    ypc = WD[cells] if WD is not None else np.full(len(cells), y_p)   # cell-centre wall distance
    utau = np.sqrt(MU * ut / ypc / RHO)              # cross-check: viscous sublayer first,
    for _ in range(30):                              # log law where y+ > 11.2
        yp = RHO * utau * ypc / MU
        log = yp > 11.225
        if not log.any():
            break
        utau[log] = ut[log] * 0.4187 / np.log(9.793 * np.maximum(yp[log], 1e-9))
    yplus_cell = RHO * utau * ypc / MU
    yf = face_field("SV_WALL_YPLUS", zi, a, len(nn))
    if yf is not None:
        yplus, source = np.ravel(yf).astype(float), "SV_WALL_YPLUS"
    else:
        yplus, source = yplus_cell, "cell velocity"
    i = int(np.argmax(yplus))
    order = np.argsort(yplus)
    cw = np.cumsum(area[order]) / area.sum()
    rows.append(dict(zone=name.split(":")[0], faces=len(nn), area_m2=f"{area.sum():.4e}", source=source,
                     yplus_mean=f"{np.average(yplus, weights=area):.3f}",
                     yplus_p95=f"{yplus[order][np.searchsorted(cw, 0.95)]:.3f}",
                     yplus_max=f"{yplus.max():.3f}",
                     max_at=" ".join(f"{v:.4f}" for v in cen[i]),
                     yplus_cell_mean=f"{np.average(yplus_cell, weights=area):.3f}",
                     yplus_cell_max=f"{yplus_cell.max():.3f}",
                     wall_dist_um_mean=f"{1e6 * np.average(ypc, weights=area):.2f}"))
    print(f"{name:36s} faces {len(nn):7d}  y+ mean {float(rows[-1]['yplus_mean']):7.3f}  p95 "
          f"{float(rows[-1]['yplus_p95']):7.3f}  max {yplus.max():7.3f} at {rows[-1]['max_at']}  ({source}; "
          f"cell check mean {float(rows[-1]['yplus_cell_mean']):.3f}, y_p {rows[-1]['wall_dist_um_mean']} um)")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print("wrote", out)
