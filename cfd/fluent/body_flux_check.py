"""Through-flow of the sampler body from a Fluent case/data pair, on a complete cut through the
chamber at mid-length.

The reporting zone interior_body_mid reads 67-79 % of the mouth flux in some meshes (the disc is
split where it meets the cavity's near-wall region), while the two mouths balance to 1e-8-1e-5.
This script rebuilds the cut from the mesh itself: every interior face inside the chamber
(r < body bore radius) whose two cells lie on opposite sides of x_mid, oriented downstream. The
flux through each face is rho * (u0 + u1)/2 . A (area vector from the face polygon), or the
solver's face mass flux SV_FLUX where the data file stores it. Reads the HDF5 files directly: no
Fluent licence.

Usage: python body_flux_check.py CASE_NAME case.cas.h5 data.dat.h5 [summary.json]
  CASE_NAME is the tims_geometry spec (phillips2000, perks2014, elliott2017, ...); the chamber
  stations are tabulated below from tims_geometry.x_stations. Prints one line: the cut flux, the mouth flux
  from summary.json and their ratio.
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

# Chamber mid-length x and bore radius (m), from tims_geometry.x_stations (cav0, cav1) and body_id;
# tabulated so the script runs without gmsh on the HPC
CHAMBER = {
    "phillips2000": (0.5100, 0.04900),
    "perks2014": (0.5050, 0.04035),
    "perks2014_flowmeter": (0.5050, 0.04035),
    "elliott2017": (0.5357, 0.05113),
    "mcdonald2010": (0.1190, 0.03175),
    "foets2020": (0.5100, 0.04900),
}
RHO = 998.2
spec, cas, dat = sys.argv[1:4]
summ = sys.argv[4] if len(sys.argv) > 4 else None
x_mid, r_in = CHAMBER[spec.split("_full")[0].split("_yaw")[0]]

m = h5py.File(cas, "r")["meshes/1"]
dres = h5py.File(dat, "r")["results/1/phase-1"]
cells = dres["cells"]

nt = m["nodes/zoneTopology"]
nmin = nt["minId"][:].astype(np.int64)
coords = np.zeros((int(nt["maxId"][:].max()) + 1, 3))
for i, key in enumerate(sorted(m["nodes/coords"].keys(), key=int)):
    c = m["nodes/coords"][key][:]
    coords[nmin[i]:nmin[i] + len(c)] = c

ft = m["faces/zoneTopology"]
names = ft["name"][0]
names = (names.decode() if isinstance(names, bytes) else names).split(";")
fmin, fmax = ft["minId"][:].astype(np.int64), ft["maxId"][:].astype(np.int64)
zid = ft["id"][:].astype(np.int64) if "id" in ft else np.arange(1, len(names) + 1)

# Face connectivity as global arrays indexed by face id - 1. In a solver .cas.h5 the faces/nodes and
# faces/c0 data are usually one section covering every face, while faces/c1 is split into sections
# (interior zones only), each with its own minId; fill them all into full-length arrays.
NF = int(fmax.max())


def full(group, fill=0):
    out = np.full(NF, fill, dtype=np.int64)
    for k in group.keys():
        ds = group[k]
        lo = int(ds.attrs["minId"][0]) if "minId" in ds.attrs else 1
        a = ds[:]
        out[lo - 1:lo - 1 + len(a)] = a
    return out


C0all = full(m["faces/c0"])
C1all = full(m["faces/c1"]) if "c1" in m["faces"] else np.zeros(NF, dtype=np.int64)

cent, avec, FID = [], [], []
for key in m["faces/nodes"].keys():
    g = m["faces/nodes"][key]
    nn = g["nnodes"][:].astype(np.int64)
    lo = int(g.attrs["minId"][0]) if "minId" in g.attrs else 1
    nodes = g["nodes"][:].astype(np.int64)
    off = np.concatenate([[0], np.cumsum(nn)])
    P = coords[nodes]
    idx = np.repeat(np.arange(len(nn)), nn)
    cx = np.bincount(idx, P[:, 0], len(nn)) / nn          # face centroid (node mean), pre-filter
    cy = np.bincount(idx, P[:, 1], len(nn)) / nn
    cz = np.bincount(idx, P[:, 2], len(nn)) / nn
    c1 = C1all[lo - 1:lo - 1 + len(nn)]
    near = (np.abs(cx - x_mid) < 0.06) & (np.hypot(cy, cz) < r_in + 1e-4) & (c1 > 0)
    for f in np.nonzero(near)[0]:
        q = P[off[f]:off[f + 1]]
        avec.append(0.5 * np.cross(q, np.roll(q, -1, axis=0)).sum(axis=0))
        cent.append((cx[f], cy[f], cz[f]))
        FID.append(lo + f)
    del P, idx, cx, cy, cz, nodes
cent, avec, FID = np.array(cent), np.array(avec), np.array(FID)
C0, C1 = C0all[FID - 1], C1all[FID - 1]

# Cell centroid x from the faces gathered (the slab is wider than any chamber cell)
cell_ids = np.unique(np.concatenate([C0, C1]))
pos = {c: i for i, c in enumerate(cell_ids)}
i0 = np.array([pos[c] for c in C0])
i1 = np.array([pos[c] for c in C1])
sx = np.bincount(np.concatenate([i0, i1]), np.concatenate([cent[:, 0], cent[:, 0]]), len(cell_ids))
cnt = np.bincount(np.concatenate([i0, i1]), minlength=len(cell_ids))
xc = sx / np.maximum(cnt, 1)
up0, up1 = xc[i0] < x_mid, xc[i1] < x_mid
cut = up0 != up1
sign = np.where(up0[cut], 1.0, -1.0)            # +: c0 upstream, so c0 -> c1 is downstream

# Face mass flux: SV_FLUX if stored on faces (global, by face id), else interpolate cell velocity
flux = None
if "faces" in dres and "SV_FLUX" in dres["faces"]:
    sv = dres["faces/SV_FLUX"]
    fl = np.full(NF, np.nan)
    for k in sv.keys():
        ds = sv[k]
        if "minId" in ds.attrs:
            lo = int(ds.attrs["minId"][0])
        elif len(ds) == NF:
            lo = 1
        else:
            continue
        a = ds[:]
        fl[lo - 1:lo - 1 + len(a)] = a
    vals = fl[FID[cut] - 1]
    if np.isfinite(vals).all():
        flux = (sign * vals).sum()
        how = "SV_FLUX"
if flux is None:
    zone = sorted(cells["SV_U"].keys(), key=int)
    if len(zone) != 1:
        print("note: several cell zones; using zone", zone[0])
    U = np.stack([cells[sv][zone[0]][:] for sv in ("SV_U", "SV_V", "SV_W")], axis=1)
    ua = 0.5 * (U[C0[cut] - 1] + U[C1[cut] - 1])
    flux = RHO * (sign * np.einsum("ij,ij->i", ua, avec[cut])).sum()
    how = "cell-velocity interpolation"

mouth = None
if summ and Path(summ).exists():
    mouth = abs(json.load(open(summ))["mean"]["mdot_mouth_down"])
print(f"{spec} {Path(dat).name}: x_mid {x_mid:.3f} m, {cut.sum()} cut faces, body flux {flux:+.4e} kg/s "
      f"({how})" + (f", mouth {mouth:.4e}, ratio {flux / mouth:.4f}" if mouth else ""))
