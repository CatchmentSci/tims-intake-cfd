"""Flow field on the symmetry plane (y = 0) of a Fluent case/data pair, one value per face taken
from the cell behind it (faces/c0).
Usage: python field_slice.py case.cas.h5 data.dat.h5 out.npz
Writes: off, xz (polygon vertices, m), names, zone, area, and u, v, w, p, k for each face, plus
mu_t, mu_lam and gamma (intermittency) where the data file has them, for the turbulent viscosity
ratio and the flow regime in the bores.
"""
import sys
import h5py
import numpy as np

cas, dat, out = sys.argv[1:4]
m = h5py.File(cas, "r")["meshes/1"]
d = h5py.File(dat, "r")["results/1/phase-1/cells"]

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
# Face sections are one per zone in a Meshing .msh.h5 but a single section in a solver .cas.h5;
# each section carries the face-id range it covers.
secs = []
for key in m["faces/nodes"].keys():
    g = m["faces/nodes"][key]
    lo = int(g.attrs["minId"][0]) if "minId" in g.attrs else None
    nn = g["nnodes"][:].astype(np.int64)
    secs.append(dict(key=key, lo=lo, nn=nn, off=np.concatenate([[0], np.cumsum(nn)])))
if any(s["lo"] is None for s in secs):          # Meshing file: sections follow zone order
    for s, lo in zip(sorted(secs, key=lambda s: int(s["key"])), fmin):
        s["lo"] = int(lo)


def zone_faces(zi):
    a, b = int(fmin[zi]), int(fmax[zi])
    s = next(s for s in secs if s["lo"] <= a and b < s["lo"] + len(s["nn"]))
    i0, i1 = a - s["lo"], b - s["lo"] + 1
    nodes = m["faces/nodes"][s["key"]]["nodes"][s["off"][i0]:s["off"][i1]].astype(np.int64)
    c0 = m["faces/c0"][s["key"]][i0:i1].astype(np.int64)
    return s["nn"][i0:i1], nodes, c0


polys, zidx, cells, keep = [], [], [], []
for zi, name in enumerate(names):
    if not name.startswith("symmetry") or name.startswith("symmetry_far"):
        continue
    nn, nodes, c0 = zone_faces(zi)
    P = coords[nodes]
    print(f"{name:24s} faces={len(nn):8d}  max|y|={np.abs(P[:, 1]).max():.1e}")
    if np.abs(P[:, 1]).max() > 1e-6:   # not on y = 0 (e.g. the flume side wall and floor)
        print("  skipped: not on the y = 0 plane")
        continue
    keep.append(name)
    o = np.concatenate([[0], np.cumsum(nn)])
    polys += [P[a:b][:, [0, 2]] for a, b in zip(o[:-1], o[1:])]
    zidx += [len(keep) - 1] * len(nn)
    cells.append(c0)

cells = np.concatenate(cells) - 1            # cell data are stored by cell id, from 1
fields = {k: d[sv]["1"][:][cells] for k, sv in
          dict(u="SV_U", v="SV_V", w="SV_W", p="SV_P", k="SV_K", mu_t="SV_MU_T",
               mu_lam="SV_MU_LAM", gamma="SV_INTERMIT").items() if sv in d}
off = np.concatenate([[0], np.cumsum([len(p) for p in polys])])
area = np.array([0.5 * abs(np.dot(p[:, 0], np.roll(p[:, 1], 1)) - np.dot(p[:, 1], np.roll(p[:, 0], 1)))
                 for p in polys])
np.savez_compressed(out, off=off, xz=np.concatenate(polys).astype(np.float32), names=np.array(keep),
                    zone=np.array(zidx, np.int16), area=area,
                    **{k: v.astype(np.float32) for k, v in fields.items()})
print("faces", len(polys), "max|U|", float(np.sqrt(fields["u"]**2 + fields["v"]**2 + fields["w"]**2).max()),
      "->", out)
