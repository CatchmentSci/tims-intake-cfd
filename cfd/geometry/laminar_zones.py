"""Forced-laminar bore runs: which cell zones lie in each laminar-bore cylinder, from the mesh file.

The driver's --laminar-bore splits the cells inside a cylinder (x0 < x < x1, r < R on the x axis)
out of every zone that straddles it. A zone lying wholly inside has nothing to split, and the
Test runs showed that such zones (the bore bodies) stayed SST-gamma. This script reads the
.msh.h5, computes every cell's centroid and volume, and lists for each bore:
  whole  zones with every cell centroid inside the cylinder (the driver takes them whole);
  split  zones with some cells inside (the driver splits them, as before);
and the cylinder's cell volume V_cyl (cells with centroid inside), which the driver uses as the
reference for its laminar volume fraction. It also predicts the fraction an earlier driver version
(split zones only) gave, to compare with the screens' slices (9 % of each Perks bore).

Cell centroids and volumes: each face's centroid and area vector (fan about the node mean); a
provisional cell centre (area-weighted mean of its face centroids); then the pyramids on the faces:
V = sum |A.(x_f - x_c)| / 3, centroid = sum V_p (x_c + 3/4 (x_f - x_c)) / V.

Usage (from cfd/): python geometry/laminar_zones.py <case> <mesh.msh.h5> [out folder]
Writes <out>/laminar_zones_<case>.json and prints the table and the driver arguments.
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

BORES = {   # (x0, x1, R) in m, the laminar-bore cylinders
    "perks2014": [(-0.1365, 0.0240, 0.0040), (0.9860, 1.1465, 0.0040)],
    "perks2014_flowmeter": [(-0.1365, 0.0240, 0.0040), (0.0245, 0.1410, 0.0050),
                            (0.9860, 1.1465, 0.0040)],
}


def text(v):
    v = v[0] if hasattr(v, "__len__") and not isinstance(v, (str, bytes)) else v
    return v.decode() if isinstance(v, bytes) else str(v)


def main():
    case, mesh = sys.argv[1], sys.argv[2]
    out = Path(sys.argv[3] if len(sys.argv) > 3 else "results/checks")
    m = h5py.File(mesh, "r")["meshes/1"]

    nt = m["nodes/zoneTopology"]
    nmin = nt["minId"][:].astype(np.int64)
    nmax = nt["maxId"][:].astype(np.int64)
    coords = np.zeros((int(nmax.max()) + 1, 3))
    for i, key in enumerate(sorted(m["nodes/coords"].keys(), key=int)):
        c = m["nodes/coords"][key][:]
        coords[nmin[i]:nmin[i] + len(c)] = c

    ct = m["cells/zoneTopology"]
    cnames = text(ct["name"]).split(";")
    cmin = ct["minId"][:].astype(np.int64)
    cmax = ct["maxId"][:].astype(np.int64)
    ncell = int(cmax.max())
    czone = np.full(ncell + 1, -1, np.int32)
    for i, (a, b) in enumerate(zip(cmin, cmax)):
        czone[a:b + 1] = i

    # Faces, section by section: centroid, area vector, cells either side
    FC, AV, C0, C1 = [], [], [], []
    for key in m["faces/nodes"].keys():
        g = m["faces/nodes"][key]
        nn = g["nnodes"][:].astype(np.int64)
        nodes = g["nodes"][:].astype(np.int64)
        fi = np.repeat(np.arange(len(nn)), nn)
        start = np.concatenate([[0], np.cumsum(nn)[:-1]])
        pos = np.arange(len(nodes)) - start[fi]
        nxt = np.where(pos + 1 == nn[fi], start[fi], np.arange(len(nodes)) + 1)
        P = coords[nodes]
        fc = np.stack([np.bincount(fi, P[:, k], len(nn)) for k in range(3)], 1) / nn[:, None]
        cr = np.cross(P - fc[fi], P[nxt] - fc[fi])
        av = 0.5 * np.stack([np.bincount(fi, cr[:, k], len(nn)) for k in range(3)], 1)
        FC.append(fc)
        AV.append(av)
        C0.append(m["faces/c0"][key][:].astype(np.int64))
        C1.append(m["faces/c1"][key][:].astype(np.int64) if key in m["faces/c1"]
                  else np.zeros(len(nn), np.int64))
        del nodes, fi, P, cr
    fc, av, c0, c1 = (np.concatenate(a) for a in (FC, AV, C0, C1))
    del FC, AV, C0, C1
    if len(c1) != len(c0):
        raise RuntimeError("c0 and c1 differ in length")
    a = np.linalg.norm(av, axis=1)

    # Provisional centres, then pyramid volumes and centroids
    side = [(c0, fc, a)] + [(c1[c1 > 0], fc[c1 > 0], a[c1 > 0])]
    w = sum(np.bincount(c, aa, ncell + 1) for c, _, aa in side)
    xc = np.stack([sum(np.bincount(c, f[:, k] * aa, ncell + 1) for c, f, aa in side)
                   for k in range(3)], 1) / np.maximum(w, 1e-300)[:, None]
    vol = np.zeros(ncell + 1)
    mom = np.zeros((ncell + 1, 3))
    for c, f, av_ in [(c0, fc, av), (c1[c1 > 0], fc[c1 > 0], av[c1 > 0])]:
        d = f - xc[c]
        vp = np.abs(np.einsum("ij,ij->i", av_, d)) / 3.0
        vol += np.bincount(c, vp, ncell + 1)
        g = xc[c] + 0.75 * d
        for k in range(3):
            mom[:, k] += np.bincount(c, vp * g[:, k], ncell + 1)
    cen = mom / np.maximum(vol, 1e-300)[:, None]
    cells = np.arange(1, ncell + 1)
    if (czone[cells] < 0).any():
        raise RuntimeError("cells outside every cell zone")
    zvol = np.bincount(czone[cells], vol[cells], len(cnames))
    print(f"{case}: {mesh}; {ncell} cells in {len(cnames)} zones; total volume {zvol.sum():.6e} m3")

    res = dict(case=case, mesh=str(mesh), cells=ncell, bores=[])
    args = []
    for i, (x0, x1, R) in enumerate(BORES[case]):
        x, y, z = cen[cells, 0], cen[cells, 1], cen[cells, 2]
        ins = (x > x0) & (x < x1) & (y * y + z * z < R * R)
        zin = np.bincount(czone[cells][ins], minlength=len(cnames))
        ztot = np.bincount(czone[cells], minlength=len(cnames))
        vin = np.bincount(czone[cells][ins], vol[cells][ins], len(cnames))
        vcyl = vin.sum()
        whole, split, rows = [], [], []
        for zi, name in enumerate(cnames):
            if zin[zi] == 0:
                continue
            kind = "whole" if zin[zi] == ztot[zi] else "split"
            (whole if kind == "whole" else split).append(name)
            rows.append(dict(zone=name, kind=kind, cells=int(ztot[zi]), cells_inside=int(zin[zi]),
                             volume=float(zvol[zi]), volume_inside=float(vin[zi])))
        old = sum(r["volume_inside"] for r in rows if r["kind"] == "split") / vcyl
        new = (sum(r["volume"] for r in rows if r["kind"] == "whole")
               + sum(r["volume_inside"] for r in rows if r["kind"] == "split")) / vcyl
        ana = np.pi * R * R * (x1 - x0)
        res["bores"].append(dict(index=i, x0=x0, x1=x1, R=R, cells=int(ins.sum()), V_cyl=float(vcyl),
                                 V_cylinder_analytic=float(ana), zones=rows, whole=whole, split=split,
                                 fraction_split_only=float(old), fraction_whole_and_split=float(new)))
        print(f"  bore {i}: x {x0}..{x1} m, R {R * 1e3:.1f} mm: {ins.sum()} cells, V_cyl {vcyl:.6e} m3 "
              f"({100 * vcyl / ana:.2f} % of pi R^2 L)")
        for r in rows:
            print(f"    {r['kind']:5s} {r['zone']:24s} {r['cells_inside']:8d} of {r['cells']:8d} cells inside, "
                  f"V inside {r['volume_inside']:.4e} of {r['volume']:.4e} m3")
        print(f"    laminar volume fraction: split zones only (23-24 Nov driver) {100 * old:.2f} %; "
              f"whole + split (fixed driver) {100 * new:.2f} %")
        args += [f"--laminar-zone {i} {n} whole" for n in whole]
        args += [f"--laminar-zone {i} {n} split" for n in split]
        args.append(f"--bore-volume {i} {vcyl:.6e}")
    res["driver_args"] = " ".join(args)
    print("  driver arguments:", res["driver_args"])
    if "," in res["driver_args"]:
        raise RuntimeError("comma in the driver arguments")
    out.mkdir(parents=True, exist_ok=True)
    (out / f"laminar_zones_{case}.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
