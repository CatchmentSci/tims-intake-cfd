"""Connected fluid regions of a Fluent case: which cells are joined through interior faces, and
which boundary zones bound each region. A sampler whose inside is cut off from the outside shows
up as a separate region bounded only by walls.
Usage: python connectivity.py case.cas.h5
"""
import sys
from collections import Counter

import h5py
import numpy as np

m = h5py.File(sys.argv[1], "r")["meshes/1"]
ft = m["faces/zoneTopology"]
names = ft["name"][0]
names = (names.decode() if isinstance(names, bytes) else names).split(";")
fmin, fmax = ft["minId"][:].astype(np.int64), ft["maxId"][:].astype(np.int64)
ncell = int(m["cells/zoneTopology"]["maxId"][:].max())

nface = int(fmax.max())
c0 = np.zeros(nface + 1, np.int64); c1 = np.zeros(nface + 1, np.int64)   # indexed by face id
for grp, arr in (("faces/c0", c0), ("faces/c1", c1)):    # c1 sections cover interior faces only
    for key in m[grp].keys():
        ds = m[grp][key]
        lo, hi = int(ds.attrs["minId"][0]), int(ds.attrs["maxId"][0])
        arr[lo:hi + 1] = ds[:]
c0, c1 = c0[1:], c1[1:]
ids = np.arange(1, nface + 1)

inner = c1 > 0
a0, a1 = c0[inner] - 1, c1[inner] - 1
lab = np.arange(ncell)                       # hook-and-compress union-find (no scipy on the HPC)
while True:
    la, lb = lab[a0], lab[a1]
    d = la != lb
    if not d.any():
        break
    np.minimum.at(lab, np.maximum(la[d], lb[d]), np.minimum(la[d], lb[d]))
    while True:
        nl = lab[lab]
        if np.array_equal(nl, lab):
            break
        lab = nl
_, lab = np.unique(lab, return_inverse=True)
n = lab.max() + 1
size = np.bincount(lab)
print(f"cells {ncell}, interior faces {inner.sum()}, regions {n}")
zs = np.argsort(fmin)
zone_of = zs[np.searchsorted(fmax[zs], ids)]         # face id -> zone index
order = np.argsort(-size)
for r in order[:10]:
    b = (~inner) & (lab[c0 - 1] == r)
    cnt = Counter(names[z] for z in zone_of[b])
    print(f"region {r}: {size[r]} cells; boundary faces by zone: {dict(cnt.most_common())}")
