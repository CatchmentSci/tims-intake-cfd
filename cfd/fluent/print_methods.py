"""Print the solution-method settings stored in a Fluent .cas.h5:
spatial discretization of each equation, the pressure scheme, gradient and time settings and the
pressure-velocity coupling. Reads the case file's rp-variables text directly (no Fluent licence).

Fluent codes (from the case file's own option lists): convection 0 first-order upwind, 1 second-order
upwind, 2 power law, 3 central, 4 QUICK, 6 third-order MUSCL, 7 bounded central; pressure 10
standard, 11 linear, 12 second order, 13 body-force weighted, 14 PRESTO!; flow/scheme 20 SIMPLE,
21 SIMPLEC, 22 PISO, 24 coupled.

Usage: python print_methods.py case.cas.h5 [case2.cas.h5 ...]
"""
import re
import sys

import h5py

# (name value) pairs in the rp-variable text; keep the solution-method ones
PAIR = re.compile(r"\(([A-Za-z][\w/?*-]*)\s+([^()\s]+)\)")
WANT = re.compile(r"(^|/)scheme$|^flow/scheme|least-square|gradient/(scheme|method|cell|node)|"
                  r"^time/|unsteady|second-order|bounded|^solve/|pseudo-time-method/(form|method)|"
                  r"warped|high-order-term|nita/|^p-v-coupling", re.I)


def texts(group, path=""):
    for k, v in group.items():
        p = f"{path}/{k}"
        if isinstance(v, h5py.Group):
            yield from texts(v, p)
        else:
            if v.dtype.kind not in "SOU":          # skip numeric arrays (mesh, data)
                continue
            try:
                raw = v[()]
            except Exception:  # noqa: BLE001
                continue
            if isinstance(raw, bytes):
                yield p, raw.decode(errors="replace")
            elif hasattr(raw, "dtype") and raw.dtype.kind in "SO" and raw.size:
                for r in raw.ravel():
                    yield p, (r.decode(errors="replace") if isinstance(r, bytes) else str(r))


for cas in sys.argv[1:]:
    print(f"==== {cas}")
    f = h5py.File(cas, "r")
    root = f["settings"] if "settings" in f else f
    found = {}
    for p, t in texts(root):
        for name, val in PAIR.findall(t):
            if WANT.search(name) and name not in found and len(val) < 40:
                found[name] = val
    for name in sorted(found):
        print(f"  {name:55s} {found[name]}")
