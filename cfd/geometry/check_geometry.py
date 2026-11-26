"""Independent checks of a generated TIMS boundary mesh.

1. Re-parses the Fluent ASCII .msh (hex indices) without using gmsh.
2. Checks that the external boundary (all non-interior zones) is a closed 2-manifold.
3. Compares zone areas and the chamber volume with analytical values from the spec.
Usage: python check_geometry.py out/phillips2000
"""
import json
import math
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np


def read_fluent_boundary_mesh(path):
    text = Path(path).read_text()
    nodes = None
    zones, names = {}, {}
    for m in re.finditer(r"\(10 \(1 1 ([0-9a-f]+) 1 3\)\(\n(.*?)\)\)", text, re.S):
        nodes = np.array(m.group(2).split(), dtype=float).reshape(-1, 3)
    for m in re.finditer(r"\(13 \(([0-9a-f]+) ([0-9a-f]+) ([0-9a-f]+) ([0-9a-f]+) 3\)\(\n(.*?)\)\)",
                         text, re.S):
        zid = int(m.group(1), 16)
        rows = np.array([int(v, 16) for v in m.group(5).split()]).reshape(-1, 5)
        if len(rows) != int(m.group(3), 16) - int(m.group(2), 16) + 1:
            raise ValueError(f"zone {zid}: face count does not match header")
        zones[zid] = rows[:, :3] - 1
    for m in re.finditer(r"\(45 \((\d+) (\S+) (\S+)\)\(\)\)", text):
        names[int(m.group(1))] = (m.group(3), m.group(2))
    return nodes, zones, names


def tri_area(nodes, tri):
    p = nodes[tri]
    return 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1).sum()


def main(stem):
    stem = Path(stem)
    nodes, zones, names = read_fluent_boundary_mesh(stem.with_suffix(".msh"))
    meta = json.loads(stem.with_name(stem.name + "_zones.json").read_text())
    s, xs = meta["spec"], meta["stations_m"]
    print(f"nodes {len(nodes)}, zones {len(zones)}")

    internal = {n for n, v in meta["zones"].items() if v["fluent_type"] == "interior"}
    ext = np.vstack([t for z, t in zones.items() if names[z][0] not in internal])
    edges = Counter()
    for a, b, c in ext:
        for e in ((a, b), (b, c), (c, a)):
            edges[tuple(sorted(e))] += 1
    bad = sum(1 for v in edges.values() if v != 2)
    print(f"external boundary: {len(ext)} faces, {len(edges)} edges, non-manifold/open edges: {bad}")

    ri, ro, rb, rbo = s["tube_id"] / 2, s["tube_od"] / 2, s["body_id"] / 2, s["body_od"] / 2
    oi, oo = (s.get("outlet_id") or s["tube_id"]) / 2, (s.get("outlet_od") or s["tube_od"]) / 2
    exp = {
        "wall_upstream_bore": math.pi * ri * s["inlet_len"],
        "wall_downstream_bore": math.pi * oi * s["outlet_len"],
        "interior_mouth_upstream": math.pi * ri**2 / 2,
        "interior_mouth_downstream": math.pi * oi**2 / 2,
        "interior_body_mid": math.pi * rb**2 / 2,
        "wall_body": math.pi * rbo * (xs["rear"] - 0.0),
        "wall_cavity": math.pi * rb * s["body_len"] + math.pi * (2 * rb**2 - ro**2 - oo**2) / 2,
        "wall_cone": math.pi * (rbo + ro) * math.hypot(s["cone_len"], rbo - ro) / 2,
        # outer surface less the length inside the rear cap, plus the two end annuli
        "wall_downstream_tube": math.pi * oo * (s["outlet_len"] - s["cap_t"]) + math.pi * (oo**2 - oi**2),
        "wall_rear": math.pi * (rbo**2 - oo**2) / 2,
    }
    print(f"{'zone':18s} {'mesh area':>12s} {'analytic':>12s} {'error %':>8s}")
    for zid, (name, _) in names.items():
        if name in exp:
            a = tri_area(nodes, zones[zid])
            print(f"{name:18s} {a:12.5e} {exp[name]:12.5e} {100 * (a / exp[name] - 1):8.3f}")

    vol = math.pi * rb**2 * s["body_len"] - math.pi * ((ro**2 - ri**2) * s["inlet_in"]
                                                      + (oo**2 - oi**2) * s["outlet_in"])
    ratio = (s["body_id"] / s["tube_id"])**2
    print(f"chamber volume {vol * 1e3:.3f} L (published 7.54 L, Smith & Owens 8.68 L); "
          f"area ratio {ratio:.1f} (published ~600)")
    print(f"inlet tube protrudes {(xs['cone_tip'] - xs['in_tip']) * 1e3:.1f} mm beyond the cone tip; "
          f"outlet {(xs['out_tip'] - xs['rear']) * 1e3:.1f} mm beyond the rear cap")


if __name__ == "__main__":
    main(sys.argv[1])
