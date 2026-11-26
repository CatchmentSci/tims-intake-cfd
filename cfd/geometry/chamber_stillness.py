"""How still is the settling chamber? One definition for every design (Results, chamber flow).

Region: the fluid faces of the symmetry-plane slice inside the chamber cavity, x between the
inner faces of the two end caps (cav0, cav1) and |z| below the body's inner radius. Solid parts
inside the chamber (tube ends, inserts, the flow meter) have no fluid faces and drop out. Reported
for the whole chamber and for its downstream 70 % (the region the flow figures label, which leaves
out the inlet jet's first 30 %).

- still share: area-weighted share of the region with |u| below 1e-3 U;
- bulk velocity: sampled flow over the chamber cross-section, Q / (U A_body) = R (d_in / D_body)^2,
  with d_in the inlet bore at the R plane (`tube_id`) and R the medium-mesh 0.3 m/s value;
- median |u| / U over the region, area weighted.

McDonald: the production run is on the v5bore mesh (R 0.1863); its only slice is from the v5 mesh
(R 0.1825), so the shares come from that slice and the bulk velocity from the production R.

Run from cfd/:  python geometry/chamber_stillness.py
Writes results/chamber_stillness.csv.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "geometry")
from tims_geometry import SPECS, x_stations  # noqa: E402
from plot_flow import load  # noqa: E402

U = 0.3
CASES = {
    "Phillips": ("phillips2000", "phillips2000_v5_s1.00", "urans_v5/phillips2000_s1.00_sst-gamma"),
    "Perks": ("perks2014", "perks2014_v5_s1.00", "urans_v5/perks2014_s1.00_sst-gamma"),
    "Perks + meter": ("perks2014_flowmeter", "perks2014_flowmeter_v5_s1.00",
                      "urans_v5/perks2014_flowmeter_s1.00_sst-gamma"),
    "Elliott": ("elliott2017", "elliott2017_v5_s1.00", "urans_v5/elliott2017_s1.00_sst-gamma"),
    "McDonald": ("mcdonald2010", "mcdonald2010_v5_s1.00", "urans_v5bore/mcdonald2010_s1.00_sst-gamma"),
    "Foets": ("foets2020", "foets2020_v5bore_s1.00", "foets/foets2020_s1.00_sst-gamma"),
}
OUT = Path("results/chamber_stillness.csv")


def wmedian(v, w):
    i = np.argsort(v)
    c = np.cumsum(w[i])
    return float(v[i][np.searchsorted(c, 0.5 * c[-1])])


def main():
    rows = []
    print(f"{'design':14s} {'L mm':>6s} {'D mm':>6s} {'R':>7s} {'bulk/U':>9s}  "
          f"{'still all':>9s} {'still 70%':>9s}  {'med|u| all':>10s} {'med|u| 70%':>10s}")
    for name, (case, sl, run) in CASES.items():
        s = SPECS[case]
        X = {k: v * 1000 for k, v in x_stations(s).items()}
        m = load(f"results/mesh_slices/{sl}_U{U:.3f}_field.npz", U)
        R = json.load(open(f"results/{run}/U{U:.3f}_summary.json"))["mean"]["R"]
        c, a = m["cen"], m["h"] ** 2
        rb = s.body_id * 500
        L = X["cav1"] - X["cav0"]
        ch = (c[:, 0] > X["cav0"]) & (c[:, 0] < X["cav1"]) & (np.abs(c[:, 1]) < rb)
        back = ch & (c[:, 0] > X["cav0"] + 0.3 * L)
        bulk = R * (s.tube_id / s.body_id) ** 2
        r = dict(design=name, chamber_len_mm=L, body_id_mm=s.body_id * 1000, R=R, bulk_over_U=bulk)
        for tag, sel in (("all", ch), ("back70", back)):
            r[f"still_{tag}"] = float((a[sel] * (m["mag"][sel] < 1e-3)).sum() / a[sel].sum())
            r[f"median_{tag}"] = wmedian(m["mag"][sel], a[sel])
        rows.append(r)
        print(f"{name:14s} {L:6.0f} {r['body_id_mm']:6.1f} {R:7.4f} {bulk:9.2e}  "
              f"{100 * r['still_all']:8.0f}% {100 * r['still_back70']:8.0f}%  "
              f"{r['median_all']:10.2e} {r['median_back70']:10.2e}")
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
