"""Split the change in R in an inlet-turbulence run into the part carried by the mouth pressure
difference and the part from extra friction in the inlet bore (turbulent viscosity carried in from
the free stream, a model artefact for a laminar bore).

Total pressure is averaged across the inlet bore (fig_budget.section) at the mouth, 10 and 50 mm in,
and 2 mm short of the inner end. For a laminar bore the loss per unit flow should not change, so
any excess over the baseline is the artefact.
Run from cfd/: python geometry/ti_bore_check.py BASE_SLICE RUN_SLICE BASE_RESULTS RUN_RESULTS [CASE]
  (slices in results/mesh_slices/, results folders holding U0.300_summary.json/_timeseries.csv)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from fig_budget import RHO, section  # noqa: E402

CASE = sys.argv[5] if len(sys.argv) > 5 else "phillips2000"


def load(f):
    d = np.load(f)
    off, xz = d["off"], d["xz"]
    n = np.diff(off)
    return dict(x=np.add.reduceat(xz[:, 0], off[:-1]) / n, z=np.add.reduceat(xz[:, 1], off[:-1]) / n,
                a=d["area"], u=d["u"], w=d["w"], p=d["p"],
                q=0.5 * RHO * (d["u"] ** 2 + d["v"] ** 2 + d["w"] ** 2))


def run_means(folder):
    d = pd.read_csv(Path(folder) / "U0.300_timeseries.csv")
    h = d[d.t > d.t.max() / 2].mean()
    return h.R, h.p_mouth_up - h.p_mouth_down


def main():
    base_sl, run_sl, base_res, run_res = sys.argv[1:5]
    rep = json.loads(Path(f"geometry/out/v5/scale_1.00/{CASE}_zones.json").read_text())
    st, r = rep["stations_m"], rep["spec"]["tube_id"] / 2
    xs = [st["mouth_up"], st["in_tip"] + 0.010, st["in_tip"] + 0.050, st["in_inner"] - 0.002]
    loss = {}
    for k, f in (("base", base_sl), ("run", run_sl)):
        pt = [section(load(f), x, r)[0] for x in xs]
        loss[k] = -np.diff(pt)
    (Rb, dpb), (Rr, dpr) = run_means(base_res), run_means(run_res)
    q = Rr / Rb
    excess = loss["run"] / loss["base"] / q - 1
    tot = loss["run"].sum() / loss["base"].sum() / q - 1
    print(f"R {Rb:.6f} -> {Rr:.6f} ({100 * (q - 1):+.2f} %); mouth dp {dpb:.3f} -> {dpr:.3f} Pa "
          f"({100 * (dpr / dpb - 1):+.2f} %)")
    print("inlet-bore loss (Pa), mouth-10 / 10-50 / 50-end mm:",
          np.round(loss["base"], 3), "->", np.round(loss["run"], 3))
    print("excess loss per unit flow:", np.round(100 * excess, 1), f"% ; whole bore {100 * tot:+.2f} %")
    # extra resistance as a share of the drive, turned into R with the flume-derived response
    share = (loss["run"].sum() - loss["base"].sum() * q) / dpr
    print(f"extra bore loss = {100 * share:.2f} % of the mouth dp")


if __name__ == "__main__":
    main()
