"""Pressure budget of Perks and Perks + meter across velocity.
Same stations and method as fig_budget.py (0.3 m/s), applied to the production-final slices in
results/mesh_slices/regime/ (the bore_regime slices). Run from cfd/: python geometry/budget_vs_U.py
"""
import sys, numpy as np
sys.path.insert(0, "geometry")
import fig_budget as fb
from pathlib import Path
def make_loader(U):
    def load_slice(stem):
        run = stem.split("_v5")[0]
        z = np.load(Path("results/mesh_slices/regime") / f"regime_{run}_U{U:.3f}_field.npz")
        off, xz = z["off"], z["xz"]; n = np.diff(off)
        return dict(x=np.add.reduceat(xz[:, 0], off[:-1]) / n, z=np.add.reduceat(xz[:, 1], off[:-1]) / n,
                    a=z["area"], u=z["u"], w=z["w"], p=z["p"], q=0.5 * fb.RHO * (z["u"]**2 + z["v"]**2 + z["w"]**2))
    return load_slice
for U in (0.3, 0.585, 1.0, 1.25, 1.5):
    fb.U = U; fb.load_slice = make_loader(U)
    # budget() reads U0.300 files by name: patch Path base via a wrapper
    def budget(label, run, stem, folder, ver, U=U):
        import json, csv, statistics as st
        m = fb.load_slice(stem)
        base = Path("results") / folder / f"{run}_s1.00_sst-gamma" / f"U{U:.3f}"
        summ = json.load(open(f"{base}_summary.json"))
        rows = [r for r in csv.DictReader(open(f"{base}_timeseries.csv")) if float(r["t"]) >= summ["average_from"] - 1e-9]
        p_inf = st.mean(float(r["p_inlet"]) for r in rows); p_up = st.mean(float(r["p_mouth_up"]) for r in rows)
        qinf = 0.5 * fb.RHO * U**2; sts = fb.stations(run, ver)
        v = {k: fb.section(m, **kw) for k, kw in sts.items() if not kw.get("report")}
        cp = lambda p: (p - p_inf) / qinf
        ptA, pA = v["A"]; ptE, pE = v["E"]; pC = v["C"][1]
        b = dict(approach=1 - cp(ptA), exit=cp(ptE) - cp(pE), base=cp(pE))
        if "B" in v:
            b.update(inlet=cp(ptA) - cp(v["B"][0]), expansion=cp(v["B"][0]) - cp(pC), out_entry=cp(pC) - cp(v["D"][0]), out_bore=cp(v["D"][0]) - cp(ptE))
        else:
            b.update(inlet=cp(ptA) - cp(pC), expansion=0.0, out_entry=cp(pC) - cp(ptE), out_bore=0.0)
        b["R"] = summ["mean"]["R"]; b["chk"] = (pA - p_up) / qinf
        return b
    for label, run, stem, folder, ver in fb.DESIGNS:
        if not run.startswith("perks"): continue
        try:
            b = budget(label, run, stem, folder, ver)
        except Exception as e:
            print(U, label, "ERR", e); continue
        drive = 1 - b["approach"] - b["base"]
        inl = b["inlet"] + b["expansion"]; outp = b["out_entry"] + b["out_bore"]
        print(f"U {U:5.3f} {label:14s} R {b['R']:.4f} drive(mouth-to-outlet) {drive:.3f} inlet+exp {inl:.3f} outlet path {outp:.3f} exit KE {b['exit']:.3f} base {b['base']:+.3f} | inlet+exp/R^2 {inl/b['R']**2:.2f} outlet/R^2 {outp/b['R']**2:.2f} exit/R^2 {b['exit']/b['R']**2:.2f} chk {b['chk']:+.4f}")
