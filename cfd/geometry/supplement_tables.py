"""CSV tables for the Supplemental Materials, one per table, in results/supplemental/.

Tables S1-S5 of the Supplemental Materials are typeset from them:
- S1_dimensions.csv: dimensions of every design (value, published / catalogue / assumed with
  source, reasoning, plausible range), with tims_geometry.py;
- S2_bore_reynolds.csv: bore Reynolds number for every bore of every design at the six velocities;
- S3_meshes.csv: each mesh of each grid-convergence set: scale s, cells, faces across each bore, R
  at 0.3 m/s, and the medium-mesh error / GCI;
- S4_runs.csv: the run list (the forced-laminar bore runs included), from run_metadata.json, with
  each run's s.d. and drift of R (% of the mean) and mouth-flow imbalance over its averaging window.

Bore Reynolds number Re = 4 Q / (pi d nu), with Q = R U A_in the sampled flow (R from the medium
mesh) and d each bore the flow passes in series; nu = mu / rho = 1.002e-3 / 998.2.
Faces across a bore = d / (0.25 mm x bore scale x s); bore scale 0.5 for McDonald and Foets.

Run from cfd/:  python geometry/supplement_tables.py   (rerun after new runs land)
"""
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, "geometry")
from table_R import metrics  # noqa: E402

RES = Path("results")
OUT = RES / "supplemental"
NU = 1.002e-3 / 998.2
SWEEP = [0.154, 0.3, 0.585, 1.0, 1.25, 1.5]

# design label, results folder, run, inlet bore at the R plane (m), bore scale, bores in series (name, d)
DESIGNS = [
    ("Phillips", "urans_v5", "phillips2000", 0.004, 1.0, [("inlet and outlet tubes", 0.004)]),
    ("Perks", "urans_v5", "perks2014", 0.00794, 1.0, [("tubing", 0.00794), ("panel-mount insert", 0.00635)]),
    ("Perks + meter", "urans_v5", "perks2014_flowmeter", 0.00794, 1.0,
     [("tubing", 0.00794), ("insert", 0.00635), ("meter bore", 0.010)]),
    ("Elliott", "urans_v5", "elliott2017", 0.004, 1.0, [("tubes", 0.004), ("fitting passage", 0.00635)]),
    ("McDonald", "urans_v5bore", "mcdonald2010", 0.002, 0.5, [("inlet tube", 0.002), ("outlet hole", 0.003)]),
    ("Foets", "foets", "foets2020", 0.004, 0.5, [("inlet tube", 0.004), ("outlet tube", 0.002)]),
]

# grid-convergence sets and errors (0.3 m/s) from geometry/mesh_gci.py: cell counts from the meshing
# logs (results/checks/mesh_min_orthogonal_quality.csv), errors recomputed from the summaries. Until
#
# -1.43 from the data (final numbers check).
import mesh_gci  # noqa: E402
_G = {name: mesh_gci.gci(name) for name in mesh_gci.SETS}
MESHES = {name: [(lab, sc, n / 1e6) for lab, sc, n in zip(("coarse", "medium", "fine"), mesh_gci.SETS[name][3], g["cells"])]
          for name, g in _G.items()}
PHILLIPS_NOTE = "Not asymptotic. R rises 0.8 % from coarse to medium and 1.8 % from medium to fine. Halving the bore face size twice on the medium mesh (16, 32 and 64 faces across) raised R to 0.1869 and 0.1884, +1.7 % in all and more than half of the medium-mesh error. The fine mesh with the bore face size halved (50 faces) gives 0.1896. At least 0.6 % of the error (0.1896 against 0.1884) therefore lies outside the bores; the outlet tube ends in the wake, where the background mesh also changes. Phillips and Elliott have the same 16 faces and a similar share of narrow-bore loss at 0.3 m/s (75 and 78 % of the drive), so bore resolution alone does not separate them. The inlet-bore pressure gradient varies by about 1.6 % across the six meshes (114.0-115.8 Pa/m). In the two best-resolved bores (114.4 and 114.6 Pa/m) it is within 1 % of the Hagen-Poiseuille gradient at the computed flow and implies R = 0.1903-0.1906. Best estimate 0.1905. R rose with every refinement, so the finest computed value, 0.1896, is the lower end. In laminar flow along a straight bore the pressure gradient is smallest once the flow is fully developed, so the mesh-independent R cannot exceed the R implied by the converged gradient. That implied R falls as the standard meshes are refined (0.1926, 0.1920 and 0.1916), and every mesh with refined bores gives less (0.1897-0.1906). The upper end is taken as 0.1915, 0.001 above the estimate of 0.1905 and just below the fine mesh's 0.1916. The error is taken against the upper end"
_NOTES = {
    "Perks": "monotone, apparent order {p:.1f} capped at 2",
    "Perks + meter": "monotone, apparent order {p:.1f} capped at 2",
    "Elliott": "monotone, apparent order {p:.1f} capped at 2",
    "McDonald": "monotone, apparent order {p:.1f} capped at 2; s 2.20 / 1.00 / 0.45 for refinement ratio >= 1.29",
    "Foets": "not asymptotic; second-order extrapolation of the two finer meshes; bore halved again +0.65 %",
}
# (medium-mesh error %, GCI_fine % or None, R best estimate, note)
MESH_ERR = {"Phillips": (100 * _G["Phillips"]["e"], None, _G["Phillips"]["R_ext"], PHILLIPS_NOTE)}
MESH_ERR.update({name: (100 * g["e"], None if name == "Foets" else 100 * g["gci_fine"], g["R_ext"],
                        _NOTES[name].format(p=g["p"])) for name, g in _G.items() if name != "Phillips"})

# S1: design, item, value, status (published / catalogue / assumed / derived), source, reasoning, range
DIMENSIONS = [
    ("Phillips", "Chamber", "98 bore x 1000", "published", "Phillips et al. 2000 pp. 2590-2591", "", ""),
    ("Phillips", "Inlet and outlet tubes", "4 bore x 150, 20 inside the chamber", "published", "Phillips et al. 2000; Russell et al. 2000; Smith & Owens 2014",
     "150 read as total length including the 20 inside; no source says so explicitly", ""),
    ("Phillips", "Body o.d.", "110", "assumed", "", "standard UK PVC pipe for a 98 mm bore; Smith & Owens 2014 Fig. 3 draws about 100", "100-114; external flow only (drive and base pressure), outside the 1-D model"),
    ("Phillips", "End caps", "10 thick, flush with body o.d.", "assumed", "", "not reported; Smith & Owens give 8.68 L internal volume against 7.54 L for the tube, which recessed caps would explain", "5-15 thick, flush or recessed; chamber length and external flow only, no bore-path change"),
    ("Phillips", "Tube o.d.", "6", "assumed", "", "standard 6 x 4 mm nylon tubing", "o.d. 6-8, external flow at the tips only; bore 4 +- 0.1 (tolerance) changes q by -8 to +9 % at 0.3 m/s (1-D)"),
    ("Phillips", "Nose", "cone 100 long, 110 base, sealed", "assumed", "'polyethylene funnel', no size given",
     "the tube leaves 150 - 20 - 10 = 120 mm ahead of the cap for cone plus protrusion; 100 mm cone leaves a 20 mm protrusion, as on the Perks samplers (the author) and within 10 mm of a standard 60 deg funnel",
     "cone 90-110 (protrusion 10-30); inlet mouth acts as a pitot tube, so a few hundredths of the dynamic pressure at most"),
    ("Phillips", "Rear end", "blunt cap", "published", "Phillips et al. 2000 Fig. 1", "", ""),
    ("Perks", "Body", "3 in PVC Class C, 80.7 bore x 1000", "catalogue", "builders' parts list; supplier data (88.9 o.d., 4.1 wall)", "paper gives 90 mm chamber, 1000 long", ""),
    ("Perks", "Tubing", "7.94 bore, 12.7 o.d.", "catalogue", "Nalgene 980, 5/16 x 1/2 in (Fisher 10448612)", "paper gives 8 mm tubes", ""),
    ("Perks", "Insert", "6.35 bore, 48.5 long, 19.5 into chamber", "catalogue", "CPC PLC42006 (catalogue 2004 p. 21)", "position assumed: centred, tubing bore ends where insert begins", "insert bore 6.35 +- 0.15 (tolerance): q -4 to +4 % at 0.3 m/s (1-D)"),
    ("Perks", "Tube length", "161.5 tip to insert inner end (113 at 7.94, 48.5 at 6.35)", "derived", "", "", "tubing +-10 each end: q -1 to +1 % at 0.3 m/s (1-D)"),
    ("Perks", "Nose and tail", "funnels: 80 cone, 18 x 37 stem; tubing 20 beyond stem tip", "catalogue", "Nalgene 4256-0414; the author (~2 cm)", "funnel interior sealed (assumed)", ""),
    ("Perks", "End caps", "5 end, 97 o.d. collar x 40", "assumed", "", "3 in solvent-weld cap", "end 3-8, collar 90-100 o.d.; chamber and external flow only"),
    ("Perks + meter", "Meter", "Cynergy3 UF25B: 69 x 38 x 45 body, 10 bore, 3/8 in BSP ports", "catalogue", "UF25B datasheet; BS EN 10226",
     "on the axis (true flow axis ~2 mm off centre, assumed centred for the half model); 20 mm tubing from insert to inlet port (the author); discharges from its own port", ""),
    ("Elliott", "Body", "4 in Sch 40 PVC, 102.3 bore x 1000", "catalogue", "Elliott et al. 2017 App. A.3; ASTM D1785", "their text gives 98 (the Phillips value) and their Fig. 1b 105", ""),
    ("Elliott", "Caps", "Spears 447-040, domed", "catalogue", "Spears sheet; ASTM D2466", "crown and floor as half-ellipses; pipe seated at socket bottom (assumed)", ""),
    ("Elliott", "Tubes", "4 bore (3/8 in nylon bored out) x 150, 22.5 deg tip chamfer", "published", "Elliott et al. 2017 s2.1, A.3", "chamfer read as included angle from photo A.2c", ""),
    ("Elliott", "Inlet fitting", "6.35 passage; tube seated 14", "catalogue", "McMaster 50775K335", "seat depth assumed; tube ends ~30 short of the chamber", "bore 4 +- 0.1 (tolerance) changes q by -7 to +8 % at 0.3 m/s (1-D)"),
    ("Elliott", "Outlet", "elbow; tube down to the cap bottom (~31 exposed)", "assumed", "McMaster 50775K349; the author from Figs. 1b-c", "text's 150 mm outflow tube would reach below the sampler", "exposed 31-46: q 0 to -4 %; all 150 exposed: -24 % at 0.3 m/s (1-D); longer leg not run in CFD"),
    ("McDonald", "Body", "63.5 bore, 228 tube (218 chamber)", "published/derived", "McDonald et al. 2010 p. 227, Fig. 2", "bore from Fig. 2 (63.5; text 6.5 cm); chamber 218 inside assumed 10 mm outlet plug", ""),
    ("McDonald", "Body wall", "3.175 (69.85 o.d.)", "assumed", "", "2.5 x 2.75 in tube", "wall 2.5-3.5; external flow only"),
    ("McDonald", "Inlet", "2 bore, 20 into the body, flush with spout tip; 96 tip to inner end", "published/derived", "McDonald et al. 2010 p. 227, Fig. 2", "length not stated; 96 follows from the assumed funnel; o.d. 4 assumed", "bore 2 +- 0.1 (tolerance) changes q by -17 to +18 % at 0.3 m/s (1-D)"),
    ("McDonald", "Funnel", "hollow, open to body; 60 deg, 1.5 wall, 8 x 20 spout, 10 collar", "assumed", "'funnel fixed over the upstream end'; the author", "Fig. 2 not to scale", "angle 45-75 deg, spout 15-25; approach flow only, outside the 1-D model"),
    ("McDonald", "Outlet", "3 mm hole, 3/4 up the bore", "assumed", "McDonald et al. 2010 p. 228; Fig. 2 (0.71-0.75)", "published as 'at the top'; Fig. 2 not to scale", "3/4 up to top of the cap: +2.6 % in CFD (run E); plug 5-15 thick: q +-1 % (1-D)"),
    ("Foets", "Outlet tube", "2 bore, 4 o.d.; length and position as Phillips", "published/assumed", "Foets et al. 2020 p. 4711", "only the diameter is said to differ; o.d. assumed", "o.d. 3-6, external flow at the outlet tip only; length 150 +- 20: q +13 to -11 %; bore 2 +- 0.1: -17 to +19 % at 0.3 m/s (1-D)"),
    ("Foets", "Everything else", "as Phillips", "published", "Foets et al. 2020 p. 4711", "", ""),
]

RUNS = [  # (purpose, domain, folder, glob)
    ("velocity sweep", "half", "urans_v5", "*_s1.00_sst-gamma"),
    ("velocity sweep", "half", "urans_v5bore", "mcdonald2010_s1.00_sst-gamma"),
    ("velocity sweep", "half", "foets", "foets2020_s1.00_sst-gamma"),
    ("grid convergence", "half", "urans_v5", "*_s0.64_sst-gamma"),
    ("grid convergence", "half", "urans_v5", "*_s1.65_sst-gamma"),
    ("grid convergence", "half", "urans_v5bore", "mcdonald2010_s[02].*"),
    ("grid convergence", "half", "foets", "foets2020_s[01].6*"),
    ("fine mesh at 1.5 m/s", "half", "urans_v5_cont2", "*"),  # the fine runs at 1.5 m/s, 12 s
    ("bore refinement", "half", "rerun_v5bore", "*"),
    ("bore refinement", "half", "rerun_v5bore4", "*"),   # bore faces x0.25, medium
    ("bore refinement", "half", "foets_bore2", "*"),
    ("time step", "half", "rerun_dt0.01", "*"),
    ("time step", "half", "rerun_dt_highU", "*"),  # Elliott at 1.25 m/s, Table 2 note d
    ("time step", "half", "rerun_dt2x", "*"),      # Elliott at 1.25 m/s, 2 dt
    ("inner iterations", "half", "rerun_inner20_full", "*"),   # 20 s run
    ("closure", "half", "rerun_sst", "*"),
    ("closure", "half", "model_rke_v5", "*"),
    ("closure", "half", "sst_highU", "*"),     # plain SST at high U
    ("closure", "half", "rke_highU", "*"),     # realizable k-epsilon at high U
    ("closure, bores forced laminar", "half", "laminar_bore", "*"),  # forced-laminar bores
    ("domain extent", "half, boundaries moved out", "domain_big", "*"),  # domain extent
    ("outlet position", "half", "mcdonald_holetop", "*"),  # McDonald hole at the top
    ("no symmetry plane", "full", "rerun_full_v5", "*"),
    ("no symmetry plane", "full", "full_meter_v5", "*"),
    ("wake perturbation", "full", "full_kick", "*"),   # 5 deg inflow kick for 2 s, then 30 s
    ("flume", "0.3 m flume (half)", "rerun_flume_v5", "*"),
    ("calibration rig", "0.3 m flume + glass tube (half)", "rig_v5", "*"),
    ("ADV comparison (flume)", "Goharrokhi flume (half)", "goharrokhi_flume_v5", "*"),
]


ORDER = ["phillips2000", "perks2014", "perks2014_flowmeter", "elliott2017", "mcdonald2010", "foets2020"]


def design_order(case):
    """Sort key: the paper's design order, then the case name."""
    hits = [i for i, d in enumerate(ORDER) if case.name.startswith(d + "_")]
    return (max(hits, key=lambda i: len(ORDER[i])) if hits else len(ORDER), case.name)


def R_of(folder, run, U, scale="1.00"):
    f = RES / folder / f"{run}_s{scale}_sst-gamma" / f"U{U:.3f}_summary.json"
    return json.load(open(f))["mean"]["R"] if f.exists() else None


def write(name, header, rows):
    OUT.mkdir(exist_ok=True)
    with (OUT / name).open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {OUT / name} ({len(rows)} rows)")


def main():
    write("S1_dimensions.csv", ["design", "item", "value_mm", "status", "source", "reasoning", "range"], DIMENSIONS)

    rows = []
    for name, folder, run, d_in, _, bores in DESIGNS:
        Q = {U: R * U * math.pi * d_in**2 / 4 for U in SWEEP if (R := R_of(folder, run, U)) is not None}
        for bname, d in bores:
            rows.append([name, bname, f"{1000 * d:.2f}"] +
                        [f"{4 * Q[U] / (math.pi * d * NU):.0f}" if U in Q else "" for U in SWEEP])
    write("S2_bore_reynolds.csv", ["design", "bore", "d_mm"] + [f"Re_U{U:g}" for U in SWEEP], rows)

    rows = []
    for name, folder, run, _, bs, bores in DESIGNS:
        err, gci, rinf, note = MESH_ERR[name]
        prev = None
        for lab, s, n in MESHES[name]:
            # effective refinement ratio from the next coarser mesh, (N / N_coarser)^(1/3) (Celik et al. 2008)
            r = f"{(n / prev) ** (1 / 3):.2f}" if prev else ""
            prev = n
            faces = "; ".join(f"{b} {round(1000 * d / (0.25 * bs * float(s)))}" for b, d in bores)
            R = R_of(folder, run, 0.3, s)
            rows.append([name, lab, s, f"{n:.2f}", r, faces, f"{R:.5f}" if R else "",
                         f"{err:+.2f}" if lab == "medium" else "", f"{gci:.2f}" if gci and lab == "fine" else "",
                         f"{rinf:.4f}" if lab == "medium" else "", note if lab == "medium" else ""])
    write("S3_meshes.csv", ["design", "mesh", "s", "cells_M", "r_from_coarser", "faces_across_bores", "R_U0.3",
                            "medium_error_pct", "GCI_fine_pct", "R_best", "note"], rows)

    rows, seen = [], set()
    for purpose, domain, folder, pat in RUNS:
        base = RES / folder
        for case in sorted(base.glob(pat), key=design_order) if base.exists() else []:
            if folder == "urans_v5" and case.name.startswith("mcdonald"):
                continue  # McDonald on v5 (8 faces across the bore) is superseded by v5bore
            for meta in sorted(case.glob("U*_run_metadata.json")):
                summ = meta.with_name(meta.name.replace("run_metadata", "summary"))
                if (folder, case.name, meta.name) in seen or not summ.exists():
                    continue
                seen.add((folder, case.name, meta.name))
                m, s = json.load(open(meta)), json.load(open(summ))
                a = m["args"]
                try:
                    st = metrics(case, s["U"])
                except Exception:            # time series without mouth reports
                    st = None
                stats = [f"{st['sd']:.3f}", f"{st['halves']:+.3f}", f"{st['imb']:.1e}"] if st else ["", "", ""]
                # time-mean mouth imbalance above the 1e-4 acceptance criterion
                flag = "above 1e-4" if st and st["imb"] > 1e-4 else ""
                closure = a.get("turbulence", "") + (", bores laminar" if folder == "laminar_bore" else "")
                rows.append([purpose, case.name, Path(a.get("mesh") or "").name.replace(".msh.h5", ""), domain,
                             closure, f"{s['U']:g}", f"{s['dt']:g}", s.get("inner", a.get("inner", "")),
                             f"{s['sim_time']:g}", f"{s['average_from']:g}", a.get("cores", ""),
                             f"{s['mean']['R']:.5f}"] + stats + [flag, f"{folder}/{case.name}"])
    write("S4_runs.csv", ["purpose", "case", "mesh", "domain", "closure", "U", "dt_s", "inner", "t_end_s",
                          "average_from_s", "cores", "R", "sd_pct", "drift_pct", "imbalance", "imbalance_flag",
                          "results_path"], rows)


if __name__ == "__main__":
    main()
