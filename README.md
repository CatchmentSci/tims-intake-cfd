# tims-intake-cfd

Code for the CFD study of the intake hydrodynamics of time-integrated mass-flux (TIMS) suspended
sediment samplers:

> Perks, M. T. (submitted). Intake hydrodynamics of time-integrated suspended sediment samplers:
> computational fluid dynamics of six sampler configurations compared with laboratory calibrations.
> *Journal of Hydraulic Engineering* (ASCE).

The code builds the geometry of six published sampler configurations, meshes them and solves the
flow around and through each in Ansys Fluent (URANS, SST k-omega with gamma intermittency). It then
reduces the runs to the paper's numbers, tables and figure data, and draws the figures in MATLAB.
This repository is the state that produced the submitted manuscript.

The run results (summaries, time series, run metadata, symmetry-plane slices), the run list and the
supplementary information are in the data record:

> Perks, M. T. (2026). Intake hydrodynamics of time-integrated suspended sediment samplers: data and
> supplementary information. Zenodo. https://doi.org/10.5281/zenodo.24368491

Meshes and Fluent case and data files are not archived (several GB per run); the scripts here
regenerate them.

## Layout

```
cfd/geometry/            geometry, post-processing, checks and figure-data scripts (Python)
cfd/fluent/              Fluent driver (PyFluent), field extraction and Slurm templates
cfd/results/             run results: unpack the Zenodo record here (git-ignored)
paper/figures/data/      the figure data CSVs as plotted in the paper (tracked)
paper/figures/matlab/    MATLAB figure scripts and the shared figure style (tims_*.m)
paper/supplemental/      typesets the supplement's Tables S1-S4 from the CSVs
```

Run the Python scripts from `cfd/` (for example `python geometry/sensitivity_bands.py`) and the
MATLAB scripts from `paper/figures/matlab/`.

## Requirements

| Software | Version used | Needed for |
|---|---|---|
| Python | 3.11 (post-processing), 3.12 (HPC runs) | everything in `cfd/` |
| numpy, scipy, pandas, matplotlib, h5py, Pillow, pypdf | see `requirements.txt` | post-processing |
| gmsh (Python API) | 4.15 | geometry only (`tims_geometry.py`) |
| Ansys Fluent and Fluent Meshing | 2025 R1, with a licence | meshing and solving only |
| PyFluent (`ansys-fluent-core`) | 0.37.2 | meshing and solving only |
| MATLAB | R2024b | drawing the figures only |

`pip install -r requirements.txt` installs the Python packages. The post-processing needs neither
Fluent nor gmsh: with the Zenodo record unpacked at the repository root, every number, table and
figure-data file below is reproduced from the archived results.

## Reproducing the paper from the archived results

1. Unpack the Zenodo record into the repository root. It holds `cfd/results/` and the zone files
   `cfd/geometry/out/v5/scale_1.00/*_zones.json` and `cfd/geometry/out/v5bore/scale_1.00/*_zones.json`.
2. From `cfd/`, run in this order (each script's docstring says what it reads and writes):

   ```
   python geometry/mesh_gci.py              # Table 2 mesh errors and GCI
   python geometry/supplement_tables.py     # Tables S1-S4 CSVs (cfd/results/supplemental/)
   python geometry/plot_RU.py               # fig_RU.csv, fig_RU_slope.csv, fig_q.csv
   python geometry/fig_budget.py            # fig_budget.csv (Fig. 4)
   python geometry/screening_band.py        # fig_RU_model_own.csv (1-D model, Fig. 3a)
   python geometry/sensitivity_bands.py     # fig_bands.csv, fig_bands_sources.csv (Table 3, Fig. 3a band)
   python geometry/rig_reading.py           # fig_phillips_*.csv (Table 4, Fig. 5)
   python geometry/plot_adv.py              # fig_adv_*.csv, open domain (Fig. 6)
   python geometry/goharrokhi_flume_adv.py --figdata ...   # flume series (Table 5, Fig. 6); see its docstring
   python geometry/fig_designs_data.py      # Fig. 1
   python geometry/fig_flow_data.py         # Fig. 2
   python geometry/figS1_mesh_data.py       # Fig. S1
   python geometry/table_R.py               # Table S7
   python geometry/pipe_bounds.py           # pipe-loss bounds, flume drive rise
   python geometry/number_checks.py         # wake, base pressure, kick test, shoulder speed, y+, bubble
   python geometry/bore_regime.py           # bore flow regime (nu_t/nu, pressure drop)
   python geometry/geometry_sensitivity_1d.py   # geometry and bore-tolerance terms
   ```

3. Draw the figures from `paper/figures/matlab/`: `matlab -batch "make_all_figures"`.
4. Typeset the supplement tables: `python paper/supplemental/make_supp_tables.py`.

The figure data CSVs in `paper/figures/data/` are the files the submitted figures were drawn from.
Re-running the scripts on the archived results reproduces them exactly.

## Paper items and the scripts that make them

| Paper item | Data | Script (`cfd/geometry/` unless noted) | Figure drawn by |
|---|---|---|---|
| Table 1, design dimensions | `tims_geometry.py` (`SPECS`, with sources); Table S1 | `supplement_tables.py` | |
| Table 2, discretisation and sensitivity | `results/supplemental/S3_meshes.csv`, `S4_runs.csv` | `mesh_gci.py`, `supplement_tables.py` | |
| Table 3, sources of uncertainty | `fig_bands_sources.csv`, `fig_bands.csv` | `sensitivity_bands.py` (mesh, time step, closure, forced-laminar bores, domain); wake row `number_checks.py`; geometry rows `geometry_sensitivity_1d.py`; turbulence row `turbulence_bias.py` | |
| Table 4, Phillips calibration | `fig_phillips_cfd.csv`, `fig_phillips_lab.csv`, `fig_phillips_beta.csv` | `rig_reading.py`, `pipe_bounds.py` | |
| Table 5, ADV comparison | `fig_adv_flume_*.csv` | `goharrokhi_flume_adv.py`, `plot_adv.py` | |
| Fig. 1, designs | `fig_designs_*.csv` | `fig_designs_data.py` | `make_fig_designs.m` |
| Fig. 2, flow fields | `fig_flow_*.csv` | `fig_flow_data.py` | `make_fig_flow.m` |
| Fig. 3, R(U), slope, 1-D model, uncertainty band | `fig_RU.csv`, `fig_RU_slope.csv`, `fig_RU_model_own.csv`, `fig_bands.csv` | `plot_RU.py`, `screening_band.py`, `sensitivity_bands.py` | `make_fig_RU.m` |
| Fig. 4, pressure budget | `fig_budget.csv` | `fig_budget.py` | `make_fig_budget.m` |
| Fig. 5, Phillips design against the laboratory relations | `fig_phillips_cfd.csv`, `fig_phillips_lab.csv` | `rig_reading.py`, `pipe_bounds.py` | `make_fig_phillips.m` |
| Fig. 6, approach flow against the ADV | `fig_adv_*.csv` | `plot_adv.py`, `goharrokhi_flume_adv.py` | `make_fig_adv.m` |
| Tables S1-S4 (dimensions, bore Reynolds numbers, meshes, run list) | `results/supplemental/*.csv` | `supplement_tables.py`; typeset by `paper/supplemental/make_supp_tables.py` | |
| Table S6, wall y+ | `results/checks/yplus_*.csv` | `cfd/fluent/wall_yplus.py` (reads case and data files) | |
| Table S7, R of every run | run summaries and time series | `table_R.py` | |
| Fig. S1, domain and mesh | `figS1_*.csv` | `figS1_mesh_data.py` | `make_fig_S1_mesh.m` |

### Headline numbers

| Number in the paper | Where it comes from |
|---|---|
| R = 0.19-0.27 at 0.3 m/s for five configurations, 0.025 for the sixth | `fig_RU.csv` (`plot_RU.py`); Table S7 |
| Numerical and closure uncertainty -2.7 to +4.0 % at 0.3 m/s, -8.7 to +5.2 % from 1.0 m/s | `sensitivity_bands.py` (prints the totals; `fig_bands.csv`) |
| Medium-mesh errors and GCI (Table 2) | `mesh_gci.py` |
| Bore tolerance, about +-8 % in intake flow per 0.1 mm of a 4 mm bore | `geometry_sensitivity_1d.py` |
| Wake (up to about +6 %), base pressure, inflow-kick test, shoulder speed, flume bed y+, bubble rise speed | `number_checks.py` |
| Pipe-loss bounds and the 19-21 % flume drive rise | `pipe_bounds.py` |
| Forced-laminar bores: change in R and the closure high side | `sensitivity_bands.py` (reads `results/laminar_bore/`) |
| Bore flow regime (nu_t/nu, pressure drop against laminar and turbulent) | `bore_regime.py` |
| Phillips mesh-independent R (0.1905, bracket 0.1896-0.1915) | `bore_hp.py`; Table S3 note |
| Courant numbers | `courant_slices.py` |
| Ambient turbulence estimate (Supplement S1) | `turbulence_bias.py`, `ti_bore_check.py`, `arriving_ti.py` |

Other checks: `check_geometry.py` (closed boundary, zone areas, chamber volume), `chamber_stillness.py`,
`rig_slice_check.py`, `budget_vs_U.py` (budget at every velocity, Perks and meter),
`fluent/locate_fast_cells.py`, `fluent/connectivity.py`, `fluent/body_flux_check.py`,
`fluent/print_methods.py`.

The digitised laboratory data (`results/phillips2000_fig2_digitised.csv` and
`results/goharrokhi2019_fig3_digitised.csv`) are in the data record. `digitise_phillips_fig2.py` and
`digitise_goharrokhi_fig3.py` made them from the published figures, which are not included here.

## Running the CFD

1. **Geometry.** Build the boundary mesh of a design at mesh scale `s`, then check it:

   ```
   cd cfd
   python geometry/tims_geometry.py phillips2000 --size-scale 1.00 --out geometry/out/v5/scale_1.00
   python geometry/check_geometry.py geometry/out/v5/scale_1.00/phillips2000
   ```

   Designs: `phillips2000`, `perks2014`, `perks2014_flowmeter`, `elliott2017`, `mcdonald2010`,
   `foets2020`; the calibration rig `phillips2000_rig`; the Goharrokhi flume sampler `goharrokhi2019`;
   the outlet-position test `mcdonald2010_holetop`. Domains (`--domain`): `open` (default),
   `open_full` (no symmetry plane), `open_big` (domain-extent test), `flume_phillips` and
   `flume_goharrokhi`. Scales: 1.65, 1.00 and 0.64 (McDonald 2.20, 1.00 and 0.45). McDonald and
   Foets were meshed with their bores at half the face size (`TIMS_BORE_SCALE=0.5`).

2. **Site settings.** Copy `cfd/fluent/site.example.sh` to `cfd/fluent/site.sh` and set the modules,
   Python environment and Fluent paths of your cluster. Pass `--account` and `--partition` to
   `sbatch` as your site requires.

3. **Mesh and solve** (from `cfd/fluent/`). Every run follows one schedule: 1500 steady iterations,
   a ramp of 100 time steps at dt/20, then URANS with 10 inner iterations per step, with R averaged
   over the second half of the run. The time step and run length scale with 1/U (0.02 s and 20 s at
   0.3 m/s); the supplement's run list gives every run.

   ```
   sbatch --export=ALL,STAGE=mesh,CASE=phillips2000,SCALE=1.00,GEOVER=v5 tims_job.slurm
   sbatch --export=ALL,STAGE=transient,CASE=phillips2000,SCALE=1.00,GEOVER=v5,VEL=0.300,TAG=urans,ARGS="--dt 0.02 --sim-time 20" tims_job.slurm
   ```

   `run_tims_case_v3.py` is the driver of the production runs. `run_tims_case_v3g.py` is the same
   driver with the forced-laminar bore options (`--laminar-bore`, `--laminar-zone`, `--bore-volume`,
   zones from `geometry/laminar_zones.py`) and an overwrite guard; it ran the forced-laminar runs
   (set `DRIVER=run_tims_case_v3g.py`). `--kick SECONDS DEG` turns the inflow for the wake-symmetry
   test. `screen.slurm` meshes a case and runs a short screen before a production run;
   `tims_continue.slurm` continues a run from its last checkpoint; `precursor.slurm` makes the
   developed inflow for the Goharrokhi flume. `python run_tims_case_v3.py -h` lists every option.

4. **Field slices.** Extract the symmetry-plane slice of a finished run (reads the HDF5 files, no
   Fluent licence needed):

   ```
   python field_slice.py U0.300_final.cas.h5 U0.300_final.dat.h5 phillips2000_v5_s1.00_U0.300_field.npz
   ```

   `wall_yplus.py` gives the first-cell y+ of every wall zone in the same way.

## Citation

See `CITATION.cff`. Please cite the paper and the data record.

## Licence

The code is released under the MIT licence (see `LICENSE`). The data record is released under
Creative Commons Attribution 4.0 International (CC BY 4.0).
