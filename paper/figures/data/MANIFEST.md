# Figure data

The CSV files the paper's figures and tables were drawn from, as submitted. Each is written by a
script in `cfd/geometry/` (run from `cfd/`); the repository README maps every file to its figure or
table and gives the run order. Column definitions are in the data record's `README_data.md`.

| File(s) | Written by | Used in |
|---|---|---|
| `fig_designs_outline.csv`, `fig_designs_stations.csv`, `fig_designs_assumed.csv` | `fig_designs_data.py` | Fig. 1, Fig. 2, Fig. S1 |
| `fig_flow_axis.csv`, `fig_flow_inlet_<design>.csv`, `fig_flow_wake_<design>.csv` | `fig_flow_data.py` | Fig. 2 |
| `fig_RU.csv`, `fig_RU_slope.csv`, `fig_q.csv` | `plot_RU.py` | Fig. 3; R values in the text |
| `fig_RU_model_own.csv` | `screening_band.py` | Fig. 3a (1-D model) |
| `fig_bands.csv`, `fig_bands_sources.csv` | `sensitivity_bands.py` | Fig. 3a band; Table 3 |
| `fig_budget.csv` | `fig_budget.py` | Fig. 4 |
| `fig_phillips_cfd.csv`, `fig_phillips_lab.csv`, `fig_phillips_beta.csv` | `rig_reading.py` (bounds from `pipe_bounds.py`) | Fig. 5; Table 4 |
| `fig_adv_profile.csv`, `fig_adv_meas.csv`, `fig_adv_axis.csv` | `plot_adv.py` | Fig. 6 |
| `fig_adv_flume_*.csv` | `goharrokhi_flume_adv.py --figdata` | Fig. 6; Table 5 |
| `figS1_domain.csv`, `figS1_mesh_inlet.csv`, `figS1_mesh_outlet.csv` | `figS1_mesh_data.py` | Fig. S1 |
