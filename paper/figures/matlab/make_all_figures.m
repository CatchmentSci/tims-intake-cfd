function make_all_figures()
%MAKE_ALL_FIGURES Regenerate every MATLAB-drawn figure of the paper.
%   Refresh the CSV inputs first (from cfd/): python geometry/plot_RU.py; python geometry/plot_adv.py;
%   python geometry/fig_budget.py; python geometry/rig_reading.py;
%   python geometry/sensitivity_bands.py; python geometry/fig_designs_data.py;
%   python geometry/fig_flow_data.py;
%   python geometry/figS1_mesh_data.py
%   Then, from this folder: matlab -batch "make_all_figures"
make_fig_RU();
make_fig_adv();
make_fig_budget();
make_fig_phillips();
% make_fig_bands() dropped: the bands figure is in neither the paper nor the
% supplement. Its CSV now feeds the uncertainty envelope of Fig. 3a (make_fig_RU).
make_fig_designs();
make_fig_flow();
make_fig_S1_mesh();   % Supplemental Figure S1
end
