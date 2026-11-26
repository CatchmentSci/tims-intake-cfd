function make_fig_adv()
%MAKE_FIG_ADV Figure fig:adv: the approach flow ahead of the Phillips inlet against
%   the ADV measurements of Goharrokhi et al. (2019).
%   (a) CFD vertical profile 5 mm ahead of the tip in their flume at 0.15/0.30/0.60 m/s (over
%   their section-mean U) and the digitised ADV profiles; (b) flume CFD streamwise velocity
%   along the inlet axis at the same three velocities, with each bore mean RU. Colour gives the
%   flume velocity in both panels.
%   Input: ../data/fig_adv_meas.csv (cfd/geometry/plot_adv.py) and
%   ../data/fig_adv_flume_{profile,axis,F}.csv (cfd/geometry/goharrokhi_flume_adv.py).
%   Output: ../fig_adv.pdf and ../fig_adv.png.  Tested with MATLAB R2024b.

here = fileparts(mfilename("fullpath"));
dataDir = fullfile(here, "..", "data");
cfg = tims_fig_config();

M = readtable(fullfile(dataDir, "fig_adv_meas.csv"));
F = readtable(fullfile(dataDir, "fig_adv_flume_profile.csv"));
FA = readtable(fullfile(dataDir, "fig_adv_flume_axis.csv"));
FR = readtable(fullfile(dataDir, "fig_adv_flume_F.csv"));

% Three discrete, dark colours, one per flume velocity (a continuous colour bar suggested
% velocities that were not run, and the 0.3 m/s series was too pale to see)
% Okabe-Ito blue, bluish green and vermillion: distinct for colour-blind readers, and none close to
% the black of the line-style key
velCols = [0 114 178; 0 158 115; 213 94 0] ./ 255;
velList = [0.15 0.30 0.60];
colourOf = @(U) velCols(find(abs(velList - U) == min(abs(velList - U)), 1), :);
% marker shape per velocity as well, so the figure reads in black and white (ASCE):
% filled for the ADV, open on the CFD lines
velMarks = ["o" "s" "^"];
markOf = @(U) velMarks(find(abs(velList - U) == min(abs(velList - U)), 1));
someOf = @(n) unique(round(linspace(1, n, 5)));

fig = figure("Color", "w", "Units", "inches", "Visible", "off", ...
    "Position", [1 1 cfg.pageWidth_in 2.95], "PaperPositionMode", "auto");
layout = tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");

% (a) profile 5 mm ahead of the tip
axP = nexttile(layout, 1);
tims_format_axes(axP, cfg);
fill(axP, [0.7 1.1 1.1 0.7], [-5 -5 2.5 2.5], cfg.shade, "EdgeColor", "none", ...
    "HandleVisibility", "off");          % the axis cells behind the reported inlet velocity
text(axP, 1.095, 2.2, "axis cells", "HorizontalAlignment", "right", ...
    "VerticalAlignment", "top", "FontSize", cfg.annotationFontSize, "Color", cfg.mutedText);
plot(axP, [0.7 1.1], [0 0], ":", "Color", cfg.errorColour, "LineWidth", 0.6, ...
    "HandleVisibility", "off");
% CFD in the flume of Goharrokhi et al. (bed, walls, developed inflow), over their section-mean
% U, one line per flume velocity in the colour of its ADV series. The open-domain profile is
% not shown: the ADV was not measured in those conditions.
for Uf = unique(F.U)'
    f = F(abs(F.U - Uf) < 1e-9, :);
    plot(axP, f.u_cfd, f.z_mm, "-", "Color", colourOf(Uf), "LineWidth", 1.2, ...
        "Marker", markOf(Uf), "MarkerIndices", someOf(height(f)), "MarkerSize", 4.5, ...
        "MarkerFaceColor", "w", "HandleVisibility", "off");
end
% Legend: one entry per source and velocity, ADV dots on the top row and CFD
% lines below, one column per velocity
hAdvV = gobjects(1, numel(velList));
hCfdV = gobjects(1, numel(velList));
% Superscripts in 8 pt text: \fontsize{10} before the ^ makes them 8 pt as printed (MATLAB sets
% them at 0.8 of the current size; no printed text below 8 pt).
for k = 1:numel(velList)
    cmS = round(100 * velList(k));
    hAdvV(k) = scatter(axP, NaN, NaN, cfg.markerArea, velCols(k, :), "filled", ...
        "Marker", velMarks(k), "MarkerEdgeColor", cfg.markerEdge, "LineWidth", cfg.markerEdgeWidth, ...
        "DisplayName", sprintf("ADV, %d cm s\\fontsize{10}^{-1}", cmS));
    hCfdV(k) = plot(axP, NaN, NaN, "-", "Color", velCols(k, :), "LineWidth", 1.2, ...
        "Marker", velMarks(k), "MarkerSize", 4.5, "MarkerFaceColor", "w", ...
        "DisplayName", sprintf("CFD, %d cm s\\fontsize{10}^{-1}", cmS));
end
series = unique(M.flume_cm_s)';
for cm = series
    m = M(M.flume_cm_s == cm, :);
    scatter(axP, m.u_over_U, m.z_mm, cfg.markerArea, colourOf(cm / 100), "filled", ...
        "Marker", markOf(cm / 100), "MarkerEdgeColor", cfg.markerEdge, "LineWidth", cfg.markerEdgeWidth, ...
        "HandleVisibility", "off");
end
xlim(axP, [0.7 1.1]);
ylim(axP, [-12 12]);
xlabel(axP, "{\itu}/{\itU}", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
ylabel(axP, "{\itz} (mm)", "Interpreter", "tex", ...
    "FontSize", cfg.labelFontSize);
% column-major fill of 2 rows x 3 columns: ADV over CFD for each velocity
lgA = legend(axP, reshape([hAdvV; hCfdV], 1, []), "NumColumns", 3, "Box", "off", ...
    "FontSize", cfg.legendFontSize);
lgA.Layout.Tile = "south";          % below both panels, clear of the data
tims_panel_label(axP, "(a)", cfg, [0.03 0.86]);   % lowered with (b), clear of its lines

% (b) approach along the inlet axis
axL = nexttile(layout, 2);
tims_format_axes(axL, cfg);
fill(axL, [0 20 20 0], [0 0 1.15 1.15], cfg.shade, "EdgeColor", "none", ...
    "HandleVisibility", "off");
text(axL, 10, 1.1, {"inside", "the tube"}, "HorizontalAlignment", "center", ...
    "VerticalAlignment", "top", "FontSize", cfg.annotationFontSize, "Color", cfg.mutedText);
plot(axL, [-5 -5], [0 1.15], ":", "Color", cfg.errorColour, "LineWidth", 0.8, ...
    "HandleVisibility", "off");
text(axL, -6.5, 0.12, {"ADV cells", "5 mm ahead"}, "HorizontalAlignment", "right", ...
    "VerticalAlignment", "bottom", "FontSize", cfg.annotationFontSize, "Color", cfg.mutedText);
% flume CFD at the three ADV velocities, over their section-mean U as in (a);
% dashed: bore mean R U from each flume run
for U = FR.U'
    a = FA(abs(FA.U - U) < 1e-9 & FA.s_mm <= 0, :);   % approach only; inside the tube the bore mean is shown
    c = colourOf(U);
    plot(axL, a.s_mm, a.u_over_U, "-", "Color", c, "LineWidth", 1.3, "Marker", markOf(U), ...
        "MarkerIndices", someOf(height(a)), "MarkerSize", 4.5, "MarkerFaceColor", "w");
    plot(axL, [0.5 10 19.5], FR.R(abs(FR.U - U) < 1e-9) * [1 1 1], "--", "Color", c, ...
        "LineWidth", 1.3, "Marker", markOf(U), "MarkerIndices", 2, "MarkerSize", 4.5, ...
        "MarkerFaceColor", "w");
end
xlim(axL, [-60 20]);
ylim(axL, [0 1.15]);
yticks(axL, 0:0.2:1);
xlabel(axL, "{\itx} (mm)", "Interpreter", "tex", ...
    "FontSize", cfg.labelFontSize);
ylabel(axL, "{\itu}/{\itU}", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
tims_panel_label(axL, "(b)", cfg, [0.03 0.86]);   % below the approach lines near the top


tims_export(fig, fullfile(here, "..", "fig_adv"), cfg);
close(fig);
end
