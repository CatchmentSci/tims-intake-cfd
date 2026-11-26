function make_fig_phillips()
%MAKE_FIG_PHILLIPS Figure fig:phillips: intake ratio of the Phillips sampler against
%   the laboratory calibrations.
%   CFD: open domain, confined 0.3 m flume, flume with the calibration glass tube fitted (true
%   intake), and what the bubble method would read from that rig for a bubble at 1.25-1.4 x
%   the bulk speed, converted with the area ratio x3.06 (the factor applied). Laboratory: the data and fit of Phillips et al.
%   (2000) as published, the ADV relation of Goharrokhi et al. (2019), and the pipe-loss upper bound.
%   Input: ../data/fig_phillips_cfd.csv and fig_phillips_lab.csv, written by
%   cfd/geometry/rig_reading.py (run it from cfd/ first).
%   Output: ../fig_phillips.pdf and ../fig_phillips.png.  Tested with MATLAB R2024b.

here = fileparts(mfilename("fullpath"));
dataDir = fullfile(here, "..", "data");
cfg = tims_fig_config();
C = readtable(fullfile(dataDir, "fig_phillips_cfd.csv"), "TextType", "string");
L = readtable(fullfile(dataDir, "fig_phillips_lab.csv"), "TextType", "string");

% Okabe-Ito colours, as Fig. 6: CFD blue, Phillips et al.'s published data and fit
% vermillion, the rig's bubble readings bluish green (the variable names keep their old roles)
blue = [0 114 178] ./ 255;
red = [213 94 0] ./ 255;
salmon = [0 158 115] ./ 255;
grey = [0.45 0.45 0.45];

fig = figure("Color", "w", "Units", "inches", "Visible", "off", ...
    "Position", [1 1 cfg.pageWidth_in 3.0], "PaperPositionMode", "auto");
ax = axes(fig);
tims_format_axes(ax, cfg);
set(ax, "XScale", "log", "YScale", "log");
handles = gobjects(0);

% laboratory relations (lines), under the CFD
series = L(L.series == "adv_fit", :);
handles(end+1) = plot(ax, series.U, series.R, "-.", "Color", grey, "LineWidth", cfg.dataLineWidth, ...
    "DisplayName", "Goharrokhi et al. (2019), ADV approach-velocity ratio");
% The open-water pipe-loss bound (pipe_bound) lies 1-4 % below the rig bound and hid under it; only
% the rig bound is drawn, and the caption gives the offset.
% The same bound for the rig: glass bore added and the drive of the 0.3 m flume (+21 %). Values from
% geometry/pipe_bounds.py (Table 4).
rigBound = [0.154 0.150; 0.30 0.276; 0.585 0.458];
handles(end+1) = plot(ax, rigBound(:, 1), rigBound(:, 2), ":", "Color", cfg.refColour, ...
    "LineWidth", 1.2, "Marker", "x", "MarkerSize", 6, ...
    "DisplayName", "pipe-loss upper bound, rig in flume");
% Goharrokhi's area-ratio correction of Eq. 1 (x 3.06/1.75) is not drawn: it assumed 1.75 was applied.
series = L(L.series == "phillips_fit", :);
handles(end+1) = plot(ax, series.U, series.R, "-", "Color", red, "LineWidth", cfg.dataLineWidth, ...
    "DisplayName", "Phillips et al. (2000), Eq. 1");
series = L(L.series == "phillips_data", :);
handles(end+1) = scatter(ax, series.U, series.R, cfg.markerArea, red, "filled", "Marker", "d", ...
    "MarkerEdgeColor", cfg.markerEdge, "LineWidth", cfg.markerEdgeWidth, ...
    "DisplayName", "Phillips et al. (2000), data as published");

% CFD
spec = { ...
    "open", "CFD, open domain", [86 180 233] ./ 255, "o", "-", false, 0.965; ...
    "flume", "CFD, 0.3 m flume", blue, "s", "-", true, 1.0; ...
    "rig", "CFD, flume with glass tube (true intake)", [0 63 110] ./ 255, "^", "--", true, 1.035};
% three blues, light to dark: open domain, flume, flume with glass tube (analysis spec)
% The open-domain and rig markers are offset by -3.5 % and +3.5 % in U so
% that they no longer hide each other or the lowest published point (the caption says so).
% The centreline and bulk-speed readings are no longer drawn: both are cases the
% paper rules out. reading_bulk is still read below, as the base of the bore-filling band.

for k = 1:size(spec, 1)
    series = C(C.series == spec{k, 1}, :);
    if isempty(series)
        continue
    end
    face = spec{k, 3};
    if spec{k, 6}
        face = [1 1 1];
    end
    Uo = series.U * spec{k, 7};
    plot(ax, Uo, series.R, spec{k, 5}, "Color", spec{k, 3}, "LineWidth", cfg.dataLineWidth, ...
        "HandleVisibility", "off");
    scatter(ax, Uo, series.R, cfg.markerArea, face, "filled", "Marker", spec{k, 4}, ...
        "MarkerEdgeColor", cfg.markerEdge, "LineWidth", cfg.markerEdgeWidth, "HandleVisibility", "off");
    handles(end+1) = plot(ax, NaN, NaN, spec{k, 5}, "Color", spec{k, 3}, ...
        "LineWidth", cfg.dataLineWidth, "Marker", spec{k, 4}, "MarkerSize", 6, ...
        "MarkerFaceColor", face, "MarkerEdgeColor", cfg.markerEdge, ...
        "DisplayName", spec{k, 2}); %#ok<AGROW>
end

% bubble that fills or nearly fills the bore: beta = 1.25-1.4 times the bulk speed (Moran et al. 2021 mass balance
% at the rig's Bond and capillary numbers), drawn as a band with beta = 1.3 on it.
% Converted with the area ratio (7/4)^2 = 3.0625, the factor applied in the calibration: J. M. Phillips,
% pers. comm., the printed 1.75 was probably an error.
% reading_bulk in the CSV is already converted with x3.06 (rig_reading.py).
% The x1.75 band is no longer drawn; the text gives it as evidence for 3.06.
series = C(C.series == "reading_bulk", :);
fill(ax, [series.U; flipud(series.U)], [1.25 * series.R; flipud(1.4 * series.R)], salmon, ...
    "FaceAlpha", 0.25, "EdgeColor", "none", "HandleVisibility", "off");
plot(ax, series.U, 1.3 * series.R, "-", "Color", salmon, "LineWidth", cfg.dataLineWidth, ...
    "HandleVisibility", "off");
scatter(ax, series.U, 1.3 * series.R, cfg.markerArea, salmon, "filled", "Marker", "v", ...
    "MarkerEdgeColor", cfg.markerEdge, "LineWidth", cfg.markerEdgeWidth, "HandleVisibility", "off");
handles(end+1) = plot(ax, NaN, NaN, "-", "Color", salmon, "LineWidth", cfg.dataLineWidth, ...
    "Marker", "v", "MarkerSize", 6, "MarkerFaceColor", salmon, "MarkerEdgeColor", cfg.markerEdge, ...
    "DisplayName", "CFD rig reading, {\it\beta} = 1.3 (band 1.25–1.4)");
% beta = 1.0, the lower limit for a bubble that fills the bore
plot(ax, series.U, series.R, ":", "Color", salmon, "LineWidth", cfg.dataLineWidth, ...
    "HandleVisibility", "off");
scatter(ax, series.U, series.R, cfg.markerArea, [1 1 1], "filled", "Marker", "v", ...
    "MarkerEdgeColor", salmon, "LineWidth", cfg.markerEdgeWidth, "HandleVisibility", "off");
handles(end+1) = plot(ax, NaN, NaN, ":", "Color", salmon, "LineWidth", cfg.dataLineWidth, ...
    "Marker", "v", "MarkerSize", 6, "MarkerFaceColor", [1 1 1], "MarkerEdgeColor", salmon, ...
    "DisplayName", "CFD rig reading, {\it\beta} = 1.0 (lower limit)");
% published points again on top, so no CFD marker hides them
series = L(L.series == "phillips_data", :);
scatter(ax, series.U, series.R, cfg.markerArea, red, "filled", "Marker", "d", ...
    "MarkerEdgeColor", cfg.markerEdge, "LineWidth", cfg.markerEdgeWidth, "HandleVisibility", "off");

xlim(ax, [0.07 0.7]);
ylim(ax, [0.03 1.2]);
xticks(ax, [0.08 0.15 0.3 0.6]);
xticklabels(ax, ["0.08" "0.15" "0.3" "0.6"]);
yticks(ax, [0.05 0.1 0.2 0.5 1]);
yticklabels(ax, ["0.05" "0.1" "0.2" "0.5" "1"]);
xlabel(ax, "{\itU} (m s^{-1})", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
ylabel(ax, "{\itR}", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
lg = legend(ax, fliplr(handles), "NumColumns", 1, "Box", "off", ...
    "FontSize", cfg.legendFontSize, "Interpreter", "tex");
% square axes on the left and a one-column legend beside them (the
% legend below the axes made the figure 4.6 in tall for 2.5 in of plot)
ax.Units = "inches";
ax.Position = [0.62 0.48 2.4 2.4];
lg.Units = "inches";
lg.Position = [3.15 0.35 3.1 2.55];
lg.ItemTokenSize = [30 10];  % long enough to show the dash

tims_export(fig, fullfile(here, "..", "fig_phillips"), cfg);
close(fig);
end
