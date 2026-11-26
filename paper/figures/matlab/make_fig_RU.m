function make_fig_RU()
%MAKE_FIG_RU Figure fig:RU: intake ratio against free-stream velocity.
%   (a) R(U) for the six designs, with the 1-D pipe-network screening band of each;
%   (b) the local slope d ln R / d ln U between adjacent velocities, at their geometric mean.
%   Input: ../data/fig_RU.csv, fig_RU_slope.csv, fig_RU_model_own.csv, written by
%   cfd/geometry/plot_RU.py (run it from cfd/ first, after new runs land), and fig_bands.csv
%   (the uncertainty envelope in (a)), written by cfd/geometry/sensitivity_bands.py.
%   Output: ../fig_RU.pdf and ../fig_RU.png.  Tested with MATLAB R2024b.

here = fileparts(mfilename("fullpath"));
dataDir = fullfile(here, "..", "data");
cfg = tims_fig_config();

R = readtable(fullfile(dataDir, "fig_RU.csv"), "TextType", "string");
S = readtable(fullfile(dataDir, "fig_RU_slope.csv"), "TextType", "string");
B = readtable(fullfile(dataDir, "fig_RU_model_own.csv"), "TextType", "string");

% A variant shares its parent's colour, with a dashed line and an open marker; every design also
% has its own marker shape, so none is told apart by colour alone.
% Okabe-Ito colours, as Figs. 4-6: Phillips and Foets blue, the Perks builds
% vermillion, Elliott bluish green, McDonald reddish purple (variable names keep their old roles)
blue = [0 114 178] ./ 255;
red = [213 94 0] ./ 255;
salmon = [0 158 115] ./ 255;
grey = [204 121 167] ./ 255;
designs = [ ...
    struct("label", "Phillips", "colour", blue, "marker", "o", "line", "-", "open", false), ...
    struct("label", "Foets", "colour", blue, "marker", "o", "line", "--", "open", true), ...
    struct("label", "Perks", "colour", red, "marker", "s", "line", "-", "open", false), ...
    struct("label", "Perks + meter", "colour", red, "marker", "s", "line", "--", "open", true), ...
    struct("label", "Elliott", "colour", salmon, "marker", "^", "line", "-", "open", false), ...
    struct("label", "McDonald", "colour", grey, "marker", "d", "line", "-", "open", false)];

fig = figure("Color", "w", "Units", "inches", "Visible", "off", ...
    "Position", [1 1 cfg.pageWidth_in 3.25], "PaperPositionMode", "auto");
layout = tiledlayout(fig, 1, 2, "TileSpacing", "compact", "Padding", "compact");
axR = nexttile(layout, 1);
tims_format_axes(axR, cfg);
axS = nexttile(layout, 2);
tims_format_axes(axS, cfg);

% (b) limit lines first, so the data sit on top
plot(axS, [0 1.6], [1 1], "--", "Color", cfg.refColour, "LineWidth", cfg.refLineWidth, ...
    "HandleVisibility", "off");
plot(axS, [0 1.6], [0 0], "--", "Color", cfg.refColour, "LineWidth", cfg.refLineWidth, ...
    "HandleVisibility", "off");
% Superscripts in 8 pt text: \fontsize{10} before the ^ makes them 8 pt as printed (MATLAB sets
% them at 0.8 of the current size; no printed text below 8 pt).
text(axS, 1.56, 1.03, "friction only, {\itR} \propto {\itU}, {\itq} \propto {\itU}\fontsize{10}^2", "Interpreter", "tex", ...
    "HorizontalAlignment", "right", "VerticalAlignment", "bottom", ...
    "FontSize", cfg.annotationFontSize, "Color", cfg.mutedText);
text(axS, 0.05, 0.03, "inertial losses only, {\itq} \propto {\itU}", "Interpreter", "tex", ...
    "HorizontalAlignment", "left", "VerticalAlignment", "bottom", ...
    "FontSize", cfg.annotationFontSize, "Color", cfg.mutedText);

% (a) numerical and closure uncertainty of each design (Table 4 total, R_lo_total to R_hi_total from
% cfd/geometry/sensitivity_bands.py), as a light envelope under everything
Bd = readtable(fullfile(dataDir, "fig_bands.csv"), "TextType", "string");
for d = designs
    q = Bd(Bd.design == d.label, :);
    if ~isempty(q)
        fill(axR, [q.U; flipud(q.U)], [q.R_lo_total; flipud(q.R_hi_total)], d.colour, ...
            "FaceAlpha", 0.15, "EdgeColor", d.colour, "EdgeAlpha", 0.6, "LineWidth", 0.4, ...
            "HandleVisibility", "off");  % thin outline per band
    end
end

% (a) screening model of each design at its own CFD drive, as thin dotted lines under everything
%
for d = designs
    b = B(B.design == d.label, :);
    if ~isempty(b)
        plot(axR, b.U, b.R_model, ":", "Color", d.colour, "LineWidth", 1.1, ...
            "HandleVisibility", "off");
    end
end

handles = gobjects(0);
for d = designs
    r = R(R.design == d.label, :);
    s = S(S.design == d.label, :);
    face = d.colour;
    if d.open
        face = [1 1 1];
    end
    plot(axR, r.U, r.R, d.line, "Color", d.colour, "LineWidth", cfg.dataLineWidth, ...
        "HandleVisibility", "off");
    plot(axS, s.U, s.slope, d.line, "Color", d.colour, "LineWidth", cfg.dataLineWidth, ...
        "HandleVisibility", "off");
    scatter(axR, r.U, r.R, cfg.markerArea, face, "filled", "Marker", d.marker, ...
        "MarkerEdgeColor", cfg.markerEdge, "LineWidth", cfg.markerEdgeWidth, ...
        "HandleVisibility", "off");
    scatter(axS, s.U, s.slope, cfg.markerArea, face, "filled", "Marker", d.marker, ...
        "MarkerEdgeColor", cfg.markerEdge, "LineWidth", cfg.markerEdgeWidth, ...
        "HandleVisibility", "off");
    % legend entry: line and marker together
    handles(end+1) = plot(axR, NaN, NaN, d.line, "Color", d.colour, ...
        "LineWidth", cfg.dataLineWidth, "Marker", d.marker, "MarkerSize", 6, ...
        "MarkerFaceColor", face, "MarkerEdgeColor", cfg.markerEdge, ...
        "DisplayName", d.label); %#ok<AGROW>
end
handles(end+1) = plot(axR, NaN, NaN, ":", "Color", [0.4 0.4 0.4], "LineWidth", 1.1, ...
    "DisplayName", "1-D model at each design's drive");

for ax = [axR axS]
    xlim(ax, [0 1.6]);
    xticks(ax, 0:0.25:1.5);
    xticklabels(ax, ["0" "" "0.5" "" "1" "" "1.5"]);
    xtickangle(ax, 0);  % horizontal labels
    xlabel(ax, "{\itU} (m s^{-1})", "Interpreter", "tex", ...
        "FontSize", cfg.labelFontSize);
end
ylim(axR, [0 0.5]);
ylabel(axR, "{\itR}", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
ylim(axS, [-0.05 1.2]);
yticks(axS, 0:0.2:1);
ylabel(axS, "d ln {\itR} / d ln {\itU}", "Interpreter", "tex", ...
    "FontSize", cfg.labelFontSize);
tims_panel_label(axR, "(a)", cfg);
tims_panel_label(axS, "(b)", cfg);

lg = legend(axR, handles, "NumColumns", 4, "Box", "off", ...
    "FontSize", cfg.legendFontSize, "Interpreter", "tex");
lg.Layout.Tile = "south";
lg.ItemTokenSize = [22 10];

tims_export(fig, fullfile(here, "..", "fig_RU"), cfg);
close(fig);
end
