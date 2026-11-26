function make_fig_flow()
%MAKE_FIG_FLOW Figure fig:flow: the flow round four designs at 0.3 m/s.
%   (i) Phillips, Perks, Elliott and McDonald, two by two: the inlet region, |U|/U0 with
%   streamlines, one window (x -60..120 mm from the inlet tip, z +-30 mm) and one colour scale;
%   (ii) the same four designs, two by two: the outlet and wake, u/U0 on one diverging scale centred on zero with the u = 0
%   contour (the recirculation bubble), one window (x -90..360 mm from the rear outer face,
%   z +-75 mm). Bottom: u/U0 on the axis for all six designs, aligned at the inlet tip.
%   Solids from fig_designs_outline.csv are drawn over the fields.
%   Input: ../data/fig_flow_inlet_<d>.csv, fig_flow_wake_<d>.csv, fig_flow_axis.csv and
%   fig_designs_outline.csv / fig_designs_stations.csv, written by cfd/geometry/fig_flow_data.py
%   and fig_designs_data.py.  Output: ../fig_flow.pdf and .png.  Tested with MATLAB R2024b.

here = fileparts(mfilename("fullpath"));
dataDir = fullfile(here, "..", "data");
cfg = tims_fig_config();
O = readtable(fullfile(dataDir, "fig_designs_outline.csv"), "TextType", "string");
S = readtable(fullfile(dataDir, "fig_designs_stations.csv"), "TextType", "string");
X = readtable(fullfile(dataDir, "fig_flow_axis.csv"), "TextType", "string");

names = ["Phillips" "Perks" "Elliott" "McDonald"];
files = ["phillips" "perks" "elliott" "mcdonald"];

% colour maps: sequential for |U| (0-1.4), diverging about zero for u (-0.4..1.4)
blue = [59 76 192] ./ 255;      % colour-map anchors for panels (i) and (ii)
red = [180 4 38] ./ 255;
% Okabe-Ito line colours for panel (iii), as Figs. 3-6
lineBlue = [0 114 178] ./ 255;
lineRed = [213 94 0] ./ 255;
light = [0.97 0.97 0.97];
seq = interp1([0 0.5 1], [light; [129 164 251] ./ 255; blue * 0.75], linspace(0, 1, 256));
nNeg = round(256 * 0.4 / 1.8);
div = [interp1([0 1], [blue; light], linspace(0, 1, nNeg)); ...
    interp1([0 0.45 1], [light; [229 112 88] ./ 255; red], linspace(0, 1, 256 - nNeg))];
solidFace = [0.62 0.62 0.62];

% Drawn at the printed size: 6.5 in, the full text width, so font sizes are as printed (ASCE: 8 pt).
W = 6.5;
left = 0.62; mid = 0.22; right = 0.20;   % room for the 120 tick label at the right edge
axW = (W - left - mid - right) / 2;
axH = axW / 3;
rowGap = 0.10;      % between the two inlet rows (no tick labels on the upper row)
wakeGap = 0.55;     % x label of the inlet rows and the (ii) panel label
cbH = 0.09;
botH = 1.30;
top = 0.25;
yRowsBottom = 0.45 + botH + 0.55 + cbH + 0.58;
yWake2 = yRowsBottom;
yWake1 = yWake2 + axH + rowGap;
yIn2 = yWake1 + axH + wakeGap;
yIn1 = yIn2 + axH + rowGap;
H = yIn1 + axH + top;
fig = figure("Color", "w", "Units", "inches", "Visible", "off", "Position", [1 1 W H], ...
    "PaperPositionMode", "auto");

% The four inlet panels go two by two. Panel (ii) was cut to Phillips
% and McDonald to save space; in ASCE's submission layout the figure has its own page, so the Perks
% and Elliott wake rows are back, two by two in the same order as (i).
xCol = [left, left + axW + mid];
inletPos = [xCol(1) yIn1; xCol(2) yIn1; xCol(1) yIn2; xCol(2) yIn2];
for r = 1:numel(names)
    ps = design_shape(O, names(r));
    T = readtable(fullfile(dataDir, "fig_flow_inlet_" + files(r) + ".csv"));
    [xg, zg, F] = grid_of(T, "umag");
    [~, ~, Uq] = grid_of(T, "u");
    [~, ~, Wq] = grid_of(T, "w");
    ax1 = axes(fig, "Units", "inches", "Position", [inletPos(r, :) axW axH]);
    field_axes(ax1, cfg);
    imagesc(ax1, xg, zg, F, "AlphaData", ~isnan(F));
    colormap(ax1, seq);
    clim(ax1, [0 1.4]);
    Uq(isnan(Uq)) = 0;
    Wq(isnan(Wq)) = 0;
    sz = [-27:3:-6, -4:1:4, 6:3:27];
    hs = streamline(ax1, xg, zg, Uq, Wq, repmat(xg(1) + 0.5, size(sz)), sz);
    set(hs, "Color", [0.15 0.15 0.15], "LineWidth", 0.35);
    plot(ax1, ps, "FaceColor", solidFace, "FaceAlpha", 1, "EdgeColor", "k", "LineWidth", 0.3);
    xlim(ax1, [xg(1) xg(end)]);
    ylim(ax1, [zg(1) zg(end)]);
    design_label(ax1, names(r), cfg);
    xticks(ax1, -40:40:120);    % no tick label at the left edge, where the next panel's ends
    if mod(r, 2) == 1
        ylabel(ax1, "{\itz} (mm)", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
    else
        ax1.YTickLabel = [];
    end
    if r <= 2
        ax1.XTickLabel = [];
    else
        xlabel(ax1, "{\itx} from the inlet tip (mm)", "Interpreter", "tex", ...
            "FontSize", cfg.labelFontSize);
    end
    if r == 1
        tims_panel_label(ax1, "(a) inlet", cfg, [0 1.04], "bottom");
    end
end

% (ii) outlet and wake: u/U0 with the u = 0 contour, the same four designs in the same order
wakePos = [xCol(1) yWake1; xCol(2) yWake1; xCol(1) yWake2; xCol(2) yWake2];
for r = 1:numel(names)
    ps = design_shape(O, names(r));
    rear = S.x_mm(S.design == names(r) & S.station == "body rear face");
    T = readtable(fullfile(dataDir, "fig_flow_wake_" + files(r) + ".csv"));
    T.x_mm = T.x_tip_mm;    % plot in inlet-tip coordinates so the outlines overlay directly
    [xg, zg, F] = grid_of(T, "u");
    ax2 = axes(fig, "Units", "inches", "Position", [wakePos(r, :) axW axH]);
    field_axes(ax2, cfg);
    imagesc(ax2, xg - rear, zg, F, "AlphaData", ~isnan(F));
    colormap(ax2, div);
    clim(ax2, [-0.4 1.4]);
    contour(ax2, xg - rear, zg, F, [0 0], "k-", "LineWidth", 0.6);
    plot(ax2, translate(ps, -rear, 0), "FaceColor", solidFace, "FaceAlpha", 1, ...
        "EdgeColor", "k", "LineWidth", 0.3);
    xlim(ax2, [xg(1) xg(end - 1)] - rear);
    ylim(ax2, [zg(1) zg(end)]);
    design_label(ax2, names(r), cfg, "right");   % top right, clear of the body
    if mod(r, 2) == 1
        ylabel(ax2, "{\itz} (mm)", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
    else
        ax2.YTickLabel = [];
    end
    if r <= 2
        ax2.XTickLabel = [];
    else
        xlabel(ax2, "{\itx} from the rear face of the body (mm)", "Interpreter", "tex", ...
            "FontSize", cfg.labelFontSize);
    end
    if r == 1
        tims_panel_label(ax2, "(b) outlet and wake", cfg, [0 1.04], "bottom");
    end
end

% colour bars under the columns
yCb = yRowsBottom - 0.55 - cbH;
cb1 = axes(fig, "Units", "inches", "Position", [left yCb axW cbH]);
colour_bar(cb1, seq, [0 1.4], "(a) |{\bfU}|/{\itU}", cfg);
cb2 = axes(fig, "Units", "inches", "Position", [left + axW + mid yCb axW cbH]);
colour_bar(cb2, div, [-0.4 1.4], "(b) {\itu}/{\itU} (black: {\itu} = 0)", cfg);

% bottom: u/U0 on the axis, all six designs
designs = [ ...
    struct("label", "Phillips", "colour", lineBlue, "line", "-"), ...
    struct("label", "Foets", "colour", lineBlue, "line", "--"), ...
    struct("label", "Perks", "colour", lineRed, "line", "-"), ...
    struct("label", "Perks + meter", "colour", lineRed, "line", "--"), ...
    struct("label", "Elliott", "colour", [0 158 115] ./ 255, "line", "-"), ...
    struct("label", "McDonald", "colour", [204 121 167] ./ 255, "line", "-")];
axB = axes(fig, "Units", "inches", "Position", [left 0.45 W - left - right botH]);
tims_format_axes(axB, cfg, false);
plot(axB, [-0.3 1.7], [0 0], "-", "Color", cfg.refColour, "LineWidth", 0.5, ...
    "HandleVisibility", "off");
plot(axB, [-0.3 1.7], [1 1], ":", "Color", cfg.refColour, "LineWidth", 0.5, ...
    "HandleVisibility", "off");
for d = designs
    q = X(X.design == d.label, :);
    plot(axB, q.x_mm / 1000, q.u, d.line, "Color", d.colour, "LineWidth", cfg.dataLineWidth, ...
        "DisplayName", d.label);
end
xlim(axB, [-0.3 1.7]);
ylim(axB, [-0.45 1.45]);
xlabel(axB, "{\itx} from the inlet tip (m)", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
ylabel(axB, "{\itu}/{\itU} on the axis", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
tims_panel_label(axB, "(c)", cfg);
lg = legend(axB, "NumColumns", 3, "Box", "off", "FontSize", cfg.legendFontSize);
lg.ItemTokenSize = [18 8];
% place it in the empty band above u/U0 = 1
axW_B = W - left - right;
lg.Units = "inches";
lg.Position = [left + 0.25 * axW_B, 0.45 + 0.76 * botH, 0.70 * axW_B, 0.22 * botH];

tims_export(fig, fullfile(here, "..", "fig_flow"), cfg);
close(fig);
end

function design_label(ax, name, cfg, side)
if nargin < 4
    side = "left";
end
x = 0.015;
if side == "right"
    x = 0.985;
end
text(ax, x, 0.95, name, "Units", "normalized", "HorizontalAlignment", side, ...
    "VerticalAlignment", "top", "FontSize", cfg.axesFontSize, "FontWeight", "bold", ...
    "Interpreter", "none");
end

function field_axes(ax, cfg)
hold(ax, "on");
box(ax, "on");
axis(ax, "equal");
set(ax, "YDir", "normal", "FontSize", cfg.axesFontSize - 1, "LineWidth", 0.6, ...
    "TickDir", "out", "Layer", "top");
end

function [xg, zg, F] = grid_of(T, var)
xg = unique(T.x_mm);
zg = unique(T.z_mm);
[~, ix] = ismember(T.x_mm, xg);
[~, iz] = ismember(T.z_mm, zg);
F = nan(numel(zg), numel(xg));
F(sub2ind(size(F), iz, ix)) = T.(var);
xg = xg';
zg = zg';
end

function colour_bar(ax, map, lims, label, cfg)
v = linspace(lims(1), lims(2), size(map, 1));
imagesc(ax, v, [0 1], [v; v]);
colormap(ax, map);
clim(ax, lims);
set(ax, "YTick", [], "FontSize", cfg.axesFontSize - 1, "TickDir", "out", "Box", "on", ...
    "LineWidth", 0.6, "XLim", lims);
% label both ends (the (b) bar ran below 0 with no tick at its minimum)
t = 0:0.5:lims(2);
xticks(ax, unique([lims(1), t(t > lims(1) + 0.2)]));
xlabel(ax, label, "Interpreter", "tex", "FontSize", cfg.labelFontSize - 1);
end

function ps = design_shape(O, name)
o = O(O.design == name, :);
xs = {};
zs = {};
for k = unique(o.loop)'
    q = o(o.loop == k, :);
    xs{end + 1} = q.x_mm; %#ok<AGROW>
    zs{end + 1} = q.z_mm; %#ok<AGROW>
end
warning("off", "MATLAB:polyshape:repairedBySimplify");
ps = polyshape(xs, zs);
end
