function make_fig_budget()
%MAKE_FIG_BUDGET Figure fig:budget: how the free-stream dynamic pressure is spent
%   along each sampler's flow path at 0.3 m/s, as stacked shares of rho U^2/2 (they sum to 1).
%   Where the outlet sits in suction (negative base pressure, McDonald), the extra drive beyond 1
%   is hatched.
%   Input: ../data/fig_budget.csv, written by cfd/geometry/fig_budget.py.
%   Output: ../fig_budget.pdf and ../fig_budget.png.  Tested with MATLAB R2024b.

here = fileparts(mfilename("fullpath"));
cfg = tims_fig_config();
T = readtable(fullfile(here, "..", "data", "fig_budget.csv"), "TextType", "string");

% Okabe-Ito colours, as Figs. 3, 5 and 6: inlet side in blues and green, outlet
% side in orange, vermillion and reddish purple, base in grey.
segments = [ ...
    struct("key", "approach", "label", "approach to the mouth", "colour", [86 180 233] ./ 255), ...
    struct("key", "inlet", "label", "inlet bore", "colour", [0 114 178] ./ 255), ...
    struct("key", "expansion", "label", "expansion into the chamber", "colour", [0 158 115] ./ 255), ...
    struct("key", "out_entry", "label", "outlet entry", "colour", [230 159 0] ./ 255), ...
    struct("key", "out_bore", "label", "outlet bore and fittings", "colour", [213 94 0] ./ 255), ...
    struct("key", "exit", "label", "kinetic energy leaving", "colour", [204 121 167] ./ 255), ...
    struct("key", "base", "label", "outlet static pressure", "colour", [0.62 0.62 0.62])];
barHeight = 0.62;

fig = figure("Color", "w", "Units", "inches", "Visible", "off", ...
    "Position", [1 1 cfg.pageWidth_in 3.0], "PaperPositionMode", "auto");
layout = tiledlayout(fig, 1, 1, "TileSpacing", "compact", "Padding", "compact");
ax = nexttile(layout);
tims_format_axes(ax, cfg, false);

n = height(T);
yLabels = strings(n, 1);
for i = 1:n
    row = T(n + 1 - i, :);                 % first design at the top, as in fig_RU
    y = i;
    left = 0;
    merged = row.expansion == 0 && row.out_bore == 0;   % meter variant
    for sgm = segments
        val = row.(sgm.key);
        if val <= 0
            continue
        end
        isMerged = merged && any(sgm.key == ["inlet" "out_entry"]);
        face = sgm.colour;
        if isMerged     % merged segments in a lighter tint, so they do not read as the legend's
            face = 0.45 * sgm.colour + 0.55;
        end
        rectangle(ax, "Position", [left y - barHeight/2 val barHeight], ...
            "FaceColor", face, "EdgeColor", "w", "LineWidth", 0.8);
        if isMerged
            txt = "inlet bore, meter, expansion";
            if sgm.key == "out_entry"
                txt = "outlet entry and bore";
            end
            text(ax, left + val/2, y, txt, "HorizontalAlignment", "center", ...
                "VerticalAlignment", "middle", "FontSize", cfg.annotationFontSize, "Color", cfg.refColour);
        end
        left = left + val;
    end
    if row.base < 0                         % outlet suction: drive beyond rho U^2 / 2
        hatch(ax, 1.0, 1.0 - row.base, y - barHeight/2, y + barHeight/2, cfg.refColour);
        % direct label: a legend patch cannot show hatching
        text(ax, 1.0 - row.base/2, y + barHeight/2 + 0.08, "suction", ...
            "HorizontalAlignment", "center", "VerticalAlignment", "bottom", ...
            "FontSize", cfg.annotationFontSize, "Color", cfg.refColour);
    end
    yLabels(i) = row.design;
    % R in its own column right of the bars (two-line tick labels get split by MATLAB)
    text(ax, 1.215, y, sprintf("%.3f", row.R), "HorizontalAlignment", "right", ...
        "VerticalAlignment", "middle", "FontSize", cfg.axesFontSize);
end
plot(ax, [1 1], [0.4 n + 0.5], "--", "Color", cfg.refColour, "LineWidth", cfg.refLineWidth);

% "R" heads its column inside the axes (placed above the frame it was cut off at
% the top of the exported figure, since the compact layout leaves no margin there)
text(ax, 1.215, n + 0.72, "{\itR}", "HorizontalAlignment", "right", ...
    "VerticalAlignment", "middle", "FontSize", cfg.axesFontSize);
xlim(ax, [0 1.24]);
xticks(ax, 0:0.2:1);
ylim(ax, [0.4 n + 0.95]);
yticks(ax, 1:n);
yticklabels(ax, yLabels);
ax.TickLength = [0 0];
ax.XAxis.TickLength = [0.01 0.01];
ax.Box = "off";  % no mirrored, unlabelled ticks on the top axis
xlabel(ax, "share of \rho{\itU}^2/2 ({\itU} = 0.3 m s^{-1})", ...
    "Interpreter", "tex", "FontSize", cfg.labelFontSize);

% legend from dummy patches
handles = gobjects(0);
for sgm = segments
    handles(end+1) = fill(ax, NaN, NaN, sgm.colour, "EdgeColor", "w", ...
        "DisplayName", sgm.label); %#ok<AGROW>
end
lg = legend(ax, handles, "NumColumns", 3, "Box", "off", "FontSize", cfg.legendFontSize);
lg.Layout.Tile = "south";
lg.ItemTokenSize = [14 10];

tims_export(fig, fullfile(here, "..", "fig_budget"), cfg);
close(fig);
end

function hatch(ax, x0, x1, y0, y1, colour)
% Diagonal hatching over [x0 x1] x [y0 y1] (data units), lines clipped to the box.
pitch = 0.012;
slope = (y1 - y0) / 0.03;                  % rise per unit x, so the lines read as ~45 degrees
for c = (x0 - (y1 - y0) / slope):pitch:x1
    xs = [c, c + (y1 - y0) / slope];
    ys = [y0, y1];
    if xs(1) < x0
        ys(1) = y0 + (x0 - xs(1)) * slope;
        xs(1) = x0;
    end
    if xs(2) > x1
        ys(2) = y1 - (xs(2) - x1) * slope;
        xs(2) = x1;
    end
    if xs(2) > xs(1)
        plot(ax, xs, ys, "-", "Color", colour, "LineWidth", 0.6, "HandleVisibility", "off");
    end
end
rectangle(ax, "Position", [x0 y0 x1 - x0 y1 - y0], "EdgeColor", colour, "LineWidth", 0.6);
end
