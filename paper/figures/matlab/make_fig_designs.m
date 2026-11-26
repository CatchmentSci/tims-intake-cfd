function make_fig_designs()
%MAKE_FIG_DESIGNS Figure fig:designs: the modelled designs on the symmetry plane.
%   One row per design (Phillips/Foets, Perks/Perks + meter, Elliott, McDonald) and one
%   view type per column, each column at one scale: (i) the whole section, aligned at the inlet
%   tip; (ii) the inlet, one common 60 mm window; (iii) the outlet, one common window placed at
%   each body's rear face. Solid grey, water white; parts assumed in Table 1 hatched; variants
%   (Foets 2 mm outlet tube, the in-chamber meter) dashed in their partner's row.
%   Input: ../data/fig_designs_outline.csv, fig_designs_assumed.csv, fig_designs_stations.csv,
%   written by cfd/geometry/fig_designs_data.py from the production meshes.
%   Output: ../fig_designs.pdf and ../fig_designs.png.  Tested with MATLAB R2024b.

here = fileparts(mfilename("fullpath"));
dataDir = fullfile(here, "..", "data");
cfg = tims_fig_config();
O = readtable(fullfile(dataDir, "fig_designs_outline.csv"), "TextType", "string");
A = readtable(fullfile(dataDir, "fig_designs_assumed.csv"), "TextType", "string", ...
    "Delimiter", ",");
S = readtable(fullfile(dataDir, "fig_designs_stations.csv"), "TextType", "string");

rows = [ ...
    struct("design", "Phillips", "variant", "Foets", "label", ...
        "Phillips, 4 mm bores; Foets 2 mm outlet dashed"), ...
    struct("design", "Perks", "variant", "Perks + meter", "label", ...
        "Perks, 7.94 mm tubing and 6.35 mm insert; meter dashed"), ...
    struct("design", "Elliott", "variant", "", "label", "Elliott, 4 mm bores"), ...
    struct("design", "McDonald", "variant", "", "label", ...
        "McDonald, 2 mm inlet, 3 mm outlet hole")];

solidFace = [0.62 0.62 0.62];
edgeCol = [0.15 0.15 0.15];
variantCol = [180 4 38] ./ 255;
winSection = [-40 1400; -82 82];          % x, z (mm) from the inlet tip
winInlet = [-10 50; -30 30];              % common 60 mm inlet window
winOutlet = [-60 150; -75 75];            % relative to each body's rear outer face

% layout (inches): three columns, four rows
% Drawn at the printed size, 6.5 in (the ASCE text width), and included unscaled so the 8 pt text
% prints at 8 pt. Was 7.0 in and scaled to 94 %.
W = 6.5; rowHs = [0.95 0.95 0.95 0.95]; top = 0.42; bottom = 0.62; gap = 0.08;
H = top + bottom + sum(rowHs);
% Columns (ii) and (iii) are axis-equal and height-limited, so each is as wide as its drawing; the
% scale bars below use colW, so a wider column would stretch them (the inlet bar was 16 % long
% while (ii) was 0.80 in wide).
axHgt = rowHs(1) - gap - 0.18;
colW = [4.45, axHgt * diff(winInlet(1, :)) / diff(winInlet(2, :)), ...
    axHgt * diff(winOutlet(1, :)) / diff(winOutlet(2, :))];
colX = [0.05, 0.05 + colW(1) + 0.15, 0.05 + colW(1) + 0.15 + colW(2) + 0.15];
fig = figure("Color", "w", "Units", "inches", "Visible", "off", "Position", [1 1 W H], ...
    "PaperPositionMode", "auto");

heads = ["(a) section on the symmetry plane", "(b) inlet", "(c) outlet"];  % (a)-(c), ASCE style
for r = 1:numel(rows)
    rowH = rowHs(r);
    y0 = H - top - sum(rowHs(1:r)) + gap / 2;
    d = rows(r);
    hasGeom = any(O.design == d.design);
    rear = 0;
    if hasGeom
        rear = S.x_mm(S.design == d.design & S.station == "body rear face");
    end
    for c = 1:3
        wAx = colW(c);
        ax = axes(fig, "Units", "inches", "Position", [colX(c) y0 wAx rowH - gap - 0.18]);
        hold(ax, "on");
        axis(ax, "equal");
        ax.Visible = "off";
        switch c
            case 1
                win = winSection; hs = 9;
            case 2, win = winInlet; hs = 1.6;
            otherwise, win = winOutlet + [rear rear; 0 0]; hs = 4.5;
        end
        xlim(ax, win(1, :));
        ylim(ax, win(2, :));
        if c == 1
            text(ax, 0, 1.0, d.label, "Units", "normalized", "FontSize", cfg.legendFontSize, ...
                "VerticalAlignment", "bottom", "Clipping", "off");
        end
        if ~hasGeom
            continue
        end
        ps = design_shape(O, d.design);
        plot(ax, ps, "FaceColor", solidFace, "FaceAlpha", 1, "EdgeColor", edgeCol, ...
            "LineWidth", 0.4);
        % assumed parts of Table 1: hatched where they are solid
        a = A(A.design == d.design, :);
        for p = unique(a.part)'
            q = a(a.part == p, :);
            box = parts_shape(q);
            hatch(ax, intersect(ps, box), hs, win);
        end
        % variants dashed
        if d.variant == "Foets"
            q = A(A.design == "Foets", :);
            plot(ax, parts_shape(q), "FaceColor", "none", "EdgeColor", variantCol, ...
                "LineStyle", "--", "LineWidth", 0.7);
        elseif d.variant == "Perks + meter"
            q = A(A.design == "Perks + meter", :);
            plot(ax, parts_shape(q), "FaceColor", "none", "EdgeColor", variantCol, "LineStyle", "--", ...
                "LineWidth", 0.7);
        end
        if c == 1   % mark the detail windows
            rectangle(ax, "Position", [winInlet(1, 1) winInlet(2, 1) diff(winInlet(1, :)) ...
                diff(winInlet(2, :))], "EdgeColor", cfg.mutedText, "LineWidth", 0.4);
            rectangle(ax, "Position", [rear + winOutlet(1, 1) winOutlet(2, 1) ...
                diff(winOutlet(1, :)) diff(winOutlet(2, :))], "EdgeColor", cfg.mutedText, ...
                "LineWidth", 0.4);
        end
    end
end

% column headings on one baseline (they sat at different heights when
% placed in each axes' normalised units, since axis equal shrinks each plot box differently)
axHead = axes(fig, "Units", "inches", "Position", [0 0 W H], "Visible", "off");
xlim(axHead, [0 W]);
ylim(axHead, [0 H]);
for c = 1:3
    text(axHead, colX(c), H - 0.14, heads(c), "FontSize", cfg.legendFontSize, ...
        "FontWeight", "bold", "VerticalAlignment", "top", "Clipping", "off");
end

% scale bars under each column, and the legend
ybar = 0.40;
bars = [200 20 50];
wins = {winSection, winInlet, winOutlet};
for c = 1:3
    ax = axes(fig, "Units", "inches", "Position", [colX(c) ybar - 0.12 colW(c) 0.24]);
    ax.Visible = "off";
    hold(ax, "on");
    w = wins{c};
    xlim(ax, w(1, :));
    ylim(ax, [0 1]);
    x0 = mean(w(1, :)) - bars(c) / 2;   % centred in its column (ii, iii)
    if c == 1, x0 = w(1, 1) + 0.02 * diff(w(1, :)); end   % (i): at the left
    plot(ax, [x0 x0 + bars(c)], [0.6 0.6], "-", "Color", edgeCol, "LineWidth", 1.2);
    text(ax, x0 + bars(c) / 2, 0.35, sprintf("%d mm", bars(c)), "HorizontalAlignment", ...
        "center", "VerticalAlignment", "top", "FontSize", cfg.annotationFontSize);
end
% key drawn by hand: a legend patch cannot show hatching
axL = axes(fig, "Units", "inches", "Position", [colX(1) + 1.6 0.02 3.2 0.22]);
axL.Visible = "off";
hold(axL, "on");
xlim(axL, [0 3.2]);                 % axes units = inches, so the swatches are square
ylim(axL, [0 0.22]);
sw = 0.14;                          % swatch size (in)
y0 = 0.04;
keyX = [0 0.75 2.05];
fill(axL, keyX(1) + [0 sw sw 0], y0 + [0 0 sw sw], solidFace, "EdgeColor", edgeCol);
text(axL, keyX(1) + sw + 0.06, y0 + sw/2, "solid", "FontSize", cfg.legendFontSize, ...
    "VerticalAlignment", "middle");
fill(axL, keyX(2) + [0 sw sw 0], y0 + [0 0 sw sw], solidFace, "EdgeColor", edgeCol);
for c = sw * [-0.5 0 0.5]           % white 45-degree hatching (lines x - y = c), as on the drawings
    if c >= 0
        xs = [c sw]; ys = [0 sw - c];
    else
        xs = [0 sw + c]; ys = [-c sw];
    end
    plot(axL, keyX(2) + xs, y0 + ys, "-", "Color", "w", "LineWidth", 0.9);
end
text(axL, keyX(2) + sw + 0.06, y0 + sw/2, "assumed (Table 1)", ...
    "FontSize", cfg.legendFontSize, "VerticalAlignment", "middle");
plot(axL, keyX(3) + [0 0.22], y0 + sw/2 * [1 1], "--", "Color", variantCol, "LineWidth", 0.7);
text(axL, keyX(3) + 0.28, y0 + sw/2, "variant", "FontSize", cfg.legendFontSize, ...
    "VerticalAlignment", "middle");

tims_export(fig, fullfile(here, "..", "fig_designs"), cfg);
close(fig);
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

function ps = parts_shape(q)
% one or more closed polygons, each starting at vertex 0
starts = find(q.vertex == 0);
ends = [starts(2:end) - 1; height(q)];
xs = {};
zs = {};
for k = 1:numel(starts)
    xs{end + 1} = q.x_mm(starts(k):ends(k)); %#ok<AGROW>
    zs{end + 1} = q.z_mm(starts(k):ends(k)); %#ok<AGROW>
end
warning("off", "MATLAB:polyshape:repairedBySimplify");
ps = polyshape(xs, zs);
end

function hatch(ax, ps, spacing, win)
% white 45-degree hatching clipped to the polyshape ps
if ps.NumRegions == 0
    return
end
[xl, zl] = boundingbox(ps);
xl = [max(xl(1), win(1, 1)) min(xl(2), win(1, 2))];
zl = [max(zl(1), win(2, 1)) min(zl(2), win(2, 2))];
if any(diff([xl; zl], 1, 2) <= 0)
    return
end
span = diff(xl) + diff(zl);
for c = (xl(1) - diff(zl)):spacing:(xl(1) + span)
    seg = [c zl(1); c + diff(zl) zl(2)];
    in = intersect(ps, seg);
    if ~isempty(in)
        plot(ax, in(:, 1), in(:, 2), "-", "Color", [1 1 1], "LineWidth", 0.45, ...
            "HandleVisibility", "off");
    end
end
end
