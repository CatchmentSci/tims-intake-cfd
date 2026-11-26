function make_fig_S1_mesh()
%MAKE_FIG_S1_MESH Supplemental Figure S1: domain, mesh-sizing regions and mesh.
%   (a) The computational domain (half model, symmetry plane) with the mesh-sizing regions and
%   their target face sizes at s = 1, and the Phillips sampler drawn to scale. (b) and (c): the
%   medium volume mesh of the Phillips sampler cut by the symmetry plane at the inlet tip and at
%   the outlet tip, showing the prism layers and the faces across the bore.
%   Input: ../data/figS1_domain.csv, figS1_mesh_inlet.csv, figS1_mesh_outlet.csv,
%   fig_designs_outline.csv (cfd/geometry/figS1_mesh_data.py, fig_designs_data.py).
%   Output: ../fig_S1_mesh.pdf and .png.  Tested with MATLAB R2024b.

here = fileparts(mfilename("fullpath"));
dataDir = fullfile(here, "..", "data");
cfg = tims_fig_config();
D = readtable(fullfile(dataDir, "figS1_domain.csv"), "TextType", "string");
O = readtable(fullfile(dataDir, "fig_designs_outline.csv"), "TextType", "string");
ps = design_shape(O, "Phillips");
solid = [0.62 0.62 0.62];
regionCols = [0.30 0.30 0.30; [59 76 192] ./ 255; [229 112 88] ./ 255; [180 4 38] ./ 255; ...
    [0.1 0.55 0.3]; [0.1 0.55 0.3]];

W = 7.0;
fig = figure("Color", "w", "Units", "inches", "Visible", "off", "Position", [1 1 W 4.15], ...
    "PaperPositionMode", "auto");

% (a) domain and sizing regions
ax = axes(fig, "Units", "inches", "Position", [0.55 2.55 6.35 1.40]);
tims_format_axes(ax, cfg, false);
axis(ax, "equal");
for k = 1:height(D)
    r = D(k, :);
    c = regionCols(min(k, size(regionCols, 1)), :);
    if r.shape == "rect"
        rectangle(ax, "Position", [r.x0_mm, r.z0_mm, r.x1_mm - r.x0_mm, r.z1_mm - r.z0_mm], ...
            "EdgeColor", c, "LineWidth", 0.8);
        if k == 2
            text(ax, r.x0_mm, r.z1_mm + 15, sprintf("%s, %g mm", r.region, r.face_mm), ...
                "Color", c, "FontSize", cfg.annotationFontSize, "VerticalAlignment", "bottom");
        elseif k == 3
            text(ax, r.x1_mm, r.z1_mm + 15, sprintf("%s, %g mm", r.region, r.face_mm), ...
                "Color", c, "FontSize", cfg.annotationFontSize, "VerticalAlignment", "bottom", ...
                "HorizontalAlignment", "right");
        end
    else
        rectangle(ax, "Position", [r.xc_mm - r.r_mm, r.zc_mm - r.r_mm, 2 * r.r_mm, 2 * r.r_mm], ...
            "Curvature", [1 1], "EdgeColor", c, "LineWidth", 0.8);
    end
end
plot(ax, ps, "FaceColor", solid, "FaceAlpha", 1, "EdgeColor", "k", "LineWidth", 0.3);
d0 = D(1, :);
text(ax, d0.x1_mm - 40, d0.z1_mm - 40, sprintf("far field %g mm", d0.face_mm), ...
    "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
    "FontSize", cfg.annotationFontSize, "Color", regionCols(1, :));
text(ax, -950, -500, "inflow", "FontSize", cfg.annotationFontSize, "Color", cfg.mutedText, ...
    "VerticalAlignment", "bottom");
text(ax, 4200, -500, "outflow", "FontSize", cfg.annotationFontSize, "Color", cfg.mutedText, ...
    "VerticalAlignment", "bottom", "HorizontalAlignment", "right");
xlim(ax, [-1050 4310]);
ylim(ax, [-580 580]);
xlabel(ax, "{\itx} from the inlet tip (mm)", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
ylabel(ax, "{\itz} (mm)", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
tims_panel_label(ax, "(a)", cfg);

% (b), (c) the mesh on the symmetry plane
wins = {[-10 30; -8 8], [1235 1285; -10 10]};
files = ["figS1_mesh_inlet.csv", "figS1_mesh_outlet.csv"];
labs = ["(b) inlet tip", "(c) outlet tip"];
for k = 1:2
    E = readmatrix(fullfile(dataDir, files(k)));
    w = wins{k};
    keep = max(E(:, [1 3]), [], 2) >= w(1, 1) & min(E(:, [1 3]), [], 2) <= w(1, 2) & ...
        max(E(:, [2 4]), [], 2) >= w(2, 1) & min(E(:, [2 4]), [], 2) <= w(2, 2);
    E = E(keep, :);
    xs = [E(:, 1) E(:, 3) nan(size(E, 1), 1)]';
    zs = [E(:, 2) E(:, 4) nan(size(E, 1), 1)]';
    axm = axes(fig, "Units", "inches", "Position", [0.55 + (k - 1) * 3.30 0.42 3.0 1.40]);
    tims_format_axes(axm, cfg, false);
    axis(axm, "equal");
    plot(axm, ps, "FaceColor", solid, "FaceAlpha", 1, "EdgeColor", "none");
    plot(axm, xs(:), zs(:), "-", "Color", [0.15 0.15 0.15], "LineWidth", 0.2);
    xlim(axm, w(1, :));
    ylim(axm, w(2, :));
    xlabel(axm, "{\itx} from the inlet tip (mm)", "Interpreter", "tex", ...
        "FontSize", cfg.labelFontSize);
    if k == 1
        ylabel(axm, "{\itz} (mm)", "Interpreter", "tex", "FontSize", cfg.labelFontSize);
    end
    text(axm, 0.02, 1.03, labs(k), "Units", "normalized", "FontSize", cfg.panelLabelFontSize, ...
        "FontWeight", "bold", "VerticalAlignment", "bottom");
end

tims_export(fig, fullfile(here, "..", "fig_S1_mesh"), cfg);
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
