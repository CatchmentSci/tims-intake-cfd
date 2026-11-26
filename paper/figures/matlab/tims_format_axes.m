function tims_format_axes(ax, cfg, makeSquare)
%TIMS_FORMAT_AXES Apply the house axes style (as format_axes in the water-surface paper).
if nargin < 3
    makeSquare = true;
end
hold(ax, "on");
grid(ax, "off");
box(ax, "on");
if makeSquare
    axis(ax, "square");
end
set(ax, "FontSize", cfg.axesFontSize, "LineWidth", cfg.axesLineWidth, ...
    "TickDir", "out", "TickLabelInterpreter", "tex", "Layer", "top");
end
