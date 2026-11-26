function tims_panel_label(ax, labelText, cfg, position, verticalAlignment)
%TIMS_PANEL_LABEL Bold panel label, e.g. "(a)", in normalised axes units.
if nargin < 4 || isempty(position)
    position = [0.03 0.97];
end
if nargin < 5
    verticalAlignment = "top";
end
text(ax, position(1), position(2), labelText, "Units", "normalized", ...
    "HorizontalAlignment", "left", "VerticalAlignment", verticalAlignment, ...
    "FontSize", cfg.panelLabelFontSize, "FontWeight", "bold", ...
    "Interpreter", "none", "Clipping", "off");
end
