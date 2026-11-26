function cfg = tims_fig_config()
%TIMS_FIG_CONFIG House figure style.
%   Matches the figures of the water-surface geometry paper
%   (Water-surface-geometry-reconstruction/code/figures, e.g. generate_figure_08.m): 9 pt tick
%   labels, 10 pt LaTeX axis labels, bold (a)/(b) panel labels inside the top-left corner, boxed
%   axes with outward ticks, filled markers with a dark edge, the cool-warm colour map for
%   magnitudes, and vector PDF plus 600 dpi PNG output.
%
%   Tested with MATLAB R2024b.

cfg = struct;
cfg.pageWidth_in = 6.3;            % full text width of the Wiley template
cfg.resolution_dpi = 600;
cfg.axesFontSize = 9;
cfg.labelFontSize = 10;
cfg.panelLabelFontSize = 10;
cfg.legendFontSize = 8;
cfg.annotationFontSize = 8;
cfg.axesLineWidth = 0.8;
cfg.dataLineWidth = 1.0;
cfg.markerArea = 30;               % scatter area (points^2), as the water-surface paper
cfg.markerEdge = [0.15 0.15 0.15];
cfg.markerEdgeWidth = 0.5;
cfg.refColour = [0.20 0.20 0.20];  % 1:1 and limit lines
cfg.refLineWidth = 0.8;
cfg.errorColour = [0.62 0.62 0.62];
cfg.mutedText = [0.35 0.35 0.35];
cfg.shade = [0.93 0.93 0.93];
end
