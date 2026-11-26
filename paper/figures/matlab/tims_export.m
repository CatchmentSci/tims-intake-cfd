function tims_export(fig, outputBase, cfg)
%TIMS_EXPORT Write <outputBase>.pdf (vector) and <outputBase>.png (600 dpi).
%   Text not drawn by the LaTeX interpreter is set in Arial first. MATLAB leaves its default
%   Helvetica unembedded in PDF output, which ASCE production rejects; Arial (TrueType) is
%   embedded as a subset.
objs = findall(fig, "-property", "FontName");
for k = 1:numel(objs)
    o = objs(k);
    if isprop(o, "Interpreter") && o.Interpreter == "latex"
        continue
    end
    o.FontName = "Arial";
end
exportgraphics(fig, outputBase + ".pdf", "ContentType", "vector");
exportgraphics(fig, outputBase + ".png", "Resolution", cfg.resolution_dpi);
fprintf("wrote %s.pdf and .png\n", outputBase);
end
