function map = tims_coolwarm(m)
%TIMS_COOLWARM Cool-warm colour map of the water-surface paper (local_coolwarm there).
if nargin < 1 || isempty(m)
    m = 256;
end
anchors = [59 76 192; 84 112 222; 129 164 251; 180 205 251; ...
    221 221 221; 241 184 156; 229 112 88; 203 62 56; 180 4 38] ./ 255;
map = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, m), "linear");
map = max(0, min(1, map));
end
