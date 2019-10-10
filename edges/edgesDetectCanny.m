function [E, O] = edgesDetectCanny(I)
% Detect Canny edges with orientation in image.

if size(I, 3) == 3, I = rgb2gray(I); end

E = single(edge(I, 'canny'));
[~, O] = imgradient(I);
O = single(deg2rad(180 - abs(O)));