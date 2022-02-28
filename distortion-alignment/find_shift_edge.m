function [edge] = find_shift_edge(img,xshift,yshift,zshift)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dim = size(img);
if xshift == 0 && yshift == 0 && zshift > 0
    edge = imtranslate(img,[0 0 zshift-dim(3)]);
elseif xshift == 0 && yshift > 0 && zshift > 0
    edge1 = imtranslate(img,[0 yshift-dim(1) 0]);
    edge2 = imtranslate(img,[0 0 zshift-dim(3)]);
    edge = edge1 + edge2;
elseif xshift == 0 && yshift > 0 && zshift == 0
    edge = imtranslate(img,[0 yshift-dim(1) 0]);
elseif xshift > 0 && yshift == 0 && zshift == 0
    edge = imtranslate(img,[xshift-dim(2) 0 0]);
elseif xshift > 0 && yshift == 0 && zshift > 0
    edge1 = imtranslate(img,[xshift-dim(2) 0 0]);
    edge2 = imtranslate(img,[0 0 zshift-dim(3)]);
    edge = edge1 + edge2;
elseif xshift > 0 && yshift > 0 && zshift == 0
    edge1 = imtranslate(img,[xshift-dim(2) 0 0]);
    edge2 = imtranslate(img,[0 yshift-dim(1) 0]);
    edge = edge1 + edge2;
elseif xshift > 0 && yshift > 0 && zshift > 0
    edge1 = imtranslate(img,[0 0 zshift-dim(3)]);
    edge2 = imtranslate(img,[0 yshift-dim(1) 0]);
    edge3 = imtranslate(img,[xshift-dim(2) 0 0]);
    edge = edge1 + edge2 + edge3;
elseif xshift == 0 && yshift == 0 && zshift == 0
    edge = zeros(dim);
end
end

