function [N,ypos,xpos,zpos] = findcentroids3D (x)
% [y,lab,N,ypos,xpos] = findcentroids (x)
%
% Find the centroids of objects in a binary image.
%
% In:
%  x      binary image (3D)
%
% Out:
%  y      labeled centroids
%  N      number of objects
%  ypos   y-coordinates
%  xpos   x-coordinates
%  zpos   z-coordinates


S = regionprops(x,'Centroid');
tmp = cell2mat(struct2cell(S));
ypos = round(tmp(2:3:end));
xpos = round(tmp(1:3:end));
zpos = round(tmp(3:3:end));

N = length(S);
end
