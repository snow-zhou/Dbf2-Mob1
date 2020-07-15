function [N,ypos,xpos] = findcentroids (x)
% [y,lab,N,ypos,xpos] = findcentroids (x)
%
% Find the centroids of objects in a binary image.
%
% In:
%  x      binary image
%
% Out:
%  y      labeled centroids
%  N      number of objects
%  ypos   y-coordinates
%  xpos   x-coordinates


S = regionprops(x,'Centroid');
tmp = cell2mat(struct2cell(S));
ypos = round(tmp(2:2:end));
xpos = round(tmp(1:2:end));

N = length(S);
end
