function y = fillholes (x)
% y = fillbg (x)
%
% Fill holes inside objects in a binary image.
%
% In:
%  x    binary image
%
% Out:
%  y    binary image with filled holes
%
% Antti Niemistö 2006/07/19

if numel(x) == sum(x(:))
  y = x;
  return
end

[labimg,N] = bwlabel(~x);
sizes = histc(labimg(:),0:N);
y = x;
y(sizes(labimg+1)<18000) = 1;
