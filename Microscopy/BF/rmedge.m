function y = rmedge (x,conn,e)
% y = rmedge (x,conn,e)
%
% Remove objects that touch the boundaries of a binary image.
%
% In:
%  x      binary image
%  conn   connectivity
%  e      edge ffrom which to remove objects (all, top, bottom, left, right)
%
% Out:
%  y      binary image with objects touching the borders removed
%
% Antti Niemistö 2010/04/27

if nargin < 2
  conn = 8;
end

if nargin < 3
  e = 'all';
end

[m,n] = size(x);
switch e
 case 'all'
  xx = ones(m+2,n+2);
  xx(2:m+1,2:n+1) = x;
  labeled = bwlabel(xx,conn);
  xx(labeled==labeled(1,1)) = 0;
 case 'top'
  xx = zeros(m+2,n+2);
  xx(2:m+1,2:n+1) = x;
  xx(1,:) = 1;
  labeled = bwlabel(xx,conn);
  xx(labeled==labeled(1,1)) = 0;
 case 'bottom'
  xx = zeros(m+2,n+2);
  xx(2:m+1,2:n+1) = x;
  xx(end,:) = 1;
  labeled = bwlabel(xx,conn);
  xx(labeled==labeled(end,end)) = 0;
 case 'left'
  xx = zeros(m+2,n+2);
  xx(2:m+1,2:n+1) = x;
  xx(:,1) = 1;
  labeled = bwlabel(xx,conn);
  xx(labeled==labeled(1,1)) = 0;
 case 'right'
  xx = zeros(m+2,n+2);
  xx(2:m+1,2:n+1) = x;
  xx(:,end) = 1;
  labeled = bwlabel(xx,conn);
  xx(labeled==labeled(end,end)) = 0;
end

y = xx(2:end-1,2:end-1);
