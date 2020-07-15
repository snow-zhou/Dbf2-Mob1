function y = distancematrix3D(ypos_prev,xpos_prev,zpos_prev,ypos,xpos,zpos)

n_cells_prev = length(ypos_prev);
n_cells = length(ypos);

y = inf(max(n_cells_prev,n_cells));
for m = 1:n_cells_prev
  for n =1:n_cells
    y(m,n) = norm([ypos_prev(m)-ypos(n),xpos_prev(m)-xpos(n),zpos_prev(m)-zpos(n)]);
    if isnan(y(m,n))
        y(m,n) = inf;
    end
  end
end
