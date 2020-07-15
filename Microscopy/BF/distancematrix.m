function y = distancematrix(ypos_prev,xpos_prev,ypos,xpos)

n_cells_prev = length(ypos_prev);
n_cells = length(ypos);

y = inf(max(n_cells_prev,n_cells));
for m = 1:n_cells_prev
  for n =1:n_cells
    y(m,n) = norm([ypos_prev(m)-ypos(n),xpos_prev(m)-xpos(n)]);
  end
end
