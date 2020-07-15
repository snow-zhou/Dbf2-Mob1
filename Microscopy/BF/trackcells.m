function trk_cells = trackcells(seg_cells,trk_cells_pre, drift)

    [n_cells,ypos_tmp,xpos_tmp] = findcentroids(seg_cells);
    sizes_tmp = histcounts(seg_cells,1:n_cells+1); % cell sizes
  % Read values from previous frame
    [n_cells_prev,ypos_prev,xpos_prev] = findcentroids(trk_cells_pre);
    sizes_prev = histcounts(trk_cells_pre,1:n_cells_prev+1); % cell sizes

%     drift = 0;
    ypos_estim = ypos_prev + drift(1);
    xpos_estim = xpos_prev + drift(2);
    
  % Estimate that small cells will grow 20 per cent
  sizes_estim = sizes_prev; 
  medsize = median(sizes_prev);
  idx = (sizes_prev<.75*medsize);
  sizes_estim(idx) = 0.8*sizes_estim(idx);

  % Calculate the distance matrix based on estimated locations (may be equal
  % to previous locations) and new locations
  distmatrix = distancematrix(ypos_estim,xpos_estim,ypos_tmp,xpos_tmp);

  % Calculate sizematrices
  sizematrix = inf(max(length(sizes_estim),length(sizes_tmp)));
  sizeratiomatrix = sizematrix;
  for m = 1:length(sizes_estim)
    for n = 1:length(sizes_tmp)
      sizematrix(m,n) = sizes_estim(m)-sizes_tmp(n);
      sizeratiomatrix(m,n) = sizes_tmp(n)/sizes_estim(m);
    end
  end

  % Combine into a single weightmatrix
  weightmatrix = distmatrix+.25*abs(sizematrix);


  % Make unlikely matches impossible
  weightmatrix(distmatrix>20) = inf;
  weightmatrix(sizeratiomatrix>2) = inf;
  weightmatrix(sizeratiomatrix<.5) = inf;

  % Impose stricter constraints for full-grown cells
%   idx = inf(max(length(sizes_tmp),length(sizes_prev)));
%   tmp = repmat(sizes_tmp>1.25*medsize,size(weightmatrix,2),1);
%   idx(1:size(tmp,1),1:size(tmp,2)) = tmp;
%   weightmatrix(sizeratiomatrix>2&idx) = inf;
%   weightmatrix(sizeratiomatrix<.5&idx) = inf;

  % Solve the matching problem by Hungarian algorithm
  match = hungarian(weightmatrix);

  % Find the mapping; if there are new cells, they get the label zero at this
  % point
  mapping = zeros(1,n_cells);
  for k = 1:n_cells
    idx = find(match(:,k));
    if isempty(idx)
      idx = 0;
    end
    mapping(k) = idx;
  end
  
  % Make further assignments based on a weight threshold
  thresh = median(weightmatrix(logical(match)))+ ...
           3*mad(weightmatrix(logical(match)));
  unassigned = find(mapping==0);
  for k = unassigned
    if min(weightmatrix(:,k)) < thresh
      idx = find(weightmatrix(:,k)==min(weightmatrix(:,k)));
      if ~any(mapping==idx)
        mapping(k) = idx;
      end
    end
  end

  % Build the labeled image
  trk_cells = zeros(size(seg_cells));
  for k = 1:n_cells
      trk_cells(seg_cells == k) = mapping(k);
  end

end

function [y,med] = mad (x)

med = median(x(:));

y = median(abs(x(:)-med));
end