%% tracking
% start from the last frame and track backwards frame by frame
trk_cells = cell(nframes,1);

f = figure;
for n = 0 : nframes-1
    if n == 0
        trk_cells{nframes} = bwlabel(bw_cells{nframes});
    else
        trk_cells{nframes-n} = trackcells(bw_cells{nframes-n},trk_cells{nframes-n+1});
    end
    imagesc(trk_cells{nframes-n})
    saveas(f, ['segmentation_frame_' num2str(nframes-n)], 'png')
end