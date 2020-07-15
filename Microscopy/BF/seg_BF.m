function [] = seg_BF(filename,nslice)

%
data = bfopen([filename '.tif']);

mkdir(filename)
cd(filename)
%% segmentation
mkdir('Segmentation')
nframes = length(data{1})/nslice;

seg_cells = cell(nframes,1);
f = figure;
set(f,'visible','off');
for i = 1 : nframes    
    seg_cells{i} = segmcells2(data{1}(nslice*(i-1)+1 : nslice*i));
    f = imagesc(seg_cells{i});
    axis('equal')
    saveas(f, ['Segmentation/segmentation_frame ' num2str(i)], 'png')
end

%% tracking
% start from the last frame and track backwards frame by frame
mkdir('Tracking')
trk_cells = cell(nframes,1);
trk_cells{nframes} = seg_cells{nframes};
ncells = length(unique(seg_cells{nframes})) - 1;

f = figure;
set(f,'visible','off');
for n = 0 : nframes-1
    if n > 0
        trk_cells{nframes-n} = trackcells(seg_cells{nframes-n},trk_cells{nframes-n+1});
    end
    imagesc(trk_cells{nframes-n},[0 ncells])
    axis('equal')
    saveas(f, ['Tracking/tracking_frame_' num2str(nframes-n)], 'png')
end

%% Plot segmentation
mkdir('Overlay')
f = figure;
set(f,'visible','off');
for n = 1 : nframes
    midslice = round(nslice/2+1);
    imshow(data{1}{nslice*(n-1)+midslice}, [])
    hold on
    for i = 1 : ncells
        b = bwboundaries(trk_cells{n} == i);
        if isempty(b)
            continue
        end
        c = max(b{1});
        plot(smooth(b{1}(:,2)),smooth(b{1}(:,1)),'g','linewidth',1);
        text(c(2),c(1),num2str(i));
    end
    hold off
    saveas(f, ['Overlay/overlay_frame ' num2str(n)], 'png')
end

%%
close all
clear data;
save('BF')
cd('..')
