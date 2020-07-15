function [] = seg_nucleoli_SPB(filename, BF, nslice)

%%
load([BF '/BF'], 'trk_cells_fix','nframes');
data = bfopen([filename '.tif']);

mkdir(filename)
cd(filename)

%% max projection

I_max = cell(nframes,1);
for t = 1 : nframes     
	I_max{t} = max(cat(3, data{1}{(t-1)*nslice+1 :t*nslice}),[],3);
end

%%
mkdir('Segmentation')
bw_nucleoli_SPB = cell(nframes,1);
I_max2 = cell(nframes,1);
[~,width] = size(I_max{1});
f = figure;
set(f,'visible','off');
for i = 1 : nframes
    mask = imdilate(trk_cells_fix{i}, strel('disk',3));
    I_bg = double(I_max{i}(~mask));
	bg = median(I_bg);
    I_max2{i} = I_max{i} - bg;
%     bw_nucleoli_SPB{i} = mask&I_max2{i}>20&imbinarize(mat2gray(I_max2{i}),0.2);
    mask_mCherry = mask&I_max2{i}>30;% this threshold is determined mannually
    bw_nucleoli_SPB{i} = mask_mCherry&imbinarize(mat2gray(I_max2{i}), 'adaptive');
    bw_nucleoli_SPB{i} = bwmorph(bw_nucleoli_SPB{i}, 'clean');
    bw_nucleoli_SPB{i} = bwareaopen(bw_nucleoli_SPB{i},3);
    imshowpair(I_max2{i}, bw_nucleoli_SPB{i},'montage')
    hold on
    cellID = unique(mask);
    for j = 2 : length(cellID)
        c = bwboundaries(trk_cells_fix{i} == cellID(j));
        plot(smooth(c{1}(:,2)),smooth(c{1}(:,1)),'g','linewidth',0.5);
        plot(smooth(c{1}(:,2))+width, smooth(c{1}(:,1)),'g','linewidth',0.5);
    end
    saveas(f, ['Segmentation/segmentation_frame_' num2str(i)], 'png')
end

%%
close all
clear data
save('nucleoli')
cd('..')