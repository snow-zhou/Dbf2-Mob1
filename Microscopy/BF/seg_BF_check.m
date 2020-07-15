function [] = seg_BF_check(BF,nslice)

%
data = bfopen([BF '.tif']);

mkdir(BF)
cd(BF)
%% segmentation
mkdir('Segmentation')
nframes = length(data{1})/nslice;
midslice = round(nslice/2+1);

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
drift = zeros(nframes-1, 2);

f = figure;
set(f,'visible','off');
for n = 0 : nframes-1
    if n > 0
        % cross correlation of two frames to calculate drift
        % use image not segmentation in case of missegmentation
        [output, ~]=dftregistration(fft2(data{1}{nslice*(nframes-n-1)+midslice}),fft2(data{1}{nslice*(nframes-n)+midslice}),10);
        drift(nframes-n, 1) = output(3); %shift in row
        drift(nframes-n, 2) = output(4); %shift in column
        trk_cells{nframes-n} = trackcells(seg_cells{nframes-n},trk_cells{nframes-n+1}, drift(nframes-n,:));
        % if there are cells not assigned, go to the next frame
        mask_unassigned = seg_cells{nframes-n}&~trk_cells{nframes-n};
        dN = n;
        loop = 1;
        while sum(sum(mask_unassigned)) > 250 && dN > 1 && loop < 3
            trk_cells_backup = trackcells(seg_cells{nframes-n},trk_cells{nframes-n+1+loop}, sum(drift(nframes-n:nframes-n+loop,:)));
            cellID_new = setdiff(unique(trk_cells_backup),unique(trk_cells{nframes-n}));
            for i = 1 : length(cellID_new)
                trk_cells{nframes-n}(trk_cells_backup == cellID_new(i)) = cellID_new(i);
            end
            mask_unassigned = seg_cells{nframes-n}&~trk_cells{nframes-n};
            dN = dN - 1;
            loop = loop+1;
        end
    end
    imagesc(trk_cells{nframes-n},[0 ncells])
    axis('equal')
    saveas(f, ['Tracking/tracking_frame_' num2str(nframes-n)], 'png')
end

%% check and fix misseggmentations
% check cell area and find abnormality
% frame-centric id
cellListID = cellfun(@nonzeros, cellfun(@unique, trk_cells, 'UniformOutput',false), 'UniformOutput',false);
cellList = cell(nframes,1);
for i = 1 : nframes
    cellList{i} = cell(length(cellListID{i}),1);
    S = regionprops(trk_cells{i},'Centroid','Area');
    for n = 1 : length(cellListID{i})
        cellList{i}{n} = struct('CellID', cellListID{i}(n), ...
                            'Centroid', S(cellListID{i}(n)).Centroid,...
                            'Area', S(cellListID{i}(n)).Area);        
    end
end

% cell-centric id
ncells = max(cellListID{end});
cells = cell(ncells,1);

for i = 1 : nframes
    for n = 1 : length(cellListID{i})
        cells{cellListID{i}(n)}{i} = cellList{i}{n};
    end
end

% Calculate cell areas
Cell_Areas = nan(ncells, nframes);
for n = 1 : ncells
    for i = 1 : nframes
        if ~isempty(cells{n}{i})
            Cell_Areas(n, i) = cells{n}{i}.Area;
        end
    end
end

%%
% identify problematic frames and replace with the frame before
close all
% figure;hold on
thresh = 100; % threshold of area changes
minA_mother = 450; % minimal area for a mother 
anomaly = cell(ncells,1);
trk_cells_fix = trk_cells;
for n = 1 : ncells
    g = figure;
    set(g,'visible','off'); hold on
    plot(Cell_Areas(n,:))
    % yeast cells should only increae cell volume, especially mothers
    if Cell_Areas(n,1) > minA_mother
        [f, xi] = ksdensity(Cell_Areas(n,:));
        [~, idx] = max(f);
        A_mother = xi(idx); % assume more frames are segmented correctly
        anomaly{n} = find(abs(Cell_Areas(n,:) - A_mother) > thresh);
    else
    % find mis-segmented time points (anomaly areas)
        first_frame = find(~isnan(Cell_Areas(n,:)), 1, 'first');
        Cell_Areas(n, find(isnan(Cell_Areas(n, first_frame:end)))+first_frame-1) = 0;
        anomaly{n} = find(abs(Cell_Areas(n,:)' - smooth(Cell_Areas(n,:))) > thresh);
    end
    
	if ~isempty(anomaly{n})
        plot(anomaly{n}, Cell_Areas(n,anomaly{n}),'x')
        for i = 1 : length(anomaly{n})
            if anomaly{n}(i) == 1
                continue
            end
            trk_cells_fix{anomaly{n}(i)}(trk_cells_fix{anomaly{n}(i)}(:) == n) = 0;
            %shift the segmentation of previous frame based on drift
            [I_max, J_max] = size(trk_cells{1});
            [I,J] = ind2sub(size(trk_cells{1}),find(trk_cells_fix{anomaly{n}(i)-1}(:) == n));  
            IJ_corrected = [I-round(drift(anomaly{n}(i)-1,1)) J-round(drift(anomaly{n}(i)-1,2))];
            %in case cell drift to boarder
            if isempty(IJ_corrected)
                continue
            end
            IJ_corrected(IJ_corrected(:,1) > I_max | IJ_corrected(:,2) > J_max, :) = [];
            IJ_corrected(IJ_corrected(:,1) < 1 | IJ_corrected(:,2) < 1, :) = [];
            trk_cells_fix{anomaly{n}(i)}(sub2ind(size(trk_cells{1}),IJ_corrected(:,1),IJ_corrected(:,2))) = n;
        end
	end
    saveas(g, ['Tracking/CellArea_cell ' num2str(n)], 'png')
end

%% Plot segmentation
mkdir('Overlay2')
f = figure;
set(f,'visible','off');
for n = 1 : nframes
    imshow(data{1}{nslice*(n-1)+midslice}, [])
    hold on
    for i = 1 : ncells
        b = bwboundaries(trk_cells_fix{n} == i);
        if isempty(b)
            continue
        end
        c = max(b{1});
        plot(smooth(b{1}(:,2)),smooth(b{1}(:,1)),'g','linewidth',1);
        text(c(2),c(1),num2str(i));
    end
    hold off
    saveas(f, ['Overlay2/overlay_frame ' num2str(n)], 'png')
end

%%
close all
clear data;
save('BF')
cd('..')
