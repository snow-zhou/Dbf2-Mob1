function [] = seg_SPB_3D(SPB, BF, nslice)

%%
load([BF '/BF'], 'trk_cells_fix','nframes','drift');
trk_cells = trk_cells_fix;
data = bfopen([SPB '.tif']);

mkdir(SPB)
cd(SPB)

%% max projection
I_zstack = cell(nframes,1);
I_max = cell(nframes,1);
for t = 1 : nframes     
    I_zstack{t} = cat(3, data{1}{(t-1)*nslice+1 :t*nslice});
	I_max{t} = max(cat(3, data{1}{(t-1)*nslice+1 :t*nslice}),[],3);
end

%% segmentation
mkdir('Segmentation3D')
bw_SPB_3D = cell(nframes,1);
bw_SPB = cell(nframes,1);
seg_SPB = cell(nframes,1);
I_max2 = cell(nframes,1);
I_zstack2 = cell(nframes,1);
[height,width] = size(I_max{1});
f = figure;
set(f,'visible','off');
for i = 1 : nframes
    mask = imdilate(trk_cells{i}, strel('disk',3));
    I_bg = double(I_max{i}(~mask));
	bg = median(I_bg);
    I_max2{i} = I_max{i} - bg;
    I_zstack2{i} = I_zstack{i} - bg;
    mask_mCherry = repmat(mask,1,1,nslice)&I_zstack2{i}>45;% this threshold is determined mannually
%     bw_SPB{i} = mask_mCherry&imbinarize(mat2gray(I_max2{i}), 'adaptive');
    bw_SPB_3D{i} = mask_mCherry&loccssegm(I_zstack2{i}, 3, 0.8, 0);
    bw_SPB{i} = max(bw_SPB_3D{i},[],3);
    seg_SPB{i} = bwlabeln(bw_SPB_3D{i});
    imshowpair(I_max2{i}, max(seg_SPB{i},[],3),'montage')
    hold on
    cellID = unique(mask);
    for j = 2 : length(cellID)
        c = bwboundaries(trk_cells{i} == cellID(j));
        plot(smooth(c{1}(:,2)),smooth(c{1}(:,1)),'g','linewidth',0.5);
        plot(smooth(c{1}(:,2))+width, smooth(c{1}(:,1)),'g','linewidth',0.5);
    end
    saveas(f, ['Segmentation3D/segmentation3D_frame_' num2str(i)], 'png')
end

%% tracking
% start from the last frame and track backwards frame by frame
mkdir('Tracking3D')
trk_SPB_3D = cell(nframes,1);
trk_SPB = cell(nframes,1);
trk_SPB_3D{nframes} = seg_SPB{nframes};
nSPBs = length(unique(trk_SPB_3D{end})) - 1;

f = figure;
set(f,'visible','off');
for n = 0 : nframes-1
    if n > 0
        trk_SPB_3D{nframes-n} = trackSPB_3D(seg_SPB{nframes-n},trk_SPB_3D{nframes-n+1}, drift(nframes-n,:));
        % if there are SPBs not assigned, go to the next frame
        mask_unassigned = seg_SPB{nframes-n}&~trk_SPB_3D{nframes-n};
        if sum(sum(sum(mask_unassigned))) > 1 && n > 1
            trk_SPB_backup = trackSPB_3D(seg_SPB{nframes-n},trk_SPB_3D{nframes-n+2}, drift(nframes-n,:)+drift(nframes-n+1,:));
            spbID_new = setdiff(unique(trk_SPB_backup),unique(trk_SPB_3D{nframes-n}));
            for i = 1 : length(spbID_new)
                trk_SPB_3D{nframes-n}(trk_SPB_backup == spbID_new(i)) = spbID_new(i);
            end
        end
    end
    trk_SPB{nframes-n} = max(trk_SPB_3D{nframes-n},[],3);
    imagesc(trk_SPB{nframes-n},[0 nSPBs])
    axis('equal')
    saveas(f, ['Tracking3D/tracking3D_frame_' num2str(nframes-n)], 'png')
end

%% summarize SPB info
% frame-centric id
SPB_ListID = cellfun(@nonzeros, cellfun(@unique, trk_SPB_3D, 'UniformOutput',false), 'UniformOutput',false);
SPB_List = cell(nframes,1);
for i = 1 : nframes
    SPB_List{i} = cell(length(SPB_ListID{i}),1);
    S = regionprops(trk_SPB_3D{i},'Centroid');
    for n = 1 : length(SPB_ListID{i})
        SPB_List{i}{n} = struct('SPB_ID', SPB_ListID{i}(n), ...
                                'Centroid', S(SPB_ListID{i}(n)).Centroid,...
                                'Intensity', max(I_zstack2{i}(trk_SPB_3D{i}==SPB_ListID{i}(n))));
    end
end

% SPB-centric id
SPBs = cell(nSPBs,1);

for i = 1 : nframes
    for n = 1 : length(SPB_ListID{i})
        SPBs{SPB_ListID{i}(n)}{i} = SPB_List{i}{n};
    end
end

%remove empty SPBs from later appeared daughters and convert type
Birthframe = zeros(nSPBs,1);
for n = 1 : nSPBs
    Birthframe(n) = find(~cellfun(@isempty,SPBs{n}),1); % first nonepmty cell
    SPBs{n} = SPBs{n}(Birthframe(n):end); 
end

%
SPBsInfo = cell(nSPBs,1);

for n = 1 : nSPBs    
    Cell_ID = nan(length(SPBs{n}),1);
    % find cellID that SPB resides
    for i = 1 : length(SPBs{n})
        i_frame = nframes - length(SPBs{n}) + i;
        if isempty(SPBs{n}{i})
            SPBs{n}{i} = struct('SPB_ID', SPBs{n}{1}.SPB_ID, ...
                                'Centroid', [NaN, NaN, NaN],...
                                'Intensity', NaN);
            Cell_ID(i) = NaN;
        else
        
            Cell_ID(i) = trk_cells{i_frame}(round(SPBs{n}{i}.Centroid(2)), round(SPBs{n}{i}.Centroid(1)));
        end
    end
    SPBs{n} = cell2mat(SPBs{n});
    SPBsInfo{n} = struct(...
                          'SPB_ID', SPBs{n}(1).SPB_ID, ...
                          'Birthframe', Birthframe(n), ...
                          'Centroid', (reshape([SPBs{n}.Centroid],3,[]))', ...
                          'Intensity', [SPBs{n}.Intensity],...
                          'CellID', Cell_ID);
end

SPBsInfo = cell2mat(SPBsInfo);
%% Plot segmentation
mkdir('Overlay')
f = figure;
set(f,'visible','off');
Color_order = parula(nSPBs);
for n = 1 : nframes
    imshow(I_max2{n}, [])
    hold on
    % plot cells
    cellID = nonzeros(unique(trk_cells{n}));
    for j = 1 : length(cellID)
        b = bwboundaries(trk_cells{n} == cellID(j));
        plot(smooth(b{1}(:,2)),smooth(b{1}(:,1)),'g','linewidth',0.5);
    end
    SPBs_d = zeros(height, width); % draw SPBs
    for i = 1 : length(SPB_ListID{n})
        % make a circle around the centroid of detected SPB
        SPBs_d(round(SPB_List{n}{i}.Centroid(2)), round(SPB_List{n}{i}.Centroid(1))) = SPB_ListID{n}(i);
    end
    SPBs_d = imdilate(SPBs_d, strel('disk', 2));
    for i = 1 : length(SPB_ListID{n})
        c = bwboundaries(SPBs_d == SPB_ListID{n}(i));
        if ~isempty(c)
            plot(smooth(c{1}(:,2)),smooth(c{1}(:,1)),'Color', Color_order(SPB_ListID{n}(i),:) ,'linewidth',0.5);
        end
    end
    
    hold off
    saveas(f, ['Overlay/overlay_frame ' num2str(n)], 'png')
end
%%
close all
clear data
save('SPB')
cd('..')