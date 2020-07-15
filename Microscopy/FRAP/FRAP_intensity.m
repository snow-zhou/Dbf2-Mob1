% function [FRAP_info, cell_info] = FRAP_intensity(filename,mask)

%% 
clear
close all
% addpath('/Users/snow/lab/Matlab/FRAP')
addpath('/Users/snow/Dropbox (MIT)/Scripts/Matlab/Downloaded/bfmatlab');

%User Input
dirname='/Users/snow/Dropbox (MIT)/Microscopy Raw Data/Koch facility/OMX/20171019/TAB6-1_dbf2';
% path2mask = '/Users/snow/lab/Raw data/Microscopy/FRAP/2016-04-08 FRAP Snow';
ref_channel = 2;

cd(dirname)
% list_file = dir('*nd2');

filename = 'Mob1-eGFP_dbf2_003.dv';

%%
data = bfopen(filename);

% imagesc(data{1}{1});

%% extract timestamps for the images (OMX series)
metadata = data{1, 2};
% gather information in the time lapse
C = str2double(data{1,1}{1,2}(find(data{1,1}{1,2} == 'C')+4)); % number of channels
Z = str2double(data{1,1}{1,2}(find(data{1,1}{1,2} == 'Z')+4)); % number of z stacks
T = length(data{1})/C/Z;% number of time points

% extract timestamps
timestamp = zeros(T,1);
for t = 1 : T
	timestamp(t) = metadata.get(['Extended header Z0 W0 T' num2str(t-1) ':timeStampSeconds']);
end


%% outline cells with frame 1
minA = 100;
maxA = 4000; % rule out big clusters
    
% take the sum of the z stack
[X,Y] = size(data{1}{1});
I_sum = cell(T,C);
for t = 1 : T
    for c = 1 : C       
        I_sum{t,c} = sum(cat(3, data{1}{(t-1)*Z*C+c:2:t*C*Z}),3);
    end
end

Im = I_sum{1}{1};

Im_ad = imadjust(Im);

level=graythresh(Im_ad);
if level<0.01
    level = 0.015;
end
I_mask = imbinarize(Im_ad,level*2);
% I_mask=imcomplement(I_mask);
% [ime,~] = edge(Im,'log');
% ime(I_mask)=0;
% cc = bwconncomp(ime,8);
cc = bwconncomp(I_mask,8);
stats = regionprops(cc,'Area','Centroid','Extrema');
idx = find([stats.Area] > minA & [stats.Area] < maxA);

% break control for failure in segmentation
if isempty(idx)
    FRAP_info = []; 
    cell_info = [];
    return
end

ncells = length(idx);
cell_mask = cell(ncells,1);
for i = 1 : ncells
    cell_mask{i} = ismember(labelmatrix(cc),idx(i));
end

%% Identify FRAP area

% compare the sum intensity image for channel at time point 1-3
I_diff = (2*I_sum{2,1} - I_sum{3,1} - I_sum{1,1});
% figure; imshow(I_diff,[])

smooth_filter = fspecial('gaussian', [9 9], 2);
Im_smoothed = imfilter(Im_ad, smooth_filter);
[ime,~]=edge(Im_smoothed,'log');
I_mask = imfill(ime,'holes');


cell_FRAP = 0;
for i = 1 : ncells
    if sum(cell_mask{i}(mask>0))>0
        cell_FRAP = i;
    end
end

if cell_FRAP == 0
    FRAP_info = []; 
    cell_info = [];
    return
end

mask_FRAP = cell_mask{cell_FRAP};
mask_FRAP(imcomplement(logical(mask))) = 0;            

%% %% plot the mask with colored boudries

figure
imshow(Im, [])
hold on
for i = 1 : ncells
    b = bwboundaries(cell_mask{i});
    c = stats(idx(i)).Extrema(1,:) - 10;
    plot(b{1}(:,2),b{1}(:,1),'g','linewidth',2);
    text(c(1),c(2),num2str(i),'backgroundcolor','g');
%     x = stats(idx(i)).Extrema(:,1);
%     y = stats(idx(i)).Extrema(:,2);
%     patch(x, y, 'g', 'FaceAlpha', alpha)
end
b = bwboundaries(logical(mask));
plot(b{1}(:,2),b{1}(:,1),'r','linewidth',2);
hold off
saveas(gcf, ['analyzed matlab/' filename '_1_mask.png'], 'png')
close

figure
imshow(data{1}{2,1}, [])
hold on
for i = 1 : ncells
    b = bwboundaries(cell_mask{i});
    c = stats(idx(i)).Extrema(1,:) - 10;
    plot(b{1}(:,2),b{1}(:,1),'g','linewidth',2);
    text(c(1),c(2),num2str(i),'backgroundcolor','g');
end
b = bwboundaries(logical(mask));
plot(b{1}(:,2),b{1}(:,1),'r','linewidth',2);
hold off
saveas(gcf, ['analyzed matlab/' filename '_2_mask.png'], 'png')
close


%% estimate background
bg_pixels = double(Im(I_mask<1));
[f,xi] = ksdensity(bg_pixels);
[~, index] = max(f);
bg = xi(index);

%% read intesnity of all cells (registrate images first)
intensity = zeros(T, ncells);
Int_FRAP = zeros(T,1);
for t = 1 : T
    for i = 1 : ncells
        I = data{1}{t};
        % register images to the after bleaching 
        I_cell = I(cell_mask{i});
        intensity(t, i) = sum(I_cell - bg);
    end
    Int_FRAP(t) = sum(I(mask_FRAP) - bg);
end


%%
% summarize for all cells    
            
for i = 1 : ncells
    cell_info(i,1) = struct(...
                'cell_mask', cell_mask{i},...
                'time', timestamp,...
                'Intensity', intensity(:,i),...
                'FRAP', 0,...
                'file', filename);
    if i == cell_FRAP
        cell_info(i).FRAP = 1;
    end
end

% summarize FRAP info
%single normalization
I_FRAP_1 = Int_FRAP/Int_FRAP(1)*intensity(1,cell_FRAP)./intensity(:,cell_FRAP);
% double normalization
I_FRAP_2 = (I_FRAP_1-I_FRAP_1(2))./(I_FRAP_1(1)-I_FRAP_1(2));
% fit FRAP curve
[parameters,fit_result] = FRAP_fit([timestamp I_FRAP_2]);
FRAP_info = struct(...
                   'time', timestamp,...
                   'cell_area',cell_mask{cell_FRAP},...
                   'FRAP_area', mask_FRAP,...
                   'total_intensity', intensity(:,cell_FRAP),...
                   'FRAP_intensity', Int_FRAP,...
                   'FRAP_singleNorm', I_FRAP_1,...
                   'FRAP_doubleNorm', I_FRAP_2,...
                   'fit_result', fit_result,...
                   'parameters', parameters);
               
end