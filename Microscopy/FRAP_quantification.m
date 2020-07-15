clear
close all

addpath(genpath('/Microscopy'))

list_file = dir('*tif');

%%
n = 1;
filename = list_file(n).name(1:end-4);
nCh = 2;
dT = 1; % 1s/frame
data = bfopen([filename '.tif']);
nframes = length(data{1})/nCh;
bleach = 3;

%%
mkdir([filename '/segmentation'])
figure
bw = cell(nframes,1);
bw1 = cell(nframes,1);
I_cyto = zeros(nframes,1);
bg = zeros(nframes,1);
for i = 1 : nframes
    I1 = data{1}{nCh*(i-1)+1};
    intensity = double(I1(:));
    [f,xi] = ksdensity(intensity);
    [~,loc] = findpeaks(f,'MinPeakProminence',2e-4);
    I_cyto(i) = xi(loc(end));
    bg(i) = xi(loc(find(xi(loc)>10,1,'first')));

    bw{i} = I1>mean([I_cyto(i) bg(i)])*1.1;
    bw{i} = bwareaopen(bw{i},5);
    imshowpair(I1,bw{i},'montage')
    saveas(gcf, [filename '/segmentation/ch1_' num2str(i)],'png')
    
    bw1{i} = bw{i}&loccssegm(mat2gray(I1-bg(i)));
    bw1{i} = bwareaopen(bw1{i},2);
    bw1{i}(I1<2*I_cyto(i)) = 0;
    imshowpair(I1,bw1{i},'montage')
    saveas(gcf, [filename '/segmentation/ch1_seg1_' num2str(i)],'png')

end

%%
%use the segmentation of frame 1 as reference
CC = regionprops(bw1{1});
L = bwlabel(imdilate(bw1{1},strel('disk',5)));
I_cyto = zeros(nframes,1);
I_spb = zeros(nframes,length(CC));
figure
for i = 1 : nframes
    I1 = data{1}{nCh*(i-1)+1} - bg(i);
    %I_cyto
    intensity = double(I1(bw{1}));
    [f,xi] = ksdensity(intensity);
    [~, idx] = max(f);
    I_cyto(i) = xi(idx);
    %I_spb
    imshow(I1,[])
    hold on
    SPB = L;
    if sum(sum((bw1{i}+L)>0)) > sum(sum(L>0))
        SPB = bwlabel(imdilate(bw1{i},strel('disk',5)));
    end
    for n = 1 : min([length(CC) 2])
        intensity2 = double(I1(SPB==n));
        [f,xi] = ksdensity(intensity2);
        [~,loc] = findpeaks(f,'MinPeakProminence',2e-5);
        I_spb(i,n) = xi(loc(end));
        c = bwboundaries(SPB==n);
        plot(smooth(c{1}(:,2)),smooth(c{1}(:,1)),'linewidth',1);
    end
    hold off
    saveas(gcf, [filename '/segmentation/SPB_' num2str(i)],'png')
end

%% plot
close all
%raw intensity
figure; hold on
plot(I_spb(2:end,:))
plot(I_cyto(2:end))
saveas(gcf, [filename '/raw intensity'],'png')

%enrichment at SPB
E_spb = I_spb./I_cyto-1;
figure; hold on
plot(E_spb(2:end,:))
saveas(gcf, [filename '/enrichment at SPB'],'png')

%% fit
%identify the bleached SPB
[~, id] = max(E_spb(bleach-1,:)-E_spb(bleach,:));
I_FRAP_1 = E_spb(2:end,id);
% double normalization
I_FRAP_2 = (I_FRAP_1-I_FRAP_1(2))./(I_FRAP_1(1)-I_FRAP_1(2));
% fit FRAP curve
timestamp = (1:1:length(I_FRAP_2))'*dT;
[parameters,fit_result] = FRAP_fit([timestamp I_FRAP_2]);
figure; hold on
plot(timestamp, I_FRAP_2, 'ko');
plot(fit_result(:,1), fit_result(:,2));
saveas(gcf, [filename '/Fit'],'png')

FRAP_info = struct(...
                   'time', timestamp,...
                   'cell_intensity', I_cyto,...
                   'SPB_intensity', I_spb,...
                   'FRAP_singleNorm', I_FRAP_1,...
                   'FRAP_doubleNorm', I_FRAP_2,...
                   'fit_result', fit_result,...
                   'HT', -reallog(2)./parameters(2,2),...
                   'F_mobile', parameters(1,2));

save([filename '/analysis'])