function [] = quantify_Cdc14_nucleoli_SPB(Cdc14, Cfi1, NLS, nslice, dT, PixelSize)

%%
data = bfopen([Cdc14 '.tif']);

load([Cfi1 '/assignment'], 'trk_cells', 'trk_nucleoli','I_max2')
load([NLS '/NLSiRFP'], 'Division')

mkdir(Cdc14)
cd(Cdc14)

I_Cfi1 = I_max2;
nframes = length(data{1})/nslice;

%% max projection

I_max = cell(nframes,1);
I_max2 = cell(nframes,1);
for t = 1 : nframes     
	I_max{t} = max(cat(3, data{1}{(t-1)*nslice+1 :t*nslice}),[],3);
    mask = imdilate(trk_cells{t}, strel('disk',3));
    I_bg = double(I_max{t}(~mask));
	bg = median(I_bg);
    I_max2{t} = I_max{t} - bg;
end

%% Quantify signal for each division
% use the nucleoli segmentation and measure the ratio of I_Cdc14 to I_Cfi1
% also measure Cdc14 intensity in the nucleolus/cytosol as well as CV
ndivisions = length(Division);
[height,width] = size(I_max{1});
mkdir('Segmentation')
f = figure; 
set(f,'visible','off');
for n = 1 : ndivisions
    T_division = length(Division{n}.Spindle);
    int_Cdc14 = nan(T_division, 5); % cyto, nucleoli, ratio to Cfi1, CV_all, CV_mother
    int_Cfi1 = nan(T_division, 5); % cyto, nucleoli, total, CV, CV_mother
    for i = 1 : T_division
        i_frame = Division{n}.Budding + i - 1;
        mask_cell = trk_cells{i_frame} == Division{n}.Mother | trk_cells{i_frame} == Division{n}.Daughter;
        mask_mother = trk_cells{i_frame} == Division{n}.Mother;
%         mask_nucleoli = trk_nucleoli{i_frame} == Division{n}.Mother | trk_nucleoli{i_frame} == Division{n}.Daughter;
        mask_nucleoli = mask_cell & imbinarize(mat2gray(I_Cfi1{i_frame}));
        % intensity for Cdc14 and Cfi1
        int_Cdc14(i, 1) = median(I_max2{i_frame}(imerode(mask_cell, strel('disk',2))&~mask_nucleoli));% cytosolic intensity
        int_Cdc14(i, 2) = median(I_max2{i_frame}(mask_nucleoli)); % nucleoli intensity
        int_Cdc14(i, 3) = mean(I_max2{i_frame}(mask_nucleoli)./I_Cfi1{i_frame}(mask_nucleoli)); % nucleoli intensity
        
        int_Cfi1(i, 1) = median(I_Cfi1{i_frame}(imerode(mask_cell, strel('disk',2))&~mask_nucleoli));% cytosolic intensity       
        int_Cfi1(i, 2) = median(I_Cfi1{i_frame}(mask_nucleoli)); % nucleoli intensity
        int_Cfi1(i, 3) = mean(I_Cfi1{i_frame}(imerode(mask_cell, strel('disk',2)))); % mean intensity in cell
        
        I2_cell = double(I_max2{i_frame}(imerode(mask_cell, strel('disk',2))));
        I2_mother = double(I_max2{i_frame}(imerode(mask_mother, strel('disk',2))));
        int_Cdc14(i, 4) = std(I2_cell)/mean(I2_cell);
        int_Cdc14(i, 5) = std(I2_mother)/mean(I2_mother);
        
        I1_cell = double(I_Cfi1{i_frame}(imerode(mask_cell, strel('disk',2))));
        I1_mother = double(I_Cfi1{i_frame}(imerode(mask_mother, strel('disk',2))));
        int_Cfi1(i, 4) = std(I1_cell)/mean(I1_cell);
        int_Cfi1(i, 5) = std(I1_mother)/mean(I1_mother);        
        
        % show the segmentation & quantify Mob1 intensity
        imshowpair(I_Cfi1{i_frame}, I_max2{i_frame},'montage')
        truesize([2*height, 2*width])
        hold on
        b = bwboundaries(imclose(mask_cell, strel('square',3)));
        if ~isempty(b)
            plot(smooth(b{1}(:,2)),smooth(b{1}(:,1)),'g','linewidth',0.5);
            plot(smooth(b{1}(:,2))+width, smooth(b{1}(:,1)),'g','linewidth',0.5);
        end
        d = bwboundaries(mask_nucleoli);
        for j = 1 : length(d)
            plot(smooth(d{j}(:,2)),smooth(d{j}(:,1)),'m','linewidth',0.5);
            plot(smooth(d{j}(:,2))+width, smooth(d{j}(:,1)),'m','linewidth',0.5);
        end
        hold off
        saveas(f, ['Segmentation/Division ' num2str(n) '_frame ' num2str(i_frame)], 'png')
    end
    Division{n}.Cdc14 = int_Cdc14;
    Division{n}.Cfi1 = int_Cfi1;
    
end

%% plot
mkdir('plots')
close all
for n = 1 : ndivisions
    if length(Division{n}.Cdc14) < 9
        Division{n} = [];
        continue;
    end
    if max(Division{n}.Spindle) < 35
        Division{n} = [];
        continue;
    end
    T = dT*((1:1:length(Division{n}.Cdc14)) - Division{n}.SPBinbud + Division{n}.Budding - 1);
    figure; hold on
    yyaxis left
    plot(T, Division{n}.Cdc14(:,2),'o-')
    plot(T, Division{n}.Cdc14(:,1),'k-')

    xlabel('Time since nucleolar segmentation onset (min)')
    ylabel('Cdc14 Intensity (au)')
    
    yyaxis right
    plot(T, Division{n}.Cdc14(:,3),'*-')
    ylabel('Cdc14 to Cfi1 ratio')

    legend({'Cdc14 in nucleolus' 'Cdc14 in cytosol' 'Cdc14 to Cfi1 ratio'}, 'Location', 'northeast')

    
    saveas(gcf, ['plots/Cdc14_intensity' num2str(n)], 'png')
    
    figure; hold on
    yyaxis left
    plot(T, Division{n}.Cdc14(:,4),'o-')
    plot(T, Division{n}.Cdc14(:,5),'x-')

    xlabel('Time since nucleolar segmentation onset (min)')
    ylabel('Cdc14 release (CV)')
    
    yyaxis right
    plot(T, Division{n}.Spindle*PixelSize)
    ylabel('Nucleolar segmentation progress (um)')
    
    saveas(gcf, ['plots/Cdc14_CV' num2str(n)], 'png')
end

%%
close all
clear data
save('Cdc14')
cd('..')