function [] = quantify_Cfi1(Cfi1, SPB, nslice, PixelSize)

%%
data = bfopen([Cfi1 '.tif']);

load([SPB '/assignment'], 'trk_cells', 'trk_SPB', 'Division')

mkdir(Cfi1)
cd(Cfi1)

nframes = length(data{1})/nslice;

%% max projection & Segment Cfi1 signal
mkdir('Segmentation')
bw_nucleoli = cell(nframes,1);
I_max = cell(nframes,1);
I_max2 = cell(nframes,1);
f = figure;
set(f,'visible','off');
for t = 1 : nframes     
	I_max{t} = max(cat(3, data{1}{(t-1)*nslice+1 :t*nslice}),[],3);
    mask = imdilate(trk_cells{t}, strel('disk',3));
    I_bg = double(I_max{t}(~mask));
	bg = median(I_bg);
    I_max2{t} = I_max{t} - bg;
    
    bw_nucleoli{t} = mask&imbinarize(mat2gray(I_max2{t}));
    bw_nucleoli{t} = bwmorph(bw_nucleoli{t}, 'clean');
%     bw_nucleoli{t} = bwareaopen(bw_nucleoli{t},3);
    imshowpair(I_max2{t}, bw_nucleoli{t},'montage')
    hold on
    cellID = unique(mask);
    [~,width] = size(I_max{1});
    for j = 2 : length(cellID)
        c = bwboundaries(trk_cells{t} == cellID(j));
        plot(smooth(c{1}(:,2)),smooth(c{1}(:,1)),'g','linewidth',0.5);
        plot(smooth(c{1}(:,2))+width, smooth(c{1}(:,1)),'g','linewidth',0.5);
    end
    saveas(f, ['Segmentation/segmentation_frame_' num2str(t)], 'png')
end

trk_nucleoli = cellfun(@(x,y) x.*y, trk_cells, bw_nucleoli, 'UniformOutput',false);
nucleoliListID = cellfun(@nonzeros, cellfun(@unique, trk_nucleoli, 'UniformOutput',false), 'UniformOutput',false);

%% quantify nucleoli seggregation (shape elongation & move to bud)

ndivisions = length(Division);
[height,width] = size(I_max{1});
mkdir('Division')
f = figure; 
set(f,'visible','off');
for n = 1 : ndivisions
    T_division = length(Division{n}.Spindle);
    nucleoli_length = nan(T_division, 2); 
    nucleoli_number = nan(T_division, 2); 
    for i = 1 : T_division
        i_frame = Division{n}.Metaphase + i - 1;
        if i_frame > nframes
            continue
        end
        mask_cell = trk_cells{i_frame} == Division{n}.Mother | trk_cells{i_frame} == Division{n}.Daughter;
        mask_nucleoli = trk_nucleoli{i_frame} == Division{n}.Mother | trk_nucleoli{i_frame} == Division{n}.Daughter;
        
        nucleoli = regionprops(double(mask_nucleoli), 'MajorAxisLength', 'BoundingBox');
        if ~isempty(nucleoli)
            nucleoli_length(i,1) = nucleoli.MajorAxisLength;
            nucleoli_length(i,2) = sqrt(nucleoli.BoundingBox(3)^2 + nucleoli.BoundingBox(4)^2);
        end
        
        nucleoli2 = regionprops(imdilate(mask_nucleoli, strel('disk', 1)), 'Centroid', 'Area');
        nucleoli_number(i,1) = length(nucleoli2);
        if length(nucleoli2) == 2
            if abs(nucleoli2(1).Area/nucleoli2(2).Area-1) < 0.5 % close to equal size
                nucleoli_number(i,2) = 2;
            else
                nucleoli_number(i,2) = 1;
            end
        else
            nucleoli_number(i,2) = 1;
        end

        
        % show the segmentation
        imshow(I_max2{i_frame},[])
        truesize([1*height, 1*width])
        hold on
        b = bwboundaries(imclose(mask_cell, strel('square',3)));
        if ~isempty(b)
            plot(smooth(b{1}(:,2)),smooth(b{1}(:,1)),'g','linewidth',0.5);
        end
        d = bwboundaries(mask_nucleoli);
        for j = 1 : length(d)
            plot(smooth(d{j}(:,2)),smooth(d{j}(:,1)),'y','linewidth',0.5);
        end
        
        % draw SPB
        SPBs_d = zeros(height, width); % draw SPBs       
        S = regionprops(trk_SPB{i_frame} == Division{n}.dSPB, 'Centroid');
        if ~isempty(S)
            SPBs_d(round(S.Centroid(2)), round(S.Centroid(1))) = 1; % dSPB
        end
        S = regionprops(trk_SPB{i_frame} == Division{n}.mSPB, 'Centroid');
        if ~isempty(S)
            SPBs_d(round(S.Centroid(2)), round(S.Centroid(1))) = 2; % mSPB
        end
        SPBs_d = imdilate(SPBs_d, strel('disk', 2));     
        
        SPBsID = nonzeros(unique(SPBs_d));
        Color_order = get(gca, 'colororder');
        for j = 1 : length(SPBsID)
            c = bwboundaries(SPBs_d == SPBsID(j));
            plot(smooth(c{1}(:,2)),smooth(c{1}(:,1)),'color', Color_order(SPBsID(j),:),'linewidth',0.5);
        end
        
        hold off
        saveas(f, ['Division/Division ' num2str(n) '_frame ' num2str(i_frame)], 'png')
        
    end
    Division{n}.nucleoli = nucleoli_length*PixelSize;
    Division{n}.n_nucleoli = nucleoli_number;
    
    %find the first frame nucleoli enters bud (continuesly)
    Division{n}.NucleoliInBud = find(~cellfun(@isempty, cellfun(@(x) find(x==Division{n}.Daughter), nucleoliListID, 'UniformOutput', 0)),1);
%     Division{n}.NucleoliInBud = nucframes(find(diff(nucframes)==1,1));
    
    %find the first frame nucleoli segregated
    Division{n}.NucleoliSeg = [find(nucleoli_number(:,1)>1,1) find(nucleoli_number(:,2)>1,1)] + Division{n}.Metaphase - 1;
    
end

%% plot
mkdir('plots')
close all
for n = 1 : ndivisions
    if max(Division{n}.Spindle(:,2)) < 4
        Division{n} = [];
        continue;
    end
    if isempty(Division{n}.Ana_onset)
        Division{n} = [];
        continue;
    end
    figure; hold on
    T = (1:1:length(Division{n}.nucleoli)) + Division{n}.Metaphase - 1;
    plot(T, Division{n}.Spindle(:,2),'o-') % Spindle
    plot(T, Division{n}.nucleoli(:,1),'o-') % nucleoli, majoraxis
    plot(T, Division{n}.nucleoli(:,2),'o-') % nucleoli, boundingbox
    if ~isempty(Division{n}.NucleoliInBud)
        plot([1 1]*Division{n}.NucleoliInBud, [0 10],'k') % Nucleoli in Bud
        plot([1 1]*Division{n}.Ana_onset, [0 10],'k:') % Anaphase onset
    end
    if ~isempty(Division{n}.NucleoliSeg)
        plot([1 1]*Division{n}.NucleoliSeg(1), [0 10],'g--') % Nucleli seg loose
        if length(Division{n}.NucleoliSeg) > 1
            plot([1 1]*Division{n}.NucleoliSeg(2), [0 10],'g-.') % Nuclei seg strict
        end
    end
    xlabel('Frame #')
    ylabel('Length (um)')
    legend({'Spindle' 'Nucleoli-MajorAxisLength', 'Nucleoli-BoundingBox', 'Nucleoli In Bud', 'Anaphase onset', 'Nucleoli seg' 'Nucleoli seg stict'}, 'Location', 'northwest')
    
    saveas(gcf, ['plots/Nucleoli_seggregation new ' num2str(n)], 'png')
    close
end

%%
close all
clear data
save('Cfi1')
cd('..')