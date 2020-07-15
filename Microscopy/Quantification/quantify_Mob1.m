function [] = quantify_Mob1(Mob1, SPB_nucleoli, NLS, nslice, dT, PixelSize)

%%
data = bfopen([Mob1 '.tif']);

load([SPB_nucleoli '/assignment'], 'trk_cells', 'trk_nucleoli','I_max2')
load([NLS '/NLSiRFP'], 'Division')

mkdir(Mob1)
cd(Mob1)

nframes = length(data{1})/nslice;
I_mCherry = I_max2;
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
ndivisions = length(Division);
% n = 4;
[height,width] = size(I_max{1});
mkdir('Segmentation')
f = figure; 
set(f,'visible','off');
for n = 1 : ndivisions
    if max(Division{n}.Spindle) < 35
        continue;
    end
    T_division = min(Division{n}.Cytokinesis+1,nframes) - Division{n}.Budding + 1;
    [~, t_exit] = max(Division{n}.Spindle);
    int_Mob1 = nan(T_division,4); % cyto, mSPB, dSPB, nucleoli
    spindle_length = nan(T_division,1); % remeasure spindle length
    mask_SPBs = cell(T_division,1);
    mask_cell_daughter = cell(T_division, 1);
    for i = 1 : T_division
        i_frame = Division{n}.Budding + i - 1;
        mask_cell = trk_cells{i_frame} == Division{n}.Mother | trk_cells{i_frame} == Division{n}.Daughter;
        mask_cell_daughter{i} = imerode(trk_cells{i_frame} == Division{n}.Daughter, strel('disk',2));
        mask_mCherry1 = mask_cell&I_mCherry{i_frame}>20;% this threshold is determined mannually
        mask_mCherry = mask_mCherry1&imbinarize(mat2gray(I_mCherry{i_frame}), 'adaptive');
        mask_mCherry = bwareaopen(mask_mCherry,2);
        mask_mCherry = bwmorph(mask_mCherry, 'clean');
        mask_Mob1b = mask_cell&I_max2{i_frame}>30;
        mask_Mob1 = mask_Mob1b&loccssegm(I_max2{i_frame}, 3, 0.7);

        % if there is Mob1 signal, use it for SPBs
        label_Mob1 = labelmatrix(bwconncomp(mask_Mob1));
        if length(unique(label_Mob1)) > 2
            % remoce bud neck
            if abs(i - t_exit) < 3 % Mob1 @bud neck only appear around the time of exit
                radius_dilate = 2;
                mask_mother = imdilate(trk_cells{i_frame} == Division{n}.Mother, strel('disk', radius_dilate));
                mask_daughter = imdilate(trk_cells{i_frame} == Division{n}.Daughter, strel('disk', radius_dilate));
                mask_budneck = mask_mother&mask_daughter&mask_cell;
                mask_Mob1(mask_budneck>0) = 0;
            end
            label_Mob1 = labelmatrix(bwconncomp(mask_Mob1));
%             S_Mob1 = regionprops(label_Mob1, 'Area', 'Centroid', 'Eccentricity');
            % remove nucleoli or random bright spot
            nSPBs = length(unique(label_Mob1)) - 1;
            if nSPBs > 2
                int_SPB = zeros(nSPBs,1);
                    for j = 1 : nSPBs
                        int_SPB(j) = max(double(I_max2{i_frame}(label_Mob1==j)));
                    end
                    % choose the brightest two
                    [~,idx] = sort(int_SPB, 'descend');
                    for j = 3 : nSPBs
                        mask_Mob1(label_Mob1 == idx(j)) = 0; 
                    end
            end
            mask_SPBs{i} = imdilate(mask_Mob1, strel('disk', 1));
            
        else
            % there is either no Mob1 signal or Mob1 @ one SPB
            mask_SPBs{i} = zeros(height, width);
            if i_frame > Division{n}.SPBinbud
%                 mask_mCherry = imopen(mask_mCherry,strel('square',2));
            end
            nucleoli_SPB = labelmatrix(bwconncomp(mask_mCherry));
            if length(unique(label_Mob1)) > 1
                mask_SPBs{i} = mask_Mob1;                                
            end
%             if length(unique(nucleoli_SPB)) < 3
%                 mask_mCherry = mask_mCherry1&imbinarize(mat2gray(I_mCherry{i_frame}), 0.3);
%                 nucleoli_SPB = labelmatrix(bwconncomp(mask_mCherry));
%             end
            S = regionprops(nucleoli_SPB, 'Area', 'Centroid', 'Eccentricity');
            idx_SPB = intersect(find([S.Area] < 11), find([S.Eccentricity] < 0.9));
            
            nSPBs = length(idx_SPB);
            if nSPBs > 0
                for j = 1 : nSPBs
                    mask_SPBs{i}(nucleoli_SPB==idx_SPB(j)) = 1;
                end

                if nSPBs > 2 - (length(unique(label_Mob1)) - 1)
                    int_SPB = zeros(nSPBs,1);
                    for j = 1 : nSPBs
                        int_SPB(j) = max(double(I_mCherry{i_frame}(nucleoli_SPB==idx_SPB(j))));
                    end
                    % choose the brightest two/one
                    [~,idx] = sort(int_SPB, 'descend');
                    for j = 3 - (length(unique(label_Mob1)) - 1) : nSPBs
                        mask_SPBs{i}(nucleoli_SPB == idx_SPB(idx(j))) = 0; 
                    end
                end                                
            end
        end
        S_spb = regionprops(mask_SPBs{i}>0, 'Centroid');
        nSPBs = length(S_spb);
        
        nucleoli = mask_mCherry;
        nucleoli(mask_SPBs{i}==1) = 0;

        % show the segmentation & quantify Mob1 intensity
        imshowpair(I_mCherry{i_frame}, I_max2{i_frame},'montage')
%         imshowpair(I_mCherry{i_frame}, I_max2{i_frame},'montage', 'Scaling','joint')
        truesize([2*height, 2*width])
        hold on
        b = bwboundaries(imclose(mask_cell, strel('square',3)));
        if ~isempty(b)
            plot(smooth(b{1}(:,2)),smooth(b{1}(:,1)),'g','linewidth',0.5);
            plot(smooth(b{1}(:,2))+width, smooth(b{1}(:,1)),'g','linewidth',0.5);
        end
        int_Mob1(i, 1) = median(I_max2{i_frame}(imerode(mask_cell, strel('disk',2))&~mask_mCherry));% cytosolic intensity
        if nSPBs > 0
            SPBs_d = zeros(height, width); % draw SPBs
            for j = 1 : nSPBs
                % make a circle around the centroid of detected SPB
                SPBs_d(round(S_spb(j).Centroid(2)), round(S_spb(j).Centroid(1))) = j;
            end
            SPBs_d = imdilate(SPBs_d, strel('disk', 2));
            for j = 1 : nSPBs
                c = bwboundaries(SPBs_d == j);
                plot(smooth(c{1}(:,2)),smooth(c{1}(:,1)),'y','linewidth',0.5);
                plot(smooth(c{1}(:,2))+width, smooth(c{1}(:,1)),'y','linewidth',0.5);
                % find out which SPB (based on which cell it is in)
            end
            SPBs_id = (SPBs_d > 0) .* trk_cells{i_frame};
            if ismember(Division{n}.Mother, SPBs_id)
                int_Mob1(i, 2) = max(I_max2{i_frame}(SPBs_id == Division{n}.Mother));
            end
            if ismember(Division{n}.Daughter, SPBs_id)
                int_Mob1(i, 3) = max(I_max2{i_frame}(SPBs_id == Division{n}.Daughter));
            end
        end
        d = bwboundaries(nucleoli);
        for j = 1 : length(d)
            plot(smooth(d{j}(:,2)),smooth(d{j}(:,1)),'m','linewidth',0.5);
            plot(smooth(d{j}(:,2))+width, smooth(d{j}(:,1)),'m','linewidth',0.5);
        end
        int_Mob1(i, 4) = median(I_max2{i_frame}(nucleoli));
               
        % remeasure the spindle length
        trk_division = double(mask_mCherry|mask_Mob1);
        spindle = regionprops(trk_division, 'MajorAxisLength');
        if ~isempty(spindle)
            spindle_length(i) = spindle.MajorAxisLength;
        end
        
        hold off
        saveas(f, ['Segmentation/Division ' num2str(n) '_frame ' num2str(i_frame)], 'png')
    end
    Division{n}.Mob1 = int_Mob1;
    Division{n}.Spindle2 = spindle_length;
    trk_SPBs = cellfun(@(x,y) x.*y, mask_cell_daughter, mask_SPBs, 'UniformOutput',false);
    Division{n}.SPBinbud = find((cellfun(@nnz, trk_SPBs))>0,1) + Division{n}.Budding - 1;
end

%% plot
mkdir('plots')
close all
for n = 1 : ndivisions
    if max(Division{n}.Spindle) < 35
        Division{n} = [];
        continue;
    end
    if isempty(Division{n}.SPBinbud)
        Division{n} = [];
        continue;
    end
    figure; hold on
    T = dT*((1:1:length(Division{n}.Mob1)) - Division{n}.SPBinbud + Division{n}.Budding - 1);
    plot(T, Division{n}.Mob1(:,2),'o-')
    plot(T, Division{n}.Mob1(:,3),'o-')
    plot(T, Division{n}.Mob1(:,4),'*-')
    plot(T, Division{n}.Mob1(:,1),'k-')
    xlabel('Time since SPB enters bud (min)')
    ylabel('Mob1 Intensity (au)')
    legend({'Mob1 at mSPB' 'Mob1 at dSPB' 'Mob1 in nucleolus' 'Mob1 in cytosol'}, 'Location', 'northwest')
    
    yyaxis right
    plot(T, -Division{n}.NLS(:,1))
    axis([T(1) T(end) -1.05 -0.35])
    ylabel('NLS-TR release (-CV)')
    
    saveas(gcf, ['plots/Mob1_Division ' num2str(n)], 'png')
    close
    
    figure; hold on
    plot(T, Division{n}.Spindle2*PixelSize)
%     plot(T, Division{n}.Spindle(1:length(Division{n}.Spindle2)))
    xlabel('Time since SPB enters bud (min)')
    ylabel('Estimated spindle length (um)')
    saveas(gcf, ['plots/Spindle_Division ' num2str(n)], 'png')
    close
%     legend({'Mob1 at mSPB' 'Mob1 at dSPB' 'Mob1 in nucleolus'})
end

% Division = Division(~cellfun(@isempty,Division));

%%
close all
clear data
save('Mob1')
cd('..')