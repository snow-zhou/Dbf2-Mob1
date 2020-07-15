function [] = seg_assign_SPB(SPB, dT, PixelSize, dZ)
% use cell centric view to identify mitosis

%%
load([SPB '/SPB'], 'trk_cells', 'bw_SPB', 'trk_SPB', 'SPBsInfo', 'I_max2')
cd(SPB)
T_ana = 90/dT; % set a reasonable time frame for anaphase: 90mins since SPB enter buds
nframes = length(I_max2);

%% assign SPB to cells
trk_SPB_cell = cellfun(@(x,y) x.*y, trk_cells, bw_SPB, 'UniformOutput',false);

% frame-centric id
cellListID = cellfun(@nonzeros, cellfun(@unique, trk_cells, 'UniformOutput',false), 'UniformOutput',false);
SPBListID = cellfun(@nonzeros, cellfun(@unique, trk_SPB_cell, 'UniformOutput',false), 'UniformOutput',false);

ncells = length(cellListID{end});

%% identify daughters to analyze
mothers = SPBListID{1}; % cells present with SPB(s) in the first frames
daughters = setdiff((1:1:ncells),mothers); % cells didn't have SPB or exit in the first frame

cells_count = SPBListID{end}; % all the cells have SPB in the last frame
daughters_count = intersect(daughters,cells_count);

%% find mothers & anaphase duration
ndivisions = length(daughters_count);
Division = cell(ndivisions,1);
for i = 1 : ndivisions
    first_frame = find(~cellfun(@isempty, cellfun(@(x) find(x==daughters_count(i)), cellListID, 'UniformOutput', 0)),1);
    %find the first frame SPB enters bud (continuesly)
    nucframes = find(~cellfun(@isempty, cellfun(@(x) find(x==daughters_count(i)), SPBListID, 'UniformOutput', 0)));
    first_nucframe = nucframes(find(diff(nucframes,2)==0,1));
    if isempty(first_nucframe)
        continue
    end
    if (nframes - first_nucframe) < 2
        continue
    end
    % remove false daughters due to mis-segmentation
    if first_nucframe == first_frame
        continue
    end
    
    % find mother by matching the dSPB that travels from mother to bud
    dSPB_ID = nonzeros(unique((trk_SPB_cell{first_nucframe}==daughters_count(i)).*trk_SPB{first_nucframe}));
    if isempty(dSPB_ID)
        continue;
    end
    if length(dSPB_ID) > 1
        continue;
    end
    mother_ID = nonzeros(unique((trk_SPB{first_nucframe-1}==dSPB_ID).*trk_SPB_cell{first_nucframe-1}));
    mother_ID = setdiff(mother_ID, daughters_count(i));
    if isempty(mother_ID)
        continue;
    end
    if length(mother_ID) > 1
        continue;
    end
    mSPB_ID = nonzeros(unique((trk_SPB_cell{first_nucframe}==mother_ID).*trk_SPB{first_nucframe}));
    mSPB_ID = setdiff(mSPB_ID, dSPB_ID);
    if isempty(mSPB_ID)
        continue;
    end
    if length(mSPB_ID) > 1
        continue;
    end
    
    % redefine first_nucframe with the dSPB info
    nucframes = find(SPBsInfo(dSPB_ID).CellID == daughters_count(i)) + SPBsInfo(dSPB_ID).Birthframe - 1;
    first_nucframe = nucframes(find(diff(nucframes)==1,1));
    
    % find mitotic exit based on spindle length (distance of SPBs)
    metaphase = max([SPBsInfo(dSPB_ID).Birthframe SPBsInfo(mSPB_ID).Birthframe]);
    if metaphase < find(~cellfun(@isempty, cellfun(@(x) find(x==mother_ID), cellListID, 'UniformOutput', 0)),1)
        continue;
    end
    spindle_length = nan(min([T_ana, nframes - first_nucframe]) + (first_nucframe - metaphase), 2); % 3D/2D length
    for t = metaphase : min([first_nucframe + T_ana, nframes])
        i_daughter = t - SPBsInfo(dSPB_ID).Birthframe + 1;
        i_mother = t - SPBsInfo(mSPB_ID).Birthframe + 1;
        spindle_length(t - metaphase + 1,1) = norm([PixelSize PixelSize dZ].*(SPBsInfo(dSPB_ID).Centroid(i_daughter,:) - SPBsInfo(mSPB_ID).Centroid(i_mother,:)));
        spindle_length(t - metaphase + 1,2) = PixelSize*norm(SPBsInfo(dSPB_ID).Centroid(i_daughter,1:2) - SPBsInfo(mSPB_ID).Centroid(i_mother,1:2));
    end    
    
    % find the frame of cytokinesis (mitotic exit + 6 frames)
    if max(spindle_length(:,1)) < 5 % if spindle elongation occurred
        continue;
    end
    [~, locs, ~, p]=findpeaks(spindle_length(:,1), 'MinPeakProminence',1);
    if ~isempty(p)
        [~, idx] = max(p);
        idx = locs(idx);
    else 
        [~, idx] = max(spindle_length(:,1));
    end
    cytokinesis = idx + metaphase + 6; %
    
    % find anaphase onset (spindle > 3um in the elongation phase)
    ana_onset = find(spindle_length(1:idx,2) < 3, 1, 'last') + metaphase;    
    if isempty(ana_onset)
        continue;
    end
        
    Division{i} = struct(...
                        'Daughter', daughters_count(i),...
                        'Mother', mother_ID,... 
                        'mSPB', mSPB_ID,...
                        'dSPB', dSPB_ID,...
                        'Budding', first_frame ,...
                        'Metaphase', metaphase,...
                        'SPBinbud', first_nucframe,...
                        'Ana_onset', ana_onset,...
                        'Cytokinesis', cytokinesis,...
                        'Spindle', spindle_length);
end

% remove empty Division cell
Division = Division(~cellfun(@isempty,Division));
ndivisions = length(Division);

%%
mkdir('Division')
f = figure;
%  set(f,'visible','off');
for i = 1 : nframes
    imshow(I_max2{i},[])
    hold on
    for n = 1 : ndivisions
        if i >= Division{n}.Budding && i <= Division{n}.Cytokinesis+1
            b = bwboundaries(imclose(trk_cells{i} == Division{n}.Mother | trk_cells{i} == Division{n}.Daughter, strel('square',3)));
            if ~isempty(b)
                plot(smooth(b{1}(:,2)),smooth(b{1}(:,1)),'g','linewidth',1);
                d = [min(b{1}(:,2))-5, max(b{1}(:,1)+5)];
                text(d(1),d(2),num2str(n),'Color','green','FontSize',12);
            end
        end
    end
    hold off
    saveas(f, ['Division/overlay_frame ' num2str(i)], 'png')
end

%%
mkdir('plots')
g = figure;
set(g,'visible','off');
for i = 1 : ndivisions
    plot(((1:1:length(Division{i}.Spindle)))+Division{i}.Metaphase-1, Division{i}.Spindle, 'o-')
    hold on
    plot([1 1]*Division{i}.Ana_onset, [min(Division{i}.Spindle(:,2)) max(Division{i}.Spindle(:,2))],'--')
    plot([1 1]*Division{i}.Cytokinesis, [min(Division{i}.Spindle(:,2)) max(Division{i}.Spindle(:,2))])
    axis([0 nframes 0 10])
    xlabel('Time (frame #)')
    ylabel('Spindle length (um)')
    legend({'Spindle length 3D' 'Spindle length 2D' 'Anaphase onset (> 3um)' 'Cytokinesis (mitotic exit + 5 frames)'}, 'location', 'northwest')
    hold off
    saveas(g, ['plots/spindle_division ' num2str(i)], 'png')
end

%%
close all
save('assignment')
cd('..')