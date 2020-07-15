function [] = seg_assign_nucleoli(BF, SPB_nucleoli, dT)


%%
load([BF '/BF'], 'trk_cells_fix','nframes','ncells');
load([SPB_nucleoli '/nucleoli'], 'bw_nucleoli_SPB','I_max2')
cd(SPB_nucleoli)
radius = 16; % #pixels (~3.5 um) as threshold for mother-daughter nucleoli distance upon anaphase onset
T_ana = 75/dT;% set a reasonable time frame for anaphase: 75mins
%% assign nucleoli
trk_cells = trk_cells_fix;
trk_nucleoli = cellfun(@(x,y) x.*y, trk_cells, bw_nucleoli_SPB, 'UniformOutput',false);

% frame-centric id
cellListID = cellfun(@nonzeros, cellfun(@unique, trk_cells, 'UniformOutput',false), 'UniformOutput',false);
nucleoliListID = cellfun(@nonzeros, cellfun(@unique, trk_nucleoli, 'UniformOutput',false), 'UniformOutput',false);

cellList = cell(nframes,1);
for i = 1 : nframes
    cellList{i} = cell(length(cellListID{i}),1);
    S = regionprops(trk_cells{i},'Centroid','Area');
    Sn = regionprops(trk_nucleoli{i}, 'Centroid','Area');
    for n = 1 : length(cellListID{i})
        if ~ismember(cellListID{i}(n),nucleoliListID{i})
            cellList{i}{n} = struct('CellID', cellListID{i}(n), ...
                                    'Centroid', S(cellListID{i}(n)).Centroid,...
                                    'Area', S(cellListID{i}(n)).Area,...
                                    'Nucleoli_SPB', struct('Area', 0, 'Centroid', [NaN, NaN]),...
                                    'SPBs', []);
        else
            Sspb = regionprops(imerode(trk_nucleoli{i} == cellListID{i}(n), strel('disk',0)), 'Centroid','Area');
            cellList{i}{n} = struct('CellID', cellListID{i}(n), ...
                                    'Centroid', S(cellListID{i}(n)).Centroid,...
                                    'Area', S(cellListID{i}(n)).Area,...
                                    'Nucleoli_SPB', Sn(cellListID{i}(n)),...
                                    'SPBs', Sspb);
        end
    end
end

%% identify daughters to analyze
mothers = nucleoliListID{1}; % cells present with nucleus in the first frames
daughters = setdiff((1:1:ncells),mothers); % cells didn't have nucleus or exit in the first frame

cells_count = nucleoliListID{end}; % all the cells have nucleus in the last frame
daughters_count = intersect(daughters,cells_count);

%% find mothers & anaphase duration
mkdir('plots')
ndivisions = length(daughters_count);
Division = cell(ndivisions,1);
for i = 1 : ndivisions
    first_frame = find(~cellfun(@isempty, cellfun(@(x) find(x==daughters_count(i)), cellListID, 'UniformOutput', 0)),1);
    %find the first frame SPB enters bud (continuesly)
    nucframes = find(~cellfun(@isempty, cellfun(@(x) find(x==daughters_count(i)), nucleoliListID, 'UniformOutput', 0)));
    first_nucframe = nucframes(find(diff(nucframes)==1,1));
    if isempty(first_nucframe)
        continue
    end
    % remove false daughters due to mis-segmentation
    if first_nucframe == first_frame
        continue
    end
    candidates = nucleoliListID{first_nucframe-1}; %all cells with nucleoli at the frame before are candidate mothers
    % eliminate unfit mothers
    qualification = ones(length(candidates),1);
    qualification2 = ones(length(candidates),1);
    for n = 1 : length(candidates)
        qualification(n) = cellList{first_nucframe - 1}{cellListID{first_nucframe - 1}==candidates(n)}.Nucleoli_SPB.Area;
        if ismember(candidates(n), daughters_count)
            first_nucframe_candidate = find(~cellfun(@isempty, cellfun(@(x) find(x==candidates(n)), nucleoliListID, 'UniformOutput', 0)),1);
            if (first_nucframe - first_nucframe_candidate) < 20
                qualification2(n) = 0;
            end
        end
    end
    candidates = candidates(qualification>5 & qualification2>0);
    % find cells whose nucleoli+SPB are within proximity
    I_center = (trk_nucleoli{first_nucframe} == daughters_count(i));
    I_range = imdilate(I_center, strel('disk', radius));
    hits = intersect(candidates,unique(I_range.*trk_nucleoli{first_nucframe - 1}));
    
    if isempty(hits)
        continue; % maybe try bigger radius?
    end
    
    if length(hits)>1
        %find the cells that recently have a SPB moved to the bud
        %find the cloest SPB to the dSPB in bud in the previous frame
        SPBs = cell(length(hits),1);
        dist_SPB = inf(length(hits),1);
        for n = 1 : length(hits)     
            SPBs{n} = cellList{first_nucframe - 1}{cellListID{first_nucframe - 1}==hits(n)}.SPBs;
            target = cellList{first_nucframe}{cellListID{first_nucframe}==daughters_count(i)}.Nucleoli_SPB;
            dist_SPB(n) = min(vecnorm(reshape([SPBs{n}.Centroid]', 2, length(SPBs{n})) - target.Centroid'));
%             nSPBs(n)
        end
        [~,idx] = min(dist_SPB);
        match = hits(idx);
    else
        match = hits;
    end
    
    % find mitotic exit based on spindle length (distance of SPBs)
    Budding = find(~cellfun(@isempty, cellfun(@(x) find(x==daughters_count(i)), cellListID, 'UniformOutput', 0)),1);
    spindle_length = nan(min([T_ana, nframes - first_nucframe]) + (first_nucframe - Budding) + 1, 1);
    for t = Budding : min([first_nucframe + T_ana, nframes])
        trk_division = double((trk_nucleoli{t}==daughters_count(i)|trk_nucleoli{t}==match));
        spindle = regionprops(trk_division, 'MajorAxisLength');
        if ~isempty(spindle)
            spindle_length(t - Budding + 1) = spindle.MajorAxisLength;
        end
    end
    
    % find the frame of cytokinesis (mitotic exit + 20 mins)
    cytokinesis = nframes;
    if max(spindle_length) > 35 % if spindle elongation occurred
        [~, idx] = max(spindle_length);
        cytokinesis = idx + Budding + 15/dT; %
    end
    
    g = figure;
    set(g,'visible','off'); hold on
    plot((Budding:1:min([first_nucframe + T_ana, nframes])), spindle_length)
    plot([1 1]*cytokinesis, [min(spindle_length) max(spindle_length)])
    axis([0 nframes 0 80])
    xlabel('Time (frame #)')
    ylabel('Spindle length (pixel)')
    hold off
    saveas(g, ['plots/spindle_division ' num2str(i)], 'png')
        
    Division{i} = struct(...
                        'Daughter', daughters_count(i),...
                        'Mother', match,...
                        'SPBinbud', first_nucframe,...
                        'Budding', first_frame ,...
                        'Cytokinesis', cytokinesis,...
                        'Spindle', spindle_length);
end

% remove empty Division cell
Division = Division(~cellfun(@isempty,Division));
ndivisions = length(Division);

%% Extract info
% cell-centric id
cells = cell(ncells,1);

for i = 1 : nframes
    for n = 1 : length(cellListID{i})
        cells{cellListID{i}(n)}{i} = cellList{i}{n};
    end
end

%remove empty cells from later appeared daughters and convert type
for n = 1 : ncells
    if isempty(cells{n})
        continue
    end
    cells{n} = cells{n}(~cellfun(@isempty,cells{n}));
    cells{n} = cell2mat(cells{n});
end

%%
cellsInfo = cell(ncells,1);

for n = 1 : ncells
    if isempty(cells{n})
        continue
    end
    cellsInfo{n} = struct(...
                          'CellID', cells{n}(1).CellID, ...
                          'Birthframe', nframes - length(cells{n}) + 1, ...
                          'Centroid', (reshape([cells{n}.Centroid],2,[]))', ...
                          'Area', [cells{n}.Area]);
end

% cellsInfo = cell2mat(cellsInfo);

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
close all
save('assignment')
cd('..')