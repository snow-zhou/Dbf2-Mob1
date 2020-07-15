function [] = quantify_NLSiRFP2(NLS, Nuc, nslice, dT)

%%
data = bfopen([NLS '.tif']);

load([Nuc '/assignment'], 'Division','trk_cells', 'trk_nucleoli','nframes')

mkdir(NLS)
cd(NLS)

% nframes = length(data{1})/nslice;

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

for n = 1 : ndivisions
    T_division = min(Division{n}.Cytokinesis+1, nframes) - Division{n}.Budding + 1;
    NLS_nuc = nan(T_division,1);
%     NLS_nuc2 = zeros(T_division,1);
    NLS_cyto = nan(T_division,1);
    NLS_nucleoli = nan(T_division,1);
    NLS_CV = nan(T_division,1);
    NLS_CV_mother = nan(T_division,1);
    for i = 1 : T_division
        i_frame = Division{n}.Budding + i - 1;
        mask_cell = trk_cells{i_frame} == Division{n}.Mother | trk_cells{i_frame} == Division{n}.Daughter;
        mask_cell = imerode(mask_cell, strel('disk',2));
        mask_nucleoli = trk_nucleoli{i_frame} == Division{n}.Mother | trk_nucleoli{i_frame} == Division{n}.Daughter;
        mask_mother = trk_cells{i_frame} == Division{n}.Mother;
        mask_mother = imerode(mask_mother, strel('disk',2));
        % intensity for NLS
        I_NLS = double(I_max2{i_frame}(mask_cell));
        I_NLS = sort(I_NLS,'descend');
        I_NLS_nucleolus = double(I_max2{i_frame}(mask_nucleoli));
        I_NLS_mother = double(I_max2{i_frame}(mask_mother));
        if isempty(I_NLS)
            continue
        end
        [f,xi] = ksdensity(I_NLS);
        [~, idx] = findpeaks(f);
        NLS_nuc(i) = xi(idx(end));
%         NLS_nuc2(i) = mean(I_NLS(1:round(end/50)));
        NLS_cyto(i) = xi(idx(1));
        NLS_nucleoli(i) = median(I_NLS_nucleolus);
        NLS_CV(i) = std(I_NLS)/mean(I_NLS);
        NLS_CV_mother(i) = std(I_NLS_mother)/mean(I_NLS_mother);
        
    end
    Division{n}.NLS = [NLS_CV NLS_nuc NLS_cyto NLS_nucleoli NLS_CV_mother];  
    
end

%% plot
mkdir('plots')
close all
for n = 1 : ndivisions
    figure; hold on
    yyaxis left
    [T_division, ~] = size(Division{n}.NLS);
    T = dT*(1:1:T_division);
    plot(T, smooth(Division{n}.NLS(:,3)./Division{n}.NLS(:,2)))
    plot(T, smooth(Division{n}.NLS(:,3)./Division{n}.NLS(:,4)))
    axis([0 dT*(T_division) 0 1.1])
    ylabel('NLS-iRFP release (I_c_y_t_o/I_n_u_c)')
%     legend({'Estimated nuclear' 'Nucleoli'}, 'Location', 'northwest')
    xlabel('Time since detected budding (min)')
    yyaxis right
    plot(T, -smooth(Division{n}.NLS(:,1)))
    plot(T, -smooth(Division{n}.NLS(:,5)))
    axis([0 dT*(T_division+1) -1 -0.4])
    ylabel('NLS-TR release (-CV)')
    legend({'Estimated nuclear (peak)' 'Nucleoli' 'mother+daughter' 'mother only'}, 'Location', 'northwest')
    
    saveas(gcf, ['plots/NLSiRFP_division' num2str(n)], 'png')
end

%%
if ~isempty(Division)
    mkdir('Division')
    close all
    f = figure;
    set(f,'visible','off');
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
end

%%
close all
clear data
save('NLSiRFP')
cd('..')