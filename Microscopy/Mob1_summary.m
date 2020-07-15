%% summarize results for comparison
clear
close all
%input sample info
addpath('/plotSpread');
dirname = '';
cd(dirname)

folder = {};
label = {};

%% combine all the data from different fields for each sample
n_sample = length(folder);
samples = cell(n_sample,1);
for n = 1 : n_sample 
    list_file = dir([folder{n} '/C2*.tif']);
    for i = 1 : length(list_file)
        load([folder{n} '/' list_file(i).name(1:end-4) '/Mob1.mat'], 'Division', 'dT', 'PixelSize', 'nframes')
        for j = 1 : length(Division) 
            if isempty(Division{j})
                continue
            end
            if max(Division{j}.Spindle) < 35
                Division{j} = [];
                continue
            end
            [~, mitotic_exit] = max(Division{j}.Spindle2);
            if (mitotic_exit + Division{j}.Budding > nframes)
                Division{j} = [];
                continue
            end
            Division{j}.filename = list_file(i).name(1:end-4);
            Division{j}.ID = j;
            Division{j}.ana_onset = find(smooth(Division{j}.Spindle2)*PixelSize > 4, 1, 'first') + Division{j}.Budding - 1; 
            if isempty(Division{j}.ana_onset)
                Division{j} = [];
                continue
            end
            Division{j}.mitotic_exit = mitotic_exit + Division{j}.Budding;
            Division{j}.max_enrichment = max(Division{j}.Mob1(:,2:4)./Division{j}.Mob1(:,1));
            Division{j}.min_NLS_CV = min(Division{j}.NLS(:,1));
        end
        Division = Division(~cellfun(@isempty,Division));
        samples{n} = [samples{n}; cell2mat(Division)]; 

    end
end

mkdir('summary')
cd('summary')

%% plot dynamics of Mob1 localization
Mob1 = cell(n_sample,1);
NLS = cell(n_sample,1);
Spindle = cell(n_sample,1);
np_left = zeros(n_sample,1);
close all
f1 = figure('Position', [10 10 400 1200]); hold on
f2 = figure; hold on
f3 = figure; hold on
f4 = figure; hold on
Color_order = get(gca,'colororder');

for i = 1 : n_sample
    np_left(i) = max([samples{i}.ana_onset] - [samples{i}.Budding]) + 1;
    np_right = max([samples{i}.Cytokinesis] - [samples{i}.ana_onset]) + 1;
    Mob1{i} = nan(length(samples{i}), np_left(i)+np_right + 1, 3);
    NLS{i} = nan(length(samples{i}), np_left(i)+np_right + 1, 2);
    Spindle{i} = nan(length(samples{i}), np_left(i)+np_right + 1);
    for n = 1 : length(samples{i})
        shift = np_left(i) - (samples{i}(n).ana_onset - samples{i}(n).Budding);
        Mob1{i}(n, shift+1:shift+length(samples{i}(n).Mob1), :) = samples{i}(n).Mob1(:,2:4)./samples{i}(n).Mob1(:,1);
        NLS{i}(n, shift+1:shift+length(samples{i}(n).Mob1), 1) = samples{i}(n).NLS(:,1);
        NLS{i}(n, shift+1:shift+length(samples{i}(n).Mob1), 2) = samples{i}(n).NLS(:,3)./samples{i}(n).NLS(:,2);
        Spindle{i}(n, shift+1:shift+length(samples{i}(n).Mob1)) = samples{i}(n).Spindle2;
    end
    
    figure(f1)
    subplot(n_sample, 1, i)
    hold on
    plot(dT*(-np_left(i):1:np_right), Mob1{i}(:,:,1)', ':', 'LineWidth', 0.1, 'Color', Color_order(1,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), Mob1{i}(:,:,2)', ':', 'LineWidth', 0.1, 'Color', Color_order(2,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Mob1{i}(:,:,1)), 'LineWidth', 2, 'Color', Color_order(1,:))
    plot(dT*(-np_left(i):1:np_right), nanmean(Mob1{i}(:,:,2)), 'LineWidth', 2, 'Color', Color_order(2,:))
    axis([-40 40 0 20])
    ylabel('Enrichment of Mob1 at SPBs')
    title(label{i})
    legend({'mSPB' 'dSPB'},'location', 'northwest')
    
    figure(f2)
    plot(dT*(-np_left(i):1:np_right), Mob1{i}(:,:,3)', ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Mob1{i}(:,:,3)), 'LineWidth', 3, 'Color', Color_order(i,:))
    
    figure(f3)
    plot(dT*(-np_left(i):1:np_right), NLS{i}(:,:,1)', ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(NLS{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
    
    figure(f4)
    plot(dT*(-np_left(i):1:np_right), Spindle{i}'*PixelSize, ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Spindle{i})*PixelSize, 'LineWidth', 3, 'Color', Color_order(i,:))
end

figure(f1)
xlabel('Time since anaphase onset (min)')
saveas(f1, 'Aligned Mob1 at SPBs','fig')
saveas(f1, 'Aligned Mob1 at SPBs','png')

figure(f2)
axis([-40 40 1 3.5])
xlabel('Time since anaphase onset (min)')
ylabel('Enrichment of Mob1 in nucleolus')
legend(label,'location','northwest')
saveas(f2, 'Aligned Mob1 in nucleolus','fig')
saveas(f2, 'Aligned Mob1 in nucleolus','png')

figure(f3)
axis([-40 40 0.4 1.2])
xlabel('Time since anaphase onset (min)')
ylabel('NLS_C_d_c_1_4 release (CV)')
legend(label,'location','northeast')
saveas(f3, 'Aligned NLS release','fig')
saveas(f3, 'Aligned NLS release','png')

figure(f4)
axis([-40 40 0 15])
xlabel('Time since anaphase onset (min)')
ylabel('Estimated spinle length (um)')
legend(label,'location','northwest')
saveas(f4, 'Aligned Spinle length','fig')
saveas(f4, 'Aligned Spinle length','png')

%% plot CI
close all
figure('Position', [10 10 600 500]); hold on

for i = 1 : n_sample
    np_right = max([samples{i}.Cytokinesis] - [samples{i}.ana_onset]) + 1;
    x = dT*(-np_left(i):1:np_right);
    
    subplot(2,3,1); hold on
    upper = nanmean(Mob1{i}(:,:,1))+1.96*nanstd(Mob1{i}(:,:,1))./sqrt(sum(~isnan(Mob1{i}(:,:,1))));
    lower = nanmean(Mob1{i}(:,:,1))-1.96*nanstd(Mob1{i}(:,:,1))./sqrt(sum(~isnan(Mob1{i}(:,:,1))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Mob1{i}(:,:,1)), 'LineWidth', 2, 'Color', Color_order(i,:))
    axis([-40 60 0 20])
    xlabel('Time since anaphase onset (min)')
    ylabel('Enrichment of Mob1 at mSPB')
    
    subplot(2,3,4); hold on
    upper = nanmean(Mob1{i}(:,:,2))+1.96*nanstd(Mob1{i}(:,:,2))./sqrt(sum(~isnan(Mob1{i}(:,:,2))));
    lower = nanmean(Mob1{i}(:,:,2))-1.96*nanstd(Mob1{i}(:,:,2))./sqrt(sum(~isnan(Mob1{i}(:,:,2))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Mob1{i}(:,:,2)), 'LineWidth', 2, 'Color', Color_order(i,:))
    axis([-40 60 0 20])
    xlabel('Time since anaphase onset (min)')
    ylabel('Enrichment of Mob1 at dSPB')
    
    subplot(2,3,2); hold on
    upper = nanmean(NLS{i}(:,:,1))+1.96*nanstd(NLS{i}(:,:,1))./sqrt(sum(~isnan(NLS{i}(:,:,1))));
    lower = nanmean(NLS{i}(:,:,1))-1.96*nanstd(NLS{i}(:,:,1))./sqrt(sum(~isnan(NLS{i}(:,:,1))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(NLS{i}(:,:,1)), 'LineWidth', 2, 'Color', Color_order(i,:))
    axis([-40 60 0.4 1.2])
    xlabel('Time since anaphase onset (min)')
    ylabel('NLS_C_d_c_1_4 release (CV)')
    
    subplot(2,3,3); hold on
    upper = nanmean(Mob1{i}(:,:,3))+1.96*nanstd(Mob1{i}(:,:,3))./sqrt(sum(~isnan(Mob1{i}(:,:,3))));
    lower = nanmean(Mob1{i}(:,:,3))-1.96*nanstd(Mob1{i}(:,:,3))./sqrt(sum(~isnan(Mob1{i}(:,:,3))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Mob1{i}(:,:,3)), 'LineWidth', 2, 'Color', Color_order(i,:))
    axis([-40 60 1 3])
    xlabel('Time since anaphase onset (min)')
    ylabel('Enrichment of Mob1 in nucleolus')
end

legend(label)
saveas(gcf, 'Aligned Mob1_CI','fig')
saveas(gcf, 'Aligned Mob1_CI','pdf')

%% plot Mob1 localization/enrichment
close all
figure('Position', [10 10 600 500]);
Color_order = get(gca,'colororder');
data = cell(n_sample,5);
for i = 1 : n_sample
    data{i,1} = max(Mob1{i}(:,np_left(i)+1:np_left(i)+11,3),[],2); % nucleolus
    data{i,2} = max(Mob1{i}(:,np_left(i)+1:np_left(i)+11,1),[],2); % mSPB
    data{i,5} = max(Mob1{i}(:,np_left(i)+1:np_left(i)+11,2),[],2); % dSPB
    data{i,3} = Color_order(i, :);
    data{i,4} = [samples{i,1}.min_NLS_CV];
end

subplot(2, 2, 1); hold on
plotSpread(data(:,1), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,1}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 1 3.5])
ylabel('Maximum enrichment in nucleolus')

subplot(2, 2, 2); hold on
plotSpread(data(:,2), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,2}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 0 32])
ylabel('Maximum enrichment at mSPB')

subplot(2, 2, 3); hold on
plotSpread(data(:,4), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,4}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 0.2 1.2])
ylabel('Minimal CV of NLS_C_d_c_1_4')

subplot(2, 2, 4); hold on
plotSpread(data(:,5), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,5}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 0 32])
ylabel('Maximum enrichment at dSPB')

saveas(gcf, 'Mob1 localization', 'fig')
saveas(gcf, 'Mob1 localization', 'pdf')

%% plot Mob1 enrichment (spread)
close all
figure('Position', [10 10 300 300])
plotSpread(data(:,1), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
hold on
for i = 1 : n_sample
    plot([-0.3 0.3]+i, [1 1]*nanmedian(data{i,1}), 'Color', Color_order(i, :),'linewidth', 3)
end

ax = gca;
axis([0.25 n_sample+0.75 1 3.5])
set(ax,'YTick',(1:1:5.5))
ylabel('Maximum enrichment in nucleolus')
saveas(gcf, 'Maximum enrichment in nucleolus3', 'fig')
saveas(gcf, 'Maximum enrichment in nucleolus3', 'pdf')

%% plot anaphase duration
close all
figure('Position', [10 10 400 1200]);
Color_order = get(gca,'colororder');
for i = 1 : n_sample
    subplot(n_sample, 1, i)
    h = histogram(dT*([samples{i}.mitotic_exit] - [samples{i}.ana_onset]), (2*dT+dT/2:dT:dT*22), 'FaceColor', Color_order(i,:));
    hold on
    plot([1 1]*median([samples{i}.mitotic_exit] - [samples{i}.ana_onset])*dT,[0 max(h.Values)],'k:', 'LineWidth', 3)
    legend(label{i})
    ylabel('Cell count')
    ax = gca;
    set(ax,'XTick',(2*dT:dT:dT*22))
end
xlabel('Anaphase duration (min)')
saveas(gcf, 'Anaphase duration distribution', 'fig')
saveas(gcf, 'Anaphase duration distribution', 'png')

%% plot cummulative anaphase duration
close all
figure; hold on
Color_order = get(gca,'colororder');
for i = 1 : n_sample
    bin = (2*dT+dT/2:dT:dT*22);
    [N,EDGES] = histcounts(dT*([samples{i}.mitotic_exit] - [samples{i}.ana_onset]), bin);
    ndist = N / sum(N);
    cdist = cumsum(ndist);
    plot(0.5*(bin(1:end-1)+bin(2:end)), cdist, 'Color', Color_order(i,:), 'LineWidth', 2);   
end
axis([0 60 0 1])
ax = gca;
set(ax,'XTick',(2*dT:dT:dT*22))
xlabel('Anaphase duration (min)')
ylabel('Fraction completed anaphase')
saveas(gcf, 'Anaphase duration distribution2', 'fig')
saveas(gcf, 'Anaphase duration distribution2', 'png')


%%
close all
save('summary')

