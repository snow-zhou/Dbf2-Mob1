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
        load([folder{n} '/' list_file(i).name(1:end-4) '/Cdc14.mat'], 'Division', 'dT', 'PixelSize', 'nframes')
        for j = 1 : length(Division) 
            if isempty(Division{j})
                continue
            end

            [~, mitotic_exit] = min(Division{j}.Cdc14(Division{j}.SPBinbud-Division{j}.Budding:end,3));
            mitotic_exit = mitotic_exit + Division{j}.SPBinbud-Division{j}.Budding -1;
            if mitotic_exit + Division{j}.Budding > nframes
                Division{j} = [];
                continue
            end
            Division{j}.filename = list_file(i).name(1:end-4);
            Division{j}.ID = j;
            Division{j}.ana_onset = find(smooth(Division{j}.Spindle)*PixelSize > 6, 1, 'first') + Division{j}.Budding - 3; 
            Division{j}.mitotic_exit = mitotic_exit + Division{j}.Budding;
            Division{j}.min_Cdc14_ratio = min(Division{j}.Cdc14(Division{j}.SPBinbud-Division{j}.Budding:end,3));
            Division{j}.min_Cdc14_CV = min(Division{j}.Cdc14(Division{j}.SPBinbud-Division{j}.Budding:end,4));
            Division{j}.min_NLS_CV = min(Division{j}.NLS(Division{j}.SPBinbud-Division{j}.Budding:end,1));
        end
        Division = Division(~cellfun(@isempty,Division));
        samples{n,1} = [samples{n,1}; cell2mat(Division)];         
    end
end

mkdir('summary')
cd('summary')

%% compare Cdc14 release (spread)
close all
figure
Color_order = get(gca,'colororder');
data = cell(n_sample,4);
for i = 1 : n_sample
    data{i,1} = [samples{i,1}.min_Cdc14_ratio];
    data{i,2} = [samples{i,1}.min_Cdc14_CV];
    data{i,3} = Color_order(i, :);
    data{i,4} = [samples{i,1}.min_NLS_CV];
end

plotSpread(data(:,1), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,1}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 -0.19 1.3])
ylabel('Minimal ratio of Cdc14 to Cfi1 in nucleolus')
saveas(gcf, 'Min ratio', 'fig')
saveas(gcf, 'Min ratio', 'pdf')

figure
plotSpread(data(:,2), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,2}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 0.2 1.5])
ylabel('Minimal CV of Cdc14')
saveas(gcf, 'Min CV', 'fig')
saveas(gcf, 'Min CV', 'pdf')

figure
plotSpread(data(:,4), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,4}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 0.2 1.2])
ylabel('Minimal CV of NLS_C_d_c_1_4')
saveas(gcf, 'Min NLS CV', 'fig')
saveas(gcf, 'Min NLS CV', 'pdf')

%% plot dynamics of Cdc14 release
Cdc14_ratio = cell(n_sample,1);
Cdc14_CV = cell(n_sample,1);
Cfi1_CV = cell(n_sample,1);
Cfi1_cell = cell(n_sample,1);
Spindle = cell(n_sample,1);
NLS_CV = cell(n_sample,1);
np_left = zeros(n_sample,1);
close all
f1 = figure; hold on
f2 = figure; hold on
f3 = figure; hold on
f4 = figure; hold on
f5 = figure; hold on
f6 = figure; hold on
for i = 1 : n_sample
    np_left(i) = max([samples{i}.SPBinbud] - [samples{i}.Budding]) + 1;
    np_right = max([samples{i}.Cytokinesis] - [samples{i}.SPBinbud]) + 1;
    Cdc14_ratio{i} = nan(length(samples{i}), np_left(i)+np_right + 1, 1);
    Cdc14_CV{i} = nan(length(samples{i}), np_left(i)+np_right + 1, 1);
    Cfi1_CV{i} = nan(length(samples{i}), np_left(i)+np_right + 1, 1);
    Cfi1_cell{i} = nan(length(samples{i}), np_left(i)+np_right + 1, 2);
    Spindle{i} = nan(length(samples{i}), np_left(i)+np_right + 1);
    NLS_CV{i} = nan(length(samples{i}), np_left(i)+np_right + 1, 1);
    for n = 1 : length(samples{i})
        shift = np_left(i) - (samples{i}(n).SPBinbud - samples{i}(n).Budding);
        Cdc14_ratio{i}(n, shift+1:shift+length(samples{i}(n).Cdc14), 1) = samples{i}(n).Cdc14(:,3);
        Cdc14_CV{i}(n, shift+1:shift+length(samples{i}(n).Cdc14), 1) = samples{i}(n).Cdc14(:,4);
        Cfi1_CV{i}(n, shift+1:shift+length(samples{i}(n).Cdc14), 1) = samples{i}(n).Cfi1(:,4);
        Cfi1_cell{i}(n, shift+1:shift+length(samples{i}(n).Cdc14), :) = samples{i}(n).Cfi1(:,2:3);% nucleoli and cell
        Spindle{i}(n, shift+1:shift+length(samples{i}(n).Cdc14)) = samples{i}(n).Spindle;
        NLS_CV{i}(n, shift+1:shift+length(samples{i}(n).NLS), 1) = samples{i}(n).NLS(:,1);
    end
    
    figure(f1)
    plot(dT*(-np_left(i):1:np_right), Cdc14_ratio{i}(:,:,1), ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Cdc14_ratio{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
    
    figure(f2)
    plot(dT*(-np_left(i):1:np_right), Cdc14_CV{i}(:,:,1), ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Cdc14_CV{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
    
    figure(f5)
    plot(dT*(-np_left(i):1:np_right), Cfi1_CV{i}(:,:,1), ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Cfi1_CV{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
    
    figure(f6)
    plot(dT*(-np_left(i):1:np_right), Cfi1_cell{i}(:,:,1), ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Cfi1_cell{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
    
    figure(f3)
    plot(dT*(-np_left(i):1:np_right), Spindle{i}*PixelSize, ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Spindle{i})*PixelSize, 'LineWidth', 3, 'Color', Color_order(i,:))
    
    figure(f4)
    plot(dT*(-np_left(i):1:np_right), NLS_CV{i}(:,:,1), ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(NLS_CV{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
end

figure(f1)
axis([-60 60 -0.1 2.5])
xlabel('Time since SPB in bud (min)')
ylabel('Ratio of Cdc14 to Cfi1 in nucleolus')
legend(label)
saveas(f1, 'Aligned Cdc14 in nucleolus','fig')
saveas(f1, 'Aligned Cdc14 in nucleolus','png')

figure(f2)
axis([-60 60 0.2 2.4])
xlabel('Time since SPB in bud (min)')
ylabel('Cdc14 release (CV)')
legend(label)
saveas(f2, 'Aligned Cdc14 release_CV','fig')
saveas(f2, 'Aligned Cdc14 release_CV','png')

figure(f3)
axis([-60 60 0 15])
xlabel('Time since SPB in bud (min)')
ylabel('Estimated nucleolar seggregation progress (um)')
legend(label,'location','northwest')
saveas(f3, 'Aligned Spinle length','fig')
saveas(f3, 'Aligned Spinle length','png')

figure(f4)
axis([-60 60 0.3 1.1])
xlabel('Time since SPB in bud (min)')
ylabel('NLS_C_d_c_1_4 release (CV)')
legend(label)
saveas(f4, 'Aligned NLS-Cdc14 release_CV','fig')
saveas(f4, 'Aligned NLS-Cdc14 release_CV','png')

figure(f5)
axis([-60 60 0.2 2.4])
xlabel('Time since SPB in bud (min)')
ylabel('Cfi1 (CV)')
legend(label)
saveas(f2, 'Aligned Cfi1_CV','fig')
saveas(f2, 'Aligned Cfi1_CV','png')

figure(f6)
axis([-60 60 0 200])
xlabel('Time since SPB in bud (min)')
ylabel('Cfi1 (CV)')
legend(label)
saveas(f2, 'Aligned Cfi1','fig')
saveas(f2, 'Aligned Cfi1','png')

%% plot for comparison with 95% CI
close all
f1 = figure; hold on
f2 = figure; hold on
f3 = figure; hold on
f4 = figure; hold on
for i = 1 : n_sample
    np_right = max([samples{i}.Cytokinesis] - [samples{i}.SPBinbud]) + 1;
    
    figure(f1)
    x = dT*(-np_left(i):1:np_right);
    upper = nanmean(Cdc14_ratio{i}(:,:,1))+1.96*nanstd(Cdc14_ratio{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_ratio{i}(:,:,1))));
    lower = nanmean(Cdc14_ratio{i}(:,:,1))-1.96*nanstd(Cdc14_ratio{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_ratio{i}(:,:,1))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Cdc14_ratio{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))

    figure(f2)
    upper = nanmean(Cdc14_CV{i}(:,:,1))+1.96*nanstd(Cdc14_CV{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_CV{i}(:,:,1))));
    lower = nanmean(Cdc14_CV{i}(:,:,1))-1.96*nanstd(Cdc14_CV{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_CV{i}(:,:,1))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Cdc14_CV{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))

    figure(f3)
    upper2 = nanmean(Spindle{i})+1.96*nanstd(Spindle{i})./sqrt(sum(~isnan(Spindle{i})));
    lower2 = nanmean(Spindle{i})-1.96*nanstd(Spindle{i})./sqrt(sum(~isnan(Spindle{i})));
    ha2 = area(x', [lower2*PixelSize; (upper2-lower2)*PixelSize]', 'HandleVisibility','off');
    set(ha2(1), 'FaceColor', 'none')
    set(ha2, 'LineStyle', 'none')
    set(ha2(2), 'FaceColor', Color_order(i,:))
    set(ha2(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Spindle{i})*PixelSize, 'LineWidth', 3, 'Color', Color_order(i,:))
    
    figure(f4)
    upper = nanmean(NLS_CV{i}(:,:,1))+1.96*nanstd(NLS_CV{i}(:,:,1))./sqrt(sum(~isnan(NLS_CV{i}(:,:,1))));
    lower = nanmean(NLS_CV{i}(:,:,1))-1.96*nanstd(NLS_CV{i}(:,:,1))./sqrt(sum(~isnan(NLS_CV{i}(:,:,1))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(NLS_CV{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
end

figure(f1)
axis([-60 60 0 2])
xlabel('Time since SPB in bud (min)')
ylabel('Ratio of Cdc14 to Cfi1 in nucleolus')
legend(label)
saveas(f1, 'Aligned Cdc14 in nucleolus2','fig')
saveas(f1, 'Aligned Cdc14 in nucleolus2','pdf')

figure(f2)
axis([-60 60 0.2 2.5])
xlabel('Time since SPB in bud (min)')
ylabel('Cdc14 release (CV)')
legend(label)
saveas(f2, 'Aligned Cdc14 release_CV2','fig')
saveas(f2, 'Aligned Cdc14 release_CV2','pdf')

figure(f3)
axis([-60 60 0 15])
xlabel('Time since SPB in bud (min)')
ylabel('Estimated nucleolar seggregation progress (um)')
legend(label,'location','northwest')
saveas(f3, 'Aligned Spinle length2','fig')
saveas(f3, 'Aligned Spinle length2','pdf')

figure(f4)
axis([-60 60 0.3 1.1])
xlabel('Time since SPB in bud (min)')
ylabel('NLS_C_d_c_1_4 release (CV)')
legend(label)
saveas(f4, 'Aligned NLS-Cdc14 release_CV2','fig')
saveas(f4, 'Aligned NLS-Cdc14 release_CV2','pdf')

%% normalize to -20min
close all
f1 = figure; hold on
f2 = figure; hold on
f3 = figure; hold on
Cdc14_ratio_norm = cell(n_sample,1);
Cdc14_CV_norm = cell(n_sample,1);
NLS_CV_norm = cell(n_sample,1);
for i = 1 : n_sample
    np_right = max([samples{i}.Cytokinesis] - [samples{i}.SPBinbud]) + 1;
    x = dT*(-np_left(i):1:np_right);
    Cdc14_ratio_norm{i} = Cdc14_ratio{i}/nanmean(Cdc14_ratio{i}(:, np_left(i)-3));
    upper = nanmean(Cdc14_ratio_norm{i})+1.96*nanstd(Cdc14_ratio_norm{i})./sqrt(sum(~isnan(Cdc14_ratio_norm{i})));
    lower = nanmean(Cdc14_ratio_norm{i})-1.96*nanstd(Cdc14_ratio_norm{i})./sqrt(sum(~isnan(Cdc14_ratio_norm{i})));
    figure(f1)
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Cdc14_ratio_norm{i}(:,:,1)), 'LineWidth', 2, 'Color', Color_order(i,:))
    
    Cdc14_CV_norm{i} = Cdc14_CV{i}/nanmean(Cdc14_CV{i}(:, np_left(i)-3));
    upper = nanmean(Cdc14_CV_norm{i})+1.96*nanstd(Cdc14_CV_norm{i})./sqrt(sum(~isnan(Cdc14_CV_norm{i})));
    lower = nanmean(Cdc14_CV_norm{i})-1.96*nanstd(Cdc14_CV_norm{i})./sqrt(sum(~isnan(Cdc14_CV_norm{i})));
    figure(f2)
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Cdc14_CV_norm{i}(:,:,1)), 'LineWidth', 2, 'Color', Color_order(i,:))
    
    NLS_CV_norm{i} = NLS_CV{i}/nanmean(NLS_CV{i}(:, np_left(i)-3));
    upper = nanmean(NLS_CV_norm{i})+1.96*nanstd(NLS_CV_norm{i})./sqrt(sum(~isnan(NLS_CV_norm{i})));
    lower = nanmean(NLS_CV_norm{i})-1.96*nanstd(NLS_CV_norm{i})./sqrt(sum(~isnan(NLS_CV_norm{i})));
    figure(f3)
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(NLS_CV_norm{i}(:,:,1)), 'LineWidth', 2, 'Color', Color_order(i,:))
end

figure(f1)
axis([-60 60 0 1.2])
xlabel('Time since SPB in bud (min)')
ylabel('Ratio of Cdc14 to Cfi1 in nucleolus-normalized')
legend(label,'location','southwest')
saveas(f1, 'Aligned Cdc14 in nucleolus_norm','fig')
saveas(f1, 'Aligned Cdc14 in nucleolus_norm','pdf')

figure(f2)
axis([-60 60 0 1.2])
xlabel('Time since SPB in bud (min)')
ylabel('Cdc14 release (CV)-normalized')
legend(label,'location','southwest')
saveas(f2, 'Aligned Cdc14 release_CV_norm','fig')
saveas(f2, 'Aligned Cdc14 release_CV_norm','pdf')

figure(f3)
axis([-60 60 0.4 1.2])
xlabel('Time since SPB in bud (min)')
ylabel('NLS_C_d_c_1_4 release (CV)-normalized')
legend(label,'location','southwest')
saveas(f3, 'Aligned NLS-Cdc14 release_CV_norm','fig')
saveas(f3, 'Aligned NLS-Cdc14 release_CV_norm','pdf')

%% Cdc14 release only
mkdir('Cdc14')
close all
figure('Position', [10 10 900 500])


for i = 1 : n_sample
    subplot(2, n_sample, i); hold on
    np_right = max([samples{i}.Cytokinesis] - [samples{i}.SPBinbud]) + 1;
    plot(dT*(-np_left(i):1:np_right), Cdc14_ratio{i}(:,:,1), ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Cdc14_ratio{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))  
    axis([-40 60 0 2])
    xlabel('Time since SPB in bud (min)')
    ylabel('Cdc14 release (ratio)')
    legend(label{i})
   
    subplot(2, n_sample, i+n_sample); hold on
    plot(dT*(-np_left(i):1:np_right), Cdc14_CV{i}(:,:,1), ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right), nanmean(Cdc14_CV{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))  
    axis([-40 60 0.2 2.4])
    xlabel('Time since SPB in bud (min)')
    ylabel('Cdc14 release (CV)')
%     legend(label{i})
end


saveas(gcf, 'Cdc14/Aligned Cdc14 release_ratio','fig')
saveas(gcf, 'Cdc14/Aligned Cdc14 release_ratio','png')

%%
close all
figure('Position', [10 10 400 500])

for i = 1 : n_sample
    subplot(2, 2, 1); hold on
    np_right = max([samples{i}.Cytokinesis] - [samples{i}.SPBinbud]) + 1;
    x = dT*(-np_left(i):1:np_right);
    upper = nanmean(Cdc14_ratio{i}(:,:,1))+1.96*nanstd(Cdc14_ratio{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_ratio{i}(:,:,1))));
    lower = nanmean(Cdc14_ratio{i}(:,:,1))-1.96*nanstd(Cdc14_ratio{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_ratio{i}(:,:,1))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Cdc14_ratio{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
    axis([-40 60 0 2])
    xlabel('Time since SPB in bud (min)')
    ylabel('Cdc14 release (ratio)')
    
    subplot(2, 2, 2); hold on
    upper2 = nanmean(Cdc14_CV{i}(:,:,1))+1.96*nanstd(Cdc14_CV{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_CV{i}(:,:,1))));
    lower2 = nanmean(Cdc14_CV{i}(:,:,1))-1.96*nanstd(Cdc14_CV{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_CV{i}(:,:,1))));
    ha2 = area(x', [lower2; (upper2-lower2)]', 'HandleVisibility','off');
    set(ha2(1), 'FaceColor', 'none')
    set(ha2, 'LineStyle', 'none')
    set(ha2(2), 'FaceColor', Color_order(i,:))
    set(ha2(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Cdc14_CV{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
    axis([-40 60 0.2 2.4])
    xlabel('Time since SPB in bud (min)')
    ylabel('Cdc14 release (CV)')
    
    subplot(2, 2, 3); hold on
    np_right = max([samples{i}.Cytokinesis] - [samples{i}.SPBinbud]) + 1;
    x = dT*(-np_left(i):1:np_right);
    upper = nanmean(Cdc14_ratio_norm{i}(:,:,1))+1.96*nanstd(Cdc14_ratio_norm{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_ratio_norm{i}(:,:,1))));
    lower = nanmean(Cdc14_ratio_norm{i}(:,:,1))-1.96*nanstd(Cdc14_ratio_norm{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_ratio_norm{i}(:,:,1))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Cdc14_ratio_norm{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
    axis([-40 60 0 1.2])
    xlabel('Time since SPB in bud (min)')
    ylabel('Cdc14 release (ratio)-normalized')
    
    subplot(2, 2, 4); hold on
    upper2 = nanmean(Cdc14_CV_norm{i}(:,:,1))+1.96*nanstd(Cdc14_CV_norm{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_CV_norm{i}(:,:,1))));
    lower2 = nanmean(Cdc14_CV_norm{i}(:,:,1))-1.96*nanstd(Cdc14_CV_norm{i}(:,:,1))./sqrt(sum(~isnan(Cdc14_CV_norm{i}(:,:,1))));
    ha2 = area(x', [lower2; (upper2-lower2)]', 'HandleVisibility','off');
    set(ha2(1), 'FaceColor', 'none')
    set(ha2, 'LineStyle', 'none')
    set(ha2(2), 'FaceColor', Color_order(i,:))
    set(ha2(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Cdc14_CV_norm{i}(:,:,1)), 'LineWidth', 3, 'Color', Color_order(i,:))
    axis([-40 60 0 1.2])
    xlabel('Time since SPB in bud (min)')
    ylabel('Cdc14 release (CV)-normalized')
end

legend(label,'location','southwest')
saveas(gcf, 'Cdc14/Aligned Cdc14 release_CI','fig')
saveas(gcf, 'Cdc14/Aligned Cdc14 release_CI','pdf')

%% Plot relative degree of Cdc14 release
close all
figure('Position', [10 10 500 800])
data_norm = cell(n_sample,3);
Cdc14_ratio_norm2 = cell(n_sample,1);
Cdc14_CV_norm2 = cell(n_sample,1);
for i = 1 : n_sample
    Cdc14_ratio_norm2{i} = Cdc14_ratio{i}./Cdc14_ratio{i}(:, np_left(i)-3);
    Cdc14_CV_norm2{i} = Cdc14_CV{i}./Cdc14_CV{i}(:, np_left(i)-3);
    data_norm{i,1} = 1- min(Cdc14_ratio_norm2{i}(:,np_left+1:end), [],2);
    data_norm{i,2} = Color_order(i, :);
    data_norm{i,3} = 1- min(Cdc14_CV_norm2{i}(:,np_left+1:end), [],2);
end
subplot(2, 1, 1); hold on
plotSpread(data_norm(:,1), ...
    'xNames', label, ...
    'distributionColors', data_norm(:,2), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*nanmedian(data_norm{i,1}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 -0.19 1.1])
ylabel('Relative change in ratio of Cdc14 to Cfi1')

subplot(2, 1, 2); hold on
plotSpread(data_norm(:,3), ...
    'xNames', label, ...
    'distributionColors', data_norm(:,2), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*nanmedian(data_norm{i,3}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 -0 1.1])
ylabel('Relative change in CV of Cdc14')

saveas(gcf, 'Cdc14/Cdc14 release', 'fig')
saveas(gcf, 'Cdc14/Cdc14 release', 'pdf')

%% check Cfi1 expression
close all
figure('Position', [10 10 500 800])
data = cell(n_sample,3);
for i = 1 : n_sample
    data{i,1} = nanmean(Cfi1_cell{i}(:,:,1), 2); % nucleoli
    data{i,2} = nanmean(Cfi1_cell{i}(:,:,2), 2); % cell
    data{i,3} = Color_order(i, :);
end

subplot(2, 1, 1); hold on
plotSpread(data(:,1), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,1}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 0 200])
ylabel('Median intensity in nucleoli')

subplot(2, 1, 2); hold on
plotSpread(data(:,2), ...
    'xNames', label, ...
    'distributionColors', data(:,3), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,2}), 'Color', Color_order(i, :),'linewidth', 3)
end
axis([0.25 n_sample+0.75 0 40])
ylabel('Mean intensity in cell')

saveas(gcf, 'Cdc14/Cfi1', 'fig')
saveas(gcf, 'Cdc14/Cfi1', 'pdf')

%% p values
p_ttest = nan(n_sample,n_sample);
p_ranksum = nan(n_sample,n_sample);
for i = 1 : n_sample
    for j = i+1 : n_sample
        [~, p_ttest(i,j)] = ttest2(data_norm{i,1},data_norm{j,1});
        p_ranksum(i,j) = ranksum(data_norm{i,1},data_norm{j,1});
    end
end

%%
save('summary')

