%% summarize results for comparison
clear
close all
addpath('/plotSpread');
%input sample info
dirname = '';
cd(dirname)

folder = {};
label = {};

%% combine all the data from different fields for each sample
n_sample = length(folder);
samples = cell(n_sample,1);
f = figure; 
set(f,'visible','off');
for n = 1 : n_sample 
    list_file = dir([folder{n} '/C3*.tif']);
    for i = 1 : length(list_file)
        load([folder{n} '/' list_file(i).name(1:end-4) '/Cfi1.mat'], 'Division', 'dT', 'PixelSize', 'nframes')
        for j = 1 : length(Division) 
            if isempty(Division{j})
                continue
            end
            if max(Division{j}.Spindle(:,2)) < 6
                Division{j} = [];
                continue
            end 
            if isempty(Division{j}.SPBinbud)
                Division{j} = [];
                continue
            end
            [~, mitotic_exit] = max(Division{j}.Spindle(:,2));
            if mitotic_exit + Division{j}.Metaphase > nframes
                Division{j} = [];
                continue
            end
            Division{j}.filename = list_file(i).name(1:end-4);
            Division{j}.ID = j;
            
            Division{j}.nuc_onset = find(Division{j}.nucleoli(1:end,2) > 3, 1, 'first') + Division{j}.Metaphase -1;
            if Division{j}.nuc_onset < Division{j}.Ana_onset
                Division{j}.nuc_onset = NaN;
            end
            if isempty(Division{j}.nuc_onset)
                Division{j}.nuc_onset = NaN;
            end
            Division{j}.mitotic_exit = mitotic_exit + Division{j}.Metaphase;
            
            T = (1:1:length(Division{j}.nucleoli)) + Division{j}.Metaphase - 1;
            plot(T, Division{j}.Spindle(:,2),'o-') % Spindle
            hold on            
            plot(T, Division{j}.nucleoli(:,1),'o-') % nucleoli, majoraxis
            plot(T, Division{j}.nucleoli(:,2),'o-') % nucleoli, boundingbox
            if ~isempty(Division{j}.NucleoliInBud)
                plot([1 1]*Division{j}.NucleoliInBud, [0 10],'k--') % Nucleoli in Bud
                plot([1 1]*Division{j}.SPBinbud, [0 10],'k:') % SPB in Bud
                plot([1 1]*Division{j}.nuc_onset, [0 10],'r--') % Nucleoli seggregation
                plot([1 1]*Division{j}.Ana_onset, [0 10],'r:') % Anaphase onset
            end
            xlabel('Frame #')
            ylabel('Length (um)')
            legend({'Spindle' 'Nucleoli-MajorAxisLength', 'Nucleoli-BoundingBox', 'Nucleoli In Bud', 'SPB in Bud', 'Nucleoli seggregation', 'Anaphase onset'}, 'Location', 'northwest')

            saveas(f, [folder{n} '/' list_file(i).name(1:end-4) '/plots/Nucleoli_seggregation ' num2str(j)], 'png')
            hold off
    
        end
        Division = Division(~cellfun(@isempty,Division));
        samples{n} = [samples{n}; cell2mat(Division)]; 

    end
end

mkdir('summary')
cd('summary')

%% plot anaphase duration
close all
figure('Position', [10 10 400 1200]);
Color_order = get(gca,'colororder');
ana_duration = cell(n_sample,1);
nuc_delay = cell(n_sample,2);
filter = cell(n_sample,1);
for i = 1 : n_sample
    subplot(n_sample, 1, i)
    ana_duration{i} = [samples{i}.mitotic_exit] - [samples{i}.Ana_onset];
    nuc_delay{i,1} = [samples{i}.nuc_onset] - [samples{i}.Ana_onset];
    nuc_delay{i,2} = [samples{i}.NucleoliInBud] - [samples{i}.Ana_onset];
    filter{i} = ([samples{i}.mitotic_exit] - [samples{i}.Metaphase])' < cell2mat(arrayfun(@(x) length(x.Spindle), samples{i}, 'UniformOutput', false));
    h = histogram(dT*(ana_duration{i}(filter{i})), (2*dT+dT/2:dT:dT*22), 'FaceColor', Color_order(i,:));
    hold on
    plot([1 1]*median(ana_duration{i}(filter{i}))*dT,[0 max(h.Values)],'k:', 'LineWidth', 3)
    title(label{i})
    ylabel('Cell count')
    ax = gca;
    set(ax,'XTick',(dT:dT:dT*22))
end
xlabel('Anaphase duration (min)')
saveas(gcf, 'Anaphase duration distribution', 'fig')
saveas(gcf, 'Anaphase duration distribution', 'png')
%% compare SPB in bud with anaphase onset
close all
figure('Position', [10 10 400 800])
ana_duration2 = cell(n_sample,1);
for i = 1 : n_sample
    subplot(n_sample, 1, i)
    ana_duration2{i} = [samples{i}.mitotic_exit] - [samples{i}.SPBinbud];
    h = histogram(dT*(ana_duration2{i}(filter{i}) - ana_duration{i}(filter{i})), (-4*dT+dT/2:dT:dT*20));
    legend(label{i})
    ylabel('Cell count')
end
xlabel('Time from SPB in bud to anaphase onset (min)')
saveas(gcf, 'SPBinbud to ana_onset distribution', 'fig')
saveas(gcf, 'SPBinbud to ana_onset distribution', 'png')

figure('Position', [10 10 300 800])
for i = 1 : n_sample
    subplot(n_sample, 1, i)
    scatter(dT*(ana_duration2{i}(filter{i}) - ana_duration{i}(filter{i})), ana_duration{i}(filter{i})*dT, [], 'filled', 'MarkerfaceAlpha', 5/sum(filter{i}));
    axis([-10 50 0 60])
    legend(label{i})
    xlabel('Time from SPB in bud to anaphase onset (min)')
    ylabel('Anaphase duration (min)')
end
saveas(gcf, 'SPBinbud to mitotic exit', 'fig')
saveas(gcf, 'SPBinbud to mitotic exit', 'png')

%% replot the anaphase duration after filter out cells with SAC activated
close all
figure('Position', [10 10 400 1000])
filter2 = cell(n_sample,1); % filter out cells have activated SAC (long metaphase where SPB wonders to bud)
for i = 1 : n_sample
    subplot(n_sample, 1, i)
    filter2{i} = ([samples{i}.Ana_onset] - [samples{i}.SPBinbud])' < 10/dT; % 10min arrest in metaphase
    h = histogram(dT*(ana_duration{i}(filter2{i})), (2*dT+dT/2:dT:dT*22), 'FaceColor', Color_order(i,:));
    hold on
    plot([1 1]*median(ana_duration{i}(filter2{i}))*dT,[0 max(h.Values)],'k:', 'LineWidth', 3)
    title(label{i})
    ylabel('Cell count')
    ax = gca;
    set(ax,'XTick',(2*dT:dT:dT*22))
end
xlabel('Anaphase duration (min)')
saveas(gcf, 'Anaphase duration distribution_clean', 'fig')
saveas(gcf, 'Anaphase duration distribution_clean', 'png')

%% align anaphase
close all
Nucleoli = cell(n_sample,1);
Spindle = cell(n_sample,1);
np_left = zeros(n_sample,1);
np_right = zeros(n_sample,1);
f1 = figure; hold on
f2 = figure; hold on
for i = 1 : n_sample
    np_left(i) = max([samples{i}.Ana_onset] - [samples{i}.Metaphase]) + 1;
    np_right(i) = max([samples{i}.Cytokinesis] - [samples{i}.Ana_onset]) + 25;
    Nucleoli{i} = nan(length(samples{i}), np_left(i)+np_right(i) + 1, 2);
    Spindle{i} = nan(length(samples{i}), np_left(i)+np_right(i) + 1, 2);
    for n = 1 : length(samples{i})
        shift = np_left(i) - (samples{i}(n).Ana_onset - samples{i}(n).Metaphase);
        Nucleoli{i}(n, shift+1:shift+length(samples{i}(n).nucleoli), :) = samples{i}(n).nucleoli;
        Spindle{i}(n, shift+1:shift+length(samples{i}(n).Spindle), :) = samples{i}(n).Spindle;
    end
    figure(f1)
    l = np_left(i)+np_right(i)+1;
    plot(dT*(-np_left(i):1:np_right(i)), Spindle{i}(:,1:l,2), ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right(i)), nanmean(Spindle{i}(:,1:l,2)), 'LineWidth', 3, 'Color', Color_order(i,:))

    figure(f2)
    plot(dT*(-np_left(i):1:np_right(i)), Nucleoli{i}(:,1:l,2), ':', 'LineWidth', 0.1, 'Color', Color_order(i,:), 'HandleVisibility','off')
    plot(dT*(-np_left(i):1:np_right(i)), nanmean(Nucleoli{i}(:,1:l,2)), 'LineWidth', 3, 'Color', Color_order(i,:))

end

figure(f1)
axis([-40 60 0 10])
xlabel('Time since anaphase onset (min)')
ylabel('Spinle length 3D (um)')
legend(label, 'location', 'northwest')
saveas(f1, 'Aligned Spinle length','png')

figure(f2)
axis([-40 60 0 10])
xlabel('Time since anaphase onset (min)')
ylabel('Nucleoli length (um)')
legend(label, 'location', 'northwest')
saveas(f2, 'Aligned Nucleoli length','png')

%% plot CI
close all
figure('Position', [10 10 600 500]); hold on

for i = 1 : n_sample
    x = dT*(-np_left(i):1:np_right(i));
    
    subplot(2,2,1); hold on
    l = np_left(i)+np_right(i)+1;
    upper = nanmean(Spindle{i}(:,1:l,2))+1.96*nanstd(Spindle{i}(:,1:l,2))./sqrt(sum(~isnan(Spindle{i}(:,1:l,2))));
    lower = nanmean(Spindle{i}(:,1:l,2))-1.96*nanstd(Spindle{i}(:,1:l,2))./sqrt(sum(~isnan(Spindle{i}(:,1:l,2))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Spindle{i}(:,1:l,2)), 'LineWidth', 2, 'Color', Color_order(i,:))
    axis([-40 60 0 10])
    xlabel('Time since anaphase onset (min)')
    ylabel('Spindle length (um)')
    
    subplot(2,2,2); hold on
    upper = nanmean(Nucleoli{i}(:,1:l,2))+1.96*nanstd(Nucleoli{i}(:,1:l,2))./sqrt(sum(~isnan(Nucleoli{i}(:,1:l,2))));
    lower = nanmean(Nucleoli{i}(:,1:l,2))-1.96*nanstd(Nucleoli{i}(:,1:l,2))./sqrt(sum(~isnan(Nucleoli{i}(:,1:l,2))));
    ha = area(x', [lower; (upper-lower)]', 'HandleVisibility','off');
    set(ha(1), 'FaceColor', 'none')
    set(ha, 'LineStyle', 'none')
    set(ha(2), 'FaceColor', Color_order(i,:))
    set(ha(2), 'FaceAlpha', 0.2)    
    plot(x, nanmean(Nucleoli{i}(:,1:l,2)), 'LineWidth', 2, 'Color', Color_order(i,:))
    axis([-40 60 0 10])
    xlabel('Time since anaphase onset (min)')
    ylabel('Nucleoli length (um)')
    
    subplot(2,2,3); hold on
    histogram(dT*(ana_duration{i}(filter{i}&filter2{i})), (-dT*1.5:dT:dT*15.5), 'Normalization','cdf', 'DisplayStyle', 'stairs','linewidth', 1.5);
    axis([0 45 0 1.1])
    xlabel('Anaphase duration (min)')
    ylabel('Cumulative density')
    
    subplot(2,2,4); hold on
    histogram(dT*(nuc_delay{i,1}(filter{i}&filter2{i})), (-dT*1.5:dT:dT*15.5), 'Normalization','cdf', 'DisplayStyle', 'stairs','linewidth', 1.5);
    axis([0 30 0 1.1])
    xlabel('Delay of nucleolar seggregation (min)')
    ylabel('Cumulative density')
end

legend(label,'location','southeast')
saveas(gcf, 'Aligned_CI','fig')
saveas(gcf, 'Aligned_CI','pdf')

%%
close all
figure('Position', [10 10 250 700]);
data = cell(n_sample,4);
for i = 1 : n_sample
    data{i,1} = dT*ana_duration{i}(filter{i}&filter2{i}); % nucleolus
    data{i,2} = dT*nuc_delay{i,1}(filter{i}&filter2{i}); % onset
    data{i,3} = dT*nuc_delay{i,2}(filter{i}&filter2{i}); % in bud
    data{i,4} = Color_order(i, :);
end

subplot(3, 1, 1); hold on
plotSpread(data(:,1), ...
    'xNames', label, ...
    'distributionColors', data(:,4), 'distributionMarkers', repmat({'o'},n_sample, 1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,1}), 'k','linewidth', 3)
end
axis([0.25 n_sample+0.75 0 60])
ylabel('Anaphase duration (min)')

subplot(3, 1, 2); hold on
plotSpread(data(:,2), ...
    'xNames', label, ...
    'distributionColors', data(:,4), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*nanmedian(data{i,2}), 'k','linewidth', 3)
end
axis([0.25 n_sample+0.75 0 40])
ylabel('Delay of nucleolar seggregation-onset (min)')

subplot(3, 1, 3); hold on
plotSpread(data(:,3), ...
    'xNames', label, ...
    'distributionColors', data(:,4), 'distributionMarkers', repmat({'o'},n_sample,1));
for i = 1 : n_sample
    plot([-0.4 0.4]+i, [1 1]*median(data{i,3}), 'k','linewidth', 3)
end
axis([0.25 n_sample+0.75 0 40])
ylabel('Delay of nucleolar seggregation-InBud (min)')

saveas(gcf, 'Spread', 'fig')
saveas(gcf, 'Spread', 'pdf')

%% compare for each cell
figure;
hold on
for i = 1 : n_sample
    scatter(dT*ana_duration{i}(filter{i}&filter2{i}), dT*nuc_delay{i,1}(filter{i}&filter2{i}),[],Color_order(i,:))
    scatter(dT*ana_duration{i}(filter{i}&filter2{i}), dT*nuc_delay{i,2}(filter{i}&filter2{i}),[],Color_order(i,:),'filled','MarkerFaceAlpha',0.2)
    
end


%%
close all
save('summary')

