%% summarize results
clear
close all
%input sample info
addpath('/plotSpread');

dirname = '';
cd(dirname)

folder = {};
label = {};

%% combine all the data for each sample
n_sample = length(folder);
samples = [];
for n = 1 : n_sample 
    list_file = dir([folder{n} '/*.tif']);
    for i = 1 : length(list_file)
        load([folder{n} '/' list_file(i).name(1:end-4) '/analysis'], 'FRAP_info')
        FRAP_info.ID = n;
        FRAP_info.file = [folder{n} '/' list_file(i).name(1:end-4)];
        samples = [samples; FRAP_info];                       
    end
end

mkdir('summary')
cd('summary')

%%
close all
figure
hold on
Color_order = get(gca,'colororder');
bleach = 3;
I_FRAP = nan(length(samples), length(samples(1).cell_intensity));
for n = 1:length(samples)
    E_spb = samples(n).SPB_intensity./samples(n).cell_intensity-1;
    % find the SPB that is bleached
    [~, id1] = max(E_spb(bleach-1,:)-E_spb(bleach,:));
    id2 = setdiff([1 2],id1); % the other SPB
    % double normalization
    I_FRAP(n,:) = (E_spb(:,id1)-E_spb(bleach,id1))./(mean(E_spb(1:bleach-1,id1))-E_spb(bleach,id1));

    subplot(2,2,1)
    hold on
    plot([0; samples(n).time]-2, samples(n).cell_intensity, 'k')        
    plot([0; samples(n).time]-2, samples(n).SPB_intensity(:,id1), 'Color', Color_order(samples(n).ID,:))
    if ~(samples(n).ID==2)
    plot([0; samples(n).time]-2, samples(n).SPB_intensity(:,id2),':', 'Color', Color_order(samples(n).ID,:),'LineWidth', 1)
    end

    subplot(2,2,2)
    hold on
    plot([0; samples(n).time]-2, E_spb(:,id1), 'Color', Color_order(samples(n).ID,:))
    if ~(samples(n).ID==2)
    plot([0; samples(n).time]-2, E_spb(:,id2),':', 'Color', Color_order(samples(n).ID,:),'LineWidth', 1)
    end

    subplot(2,2,3)
    hold on
    plot([0; samples(n).time]-2, I_FRAP(n,:), 'Color', Color_order(samples(n).ID,:))
    
    subplot(2,2,4)
    hold on
    plot(samples(n).fit_result(:,1) -2, samples(n).fit_result(:,2), 'Color', Color_order(samples(n).ID,:))
end

subplot(2,2,1)
xlabel('Time (s)')
ylabel('Raw intensity')
legend({'Cytosol' 'Bleached SPB' 'Ref SPB'})
subplot(2,2,2)
xlabel('Time (s)')
ylabel('Enrichment at SPB')
legend({'Bleached SPB' 'Ref SPB'})
subplot(2,2,3)
xlabel('Time (s)')
ylabel('Double normalized I')
subplot(2,2,4)
xlabel('Time (s)')
ylabel('Fit')

%% plot average for dSPB
close all
figure('Position',[1000 200 250 500])
T = samples(n).time-2;
I_FRAP1 = [samples.FRAP_doubleNorm]';
subplot(2,1,1)
errorbar(T, mean(I_FRAP1(1:5,:)), std(I_FRAP1(1:5,:)),'ko')
hold on
I_fit = [samples.fit_result]';
I_fit(1:2:end,:) = [];
upper = mean(I_fit(1:5,:))+std(I_fit(1:5,:));
lower = mean(I_fit(1:5,:))-std(I_fit(1:5,:));
ha = area(T, [lower; (upper-lower)]', 'HandleVisibility','off');
set(ha(1), 'FaceColor', 'none')
set(ha, 'LineStyle', 'none')
set(ha(2), 'FaceColor', Color_order(2,:))
set(ha(2), 'FaceAlpha', 0.3)  
plot(T, mean(I_fit(1:5,:)),'Color', Color_order(2,:),'LineWidth',2)
xlabel('time (s)')
ylabel('Fraction of initial intensity')
axis([-2 30 0 1.3])

subplot(2,1,2)
plot(T, mean(I_FRAP1(1:5,:)),'ko')
hold on
ha = area(T, [lower; (upper-lower)]', 'HandleVisibility','off');
set(ha(1), 'FaceColor', 'none')
set(ha, 'LineStyle', 'none')
set(ha(2), 'FaceColor', Color_order(2,:))
set(ha(2), 'FaceAlpha', 0.2)  
plot(T, mean(I_fit(1:5,:)),'Color', Color_order(2,:),'LineWidth',2)
xlabel('time (s)')
ylabel('Fraction of initial intensity')
axis([-2 30 0 1.3])
saveas(gcf, 'FRAP_CI_dSPB','pdf')

% average for dSPB
HT = [samples.HT];
stats_HT_dSPB = [mean(HT(1:5)) std(HT(1:5))];

%% plot average for SPB
close all
figure('Position',[1000 200 250 500])
subplot(2,1,1)
errorbar(T, mean(I_FRAP1(1:6,:)), std(I_FRAP1(1:6,:)),'ko')
hold on
upper = mean(I_fit(1:6,:))+std(I_fit(1:6,:));
lower = mean(I_fit(1:6,:))-std(I_fit(1:6,:));
ha = area(T, [lower; (upper-lower)]', 'HandleVisibility','off');
set(ha(1), 'FaceColor', 'none')
set(ha, 'LineStyle', 'none')
set(ha(2), 'FaceColor', Color_order(2,:))
set(ha(2), 'FaceAlpha', 0.3)  
plot(T, mean(I_fit(1:6,:)),'Color', Color_order(2,:),'LineWidth',2)
xlabel('time (s)')
ylabel('Fraction of initial intensity')
axis([-2 30 0 1.3])

subplot(2,1,2)
plot(T, mean(I_FRAP1(1:6,:)),'ko')
hold on
ha = area(T, [lower; (upper-lower)]', 'HandleVisibility','off');
set(ha(1), 'FaceColor', 'none')
set(ha, 'LineStyle', 'none')
set(ha(2), 'FaceColor', Color_order(2,:))
set(ha(2), 'FaceAlpha', 0.2)  
plot(T, mean(I_fit(1:6,:)),'Color', Color_order(2,:),'LineWidth',2)
xlabel('time (s)')
ylabel('Fraction of initial intensity')
axis([-2 30 0 1.3])
saveas(gcf, 'FRAP_CI','pdf')

% average for all SPB
stats_HT = [mean(HT(1:6)) std(HT(1:6))];

%%
close all
save('summary')