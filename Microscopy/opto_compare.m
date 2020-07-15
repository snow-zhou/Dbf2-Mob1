%% summarize results for comparison
clear
close all
%input sample info
dirname = '';
cd(dirname)

load('intensity')
label_stages = {'G1-S' 'meta' 'ana'};
%% plot CI
close all
figure('Position', [10 10 600 500]); hold on
Color_order = get(gca,'colororder');
for i = 1 : n_sample/2
    subplot(2,3,i); hold on
    for j = 1 : n_stage
        slope = [samples{i,j}.data.slope];
        slope = slope(:, 1:2:end);
    
        upper = nanmean(slope,2)+1.96*nanstd(slope,[],2)./sqrt(sum(~isnan(slope),2));
        lower = nanmean(slope,2)-1.96*nanstd(slope,[],2)./sqrt(sum(~isnan(slope),2));
        ha = area(T', [lower (upper-lower)], 'HandleVisibility','off');
        set(ha(1), 'FaceColor', 'none')
        set(ha, 'LineStyle', 'none')
        set(ha(2), 'FaceColor', Color_order(j,:))
        set(ha(2), 'FaceAlpha', 0.2)    
        plot(T, nanmean(slope,2),'LineWidth', 2, 'Color', Color_order(j,:))
    end
    axis([0 5 0 1])
    ax = gca;
    set(ax,'XTick',[0 0.5 1.5 2.5 3.5 4.5])
    set(ax,'YTick',(0:0.5:1))
    ylabel('Relative enrichment of PIF in anchored region')
    xlabel('Time of 650 on (min)')
    legend(label_stages,'location','northwest')
    title(label{i})
end

for i = 3 : n_sample
    subplot(2,3,i+1); hold on
    for j = 1 : n_stage
        slope = [samples{i,j}.data.slope];
        slope = slope(:, 1:2:end);
    
        upper = nanmean(slope,2)+1.96*nanstd(slope,[],2)./sqrt(sum(~isnan(slope),2));
        lower = nanmean(slope,2)-1.96*nanstd(slope,[],2)./sqrt(sum(~isnan(slope),2));
        ha = area(T', [lower (upper-lower)], 'HandleVisibility','off');
        set(ha(1), 'FaceColor', 'none')
        set(ha, 'LineStyle', 'none')
        set(ha(2), 'FaceColor', Color_order(j,:))
        set(ha(2), 'FaceAlpha', 0.2)    
        plot(T, nanmean(slope,2),'LineWidth', 2, 'Color', Color_order(j,:))
    end
    axis([0 5 0 2])
    ax = gca;
    set(ax,'XTick',[0 0.5 1.5 2.5 3.5 4.5])
    set(ax,'YTick',(0:0.5:2))
    ylabel('Relative enrichment of PIF in anchored region')
    xlabel('Time of 650 on (min)')
    legend(label_stages,'location','northwest')
    title(label{i})
end
saveas(gcf, 'Compare_slope','fig')
saveas(gcf, 'Compare_slope','pdf')


%%



