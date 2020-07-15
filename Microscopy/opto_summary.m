%% summarize results for comparison
clear
close all
%input sample info
dirname = '';
cd(dirname)

folder = {};
label = {};
subfolder = {'G1-S' 'M' 'A-T'};
sublabel = {'G1S' 'M' 'A'};
T = [0 0.5 1.5 2.5 3.5 4.5]; % minute of 650nm on
tol = 0.8; % tolerence for bad fit

%%
n_sample = length(folder);
n_stage = length(subfolder);
linecolor = ['k','b','m'];
samples = cell(n_sample,n_stage);
for i = 1 : n_sample
    f0 = figure; hold on
    f1 = figure; hold on
    f2 = figure; hold on
    f3 = figure; hold on
    f4 = figure; hold on
    for j = 1 : n_stage
        A = exist([folder{i} '/' subfolder{j}],'dir');
        if A == 0
            continue
        end
        load([folder{i} '/' subfolder{j} '/analyzed matlab/intensity.mat'], 'summary','fit_rho', 'fit_int', 'fit_slope');
        samples{i,j} = struct(...
                        'name', folder{i},...
                        'stage', subfolder{j},...
                        'label', label{i},...
                        'data', summary,...
                        'fit_rho', fit_rho,...
                        'fit_int', fit_int,...
                        'fit_slope', fit_slope);
                    
        int_PIF = [summary.int_PIF];
        int_PIF_anchor = cellfun(@median,int_PIF(:,2:3:end));
        int_PIF_cyto = cellfun(@median,int_PIF(:,3:3:end));
        slope = [summary.slope];
        slope = slope(:, 1:2:end);
        rho = [summary.rho];
        rho = rho(:, 1:2:end);
        %calculate half time (rate) and filter out bad fits
        fit_rho = fit_rho(:,fit_rho(4,:) > tol);
        fit_int = fit_int(:,fit_rho(4,:) > tol);
        fit_slope = fit_slope(:,fit_rho(4,:) > tol);
        HT_rho = -reallog(2)./fit_rho(2,:) * 60; % convert to s
        HT_int = -reallog(2)./fit_int(2,:) * 60;
        HT_slope = -reallog(2)./fit_slope(2,:) * 60;
        
        figure(f0)
        plot(T, int_PIF_anchor./int_PIF_cyto - 1, [linecolor(j) ':'])
        plot(T, nanmean(int_PIF_anchor./int_PIF_cyto - 1,2),[linecolor(j) 'o-'])
        
        figure(f1)
        plot(T, slope, [linecolor(j) ':']);
        plot(T, nanmean(slope,2),[linecolor(j) 'o-'])      
        
        figure(f2)
        plot(T, rho, [linecolor(j) ':']);
        plot(T, nanmean(rho,2),[linecolor(j) 'o-'])      
        
        figure(f3)
        scatter(HT_slope, HT_int, [], linecolor(j))
        scatter(HT_rho, HT_int, [], linecolor(j),'filled','MarkerFaceAlpha', .3)
        
        figure(f4)
        scatter(fit_slope(3,:), fit_int(3,:), [], linecolor(j))
        scatter(fit_rho(3,:), fit_int(3,:), [], linecolor(j),'filled','MarkerFaceAlpha', .3)
    end
    
    figure(f0)    
    axis([0 5 -0.1 2.5])
    ax = gca;
    set(ax,'XTick',[0 0.5 1.5 2.5 3.5 4.5])
    set(ax,'YTick',(0:0.5:2.5))
    ylabel('Relative enrichment of PIF in anchored region')
    xlabel('Time of 650 on (min)')
    saveas(gcf, [folder{i} '/compare_cell_cycle_enrichment'],'fig')
    saveas(gcf, [folder{i} '/compare_cell_cycle_enrichment'],'png')
    close
    
    figure(f1)    
    axis([0 5 -0.1 1.8])
    ax = gca;
    set(ax,'XTick',[0 0.5 1.5 2.5 3.5 4.5])
    set(ax,'YTick',(0:0.2:2))
    ylabel('Slope of PIF relative to PhyB')
    xlabel('Time of 650 on (min)')
    saveas(gcf, [folder{i} '/compare_cell_cycle_slope'],'fig')
    saveas(gcf, [folder{i} '/compare_cell_cycle_slope'],'png')
    close

    figure(f2)    
    axis([0 5 -0.1 1.1])
    ax = gca;
    set(ax,'XTick',[0 0.5 1.5 2.5 3.5 4.5])
    set(ax,'YTick',(0:0.2:2))
    ylabel('Colocalization/Correlation bewteen PIF and PhyB')
    xlabel('Time of 650 on (min)')
    saveas(gcf, [folder{i} '/compare_cell_cycle_rho'],'fig')
    saveas(gcf, [folder{i} '/compare_cell_cycle_rho'],'png')
    close
    
    figure(f3)    
    axis([0 100 0 100])
    xlabel('Half time of correlation/slope (s)')
    ylabel('Half time of median intensity (s)')
    legend({'Slope' 'Colocalization/correlation coefficient'})
    saveas(gcf, [folder{i} '/compare_cell_cycle_HT'],'fig')
    saveas(gcf, [folder{i} '/compare_cell_cycle_HT'],'png')
    close
    
    figure(f4)    
    xlabel('Range of change for correlation/slope (s)')
    ylabel('Range of change for median intensity (s)')
    legend({'Slope' 'Colocalization/correlation coefficient'})
    saveas(gcf, [folder{i} '/compare_cell_cycle_range'],'fig')
    saveas(gcf, [folder{i} '/compare_cell_cycle_range'],'png')
    close
end

mkdir('summary')

%%
save('summary/intensity')



