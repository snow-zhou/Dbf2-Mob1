clear
close all

addpath(genpath('/Microscopy'))

dirname = '';
Ch_PIF = 2;
Ch_PhyB = 3;
nCh = 3;
nframe = 6;
bg_PIF = 150;
bg_PhyB = 120;
max_PhyB = 350; % limit to filter out bright background spot
T = [0 0.5 1.5 2.5 3.5 4.5]; % minute of 650nm on
%%
cd(dirname)
list_file = dir('*tif');
mkdir('analyzed matlab')
%% extract the intensity 
nfile = length(list_file);
intensity = cell(nfile,1);
f = figure;
set(f,'visible','off');
for n = 1 : nfile
    filename = list_file(n).name(1:end-4);
    mkdir([filename '/segmentation'])
    data = bfopen(list_file(n).name);
    rho = zeros(nframe,2);
    slope = zeros(nframe,2);
    int_PIF = cell(nframe,3);
    int_PhyB = cell(nframe,2);

    % use the first frame PIF image as mask
    I_PIF = data{1}{Ch_PIF};
    mask = imopen(imbinarize(mat2gray(I_PIF)),strel('square',5));
    c = bwboundaries(mask);
    imshowpair(I_PIF, mask, 'montage')
	saveas(f,[filename '/segmentation/cell_mask'],'png')
    for i = 1 : nframe
        I_PhyB = data{1}{nCh*(i-1)+Ch_PhyB};
        I_PIF = data{1}{nCh*(i-1)+Ch_PIF};        
        I_PIF2 = I_PIF - bg_PIF;
        I_PhyB2 = I_PhyB - bg_PhyB;        

        % segment PhyB anchor
%         bw_PhyB = mask&loccssegm(mat2gray(I_PhyB2),8);
        I_PhyB2(I_PhyB2 > max_PhyB) = 0;
        I_PhyB2(~mask) = 0;
        bw_PhyB = imbinarize(mat2gray(I_PhyB2),0.5);
%         bw_PhyB = imerode(bw_PhyB, strel('square',2));
        bw_PhyB = bwareaopen(bw_PhyB, 5);
        imshowpair(I_PhyB2, bw_PhyB, 'montage')
        saveas(f,[filename '/segmentation/PhyB_frame-' num2str(i)],'png')
        
        % show segmentation
        g = figure;
        set(g,'visible','off');
        imshowpair(I_PIF2, I_PhyB2, 'montage')
        hold on
        plot(smooth(c{1}(:,2)),smooth(c{1}(:,1)),'y','linewidth',1);
        [~,width] = size(I_PIF); 
        plot(smooth(c{1}(:,2))+width,smooth(c{1}(:,1)),'y','linewidth',1);
        d = bwboundaries(imdilate(bw_PhyB,strel('disk',1)));
        nobj = length(d);
        for j = 1 : nobj
            plot(smooth(d{j}(:,2)),smooth(d{j}(:,1)),'m','linewidth',1);
            plot(smooth(d{j}(:,2))+width,smooth(d{j}(:,1)),'m','linewidth',1);
        end
        saveas(g,[filename '/segmentation/seg_frame-' num2str(i)],'png')
        close
        
        % intensity
        int_PIF{i,1} = double(I_PIF2(mask)); % all pixel for PIF
        int_PhyB{i,1} = double(I_PhyB2(mask)); % all pixel for PhyB
        int_PIF{i,2} = double(I_PIF2(bw_PhyB)); % PIF in anchored region
        int_PIF{i,3} = double(I_PIF2(mask & ~imdilate(bw_PhyB,strel('disk',1)))); % cytoplasmic 
        int_PhyB{i,2} = double(I_PhyB2(bw_PhyB)); % PhyB in anchored region
        
        [rho(i,1), rho(i,2)] = corr(int_PIF{i,1}, int_PhyB{i,1});
        slope(i,:) = polyfit(int_PhyB{i,1}, int_PIF{i,1}, 1);
        
        h = figure;
        set(h,'visible','off');
        scatter(int_PhyB{i,1},int_PIF{i,1},'filled', 'MarkerFaceAlpha',.2)
        hold on
        scatter(int_PhyB{i,2},int_PIF{i,2})
        plot(int_PhyB{i,1}, slope(i,1)*int_PhyB{i,1}+slope(i,2))
%         axis([0 120 0 80])
        xlabel('PhyB intensity in cell')
        ylabel('PIF intensity in cell')
        saveas(h,[filename '/segmentation/intensity_frame-' num2str(i)],'png')
        close
        
    end
    
    % fit rate
    timestamp = T';
    %fit rate of increase for correlation coefficient
    [par_rho,fit_rho] = fit_rate([timestamp rho(:,1)]);
    %fit rate of increase for intensity (relative enrichment in anchored region)
    if isnan(mean(cellfun(@median,int_PIF(:,2))./cellfun(@median,int_PIF(:,3))))
        par_int = NaN;
    else
        [par_int,fit_int] = fit_rate([timestamp cellfun(@median,int_PIF(:,2))./cellfun(@median,int_PIF(:,3))]);
    end
    %fit rate of increase for slope of correlation
    [par_slope,fit_slope] = fit_rate([timestamp slope(:,1)]);
    
    % summary
    h2= figure;
    set(h2,'visible','off');
    scatter(timestamp, rho(:,1))
    hold on
    plot(fit_rho(:,1), fit_rho(:,2))
    xlabel('Time (min)')
    ylabel('colocalization/correlation coefficient')
    saveas(gcf,[filename '/colocalization'],'png')
    clf
    
    plot(timestamp, cellfun(@median,int_PIF(:,3)),'o-')
    hold on
    plot(timestamp, cellfun(@median,int_PIF(:,2)),'o-')
    xlabel('Time (min)')
    ylabel('Intensity')
    yyaxis right
    scatter(timestamp, cellfun(@median,int_PIF(:,2))./cellfun(@median,int_PIF(:,3)),'LineWidth',1.5)
    plot(fit_int(:,1), fit_int(:,2), 'LineWidth',2)
    ylabel('Relative enrichment in anchored region')
    legend({'Median intensity in cytosol' 'Median intensity in anchored region' 'Relative enrichment in anchored region' 'Exponential Fit'}, 'location', 'best')
    saveas(gcf,[filename '/median intensity'],'png')
    clf
    
    plot(timestamp, cellfun(@sum,int_PIF(:,1)))
    hold on
    plot(timestamp, cellfun(@sum,int_PIF(:,2)))
    xlabel('Time (min)')
    ylabel('Total Intensity')
    yyaxis right
    plot(timestamp, cellfun(@sum,int_PIF(:,2))./cellfun(@sum,int_PIF(:,1)))
    ylabel('Fraction of signal in anchored region')
    legend({'Total intensity in cell' 'Total intensity in anchored region' 'Fraction in anchored region'}, 'location', 'best')
    saveas(gcf,[filename '/sum intensity'],'png')
    clf
    
    scatter(timestamp, slope(:,1))
    hold on
    plot(fit_slope(:,1), fit_slope(:,2))
    xlabel('Time (min)')
    ylabel('Slope of correlation between PhyB and PIF')
    saveas(gcf,[filename '/slope'],'png')
    close
    
    intensity{n} = struct(...
                    'file', filename,...
                    'timestamp', timestamp,...
                    'int_PIF', {int_PIF},...
                    'int_PhyB', {int_PhyB},...
                    'rho', rho,...
                    'slope', slope,...
                    'fit_rho', par_rho,...
                    'fit_enrichment', par_int,...
                    'fit_slope', par_slope);
end
    
%% summary
close all

summary = cell2mat(intensity);

% plot all the curves together
fit_rho = [summary.fit_rho];
fit_int = [summary.fit_enrichment];
fit_slope = [summary.fit_slope];
sampling = linspace(timestamp(1), timestamp(end), 30);
linecolor = colormap;
f1 = figure; hold on
f2 = figure; hold on
f3 = figure; hold on
for n = 1 : nfile
    figure(f1)
    scatter(timestamp, summary(n).rho(:,1),[],'MarkerEdgeColor',linecolor(n,:));
    plot(sampling, fit_rho(1,n)*(1-exp(fit_rho(2,n)*(sampling)))+fit_rho(3,n), 'Color', linecolor(n,:));
    figure(f2)
    scatter(timestamp, cellfun(@median,summary(n).int_PIF(:,2))./cellfun(@median,summary(n).int_PIF(:,3)),[],'MarkerEdgeColor',linecolor(n,:));
    plot(sampling, fit_int(1,n)*(1-exp(fit_int(2,n)*(sampling)))+fit_int(3,n), 'Color', linecolor(n,:));
    figure(f3)
    scatter(timestamp, summary(n).slope(:,1),[],'MarkerEdgeColor',linecolor(n,:));
    plot(sampling, fit_slope(1,n)*(1-exp(fit_slope(2,n)*(sampling)))+fit_slope(3,n), 'Color', linecolor(n,:));
end
figure(f1)
xlabel('Time (min)')
ylabel('colocalization/correlation coefficient')
saveas(gcf,'analyzed matlab/colocalization','png')

figure(f2)
xlabel('Time (min)')
ylabel('Relative enrichment in anchored region')
saveas(gcf,'analyzed matlab/enrichment','png')

figure(f3)
xlabel('Time (min)')
ylabel('slope of correlation between PhyB and PIF')
saveas(gcf,'analyzed matlab/slope','png')

% plot the rate as half time
%filter out bad fits (rsquare >0.8)
HT_rho = -reallog(2)./fit_rho(2,fit_rho(4,:)>0.8) * 60; % convert to s
HT_int = -reallog(2)./fit_int(2,fit_rho(4,:)>0.8) * 60;
HT_slope = -reallog(2)./fit_slope(2,fit_rho(4,:)>0.8) * 60;
figure; hold on
scatter(HT_rho, HT_int)
scatter(HT_slope, HT_int)
xlabel('Half time of recruitment by correlation/slope (s)')
ylabel('Half time of recruitment by intensity enrichment (s)')
legend({'Colocalization/correlation coefficient' 'Slope'}, 'location', 'best')
saveas(gcf,'analyzed matlab/rate','png')

% plot the range of change relative to initial
figure
hold on
scatter(fit_rho(3,fit_rho(4,:)>0.8), fit_rho(1,fit_rho(4,:)>0.8))
scatter(fit_int(3,fit_rho(4,:)>0.8)-1, fit_int(1,fit_rho(4,:)>0.8))
scatter(fit_slope(3,fit_rho(4,:)>0.8), fit_slope(1,fit_rho(4,:)>0.8))
xlabel('Intial value')
ylabel('Range of Change')
legend({'Colocalization/correlation coefficient' 'Relative enrichment' 'Slope'}, 'location', 'best')
saveas(gcf,'analyzed matlab/range2intial','png')

% plot the range of change relative to PhyB intensity
int_PhyB = [summary.int_PhyB];
int_PhyB_ini = cellfun(@median, int_PhyB(1,2:2:end));
figure
hold on
scatter(int_PhyB_ini(fit_rho(4,:)>0.8), fit_rho(1,fit_rho(4,:)>0.8))
scatter(int_PhyB_ini(fit_rho(4,:)>0.8), fit_int(1,fit_rho(4,:)>0.8))
scatter(int_PhyB_ini(fit_rho(4,:)>0.8), fit_slope(1,fit_rho(4,:)>0.8))
xlabel('Median PhyB intensity')
ylabel('Range of Change')
legend({'Colocalization/correlation coefficient' 'Relative enrichment' 'Slope'}, 'location', 'best')
saveas(gcf,'analyzed matlab/range2PhyB','png')


%%
save('analyzed matlab/intensity')

%%
function [parameters,fit_result] = fit_rate(traces)
%assume the frist column of traces is time stamp
%fit as a*exp(b*x)+c, return parameters of a, b. 
    s = 30; % sampling points for the fit results
    [~,n] = size(traces);
    fit_result = zeros(s,n);
    fit_result(:,1) = linspace(traces(1,1),traces(end,1),s)';
    parameters = zeros(4,n-1);
    
    g = fittype('a*(1-exp(b*x))');
    
    for i = 2 : n
        %estimate parameter
        Ea = max(traces(:,i))-min(traces(:,i));
        Ec = traces(1,i);
        Eb = log(1-(traces(2,i) - Ec)/Ea)/(traces(2,1));
%         Eb = -0.02;
        [f, g]= fit(traces(:,1), traces(:,i) - Ec,g,'StartPoint',[Ea, Eb]);
        parameters(:,i-1) = [f.a; f.b; Ec; g.rsquare];
        fit_result(:,i) = f.a*(1-exp(f.b*(fit_result(:,1))))+Ec;
    end

end