%% analyze MS data
clear
close all
dirname = '/Users/snow/AmonLab/Projects/MEN/Proteomics/Proximity based/TurboID/MS/20190207_MEN proteins/';
file_Mob1 = '190215_626_SZ_953_edit.xlsx';
file_Tem1 = '190215_626_SZ_997_edit.xlsx';
file_Cdc15 = '190215_626_SZ_999_edit.xlsx';
file_Dbf2 = '190215_626_SZ_1001_edit.xlsx';
file_control = '190215_626_SZ_2588_edit.xlsx';

cd(dirname)
opts = detectImportOptions(file_control);
opts.Sheet = 'summary';
M_control = readtable(file_control,'Sheet', opts.Sheet);

opts = detectImportOptions(file_Mob1);
opts.Sheet = 'summary';
M_Mob1 = readtable(file_Mob1,'Sheet', opts.Sheet);

opts = detectImportOptions(file_Tem1);
opts.Sheet = 'summary';
M_Tem1 = readtable(file_Tem1,'Sheet', opts.Sheet);

opts = detectImportOptions(file_Cdc15);
opts.Sheet = 'summary';
M_Cdc15 = readtable(file_Cdc15,'Sheet', opts.Sheet);

opts = detectImportOptions(file_Dbf2);
opts.Sheet = 'summary';
M_Dbf2 = readtable(file_Dbf2,'Sheet', opts.Sheet);
mkdir('analysis')
cd('analysis')

M_control = removevars(M_control, 'SpectrumFile');
M_Mob1 = removevars(M_Mob1, 'SpectrumFile');
M_Tem1 = removevars(M_Tem1, 'SpectrumFile');
M_Cdc15 = removevars(M_Cdc15, 'SpectrumFile');
M_Dbf2 = removevars(M_Dbf2, 'SpectrumFile');
%% combine and compare
ProteinAccessions = unique([M_control.ProteinAccessions; M_Mob1.ProteinAccessions; M_Tem1.ProteinAccessions; M_Cdc15.ProteinAccessions; M_Dbf2.ProteinAccessions]);
T = table(ProteinAccessions);
T = outerjoin(T, M_control,'keys',1);
T = removevars(T, 'ProteinAccessions_M_control');
T.Properties.VariableNames{1} = T.Properties.VariableNames{1}(1:end-2);
T = AppendTable(T, M_Mob1, 1, [1 2], 'Mob1');
T = AppendTable(T, M_Tem1, 1, [1 2], 'Tem1');
T = AppendTable(T, M_Cdc15, 1, [1 2], 'Cdc15');
T = AppendTable(T, M_Dbf2, 1, [1 2], 'Dbf2');

%%
TotPep_control = T.TotPep;
TotPep = [T.TotPep_Mob1, T.TotPep_Tem1, T.TotPep_Cdc15, T.TotPep_Dbf2];

UniPep_control = T.UniPep;
UniPep = [T.UniPep_Mob1, T.UniPep_Tem1, T.UniPep_Cdc15, T.UniPep_Dbf2];

%% plot comparison
close all
TotPep_control(isnan(T.TotPep)) = 0.1;
TotPep(isnan(T.TotPep_Mob1),1) = 0;
TotPep(isnan(T.TotPep_Tem1),2) = 0;
TotPep(isnan(T.TotPep_Cdc15),3) = 0;
TotPep(isnan(T.TotPep_Dbf2),4) = 0;

UniPep_control(isnan(T.TotPep)) = 0.1;
UniPep(isnan(T.UniPep_Mob1),1) = 0;
UniPep(isnan(T.UniPep_Tem1),2) = 0;
UniPep(isnan(T.UniPep_Cdc15),3) = 0;
UniPep(isnan(T.UniPep_Dbf2),4) = 0;

f1 = figure('position', [0 0 800 600]); 
n_sample = 4;
label = {'Mob1' 'Tem1' 'Cdc15' 'Dbf2'};
for i = 1 : n_sample
    subplot(3,n_sample,i)
    scatter(TotPep_control, TotPep(:,i),'o')
    hold on
    scatter(UniPep_control, UniPep(:,i),'o')
    axis([0 1000 0 1000])
    xlabel('Peptides detected in Control')
    ylabel('Peptides detected in Sample')
    legend({'Total peptides' 'Unique peptides'},'location','northeast')
    title(label{i})
    subplot(3,n_sample,n_sample+i)
    loglog(TotPep_control, TotPep(:,i),'o')
    axis([0 1000 0 1000])
    xlabel('Peptides detected in Control')
    ylabel('Peptides detected in Sample')
    title('log scale')
    subplot(3,n_sample,2*n_sample+i)
    loglog(TotPep_control./UniPep_control, TotPep(:,i)./UniPep(:,i),'o')
%     axis([0 1000 0 1000])
    xlabel('Peptides detected in Control')
    ylabel('Peptides detected in Sample')
    title('normalized to UniPep')
end
saveas(gcf, 'Compare','png')
%% ratio
close all
thresh = 10; % # STD above mean ratio
TotPep_control(isnan(T.TotPep)) = 0.1;
ratio = TotPep./TotPep_control;
ratio_mean = [mean(ratio(~isnan(T.TotPep)&~isnan(T.TotPep_Mob1),1)),...
              mean(ratio(~isnan(T.TotPep)&~isnan(T.TotPep_Tem1),2)),...
              mean(ratio(~isnan(T.TotPep)&~isnan(T.TotPep_Cdc15),3)),...
              mean(ratio(~isnan(T.TotPep)&~isnan(T.TotPep_Dbf2),4))];
ratio_std = [std(ratio(~isnan(T.TotPep)&~isnan(T.TotPep_Mob1),1)),...
              std(ratio(~isnan(T.TotPep)&~isnan(T.TotPep_Tem1),2)),...
              std(ratio(~isnan(T.TotPep)&~isnan(T.TotPep_Cdc15),3)),...
              std(ratio(~isnan(T.TotPep)&~isnan(T.TotPep_Dbf2),4))];
figure('position', [0 0 800 250])
for i = 1 : n_sample
    subplot(1,n_sample,i)
    plot(sort(log2(ratio(:,i))),'o')
    hold on
    plot([0 1500],[1 1]*log2(ratio_mean(i)),'--','linewidth', 2)
    plot([0 1500],[1 1]*log2(ratio_mean(i)+ratio_std(i)*thresh),':','linewidth', 2)
    axis([0 1500 -5 10])
    title(label{i})
    xlabel('Peptides detected')
    ylabel('log2 of Ratio')
end

saveas(gcf, 'Ratio','png')

%% find hits
close all
hits_ratio = ratio > (ratio_mean+ratio_std*thresh);
filter = UniPep>1;
figure('position', [0 0 700 200])
for i = 1 : n_sample
    subplot(1,n_sample,i)
    loglog(TotPep_control, TotPep(:,i),'o')
    hold on
    scatter(TotPep_control(hits_ratio(:,i)&filter(:,i)), TotPep(hits_ratio(:,i)&filter(:,i),i),'filled','MarkerFaceAlpha',0.3)
%     loglog(TotPep_control(hits_ratio(:,i)), TotPep(hits_ratio(:,i),i),'o','MarkerFaceColor',[0.8500 0.3250 0.0980])
    plot([1 1000], [1 1000]*ratio_mean(i),'linewidth',1.5)
    plot([0.1 1000], [0.1 1000]*(ratio_mean(i)+thresh*ratio_std(i)),'--','linewidth',1)
    axis([0 1000 1 1000])
    title(label{i})
end
saveas(gcf, 'Hits','png')
saveas(gcf, 'Hits','pdf')
%% write to file
hits_sum = sum(hits_ratio,2);
T = T(:,1:12);
T = addvars(T, hits_sum);
ratio_Mob1 = ratio(:,1);
ratio_Tem1 = ratio(:,2);
ratio_Cdc15 = ratio(:,3);
ratio_Dbf2 = ratio(:,4);
T = addvars(T, ratio_Mob1);
T = addvars(T, ratio_Tem1);
T = addvars(T, ratio_Cdc15);
T = addvars(T, ratio_Dbf2);
T_hits = sortrows(T(hits_sum>0&sum(filter,2)>0,:), 'hits_sum', 'descend');
writetable(T_hits,'Compare_hits.xlsx');

%%
close all
save('analysis_20200330')
%%
function Tblo = AppendTable(A, B, key, keys_con, NameB)
    T = outerjoin(A, B,'keys', key);
    [~, n_var] = size(B);
    [~, n_var_A] = size(A);
    Tblo = A;
    for i = 1 : n_var
        if ismember(i,keys_con)
            Tblo(cellfun('isempty', T{:,i}),i) = T(cellfun('isempty', T{:,i}),i+n_var_A);
        else
            Tblo = [Tblo T(:,i+n_var_A)];
            Tblo.Properties.VariableNames{end} = [Tblo.Properties.VariableNames{end}(1:end-2) '_' NameB];            
        end
    end
end