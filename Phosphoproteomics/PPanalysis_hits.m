%%
% find hits

clear
close all
dirname = '/Users/snow/AmonLab/Projects/MEN/Proteomics/Kinase targets/Phospho-proteomics/Second run_201906_Cdc15-Cdc5_aF release/Results/analysis_20200417';
cd(dirname)
load('analysis_20200417')

%% find cutoff for ratio
% use the two controls to find the cutoff that allows 5% or 1% FDR 
% also calculate the stats for all, and non-phosphorylated peptides
Ratio_control = [CDC15(:,1)./CDC5(:,1) CDC5(:,1)./CDC15(:,1)]; % Cdc15/Cdc5; then Cdc5/Cdc15
[~,p_control] = ttest2(CDC15_raw, CDC5_raw,'Dim', 2);
p_control_adj = mafdr(p_control,'BHFDR',true);

stats_ratio_control = [nanmean(log2(Ratio_control)); nanmedian(log2(Ratio_control)); nanstd(log2(Ratio_control)); prctile(Ratio_control,95); prctile(Ratio_control,99)];

% number of peptides (false positives) in the set ratio cutoff: R = 2, 3, 4
n_ratio_control = [sum(Ratio_control>2&p_control_adj<0.05); sum(Ratio_control>3&p_control_adj<0.05); sum(Ratio_control>4&p_control_adj<0.05)];
n_ratio = [sum(Ratio_Cdc15>2&p_Cdc15_adj<0.05) sum(Ratio_Cdc5>2&p_Cdc5_adj<0.05); ...
            sum(Ratio_Cdc15>3&p_Cdc15_adj<0.05) sum(Ratio_Cdc5>3&p_Cdc5_adj<0.05);  ...
            sum(Ratio_Cdc15>4&p_Cdc15_adj<0.05) sum(Ratio_Cdc5>4&p_Cdc5_adj<0.05)];

%% find cutoff for missing data
% use the controls to find 5% FDR cutoff
n_replicates = 7;
n_replicates_control = zeros(n_replicates,2);
n_replicates_Cdc15 = zeros(n_replicates,2);
n_replicates_Cdc5 = zeros(n_replicates,2);
for i = 1 : n_replicates
    n_replicates_control(i,:) = [sum(CDC15(CDC5(:,2)==0,2)==i) sum(CDC5(CDC15(:,2)==0,2)==i)];
    n_replicates_Cdc15(i,:) = [sum(CDC15(cdc15as1(:,2)==0,2)==i) sum(cdc15as1(CDC15(:,2)==0,2)==i)];
    n_replicates_Cdc5(i,:) = [sum(CDC5(cdc5as1(:,2)==0,2)==i) sum(cdc5as1(CDC5(:,2)==0,2)==i)];
    
end
% number of peptides (false positives) if set cutoff to >= 3, 4, 5
d_replicates = [sum(n_replicates_control) sum(n_replicates_Cdc15) sum(n_replicates_Cdc5);...
                sum(n_replicates_control(3:end,:)) sum(n_replicates_Cdc15(3:end,:)) sum(n_replicates_Cdc5(3:end,:));...
                sum(n_replicates_control(4:end,:)) sum(n_replicates_Cdc15(4:end,:)) sum(n_replicates_Cdc5(4:end,:));...
                sum(n_replicates_control(5:end,:)) sum(n_replicates_Cdc15(5:end,:)) sum(n_replicates_Cdc5(5:end,:))];

close all
figure('position',[10 10 600 750])
subplot(3,2,1)
histogram(CDC15(cdc15as1(:,2)==0,2))
xlabel('# of replicates detected in wt but missing in as')
ylabel('# of peptides')
title('Cdc15')
subplot(3,2,2)
histogram(CDC5(cdc5as1(:,2)==0,2))
xlabel('# of replicates detected in wt but missing in as')
ylabel('# of peptides')
title('Cdc5')
subplot(3,2,3)
histogram(cdc15as1(CDC15(:,2)==0,2))
xlabel('# of replicates detected in as but missing in wt')
ylabel('# of peptides')
subplot(3,2,4)
histogram(cdc5as1(CDC5(:,2)==0,2))
xlabel('# of replicates detected in as but missing in wt')
ylabel('# of peptides')
subplot(3,2,5)
histogram(CDC15(CDC5(:,2)==0,2))
xlabel('# of replicates detected in CDC15 but missing in CDC5')
ylabel('# of peptides')
subplot(3,2,6)
histogram(CDC5(CDC15(:,2)==0,2))
xlabel('# of replicates detected in CDC5 but missing in CDC15')
ylabel('# of peptides')
saveas(gcf, 'Distribution of missing data','png')

figure('position',[10 10 600 500])
subplot(2,2,1)
histogram(CDC15(:,2))
xlabel('# of replicates detected in wt')
ylabel('# of peptides')
title('Cdc15')
subplot(2,2,2)
histogram(CDC5(:,2))
xlabel('# of replicates detected in wt')
ylabel('# of peptides')
title('Cdc5')
subplot(2,2,3)
histogram(cdc15as1(:,2))
xlabel('# of replicates detected in as')
ylabel('# of peptides')
subplot(2,2,4)
histogram(cdc5as1(:,2))
xlabel('# of replicates detected in as')
ylabel('# of peptides')
saveas(gcf, 'Distribution of replicates','png')
close all
%% find hits
thresh_ratio = 2; % 5% FDR
thresh_p = 0.05; % for adjusted p-value, use BHFDR 5% FDR
cutoff = 4; % >4 out of 7 (at least 5)

M_ppSites_1 = M_ppSites;

% add additional criteria (more than 1 replicate in wt)
hits_Cdc15 = (M_ppSites_1.Ratio_Cdc15 > thresh_ratio) & (M_ppSites_1.p_Cdc15_adj < thresh_p) & (M_ppSites_1.CDC15(:,2)>1);
hits_Cdc5 = (M_ppSites_1.Ratio_Cdc5 > thresh_ratio) & (M_ppSites_1.p_Cdc5_adj < thresh_p) & (M_ppSites_1.CDC5(:,2)>1);

% for no ratios due to missing data
hits_Cdc15(:,2) = isnan(M_ppSites_1.Ratio_Cdc15) & M_ppSites_1.CDC15(:,2)>cutoff;
hits_Cdc5(:,2) = isnan(M_ppSites_1.Ratio_Cdc5) & M_ppSites_1.CDC5(:,2)>cutoff;

% for bad p-values due to missing data
hits_Cdc15(:,3) = M_ppSites_1.cdc15as1(:,2)==1 & M_ppSites_1.CDC15(:,2)>cutoff & M_ppSites_1.Ratio_Cdc15 > thresh_ratio;
hits_Cdc5(:,3) = M_ppSites_1.cdc5as1(:,2)==1 & M_ppSites_1.CDC5(:,2)>cutoff & M_ppSites_1.Ratio_Cdc5 > thresh_ratio;

M_ppSites_1 = addvars(M_ppSites_1, hits_Cdc15,'after','dN_Cdc15');
M_ppSites_1 = addvars(M_ppSites_1, hits_Cdc5,'after','dN_Cdc5');

%% write all the pp sites
% M_ppSites = addvars(M_ppSites, sum(hits_Cdc15,2)>0,'after','dN_Cdc15');
M_ppSites.Properties.VariableNames{13} = 'Cdc15_dependent';
M_ppSites = addvars(M_ppSites, sum(hits_Cdc5,2)>0,'after','dN_Cdc5');
M_ppSites.Properties.VariableNames{end} = 'Cdc5_dependent';

[~, Site] = parse_ppSites(M_ppSites.PP_peptides);
M_ppSites = addvars(M_ppSites, Site,'after','Gene');
M_ppSites = removevars(M_ppSites,'PP_peptides');
M_ppSites = sortrows(M_ppSites,[2 3]);

writetable(M_ppSites, 'PPsites_all.xlsx');

%% find hits of pp_sites with at least one positive peptide
PP_sites_Cdc15 = unique(M_ppSites_1.PP_peptides(sum(hits_Cdc15,2)>0));
PP_sites_Cdc5 = unique(M_ppSites_1.PP_peptides(sum(hits_Cdc5,2)>0));
PP_sites_Cdc15n5 = intersect(PP_sites_Cdc15,PP_sites_Cdc5);

%consolidate all the peptides for each positve pp sites
[PP_sites_detected, ia2, ic2] = unique(M_ppSites_1.PP_peptides,'stable');
[~,idx_Cdc15,~] = intersect(PP_sites_detected, PP_sites_Cdc15);
[~,idx_Cdc5,~] = intersect(PP_sites_detected, PP_sites_Cdc5);

M_ppSites_Cdc15 = M_ppSites_1(ismember(ic2,idx_Cdc15),:);
M_ppSites_Cdc5 = M_ppSites_1(ismember(ic2,idx_Cdc5),:);

writetable(M_ppSites_Cdc15, 'PPsites_hits_all.xlsx', 'Sheet',1);
writetable(M_ppSites_Cdc5, 'PPsites_hits_all.xlsx', 'Sheet',2);

%% find hits of pp_sites taking all peptides into account
% single pp site peptides should be weighted more since additional site(s) might bias the outcome especially if they goes the other way
% strict rule: all unique site peptides need to be positive hits for a pp site
% and if there is only duel/multi sites peptides, if other sites are hits, call it a maybe
[~, ~, idx2] = unique(M_ppSites_Cdc15.PP_peptides,'stable');
idx_Cdc15 = ones(length(PP_sites_Cdc15),1);
n_pp_peptides_Cdc15 = cellfun('length', M_ppSites_Cdc15.Phosphosites);
thresh_ratio2 = 1.5;
for i = 1 : length(PP_sites_Cdc15)
    % check for single site peptides, all peptides need to satisfy ratio cutoff, loose the p-value threshold in case of detection issue    
    if sum(M_ppSites_Cdc15.Ratio_Cdc15(idx2==i&n_pp_peptides_Cdc15<2,:)<=thresh_ratio2) > 0
        idx_Cdc15(i) = 0;
        % find sites that didn't fit above but has at least one positive peptide to preserve more hits
        if sum(sum(M_ppSites_Cdc15.hits_Cdc15(idx2==i&n_pp_peptides_Cdc15<2,:),2)>0) > 0
            idx_Cdc15(i) = 0.3;
        end
    end
    
    % if no single site peptide, count as weak hits
    if sum(idx2==i&n_pp_peptides_Cdc15<2)<1
        idx_Cdc15(i) = 0.5;
    end
end
for i = 1 : length(PP_sites_Cdc15)
    % if no single site peptide, count as weak hits
    if sum(idx2==i&n_pp_peptides_Cdc15<2)<1
        % if the additional pp sites is a strong hit, makes it weaker
        idx = strfind(PP_sites_Cdc15{i}, '_')+1;
        target_pp = str2double(PP_sites_Cdc15{i}(idx:end));
        co_pp = unique([M_ppSites_Cdc15.Phosphosites{idx2==i}]);
        co_pp(co_pp==target_pp) = [];
        for j = 1 : length(co_pp)
            if ismember([PP_sites_Cdc15{i}(1:idx-1) num2str(co_pp(j))], PP_sites_Cdc15(idx_Cdc15>0.5))
                idx_Cdc15(i) = 0.2;
            end
        end
    end
end
PP_sites_Cdc15_filter = struct('PP_sites', PP_sites_Cdc15(idx_Cdc15>0),...
                                'Score', num2cell(idx_Cdc15(idx_Cdc15>0)));
[~,idx_Cdc15_filter,~] = intersect(PP_sites_detected, PP_sites_Cdc15(idx_Cdc15>0));
M_ppSites_Cdc15_filter = M_ppSites_1(ismember(ic2,idx_Cdc15_filter),:);

[~, ~, idx2] = unique(M_ppSites_Cdc5.PP_peptides,'stable');
idx_Cdc5 = ones(length(PP_sites_Cdc5),1);
n_pp_peptides_Cdc5 = cellfun('length', M_ppSites_Cdc5.Phosphosites);
for i = 1 : length(PP_sites_Cdc5)
    % check for single site peptide, at least half of the peptides are hits     
    if sum(M_ppSites_Cdc5.Ratio_Cdc5(idx2==i&n_pp_peptides_Cdc5<2,:)<=thresh_ratio2) > 0
        idx_Cdc5(i) = 0;
        % find sites that didn't fit above but has at least one positive peptide to preserve more hits
        if sum(sum(M_ppSites_Cdc5.hits_Cdc5(idx2==i&n_pp_peptides_Cdc5<2,:),2)>0) > 0
            idx_Cdc5(i) = 0.3;
        end
    end
    % if no single site peptide, count as weak hits
    if sum(idx2==i&n_pp_peptides_Cdc5<2)<1
        idx_Cdc5(i) = 0.5;
    end
end
for i = 1 : length(PP_sites_Cdc5)
    if sum(idx2==i&n_pp_peptides_Cdc5<2)<1
        % if the additional pp sites is a strong hit, makes it weaker
        idx = strfind(PP_sites_Cdc5{i}, '_')+1;
        target_pp = str2double(PP_sites_Cdc5{i}(idx:end));
        co_pp = unique([M_ppSites_Cdc5.Phosphosites{idx2==i}]);
        co_pp(co_pp==target_pp) = [];
        for j = 1 : length(co_pp)
            if ismember([PP_sites_Cdc5{i}(1:idx-1) num2str(co_pp(j))], PP_sites_Cdc5(idx_Cdc5>0.5))
                idx_Cdc5(i) = 0.2;
            end
        end
    end
end

PP_sites_Cdc5_filter = struct(...
                        'PP_sites', PP_sites_Cdc5(idx_Cdc5>0),...
                        'Score', num2cell(idx_Cdc5(idx_Cdc5>0)));
[~,idx_Cdc5_filter,~] = intersect(PP_sites_detected, PP_sites_Cdc5(idx_Cdc5>0));
M_ppSites_Cdc5_filter = M_ppSites_1(ismember(ic2,idx_Cdc5_filter),:);
PP_sites_Cdc15n5_filter = intersect(PP_sites_Cdc5(idx_Cdc5>0), PP_sites_Cdc15(idx_Cdc15>0));

writetable(M_ppSites_Cdc15_filter, 'PPsites_hits_filtered.xlsx', 'Sheet',1);
writetable(M_ppSites_Cdc5_filter, 'PPsites_hits_filtered.xlsx', 'Sheet',2);

%% List all the sites with scores
H_ppSites_Cdc15 = struct2table(PP_sites_Cdc15_filter);
[ProteinName, Site] = parse_ppSites(H_ppSites_Cdc15.PP_sites);
H_ppSites_Cdc15 = addvars(H_ppSites_Cdc15, ProteinName,'after',1);
H_ppSites_Cdc15 = addvars(H_ppSites_Cdc15, Site,'after','ProteinName');
H_ppSites_Cdc15 = removevars(H_ppSites_Cdc15,'PP_sites');
H_ppSites_Cdc15 = sortrows(H_ppSites_Cdc15,[1 2]);
[Protein_Cdc15, ~, ic_Cdc15] = unique(H_ppSites_Cdc15.ProteinName,'stable');

H_ppSites_Cdc5 = struct2table(PP_sites_Cdc5_filter);
[ProteinName, Site] = parse_ppSites(H_ppSites_Cdc5.PP_sites);
H_ppSites_Cdc5 = addvars(H_ppSites_Cdc5, ProteinName,'after',1);
H_ppSites_Cdc5 = addvars(H_ppSites_Cdc5, Site,'after','ProteinName');
H_ppSites_Cdc5 = removevars(H_ppSites_Cdc5,'PP_sites');
H_ppSites_Cdc5 = sortrows(H_ppSites_Cdc5,[1 2]);
[Protein_Cdc5, ~, ic_Cdc5] = unique(H_ppSites_Cdc5.ProteinName,'stable');
%% Summarize the hit sites for proteins
PG_summary_Cdc15 = struct('ProteinName', Protein_Cdc15);
Sites = cell(length(PG_summary_Cdc15),2);
for p = 1 : length(PG_summary_Cdc15)
    [~,idx,~] = intersect({Proteome.GN}, Protein_Cdc15{p});
    Sites{p,1} = convert_pp_sites(H_ppSites_Cdc15.Site(ic_Cdc15==p), Proteome(idx).Sequence);
    Sites{p,2} = convert_pp_sites(H_ppSites_Cdc15.Site(ic_Cdc15==p&H_ppSites_Cdc15.Score>=0.5), Proteome(idx).Sequence);
end
[PG_summary_Cdc15.PP_Cdc15_all] = Sites{:,1};
[PG_summary_Cdc15.PP_Cdc15_strong] = Sites{:,2};

PG_summary_Cdc5 = struct('ProteinName', Protein_Cdc5);
Sites = cell(length(PG_summary_Cdc5),2);
for p = 1 : length(PG_summary_Cdc5)
    [~,idx,~] = intersect({Proteome.GN}, Protein_Cdc5{p});
    Sites{p,1} = convert_pp_sites(H_ppSites_Cdc5.Site(ic_Cdc5==p), Proteome(idx).Sequence);
    Sites{p,2} = convert_pp_sites(H_ppSites_Cdc5.Site(ic_Cdc5==p&H_ppSites_Cdc5.Score>=0.5), Proteome(idx).Sequence);
end
[PG_summary_Cdc5.PP_Cdc5_all] = Sites{:,1};
[PG_summary_Cdc5.PP_Cdc5_strong] = Sites{:,2};

%% map hits to total phospho-protein
PG_summary_hits = sortrows(struct2table(PG_summary_pp),'ProteinName');
PG_summary_Cdc15 = sortrows(struct2table(PG_summary_Cdc15),'ProteinName');
PG_summary_Cdc5 = sortrows(struct2table(PG_summary_Cdc5),'ProteinName');
[~,pp_idx_Cdc15,~] = intersect(PG_summary_hits.ProteinName,PG_summary_Cdc15.ProteinName, 'stable');
[~,pp_idx_Cdc5,~] = intersect(PG_summary_hits.ProteinName,PG_summary_Cdc5.ProteinName, 'stable');
PP_Cdc15 = cell(length(PG_summary_pp),2);
PP_Cdc5 = cell(length(PG_summary_pp),2);
PP_Cdc15(pp_idx_Cdc15,1) = PG_summary_Cdc15.PP_Cdc15_all;
PP_Cdc15(pp_idx_Cdc15,2) = PG_summary_Cdc15.PP_Cdc15_strong;
PP_Cdc5(pp_idx_Cdc5,1) = PG_summary_Cdc5.PP_Cdc5_all;
PP_Cdc5(pp_idx_Cdc5,2) = PG_summary_Cdc5.PP_Cdc5_strong;
PG_summary_hits = addvars(PG_summary_hits,PP_Cdc15);
PG_summary_hits = addvars(PG_summary_hits,PP_Cdc5);

%%
PP_Cdc15n5 = cell(length(PG_summary_pp),1);
PP_Cdc5only = cell(length(PG_summary_pp),1);
for i = 1 : length(PG_summary_pp)
    if ~isempty(PP_Cdc5{i,1})
        if ~isempty(PP_Cdc15{i,1})
            PP_Cdc15n5{i} = char(join(intersect(split(PP_Cdc15{i,1},','), split(PP_Cdc5{i,1},','),'stable'),','));
            PP_Cdc5only{i} = char(join(setdiff(split(PP_Cdc5{i,1},','), split(PP_Cdc15{i,1},','),'stable'),','));
        else
            PP_Cdc5only{i} = PP_Cdc5{i,1};
        end
    end
end
PG_summary_hits = addvars(PG_summary_hits,PP_Cdc15n5);
PG_summary_hits = addvars(PG_summary_hits,PP_Cdc5only);

writetable(PG_summary_hits, 'PG_summary_hits.xlsx');

%% parse pp_sites
function [ProteinName, site] = parse_ppSites(PP_sites)
    ProteinName = cell(length(PP_sites),1);
    site = zeros(length(PP_sites),1);
    for i = 1 : length(PP_sites)
        temp = split(PP_sites{i},'_');
        ProteinName{i} = temp{1};
        site(i) = str2double(temp{2});
    end
end

%% convert pp sites from number to string (group)
function pp_sites = convert_pp_sites(pp, Seq)
    if ~isempty(pp)
        temp = cell(length(pp),1);
        for i = 1 : length(pp)
            temp{i} = [Seq(pp(i)) num2str(pp(i))];
        end
        pp_sites = char(join(temp,','));
    else
        pp_sites = [];
    end
end