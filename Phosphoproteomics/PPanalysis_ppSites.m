%%
% new analysis pipeline starting from scratch
% need to take into acount cases where the same pp site is detected in
% multiple peptides; change from peptide centric to pp site centric
% this script only parse the data (no hit finding)

clear
close all
dirname = '/Users/snow/AmonLab/Projects/MEN/Proteomics/Kinase targets/Phospho-proteomics/Second run_201906_Cdc15-Cdc5_aF release/Results/';
cd(dirname)

filename = 'analysis_20200331/20190717_130531_ALL_Snow_28files_directDIA_Report_nan.xlsx';
file_proteome = '/Users/snow/AmonLab/Projects/MEN/Proteomics/Kinase targets/Phospho-proteomics/uniprot-proteome_UP000002311.fasta';

M = readtable(filename);
Proteome = fastaread(file_proteome);

mkdir('analysis_20200417')
cd('analysis_20200417')

%% parse peptides
peptides_raw = M.EG_PrecursorId;
Peptides_clean = cell(length(peptides_raw),1);
Phosphosites_raw = cell(length(peptides_raw),1);
for i = 1 : length(peptides_raw)
    temp = peptides_raw{i}(2:end-3);
    
    idx_Modifications_start = find(temp=='[');
    idx_Modifications_end = find(temp==']');
    idx_Modifications = [];
    for j = 1 : length(idx_Modifications_start)
        idx_Modifications = [idx_Modifications idx_Modifications_start(j):1:idx_Modifications_end(j)];
    end
    
    idx_PP = strfind(temp, 'Phospho') - 2;
    
    temp(idx_Modifications) = [];
    Peptides_clean{i} = temp;
    if ~isempty(idx_PP)
        for p = 1 : length(idx_PP)
            Phosphosites_raw{i}(p) = idx_PP(p) - length(find(idx_Modifications < idx_PP(p)));
        end
    end
    
end

%% parse the proteome
ProteinAcessions = cell(length(Proteome),1);
GN = cell(length(Proteome),1);
for i = 1 : length(Proteome)
    idx = strfind(Proteome(i).Header, '|');
    ProteinAcessions{i} = Proteome(i).Header(idx(1)+1:idx(2)-1);
    idx1 = strfind(Proteome(i).Header, 'GN=')+3;
    idx2 = strfind(Proteome(i).Header, 'PE=')-2;
    GN{i} = Proteome(i).Header(idx1:idx2);
end
[Proteome.ProteinAcessions] = ProteinAcessions{:};
[Proteome.GN] = GN{:};
%% map peptides and phospho sites to sequence
% identify proteins that is detected 
[Protein_detected, ia, ic] = unique(M.PG_ProteinAccessions,'stable');
n_protein_detected = length(Protein_detected);
PG_coverage = cell(n_protein_detected,1);
PG_pp = cell(n_protein_detected,2);
PG_GN = cell(n_protein_detected,1);
% identify pp sites for each protein detected
Phosphosites = cell(length(peptides_raw),1);
for p = 1 : n_protein_detected
    PG_idx = find(strcmp(ProteinAcessions, Protein_detected{p})>0);
    if isempty(PG_idx)
        continue;
    end
    idx_peptides = find(ic==p);
    coverage = zeros(length(Proteome(PG_idx).Sequence),1);
    for i = idx_peptides(1) : idx_peptides(end)
        idx = strfind(Proteome(PG_idx).Sequence, Peptides_clean{i});
        if length(idx) > 1 % can't be mapped uniquely to one site due to repeating sequence
            continue
        end
        coverage(idx:(idx+length(Peptides_clean{i})-1)) = coverage(idx:(idx+length(Peptides_clean{i})-1)) + 1;
        if ~isempty(Phosphosites_raw{i})
            Phosphosites{i} = Phosphosites_raw{i} + idx - 1;
        end
    end

    PG_coverage{p} = sum(coverage>0)/length(coverage)*100;
    PG_pp{p,1} = unique([Phosphosites{idx_peptides}]);
    PG_pp{p,2} = convert_pp_sites(PG_pp{p,1}, Proteome(PG_idx).Sequence);
    PG_GN{p} = Proteome(PG_idx).GN;
end

% summary for all protein detected
PG_summary = struct(...
                    'ProteinAccessions', Protein_detected,...
                    'ProteinName', PG_GN,...
                    'Coverage', PG_coverage,...
                    'Phosphorylations', PG_pp(:,2));
                
% remove non-pp proteins
n_pp = cellfun('length',PG_pp(:,1));
PG_summary_pp = PG_summary(n_pp>0);
n_sites = sum(n_pp);

%% organize peptide centric data
CDC15_raw = [M.x_1_20190626_phos_1_f M.x_2_20190712_phos_1a M.x_3_20190712_phos_1b M.x_4_20190712_phos_1c M.x_5_20190712_phos_1d M.x_6_20190712_phos_1e M.x_7_20190712_phos_1f];
cdc15as1_raw = [M.x_8_20190626_phos_2_f M.x_9_20190712_phos_2a M.x_10_20190712_phos_2b M.x_11_20190712_phos_2c M.x_12_20190712_phos_2d M.x_13_20190712_phos_2e M.x_14_20190712_phos_2f];
CDC5_raw = [M.x_15_20190626_phos_3_f M.x_16_20190712_phos_3a M.x_17_20190712_phos_3b M.x_18_20190712_phos_3c M.x_19_20190712_phos_3d M.x_20_20190712_phos_3e M.x_21_20190712_phos_3f];
cdc5as1_raw = [M.x_22_20190626_phos_4_f M.x_23_20190712_phos_4a M.x_24_20190712_phos_4b M.x_25_20190712_phos_4c M.x_26_20190712_phos_4d M.x_27_20190712_phos_4e M.x_28_20190712_phos_4f];

% remove "1" as missing data
CDC15_raw(CDC15_raw==1) = nan;
cdc15as1_raw(cdc15as1_raw==1) = nan;
CDC5_raw(CDC5_raw==1) = nan;
cdc5as1_raw(cdc5as1_raw==1) = nan;

% calculate ratio
CDC15 = [nanmean(CDC15_raw,2) sum(~isnan(CDC15_raw),2)];
cdc15as1 = [nanmean(cdc15as1_raw,2) sum(~isnan(cdc15as1_raw),2)];
CDC5 = [nanmean(CDC5_raw,2) sum(~isnan(CDC5_raw),2)];
cdc5as1 = [nanmean(cdc5as1_raw,2) sum(~isnan(cdc5as1_raw),2)];
Ratio_Cdc15 = CDC15(:,1)./cdc15as1(:,1);
Ratio_Cdc5 = CDC5(:,1)./cdc5as1(:,1);

% calculate p-values with t test
[~,p_Cdc15] = ttest2(CDC15_raw, cdc15as1_raw,'Dim', 2);
[~,p_Cdc5] = ttest2(CDC5_raw, cdc5as1_raw,'Dim', 2);

% adjust p-values (multiple testing problem)
p_Cdc15_adj = mafdr(p_Cdc15,'BHFDR',true);
p_Cdc5_adj = mafdr(p_Cdc5,'BHFDR',true);

dN_Cdc15 = CDC15(:,2) - cdc15as1(:,2);
dN_Cdc5 = CDC5(:,2) - cdc5as1(:,2);

ProteinAccession = M.PG_ProteinAccessions;
Gene = M.PG_Genes;
Peptide = M.EG_PrecursorId;
M_processed = table(ProteinAccession, Gene, Peptide, Phosphosites, CDC15_raw, CDC15, cdc15as1_raw, cdc15as1, Ratio_Cdc15, p_Cdc15_adj, dN_Cdc15, CDC5_raw, CDC5, cdc5as1_raw, cdc5as1, Ratio_Cdc5, p_Cdc5_adj, dN_Cdc5);
writetable(M_processed, 'Peptides_mapped.xlsx');

%% make pp-sites centric data
M_processed(cellfun('isempty',M_processed.Phosphosites),:) = [];
n_pp_peptides = cellfun('length',M_processed.Phosphosites);
PP_peptides = cell(sum(n_pp_peptides),1);
idx1 = zeros(sum(n_pp_peptides),1);
for i = 1 : length(M_processed.Peptide)
    for j = 1 : n_pp_peptides(i)
        PP_peptides{sum(n_pp_peptides(1:i-1))+j} = [M_processed.Gene{i} '_' num2str(M_processed.Phosphosites{i}(j))];
    end
    idx1(sum(n_pp_peptides(1:i-1))+1:sum(n_pp_peptides(1:i))) = i;
end

M_ppSites = M_processed(idx1,:);
M_ppSites = addvars(M_ppSites, PP_peptides,'before',1);
M_ppSites = sortrows(M_ppSites,'PP_peptides');
[PP_sites_detected, ia2, ic2] = unique(PP_peptides,'stable');

%%
figure
histogram([1; diff(ia2)])
xlabel('# of peptides detected for a single pp site')
ylabel('Count (pp sites)')
title('Distribution of peptides per pp site (17123 total)')
saveas(gcf, 'Distribution of peptides per pp site','png')
%%
close all
clear M M_processed 
save('analysis_20200417')                
%% convert pp sites from number to string
function pp_sites = convert_pp_sites(pp, Seq)
    if ~isempty(pp)
        pp_sites = [Seq(pp(1)) num2str(pp(1))];
        for i = 2 : length(pp)
            pp_sites = [pp_sites ',' Seq(pp(i)) num2str(pp(i))];
        end
    else
        pp_sites = [];
    end
end