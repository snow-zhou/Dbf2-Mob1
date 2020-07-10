%
% find Dbf2 substrates in Cdc15 targets
clear
close all
dirname = '/Users/snow/AmonLab/Projects/MEN/Proteomics/Kinase targets/Phospho-proteomics/Second run_201906_Cdc15-Cdc5_aF release/Results/analysis_20200417';
cd(dirname)
load('analysis_20200417')
filename = 'PG_summary_hits.xlsx';

M = readtable(filename);

%% identify Dbf2 consenses motif RxxS
% identify pp sites for each protein detected
n_protein_detected = length(M.ProteinName);
PG_pp_Dbf2 = cell(n_protein_detected,4);
n_pp_Dbf2 = zeros(n_protein_detected,2);

for p = 1 : n_protein_detected
    PG_idx = find(strcmp(ProteinAcessions, M.ProteinAccessions{p})>0);
    if isempty(PG_idx)
        continue;
    end
    [~, PG_pp_Dbf2{p,1}] = regexp(Proteome(PG_idx).Sequence, 'R..S');
    PG_pp_Dbf2{p,2} = convert_pp_sites(PG_pp_Dbf2{p,1}, Proteome(PG_idx).Sequence);
    if isempty(PG_pp_Dbf2{p,2})
        continue
    end
    PG_pp_Dbf2{p,3} = char(join(intersect(split(PG_pp_Dbf2{p,2},','), split(M.Phosphorylations{p},','),'stable'),',')); % detected pp sites that fit Dbf2 motif
    PG_pp_Dbf2{p,4} = char(join(intersect(split(PG_pp_Dbf2{p,3},','), split(M.PP_Cdc15_1{p},','),'stable'),',')); % Cdc15 dependent pp sites that fit Dbf2 motif
    if ~isempty(PG_pp_Dbf2{p,3})
        n_pp_Dbf2(p,1) = sum(PG_pp_Dbf2{p,3}==',')+1;
    end
    if ~isempty(PG_pp_Dbf2{p,4})
        n_pp_Dbf2(p,2) = sum(PG_pp_Dbf2{p,4}==',')+1;
    end
end

M_Dbf2 = addvars(M(:,[1 2 3 4 5 7 9 10]), PG_pp_Dbf2(:,4));
M_Dbf2.Properties.VariableNames{end} = 'PP_Dbf2';
writetable(M_Dbf2, 'PG_summary_hits_Dbf2.xlsx');

% calculate # of Dbf2 sites detected as pped and Cdc15 dependent;
n_Dbf2_sites = sum(n_pp_Dbf2);
%% protein targets
% target proteins for each category (Cdc15, Cdc5only, both, Cdc15 with Cdc5only sites, Dbf2, Dbf2 with Cdc5only sites)
n_substrates = [sum(~cellfun('isempty', M_Dbf2.PP_Cdc15_1)) ...
                sum(~cellfun('isempty', M_Dbf2.PP_Cdc5only)) ...
                sum(~cellfun('isempty', M_Dbf2.PP_Cdc15n5)) ...
                sum(~cellfun('isempty', M_Dbf2.PP_Cdc15_1)&~cellfun('isempty', M_Dbf2.PP_Cdc5only)) ...
                sum(~cellfun('isempty', M_Dbf2.PP_Dbf2)) ...
                sum(~cellfun('isempty', M_Dbf2.PP_Dbf2)&~cellfun('isempty', M_Dbf2.PP_Cdc5only))];
            
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