clear
close all

addpath(genpath('/Microscopy'))

dirname = '';

nslice = 7;
dT = 3; %min/frame
PixelSize = 0.2146; % um
cd(dirname)

list_dir = dir;
list_dir = list_dir([list_dir.isdir]>0);
list_dir = list_dir(3:end);
%%
for f = 1 : length(list_dir)
    cd(list_dir(f).name)

%%
    list_file = dir('C1*.tif');
%%
    for i = 1 : length(list_file)
        BF = list_file(i).name(1:end-4);
        Mob1 = ['C2-' list_file(i).name(4:end-4)];
        SPB_nucleoli = ['C3-' list_file(i).name(4:end-4)];
        NLS = ['C4-' list_file(i).name(4:end-4)];
    %%    
        if exist(BF,'dir')
            continue
        end
        
        %% segment BF
        seg_BF_check(BF, nslice);
     
        %% segment SPB and nucleoli
        seg_nucleoli_SPB(SPB_nucleoli, BF, nslice);
        
        %% find cells undergone mitosis
        seg_assign_nucleoli(BF, SPB_nucleoli, dT);
        
        %% quantify NLS
        quantify_NLSiRFP2(NLS, SPB_nucleoli, nslice, dT)
        
        %% quantify Mob1
        quantify_Mob1(Mob1, SPB_nucleoli, NLS, nslice, dT, PixelSize)
        
        %% quantify Mob1 mutants (use if mutant arrested in anaphase)
%        quantify_Mob1_b(Mob1, SPB_nucleoli, NLS, nslice, dT, PixelSize)
        
    end

    cd(dirname)
end