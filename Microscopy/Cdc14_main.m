clear
close all

addpath(genpath('/Microscopy'))

dirname = '';

nslice = 7;
dT = 5; %min/frame
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
	    Cdc14 = ['C2-' list_file(i).name(4:end-4)];
	    Cfi1 = ['C3-' list_file(i).name(4:end-4)];
	    NLS = ['C4-' list_file(i).name(4:end-4)];
	    %% segment BF
	    seg_BF_check(BF, nslice);
	 
	    %% segment SPB and nucleoli
	    seg_nucleoli_SPB(Cfi1, BF, nslice);
	    
	    %% find cells undergone mitosis
	    seg_assign_nucleoli(BF, Cfi1, dT);
	    
	    %% quantify NLS
	    quantify_NLSiRFP2(NLS, Cfi1, nslice, dT)
	    
	    %% quantify Cdc14
	    quantify_Cdc14_nucleoli_SPB(Cdc14, Cfi1, NLS, nslice, dT, PixelSize);
	end

    cd(dirname)
end