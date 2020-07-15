clear
close all

addpath(genpath('/Microscopy'))

dirname = '';
nslice = 8;
dT = 3; %min/frame
PixelSize = 0.2146; % um
dZ = 1; % z stack spacing um
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
        SPB = ['C2-' list_file(i).name(4:end-4)];
        Cfi1 = ['C3-' list_file(i).name(4:end-4)];
        %% segment BF
        seg_BF_check(BF, nslice);
     
        %% segment SPB
        seg_SPB_3D(SPB, BF, nslice);
        
        %% find cells undergone mitosis
        seg_assign_SPB(SPB, dT, PixelSize, dZ);
        
        %% quantify nucleioli seggregation
        quantify_Cfi1(Cfi1, SPB, nslice, dT, PixelSize)

    end

    cd(dirname)
end