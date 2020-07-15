function [WS2] = segmcells2(data)
  
% parameters
min_size = 500; % minimal size for an isolated cell
min_cell = 100; % minimal size for a bud

nslice = length(data);
% read brightfield (BF) image
I1 = mat2gray(data{nslice});
I2 = mat2gray(data{1});

% subtract the two images from each other
D = mat2gray(I1-I2);

% threshold by Otsu and remove small objects
T = graythresh(D);
BW1 = D>T;
if sum(BW1(:))/numel(BW1) > .5
  T = graythresh(D(BW1));
  BW1 = D>T*0.95;
end

% closing to fill holes
BW1 = bwareaopen(BW1,30);
BW1 = imopen(BW1,strel('square',1));
BW1 = imclose(bwmorph(BW1,'thin'),strel('square',5));
BW1 = bwareaopen(BW1,30);
BW1 = imdilate(BW1,strel('disk',1));

% find crude cell area mask
h = ones(5)/5^2; % kernel
meanimg = conv2(I1,h,'valid');
tmp = conv2(I1.^2,h,'valid')-meanimg.^2; % local variance
varimg = zeros(size(I1));
varimg(3:end-2,3:end-2) = tmp;
varimg_n = (varimg-min(varimg(:)))/max(varimg(:));
BW1s = varimg_n>graythresh(varimg_n)*.8;
% BW1s(1:2,:) = 1;
% BW1s(end-1:end,:) = 1;
% BW1s(:,1:2) = 1;
% BW1s(:,end-1:end) = 1;
se = strel('disk',7);
BW1d = imdilate(BW1s,se);
BW1f = imfill(BW1d,'holes');
BW1e = imerode(BW1f,se);
BW1c = bwareaopen(BW1e,min_size); % remove small objects that are likely not cells

% subtract cell membranes from the cell area mask
BW = BW1c&~BW1;

% detect insides of objects as cells (second detection method)
BWh = BW1;
BWh(:,1) = 0;
BWh(:,end) = 0;
BWf = ~imfill(BWh,[1 1]);

% combine two methods
if sum(BW(:))/numel(BW) > .5
    BW = BWf;
else
    BW = BW|BWf;
end

% opening to smooth edges
BW = imopen(BW,strel('square',3));

% remove small objects
BW = bwareaopen(BW,min_cell);

% remove large objects
% BW = BW&~bwareaopen(BW,3000);

% remove objects that touch the edge
BW_clean = rmedge(BW);
if sum(BW_clean(:))/sum(BW(:)) < .8 
    BW_clean = BW;
    BW_clean(1:2,:) = 0;
    BW_clean(end-1:end,:) = 0;
    BW_clean(:,1:2) = 0;
    BW_clean(:,end-1:end) = 0;
end
BW = BW_clean; 

% segment by watershed
cellarea = imfill(BW,'holes');
D = bwdist(~cellarea);
WS=WS_segmentationY(D,1);

% dilate each object to smooth edges
WS = bwlabel(WS);
ncells = length(unique(WS)) - 1;
se1 = strel('disk',5);
se2 = strel('disk',2);
WS2 = zeros(size(WS));
for i = 1 : ncells
    cell_mask = zeros(size(WS));
    cell_mask(WS == i) = 1;
    cell_mask = imclose(cell_mask,se1);
    cell_mask = imdilate(cell_mask,se2);
    WS2(cell_mask == 1) = i;

end
