function img_segm = loccssegm(img, radius, alpha,lpf)
% This Matlab function is part of the code package released as a
% supplementary material for the article: 
%
% "Evaluation of methods for detection of fluorescence labeled 
% subcellular objects in microscope images" by P. Ruusuvuori et al. 
%
% We kindly request you to acknowledge the authors properly 
% (citation or request for permission from the authors) when using this
% function.
%
% Website: http://www.cs.tut.fi/sgn/csb/subcell
%
%
%LOC_CS_SEGM    Local comparison and selection segmentation algorithm
%
%   usage: img_segm = loc_cs_segm(img, radius, alpha)
%
%   where:
%       img           3D-array that contains the original image
%       radius  	  radius of the disc
%       alpha         magic parameter for the segmentation 
%       img_segm      3D-array that contains the image after the
%                     segmentation
%       
%   Please cite article "Ruusuvuori et al.: Evaluation of methods for
%   detection of fluorescence labeled subcellular objects in microscope
%   images" when using this function.
%
%   Author: Tarmo Äijö <tarmo.aijo@tut.fi>

% size of the image
img_size = size(img);

img = double(img);

% image should be 2D or 3D
if (length(img_size) ~= 2 && length(img_size) ~= 3) || nargin == 0
    error '2D or 3D image is required.'
end

% default radius is 5
if nargin < 2
    radius = 5;
end

% default alpha is 0.8
if nargin < 3
    alpha = 0.8;
end
if nargin < 4
    % perform low-pass filtering 1 = yes, 0 = no
    lpf = 1;
end
% low-pass filtering
if lpf == 1
    lpfl = 3;
    img = padarray(img,[lpfl lpfl],'both','replicate');
    img = filter2(ones(lpfl)/lpfl^2,img,'same');
    img = img(lpfl+1:end-lpfl,lpfl+1:end-lpfl);
end

% let's start with a circular averaging filter
H = fspecial('disk',radius); 

% fix the coefficents for our purposes
H(find(H)) = H(find(H))./H(find(H));
H = H/sum(sum(H(1:radius+1,radius+1:end)));

% we have to separate the quarters
H_quarters = zeros(2*radius+1, 2*radius+1, 4);
H_quarters(1:radius+1,radius+1:end,1) = H(1:radius+1,radius+1:end);
H_quarters(1:radius+1,1:radius+1,2) = H(1:radius+1,1:radius+1);
H_quarters(radius+1:end,1:radius+1,3) = H(radius+1:end,1:radius+1);
H_quarters(radius+1:end,radius+1:end,4) = H(radius+1:end,radius+1:end);

% convolute the whole image slide by slide with every quarter mask
if length(img_size) == 2 % 2D image
    img_conv = zeros(img_size(1),img_size(2), 4);
    for quarter=1:4
        img_conv(:,:,quarter) = imfilter(img,H_quarters(:,:,quarter),'same','conv'); 
    end
elseif length(img_size) == 3 % 3D image
    img_conv = zeros([img_size(1) img_size(2) img_size(3) 4]);
    for quarter=1:4
        img_conv(:,:,:,quarter) = imfilter(img,H_quarters(:,:,quarter),'same','conv');
    end
end

% decide if the pixel belongs to object {0, 1}
img_segm = zeros(img_size);
img_segm = max(img_conv,[],length(size(img_conv))) < alpha*img;
