function [label_matrix, CC, edge_image] = find_cells(im0, params)
% Implement image pre-processing, marker-controlled watershed segmentation,
% and post-processing on image im0 based on parameters params.
% Segementation results are return as a label matrix in label_matrix and
% also as a connected components structure in CC. Edge_image is a binary
% image highlighting the edges of cells.  
%
% Image analysis process
%   1. Pre-processing
%       > Import image (from "Path to image:")
%       > Adaptive histogram equalization (using "Equalization clip limit")
%       > Background subtraction using a median filter (using "background size")
%       > Median filter to smooth (using "Median size")
%       > Gaussian filter to remove noise (using "Gaussian radius")
%   2. Watershed segmentation
%       > Determine background using a conservative Otsu's threshold
%       > Watershed segmentation, with background enforced
%   3. Post-processing
%       > Remove cells that are too small (using "Minimum area")
%       > Remove cells that are too big (using "Maximum area")
%       > Remove cells that are too dim (using "Minimum signal")
%   4. Display result
%       > Display raw image
%       > Display cell outlines (using "Cell outline alpha")
%       > Display the number of cells found
%
% See watershed_cells_gui and associated documentation PDF for more
% information 
%
% Lena Bartell
%
% See also: imread adapthisteq medfilt2 imfilter fspecial multithresh
%           watershed bwconncomp bwperim regionprops labelmatrix


% read image and convert to grayscale
switch size(im0,3)
    case 1
        img = mat2gray(im0);
    case 3
        img = mat2gray(rgb2gray(im0));
end

% adaptive histogram equalization
if params.equalization
    img = adapthisteq(img, 'cliplimit', params.equalization);
end

% background subtraction
if params.background
    im_bg = medfilt2(img, [1 1]*params.background);
    img = imsubtract(img, im_bg);
end

% median filter
if params.median
    img = medfilt2(img, [1 1]*params.median);
end

% gaussian filter
if params.gaussian
    
    img = imfilter(img, ...
        fspecial('gaussian', round(params.gaussian*3), params.gaussian));
end

% re-expand intensity range
img = mat2gray(img);

% set background pixels
thresh = multithresh(img, 2);
bg = img < thresh(1);

% remove areas from the foreground that are too small
if params.minarea
    bg = ~bwareaopen(~bg, params.minarea);
end

% watershed segmentation to find cells
L = marker_watershed(img, 'background', bg);

% remove objects that are too small, too big, or too dim
stats = regionprops(L, img, {'Area','MeanIntensity'});
areas = cat(1, stats.Area);
signals = cat(1, stats.MeanIntensity);
keep = true(size(areas));
if params.minarea
    keep = keep & (areas >= params.minarea);
end
if params.maxarea
    keep = keep & (areas <= params.maxarea);
end
if params.minsignal
    keep = keep & (signals >= params.minsignal);
end

% renumber label matrix
newix = zeros(1, max(L(:)));
newix(keep) = 1:sum(keep);
newix(~keep) = 0;
newix = [0 newix];
newL = newix(L+1);

% create an image showing the cells outlined
CC = bwconncomp(newL>0);
edgeim = bwperim(newL>0);

% put segmentation result into a struct 
label_matrix = newL;
CC = CC;
edge_image = edgeim;
