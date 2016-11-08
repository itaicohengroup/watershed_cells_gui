function [label_matrix, edge_image] = find_cells(im0, params, params_on)
% Implement image pre-processing, marker-controlled watershed segmentation,
% and post-processing on image im0 based on parameters params.
% Segementation results are return as a label matrix in label_matrix.
% Edge_image is a binary image highlighting the edges of cells.  
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
if params_on.equalization_cliplim
    img = adapthisteq(img, 'cliplimit', params.equalization_cliplim);
end

% background subtraction
if params_on.background_size
    im_bg = medfilt2(img, [1 1]*params.background_size);
    img = imsubtract(img, im_bg);
end

% median filter
if params_on.median_size
    img = medfilt2(img, [1 1]*params.median_size);
end

% gaussian filter
if params_on.gaussian_sigma
    
    img = imfilter(img, ...
        fspecial('gaussian', round(params.gaussian_sigma*3), params.gaussian_sigma));
end

% re-expand intensity range
img = mat2gray(img);

% set background pixels
thresh = multithresh(img, 2);
bg = img < thresh(1);

% remove areas from the foreground that are too small
if params_on.minimum_area
    bg = ~bwareaopen(~bg, params.minimum_area);
end

% watershed segmentation to find cells
L = marker_watershed(img, 'background', bg);

% remove objects that are too small, too big, or too dim
stats = regionprops(L, img, {'Area','MeanIntensity'});
areas = cat(1, stats.Area);
signals = cat(1, stats.MeanIntensity);
keep = true(size(areas));
if params_on.minimum_area
    keep = keep & (areas >= params.minimum_area);
end
if params_on.maximum_area
    keep = keep & (areas <= params.maximum_area);
end
if params_on.minimum_signal
    keep = keep & (signals >= params.minimum_signal);
end

% renumber label matrix
newix = zeros([1 max(L(:))], 'like', L);
newix(keep) = 1:sum(keep);
newix(~keep) = 0;
newix = [0 newix];
newL = newix(L+1);

% create an image showing the cells outlined
edgeim = bwperim(newL>0);

% put segmentation result into a struct 
label_matrix = newL;
edge_image = edgeim;
