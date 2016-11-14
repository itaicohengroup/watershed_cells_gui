function [label_matrix, edge_image] = find_regions(im0, params)
% Implement image pre-processing, marker-controlled watershed segmentation,
% and post-processing on image IM0 based on parameters structrue PARAMS and
% return the resulting label matrix LABEL_MATRIX and a binary image
% highlighting region edges in EDGE_IMAGE 
%
% Note: This function does not check inputs for validity. See
% watershed_cells_gui for more documentation.
%
% Inputs:
%   im0     MxNx1 (grayscale) or MxNx3 (rgb color) image matrix 
%
%   params  Structure specifying the processing parameters. Each field of
%           params is itself a structure with the fiels 'on' and 'value';
%           'on' is a true/false boolean specifying to keep/skip that
%           step of the processing and 'value' specifies the numerical
%           value of the parameter (only used if the step is 'on'). The
%           fields of params are associated with the processing step
%           outlined below should be named: 
%
%           'equalization_cliplim'  'value' should be in the range [0 1]
%           'background_size'       'value' should be an odd positive integer
%           'median_size'           'value' should be an odd positive integer
%           'gaussian_sigma'        'value' should be a positive number
%           'minimum_area'          'value' should be a positive integer
%           'maximum_area'          'value' should be a positive integer larger than the value of 'minimum_area'
%           'maximum_area'          'value' should be a positive integer
%           'minimum_signal'        'value' should be in the range [0 1]
%
% Outputs:
%   label_matrix    2D label matrix/image the same size as im0. In this
%                   matrix, pixels labeled as '0' are the background, 
%                   pixels labeled as '1' are part of region 1, pixels 
%                   labeled as '2' are part of region 2, etc.
%
%   edge_image      2D binary matrix image that is true (1) at the edge of
%                   each region in label_matrix and false (0) elsewhere.
%
% Image analysis process:
%   1. Pre-processing
%       > Adaptive histogram equalization (using "Equalization clip limit" parameter)
%       > Background subtraction using a median filter (using "background size" parameter)
%       > Median filter to smooth (using "Median size" parameter)
%       > Gaussian filter to remove noise (using "Gaussian radius" parameter)
%   2. Watershed segmentation
%       > Determine background using the lower value of a 2-level Otsu's threshold
%       > Watershed segmentation, with background enforced
%   3. Post-processing
%       > Remove regions that are too small (using "Minimum area" parameter)
%       > Remove regions that are too big (using "Maximum area" parameter)
%       > Remove regions that are too dim (using "Minimum signal" parameter)
%   4. Return result
%       > label_matrix matrix specifying the segmentation result
%       > edge_image specifying the found regions' outlines 
%
% Lena Bartell
%
% See also: imread adapthisteq medfilt2 imfilter fspecial multithresh
%           watershed bwconncomp bwperim regionprops labelmatrix


% convert image to double grayscale image in the range [0 1]
im0 = im2double(im0);
switch size(im0,3)
    case 1
        img = im0;
    case 3
        img = rgb2gray(im0);
end

% adaptive histogram equalization
if params.equalization_cliplim.on
    img = adapthisteq(img, 'cliplimit', params.equalization_cliplim.value);
end

% background subtraction
if params.background_size.on
    im_bg = medfilt2(img, [1 1]*params.background_size.value);
    img = imsubtract(img, im_bg);
end

% median filter
if params.median_size.on
    img = medfilt2(img, [1 1]*params.median_size.value);
end

% gaussian filter
if params.gaussian_sigma.on    
    img = imfilter(img, ...
        fspecial('gaussian', round(params.gaussian_sigma.value*3), params.gaussian_sigma.value));
end

% re-expand intensity range
img = mat2gray(img);

% set background pixels
thresh = multithresh(img, 2);
bg = img < thresh(1);

% remove areas from the foreground that are too small
if params.minimum_area.on
    bg = ~bwareaopen(~bg, params.minimum_area.value);
end

% watershed segmentation to find regions
L = marker_watershed(img, 'background', bg);

% remove objects that are too small, too big, or too dim
stats = regionprops(L, img, {'Area','MeanIntensity'});
areas = cat(1, stats.Area);
signals = cat(1, stats.MeanIntensity);
keep = true(size(areas));
if params.minimum_area.on
    keep = keep & (areas >= params.minimum_area.value);
end
if params.maximum_area.on
    keep = keep & (areas <= params.maximum_area.value);
end
if params.minimum_signal.on
    keep = keep & (signals >= params.minimum_signal.value);
end

% renumber label matrix
newix = zeros([1 max(L(:))], 'like', L);
newix(keep) = 1:sum(keep);
newix(~keep) = 0;
newix = [0 newix];
newL = newix(L+1);

% create an image showing the regions outlined
edgeim = bwperim(newL>0);

% put segmentation result into a struct 
label_matrix = newL;
edge_image = edgeim;

function L = marker_watershed(im, varargin)
% 
% L = MARKER_WATERSHED(IM) computes the watershed transform of image IM
%
% L = MARKER_WATERSHED(IM, Name, Value ...) forces the foreground and/or
% background pixels to be basins, as indicated by the name, value pairs.
%   'Background'    Binary image the same size as IM where pixels labeled
%                   as true are forced to be background basins
%   'Foreground'    Binary image the same size as IM where pixels labeled
%                   as true are forced to be local basins
%  
% The process is like MATLAB's watershed, but modified to:
%   1. invert im and convert intensity image in the domain [-1 0]
%   2. set positive background pixels (indicated in the binary background
%      input) to  -2
%   3. set positive forground pixels (indicated in the binary foreground
%      input) to -1
%   4. compute watershed transformation
%   5. enforce all regions with any background pixels to have the label 0
%   6. renumber all other regions to have unique consecutive positive
%      integer labels
%   7. preserve boundaries labeled as 0
%
% Lena Bartell, June 2016
%

% parse inputs
p = inputParser;
addOptional(p, 'Background', false(size(im)), @(y)isequal(size(im), size(y)));
addOptional(p, 'Foreground', false(size(im)), @(y)isequal(size(im), size(y)));
parse(p, varargin{:})
bg = p.Results.Background;
fg = p.Results.Foreground;

% invert raw image, projecting onto the domain [-1 0]
wsim = -mat2gray(im); 

% force background to be one basin
wsim(bg) = -2; 

% force foreground to be basin
wsim(fg) = -1;

% computer watershed transformation
L0 = watershed(wsim);

% convert to signed integer format & initialize output
[L0, L] = deal(L0);

% enforce all regions with any background pixels to be labeled 0
bgix = unique(L0(bg));
for ii = 1:length(bgix)
    L(L0==bgix(ii)) = 0;
end

% renumber label matrix to enforce all regions with any background pixels
% to be labeled 0 and to renumber foreground labels with unique
% consecutive positive integers (preserving boundaries labeled as 0)
newix = zeros([1 max(L0(:))], 'like', L);

allix = unique(L0);
bgix = unique(L0(bg));
fgix = allix(~ismember(allix, [0; bgix]));

newix(bgix) = 0;
newix(fgix) = 1:length(fgix);
newix = [0 newix];

L = newix(L0+1);