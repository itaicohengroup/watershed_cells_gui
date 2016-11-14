function [bw1, bw2, n1, n2, threshold] = apply_threshold(label_matrix, values, threshold, auto)
% Apply threshold to the values associated with each region in label_matrix
% and return binary images containing the regions that have values above
% (bw1) or below (bw2) the threshold.
%
% Note: This function does not check inputs for validity. See
% watershed_cells_gui for more documentation.
%
% Inputs:
%   label_matrix	3D MxN matrix of integers specifying the segmentation.
%                   Region labels should be monotonically increasing
%                   postivie integers, where 0 indicates background/ignored
%                   regions, thus the maximum value in label_matrix (K)
%                   indicates the number of regions.
%
%   values          Kx1 vector, where values(i) specifies the scalar
%                   value associated with region i that will be thresholded
%
%   threshold       Scalar number specifying the threshold to apply to the
%                   fvalues.
%
%   auto            Boolean. If true, override the input threshold and use
%                   an automatic Otsu's thresholding.
%
% Outputs:
%   bw1             MxN binary matrix/image specifying the pixels that
%                   belong to regions with a value > threshold 
%
%   bw2             MxN binary matrix/image specifying the pixels that
%                   belong to regions with a value < threshold 
%   
%   n1              Number of regions with value > threshold
%
%   n2              Number of regions with value < threshold
%
%   threshold       Scalar threshold used during segmentation (if auto was
%                   true, return the calculated Otsu's threshold)

% get region numbers
num_regions = length(values);
regions_ix = 1:num_regions;

% if auto, override the reshold
if auto
    threshold = multithresh(values, 1);
end

% create images of regions above/below threshold
bw1 = ismember(label_matrix, regions_ix(values > threshold));
bw2 = ismember(label_matrix, regions_ix(values < threshold));
n1 = sum(values > threshold);
n2 = sum(values < threshold);

