function fvalues = apply_function(colorim, label_matrix, fun)
% For each region in label_matrix, apply fun to the RGB color values in
% that region of colorim, and return the list of scalar outputs
%
% Note: This function does not check inputs for validity. See
% watershed_cells_gui for more documentation.
%
% Inputs:
%   fun             handle to a function. FUN should accept 3 vector 
%                   inputs: a list of red, green, and blue pixel values,
%                   and return on scalar value. 
%                   For example: fun = @(R,G,B) = mean(R);
%
%   colorim         3D, MxNx3 matrix specifying the true-color image
%
%   label_matrix	3D MxN matrix of integers specifying the segmentation.
%                   Region labels should be monotonically increasing
%                   postivie integers, where 0 indicates background/ignored
%                   regions, thus the maximum value in label_matrix (K)
%                   indicates the number of regions.
%
% Outputs:
%   fvalues         Kx1 vector, where fvalues(i) specifies the scalar
%                   function output from applying fun to region i.
%

% separate out RGB images;
imR = colorim(:,:,1);
imG = colorim(:,:,2);
imB = colorim(:,:,3);

% number of regions
num_regions = max(label_matrix(:));

% initialize output
fvalues = NaN(num_regions, 1);
hit_error = false;

% apply function to pixels in each region
for rr = 1:num_regions
    pix = label_matrix==rr;
    try
        fvalues(rr) = fun(imR(pix), imG(pix), imB(pix));
    catch err
        fvalues(rr) = NaN;
        hit_error = true;
    end
end

if hit_error
    warning(sprintf('Error applying function %s: %s Setting function output to NaN instead.',...
        func2str(fun), err.message));
end

