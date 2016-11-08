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
