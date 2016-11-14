function [valid, iscolor] = check_image(rawim)
% check if valid color or grayscale image

if size(rawim, 3) == 3
    iscolor = true;
    valid = true;
    
elseif size(rawim, 3) == 1
    iscolor = false;
    valid = true;
    
else
    iscolor = false;
    valid = false;
end
