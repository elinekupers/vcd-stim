function inside = isInsideAperture(pos, center, radius)
% 
% function in = isInsideRegion(pos, center, radius)
% 
% checks if pos is inside a specified region
% 
% Input
%     pos       [x y] or two-column matrix 
%     center    the center of the ellipse [x y]
%     radius    horizontal and vertical radii of the ellipse [rx ry]
%
% Output
%   inside      1 if pos is inside the aperture, 0 otherwise
% 
%
% History:
% 2024.12.12 Adapted from RandomDots code written by Roozbeh Kiani 09/26/07
% 

n = size(pos,1);
d = ((pos-repmat(center,n,1))./repmat(radius,n,1));
inside = sqrt(sum(d.^2,2))<=1;

end