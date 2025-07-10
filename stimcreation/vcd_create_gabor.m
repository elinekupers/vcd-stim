function img = vcd_create_gabor(img_sz_pix,gauss_std_pix,sf,ori_deg, ph_deg, contrast, grayval)
% 
%  img = vcd_create_gabor(img_sz_pix,gauss_std_pix,sf,ori_deg, ph_deg, contrast)
%
% INPUTS:
%   img_sz_pix      : (double) width (or height) of square support (pixels)
%   gauss_std_pix   : (double) std of gaussian window (pixels)
%   sf              : (double) sf of grating (cycles per pixels)
%   ori_deg         : (double) orientation of grating, 0 degrees = north
%   ph_deg          : (double) phase of grating (degrees)
%   contrast        : (double) Michelson contrast of gabor (fraction 0-1)
%   grayval         : (double) background gray value (128)
%
% OUTPUT:
%   img             : (matrix) [x,y] gabor image uint8 pixels with
%                     luminance values (range = [1 255])
%
% Written by Eline Kupers @ UMN 2025/02/04

%%

% Create spatial support
x = (0:(img_sz_pix - 1));
y = (0:(img_sz_pix - 1));
x = x - x(end) / 2;
y = y - y(end) / 2;
[X, Y] = meshgrid(x, y);

% clean up
clear x y

% Convert deg to radians
ori_rad = ori_deg * (pi/180);
ph_rad  = ph_deg * (pi/180);

% Create Gaussian window
G = exp( - ((X.^2)/(2*gauss_std_pix^2) + (Y.^2)/(2*gauss_std_pix^2))); % Gaussian window
G = G ./ max(G(:)); % normalize height

% If you want to create circular mask to crop the support of the gabor gaussian window
% G_mask        = exp( - ((X.^2)/(crop_sd_pix^2) + (Y.^2)/(crop_sd_pix^2)));
% G_center      = floor(size(G_mask,1)/2)+1;
% thresh        = G_mask(G_center - ceil((crop_sd_pix)/2), G_center);
% outsideWindow = G_mask < thresh;
% mask          = double(~outsideWindow);
% mask(mask==0) = 0; 
% mask          = mask.*255;
% 
% % clean up
% clear G_mask insideWindow

% Create harmonic
H = cos(2 * pi * sf * ...
    (cos(ori_rad) * X + sin(ori_rad) * Y) + ph_rad);

% Create Gabor image
img0 = G .* H; % Gabor image = Gaussian .* harmonic

% % Contrast normalization
% maxval = abs(max(img0(:)));
% minval = abs(min(img0(:)));
% img1 = img0./max(maxval,minval); % VCD uses 128 as mid-grey level
% img2 = contrast * img1;

% Convert to image range [1 255]
img = (contrast*(img0*127))+double(grayval);

end