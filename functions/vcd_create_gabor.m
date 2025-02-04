function img = vcd_create_gabor(img_sz_pix,gauss_std_pix,sf,ori_deg, ph_deg, contrast)
% 
%  img = vcd_create_gabor(img_sz_pix,gauss_std_pix,sf,ori_deg, ph_deg, contrast)
%
% INPUTS:
%
%
% OUTPUT:
%
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

% Mask Gaussian window support to 7 stds to avoid funny border effects
G_mask   = exp( - ((X.^2)/(2*3.5*gauss_std_pix^2) + (Y.^2)/(2*3.5*gauss_std_pix^2)));
G_center = floor(size(G_mask,1)/2)+1;
thresh   = G_mask(G_center - ceil(3.5*gauss_std_pix), G_center);
mask_idx = G_mask < thresh;
G(mask_idx)  = 0;

% clean up
clear G_mask mask_idx


% Create harmonic
H = cos(2 * pi * sf * ...
    (cos(ori_rad) * X + sin(ori_rad) * Y) + ph_rad);

% Create Gabor image
img0 = G .* H; % Gabor image = Gaussian .* harmonic

% Contrast normalization
maxval = abs(max(max(img0)));
minval = abs(min(min(img0)));
img0 = img0./max(maxval,minval);

img = contrast * img0;

% Convert to image range [0 255]
img  = floor((img*255)+127);

% % alternative normalization
% maxval_c = abs(max(max(img_c)));
% minval_c = abs(min(min(img_c)));
% k        = 127/max(maxval_c,minval_c); % assume 127 is mid-grey level
% img_c    = floor(img_c*k+127);

end