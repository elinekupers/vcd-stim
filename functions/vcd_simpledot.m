function [img, coord, p] = vcd_simpledot(p, disp)
%
%  [img] = vcd_simpledot(p)
%
% Purpose:
%   Create a simple dot image for experimental display.
%   See vcd_setStimParams.m for gabor parameters.
%
% INPUTS:
%   p       : dot params   (see vcd_setStimParams.m)
%   disp    : display params (see vcd_getDisplayParams.m)
%
% OUTPUTS:
%   img     : dot image 
%   p       : updated params

% Written by Eline Kupers 2024/12
%

%% Check inputs

% Make sure the image has an uneven number of pixels, so we have center pix
if mod(p.img_sz_pix,2)==0
    p.img_sz_pix = p.img_sz_pix +1;
end

% Create spatial support
thetas = linspace(0,2*pi, p.img_sz_pix);
X = p.radius_pix * cos(thetas);
Y = p.radius_pix * sin(thetas);

coord  = [X; Y];

%% Get polar angle bins



bins


% translate
coord(1,:) = coord(1,:);% + p.x0_pix;
coord(2,:) = coord(2,:);% + p.y0_pix;

% figure(1); clf; hold all;
% r = rectangle('Position',  [-0.5*p.img_sz_pix, -0.5*p.img_sz_pix, p.img_sz_pix, p.img_sz_pix], ...
%     'FaceColor', [127 127 127]./255, 'EdgeColor', 'none');
% h = drawcircle('Center',[p.x0_pix,p.y0_pix],'Radius',p.radius_pix, 'Color',p.color, 'InteractionsAllowed', 'none', 'FaceAlpha', 1, 'LineWidth', 1);
% xlim([-0.5,0.5].*p.img_sz_pix)
% ylim([-0.5,0.5].*p.img_sz_pix)
% colormap gray; axis square; axis off;

% f = getframe(gcf);
% img = f.cdata;

return







