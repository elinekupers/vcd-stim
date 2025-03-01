function [simple_dot, p] = vcd_simpledot(p)
%
%  [simple_dot, p] = vcd_simpledot(p)
%
% Purpose:
%   Create a simple dot image for experimental display.
%   See vcd_setStimParams.m for gabor parameters.
%
% INPUTS:
%   p       : dot params   (see vcd_setStimParams.m)
%
% OUTPUTS:
%   simple_dot     : dot image 
%   p              : updated params

% Written by Eline Kupers 2024/12
%

%% Check inputs

% Make sure the image has an uneven number of pixels, so we have center pix
if mod(p.stim.dot.img_sz_pix,2)==0
    p.stim.dot.img_sz_pix = p.stim.dot.img_sz_pix +1;
end

% Create spatial support
x = (0:(p.stim.dot.img_sz_pix - 1));
y = (0:(p.stim.dot.img_sz_pix - 1));
x = x - x(end) / 2;
y = y - y(end) / 2;
[X, Y] = meshgrid(x, y);

% clean up
clear x y

% Center at zero first
centerY = 0;
centerX = 0;

simple_dot = (Y - centerY).^2 ...
    + (X - centerX).^2 <= p.stim.dot.radius_pix.^2;

% convert to uint8 
simple_dot = uint8(simple_dot);
simple_dot(simple_dot==0) = p.stim.bckgrnd_grayval;
simple_dot(simple_dot==1) = p.stim.dot.color(1);

fprintf('\nDone!')


bins = [1:length(p.stim.dot.loc_deg)];

if ~isempty(p.stim.dot.delta_from_ref)
    dot_ref_locs = [0, p.stim.dot.delta_from_ref];
else
    dot_ref_locs = 0;
end

all_angles_deg  = p.stim.dot.loc_deg + dot_ref_locs';

% Wrap around
all_angles_deg(all_angles_deg < 0) = 360+all_angles_deg(all_angles_deg < 0);

% Rerotate to make 0 deg upper vertical     
all_angles_deg = all_angles_deg - 90; % 
all_angles_rad = deg2rad(all_angles_deg);

% add conditions to table
info = table(repmat(bins',length(dot_ref_locs),1), ...
             reshape(all_angles_deg',1,[])', ...
             reshape(all_angles_rad',1,[])', ...
             repelem(dot_ref_locs,length(p.stim.dot.loc_deg))');
% add column names
info.Properties.VariableNames = {'bin','ori_deg','ori_rad','delta_deg_ref'};

% Store
if p.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(p.stim.dot.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',p.stim.dot.stimfile,datestr(now,30))),'simple_dot','info','-v7.3');
    
    saveDir = fileparts(fullfile(p.stim.dot.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',p.stim.dot.infofile,datestr(now,30))))
end





% visualize dot locations
% display = vcd_getDisplayParams;
% bckground = uint8(ones(display.h_pix,display.w_pix))*p.stim.bckgrnd_grayval;
% 
% im1 = bckground;
% 
% for ang = p.stim.dot.loc_deg
%     angle = deg2rad(ang-90);
%     [x_shift,y_shift] = pol2cart(angle,p.stim.dot.iso_eccen);
% 
%     ys = display.yc + round(y_shift*display.ppd);
%     xs = display.xc + round(x_shift*display.ppd);
%     dot_halfsz = (size(simple_dot,1)/2)-0.5;
%     dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz);
%     dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz);
%     
%     im1( dot_coords_y, dot_coords_x) = simple_dot;
% end
% 
% figure;
% imshow(im1,[1 256]);
% hold on;
% plot(display.xc+[-5:5], display.yc.*ones(1,11),'k',display.xc*ones(1,11),display.yc+[-5:5],'k')
% title('simple dot - one hemifield');
% axis image;

return







