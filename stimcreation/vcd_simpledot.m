function [simple_dot, mask, info, p] = vcd_simpledot(p)
%
%  [simple_dot, p] = vcd_simpledot(p)
%
% Purpose:
%   Create a simple dot image for experimental display.
%   See vcd_setStimParams.m for gabor parameters.
%
% INPUTS:
%   p              : params struct, should contain field "p.stim.dot" 
%                    (see vcd_setStimParams.m)
%
% OUTPUTS:
%   simple_dot     : dot image (w (pixels) by h (pixels) x 3 (rgb))
%   masks          : alpha masks to crop out image edges
%   info           : table with info about simple_dot image features
%   p              : updated params struct

% Written by Eline Kupers 2024/12
%

%% Check inputs

% Make sure the image has an uneven number of pixels, so we have center pix
if mod(p.stim.dot.img_sz_pix,2)~=0
    error('[%s]: image support size does not have an even nr of pixels!', mfilename)
end

% Create spatial support
x = (0:(p.stim.dot.img_sz_pix - 1));
y = (0:(p.stim.dot.img_sz_pix - 1));
x = x - x(end) / 2;
y = y - y(end) / 2;
[X, Y] = meshgrid(x, y);

% clean up
clear x y

% Center at zero for now (stimpresentation code will deal with x,y offset)
centerY = 0;
centerX = 0;

simple_dot = (Y - centerY).^2 ...
    + (X - centerX).^2 <= p.stim.dot.radius_pix.^2;

% convert to uint8 
simple_dot = uint8(simple_dot);
simple_dot(simple_dot==0) = p.stim.bckgrnd_grayval;
simple_dot(simple_dot==1) = p.stim.dot.color(1);

% Create alpha mask (same for all dots)
mask  = uint8(zeros(size(simple_dot,1),size(simple_dot,2)));
mask0 = (Y - centerY).^2 + (X - centerX).^2 <= (p.stim.dot.alpha_mask_diam_pix).^2;
mask(mask0) = 255;
mask        = uint8(mask);

% Get nr of angles
nr_angles = [1:length(p.stim.dot.loc_deg)];

% Add baseline location (no delta)
if ~isempty(p.stim.dot.delta_from_ref)
    dot_ref_locs = [0, p.stim.dot.delta_from_ref];
else
    dot_ref_locs = 0;
end

% Create reference angle and [x,y]-coords matrix
all_angles_deg  = cat(1,p.stim.dot.ang_deg,p.stim.dot.ang_deg_delta);
all_xpos_pix   = cat(1,p.stim.dot.x0_pix,p.stim.dot.x0_pix_delta);
all_ypos_pix   = cat(1,p.stim.dot.y0_pix,p.stim.dot.y0_pix_delta);

% Wrap around 360 
all_angles_deg(all_angles_deg < 0) = 360+all_angles_deg(all_angles_deg < 0);

% convert degrees to radians (pol2cart expects angle to be in radians)
all_angles_rad = deg2rad(all_angles_deg);

% add conditions to table
stim_loc      = repmat({'left','right'},(length(nr_angles)/2),length(dot_ref_locs)); % stim loc refers to hemifield on display. We divide nr angles by 2, because they contain both L/R
bin           = repmat(nr_angles,1,length(dot_ref_locs));
dot_angle_deg = reshape(all_angles_deg',1,[])';
dot_eccen     = repmat(p.stim.dot.iso_eccen, size(dot_angle_deg,1),1);
dot_xpos_pix  = reshape(all_xpos_pix',1,[])';
dot_ypos_pix   = reshape(all_ypos_pix',1,[])';

dot_radians   = reshape(all_angles_rad',1,[])';
dot_ref_locs  = repelem(dot_ref_locs,length(p.stim.dot.loc_deg))';
unique_im     = [1:length(all_angles_deg), NaN(1,length(all_angles_deg)*length(p.stim.dot.delta_from_ref))];

info = table(unique_im(:), ... %             repmat(nr_angles',2*length(dot_ref_locs),1), ...
             bin(:), ...
             stim_loc(:), ...
             dot_angle_deg(:), ...
             dot_eccen(:),...
             dot_radians(:), ...
             dot_ref_locs, ...
             dot_xpos_pix(:), ...
             dot_ypos_pix(:));
         
% add column names
info.Properties.VariableNames = {'unique_im','bin','stim_pos','angle_deg','eccen_deg','angle_rad','delta_deg_ref','dot_xpos_pix','dot_ypos_pix'};

% Store
if p.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(p.stim.dot.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',p.stim.dot.stimfile,datestr(now,30))),'simple_dot','mask','info','-v7.3');
    
    saveDir = fileparts(fullfile(p.stim.dot.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',p.stim.dot.infofile,datestr(now,30))))
end


% debug figure
% figure(99); clf;
% subplot(131); imagesc(simple_dot); colormap gray; axis image; set(gca, 'CLim', [1 255]); 
% title('simple dot'); xlabel('pixels'); ylabel('pixels')
% subplot(132); imagesc(mask); colormap gray; axis image;  set(gca, 'CLim', [1 255]); 
% title('alpha mask'); xlabel('pixels'); ylabel('pixels')
% subplot(133); imagesc(simple_dot, 'AlphaData',mask); colormap gray; axis image;  set(gca, 'CLim', [1 255]); 
% title('dot+alpha mask'); xlabel('pixels'); ylabel('pixels')

cmap = [0,0,0; lines(4)];
sz = 50*ones(1,5);
for ii = 1:length(unique(info.unique_im(~isnan(info.unique_im))))
    
    idx1 = find(info.unique_im==ii);
    idx2 = find(info.bin==info.bin(idx1) &  sum(info.delta_deg_ref==p.stim.dot.delta_from_ref,2));
    
    figure(1); clf;
    pax = polaraxes;
    polarscatter(info.angle_rad([idx1;idx2])', info.eccen_deg([idx1;idx2])',sz,cmap,'LineWidth',3);
    pax.LineWidth = 3;
    pax.ThetaZeroLocation = 'top';
    pax.ThetaDir = 'clockwise';
    pax.FontSize = 20;
    
    title(sprintf('Simple dot location %d + ref: [%d,%d,0,%d,%d]',ii, p.stim.dot.delta_from_ref));
    print(fullfile(vcd_rootPath,'figs',sprintf('%02d_simpledot', ii)),'-dpng');
end



% 
%     subplot(2,5,ii+5)
%     if ii == 1
%         ax = polarscatter(info.ori_rad(strcmp(info.stim_pos,{'right'})& info.delta_deg_ref==dot_ref_locs(ii)), info.eccen_deg(strcmp(info.stim_pos,{'right'})& info.delta_deg_ref==dot_ref_locs(ii)),[],'k');
%     else
%         ax = polarscatter(info.ori_rad(strcmp(info.stim_pos,{'right'})& info.delta_deg_ref==dot_ref_locs(ii)), info.eccen_deg(strcmp(info.stim_pos,{'right'})& info.delta_deg_ref==dot_ref_locs(ii)),[],cmap2);
%     end
%     ax.LineWidth = 3;
%     title(sprintf('Right dot location ref: %d', dot_ref_locs(ii)));
% 
% 
% figure(2); clf;
% ax = polarscatter(info.ori_rad, info.eccen_deg,[],'k');
% ax.LineWidth = 3;
%  title('ALL dot locations')
%     
% figure(3); clf;
% ax = polarscatter(info.ori_rad(info.delta_deg_ref==0), info.eccen_deg(info.delta_deg_ref==0),[],'k');
% ax.LineWidth = 3;
% title('unique im dot locations')




return







