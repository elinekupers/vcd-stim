function [simple_dot, mask, p] = vcd_simpledot(p)
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
%   masks          : alpha masks to crop out image edges
%   p              : updated params

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

% Create info about dots
bins = [1:length(p.stim.dot.loc_deg_L)];

% Add baseline location (no delta)
if ~isempty(p.stim.dot.delta_from_ref)
    dot_ref_locs = [0, p.stim.dot.delta_from_ref];
else
    dot_ref_locs = 0;
end

all_angles_deg_L  = p.stim.dot.loc_deg_L + dot_ref_locs';
all_angles_deg_R  = p.stim.dot.loc_deg_R + dot_ref_locs';

% Wrap around 360 
all_angles_deg_L(all_angles_deg_L < 0) = 360+all_angles_deg_L(all_angles_deg_L < 0);
all_angles_deg_R(all_angles_deg_R < 0) = 360+all_angles_deg_R(all_angles_deg_R < 0);

all_angles_rad_L = deg2rad(all_angles_deg_L);
all_angles_rad_R = deg2rad(all_angles_deg_R);

% add conditions to table
stim_loc        = repmat({'left','right'},length(dot_ref_locs)*length(all_angles_deg_L),1);
dot_angles_L    = reshape(all_angles_deg_L',1,[])';
dot_angles_R    = reshape(all_angles_deg_R',1,[])';
dot_eccen_L     = repmat(p.stim.dot.iso_eccen, size(stim_loc,1),1);
dot_eccen_R     = dot_eccen_L;
dot_radians_L   = reshape(all_angles_rad_L',1,[])';
dot_radians_R   = reshape(all_angles_rad_R',1,[])';
dot_ref_locs_L  = repelem(dot_ref_locs,length(p.stim.dot.loc_deg_L))';
dot_ref_locs_R  = repelem(dot_ref_locs,length(p.stim.dot.loc_deg_R))';
unique_im = [1:length(all_angles_deg_L), NaN(1,length(all_angles_deg_L)*length(p.stim.dot.delta_from_ref)), ...
             1:length(all_angles_deg_R), NaN(1,length(all_angles_deg_R)*length(p.stim.dot.delta_from_ref))];

info = table(unique_im(:), ...
             repmat(bins',2*length(dot_ref_locs),1), ...
             stim_loc(:), ...
             cat(1, dot_angles_L(:), dot_angles_R(:)), ...
             cat(1, dot_radians_L(:), dot_radians_R(:)), ...
             cat(1, dot_ref_locs_L,dot_ref_locs_R), ...
             cat(1, dot_eccen_L(:),dot_eccen_R(:)));
         
% add column names
info.Properties.VariableNames = {'unique_im','bin','stim_pos','ori_deg','ori_rad','delta_deg_ref','eccen_deg'};

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
% subplot(131); imagesc(simple_dot); colormap gray; axis image; set(gca, 'CLim', [1 255])
% subplot(132); imagesc(mask(:,:,2)); colormap gray; axis image;  set(gca, 'CLim', [1 255])
% subplot(133); imagesc(simple_dot, 'AlphaData',mask(:,:,2)); colormap gray; axis image;  set(gca, 'CLim', [1 255])

% cmap = [0,0,0; parula(4)];
% sz = 50*ones(1,5);
% positions = {'left','right'};
% for pp = 1:length(positions)
%     for ii = 1:length(unique(info.unique_im(~isnan(info.unique_im))))
%         
%         idx1 = find(strcmp(info.stim_pos,positions(pp)) & info.delta_deg_ref==0 & info.unique_im==ii);
%         idx2 = find(strcmp(info.stim_pos,positions(pp)) & info.bin==info.bin(idx1) &  sum(info.delta_deg_ref==dot_ref_locs(2:end),2));
%         
%         figure(1); clf;
%         pax = polaraxes;
%         polarscatter(info.ori_rad([idx1;idx2])', info.eccen_deg([idx1;idx2])',sz,cmap,'LineWidth',3);
%         pax.LineWidth = 3;
%         pax.ThetaZeroLocation = 'top';
%         pax.ThetaDir = 'clockwise';
%         pax.FontSize = 20;
%         
%         title(sprintf('%s simple dot location %d + ref: [%d,%d,%d,%d,%d]', positions{pp},ii, dot_ref_locs));
%         print(fullfile(vcd_rootPath,'figs',sprintf('%02d_simpledot_%s', ii ,positions{pp})),'-dpng');
%     end
% end


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

% visualize dot locations
display = vcd_getDisplayParams;
bckground = uint8(ones(display.h_pix,display.w_pix))*p.stim.bckgrnd_grayval;

im1 = bckground;
dot_halfsz = (size(simple_dot,1)/2)-0.5;
cmap = varysat(parula(16),4);
angles_to_plot = [all_angles_deg_L,all_angles_deg_R];

for aa = 1:length(angles_to_plot)
    for bb = 1:length(dot_ref_locs)

        angle = deg2rad(angles_to_plot(aa) + dot_ref_locs(bb));
        [x_pos,y_pos] = pol2cart(angle,p.stim.dot.iso_eccen);
        
        ys = display.yc + round(x_pos*display.ppd);
        xs = display.xc + round(y_pos*display.ppd);
        
        figure(99); clf;
        imshow(im1,[1 255]);
        hold all;
        h0 = drawcircle('Center',[display.xc,display.yc],'Radius',16,'color', [1 1 1]);
        h0.InteractionsAllowed = 'none';
        h0.HandleVisibility = 'off';
        h1 = drawcircle('Center',[xs,ys],'Radius',dot_halfsz,'color', 'w','SelectedColor',[1 1 1]);
        h1.InteractionsAllowed = 'none';
        h1.FaceAlpha=1;
        
        title(sprintf('Angle: %02d; Delta: %02d',angles_to_plot(aa),dot_ref_locs(bb)), 'FontSize',20);
        axis image;
        print(fullfile(vcd_rootPath,'figs',sprintf('%02d_%02d_simpledot', aa, bb)),'-dpng','-r150','-painters');
    end
end


return







