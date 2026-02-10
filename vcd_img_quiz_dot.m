function [img_quiz_dot, mask] = vcd_img_quiz_dot(params, verbose, store_imgs)
% VCD function
%  [img_quiz_dot, mask] = vcd_img_quiz_dot(params, verbose, store_imgs)
% 
% Purpose:
%   Create a simple white quiz dot image for experimental display.
%   See vcd_setStimParams.m for imagery quiz dot parameters.
%
%   For classic stimuli, we have 2 quiz dots per stimulus location (2 left, 
%   2 right). Images will be presented within the extent of the stimulus (+
%   some buffer).
% 
%   NOTE 1: that these dots are identical, hence we only create one dot image
%   and one transparency mask. We will log all the angles and unique image
%   numbers in the info table for each stimulus class.
%
%   NOTE 2: We smooth the circle edge (anti alias) with a Savitzky-Golay 
%   sliding polynomial filter.
% 
%   When params.stim.store_imgs = true, 
%   * We store a single matfile with all the dot uint8 images, as well as
%   corresponding transparency masks, and big info table. 
%   Single mat file is defined by params.stim.img_quiz_dot_file.stimfile:
%   fullfile(vcd_rootPath,'workspaces','stimuli',
%   params.disp.name,'img_quiz_dot_<params.disp.name>_<datestr(now,'yyyymmdd')>.mat')
%
%   When params.verbose = true, we create 2 types of figures stored in the
%   folder fullfile(vcd_rootPath, 'figs', <params.disp.name>,'simple_dot');
%
% INPUTS:
%  params                  : stim params struct (see vcd_setStimParams.m)
%    *** this function requires the following struct fields ***
%    bckgrnd_grayval        : (int) background gray value (128)
%    stim.dot.img_sz_pix    : (int) size of image support (pixels) of the
%                               dot. Must be an even number.
%    stim.img.quiz_dot_diam_deg  : (double): diameter of an individual quiz dot in degrees
%    stim.img.quiz_dot_diam_pix  : (double): diameter of an individual quiz dot in pixels (rounded to the nearest integer)
%    stim.img.img_sz_pix         = (double): spatial support for the quiz dot in pixels (diameter + 6 pixels) (ensure even nr of pixels)
%    stim.img.stimfile           = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('img_quiz_dot_%s',disp_params.name)); % mat file
%  verbose           : (logical) show debug figures
%  store_imgs        : (logical) store stimuli and debug figures as pngs 
%
% OUTPUTS:
%  dot               : (uint8) imagery quiz dot image used for VCD experiment
%                              height (pixels) by width (pixels) x 3 (rgb)
%  masks             : (uint8) alpha masks used for VCD experiment, to crop out image edges:
%                              height (pixels) by width (pixels))
% Written by Eline Kupers 2026/1

%% Check inputs

% Make sure the image has an uneven number of pixels, so we have center pix
if mod(params.stim.dot.img_sz_pix,2)~=0
    error('[%s]: image support size does not have an even nr of pixels!', mfilename)
end

% Create spatial support
x = (0:(params.stim.dot.img_sz_pix - 1) + 6); % add 6 pixels for spatial support
y = (0:(params.stim.dot.img_sz_pix - 1) + 6); % add 6 pixels for spatial support
x = x - x(end) / 2;
y = y - y(end) / 2;
[X, Y] = meshgrid(x, y);

% clean up
clear x y

% Center at zero for now (stimpresentation code will deal with x,y offset)
centerY = 0;
centerX = 0;

simple_dot_mask = (Y - centerY).^2 ...
    + (X - centerX).^2 <= params.stim.dot.radius_pix.^2;

% convert logical to double
simple_dot_mask_inv = (~simple_dot_mask);

% Smooth circle edge (anti alias) with a Savitzky-Golay sliding polynomial filter
smooth_dot = dealias(double(simple_dot_mask), ...
                     params.stim.dot.antialias.fcutoff_x, ...
                     params.stim.dot.antialias.fcutoff_y, ...
                     params.stim.dot.antialias.butterworth_order, ...
                     params.stim.dot.antialias.tapered_pad_pix);

% rescale range to [1 255]
pixelrange = [1 255]; % pixel range
scale_image = @(x) uint8((x - min(x(:))) / (max(x(:)) - min(x(:))) * diff(pixelrange) + min(pixelrange));
smooth_dot2 = scale_image(smooth_dot);

% correct colors: take everything outside the circle mask and make it gray
smooth_dot2(simple_dot_mask_inv) = params.stim.bckgrnd_grayval(1); % grayval is double

% now find all the pixels that are close to white, but not quite..
idx = smooth_dot2 > 230;
% for them to be white
smooth_dot2(idx) = params.stim.dot.color(1); % grayval is double

% add RGB copy in third dim
dot = repmat(smooth_dot2, [1 1 3]);

% Create alpha mask (same for all dots)
mask        = uint8(zeros(size(dot,1),size(dot,2)));
mask0       = (Y - centerY).^2 + (X - centerX).^2 <= (params.stim.dot.alpha_mask_diam_pix).^2;
mask(mask0) = 255;
mask        = uint8(mask);

% Get nr of angles
nr_angles = [1:length(params.stim.dot.ang_deg)];

% Add baseline location (no delta)
if ~isempty(params.stim.dot.delta_from_ref)
    dot_ref_locs = [0, params.stim.dot.delta_from_ref];
else
    dot_ref_locs = 0;
end

% Create reference angle and [x,y]-coords matrix
all_angles_deg = cat(1,params.stim.dot.ang_deg,params.stim.dot.ang_deg_delta);
all_xpos_pix   = cat(1,params.stim.dot.x0_pix,params.stim.dot.x0_pix_delta);
all_ypos_pix   = cat(1,params.stim.dot.y0_pix,params.stim.dot.y0_pix_delta);

% Wrap around 360 
all_angles_deg(all_angles_deg < 0) = 360+all_angles_deg(all_angles_deg < 0);

% convert degrees to radians (pol2cart expects angle to be in radians)
all_angles_rad = deg2rad(all_angles_deg);

% add conditions to table
stim_loc      = repmat({'left','right'},(length(nr_angles)/2),length(dot_ref_locs)); % stim loc refers to hemifield on display. We divide nr angles by 2, because they contain both L/R
stim_loc_idx  = repmat([1,2],(length(nr_angles)/2),length(dot_ref_locs)); % stim loc refers to hemifield on display. We divide nr angles by 2, because they contain both L/R
dot_ang_idx   = repmat(nr_angles,1,length(dot_ref_locs));
dot_angle_deg = reshape(all_angles_deg',1,[])';
dot_eccen     = repmat(params.stim.dot.iso_eccen, size(dot_angle_deg,1),1);
dot_xpos_pix  = reshape(all_xpos_pix',1,[])';
dot_ypos_pix  = reshape(all_ypos_pix',1,[])';

dot_radians     = reshape(all_angles_rad',1,[])';
dot_ref_locs    = repelem(dot_ref_locs,length(params.stim.dot.ang_deg))';
dot_ref_loc_idx = repelem([0:length(params.stim.dot.delta_from_ref)],length(params.stim.dot.ang_deg))';
unique_ref_im   = reshape(params.stim.dot.unique_im_nrs_wm_test,4,[])';
unique_im       = [params.stim.dot.unique_im_nrs_core, unique_ref_im(:)'];
specialcore_bool  = ismember(unique_im,params.stim.dot.unique_im_nrs_specialcore)';
assert(sum(specialcore_bool)==8);

info = table(unique_im(:), ... 
             stim_loc(:), ...
             stim_loc_idx(:),...
             dot_angle_deg(:), ...
             dot_ang_idx(:), ...
             dot_eccen(:),...
             dot_radians(:), ...
             dot_ref_locs(:), ...
             dot_ref_loc_idx(:),...
             dot_xpos_pix(:), ...
             dot_ypos_pix(:), ...
             specialcore_bool(:));
         
% add column names
info.Properties.VariableNames = {'unique_im','stim_pos','stim_pos_i','angle_deg','angle_i',...
    'eccen_deg','angle_rad','delta_deg','delta_i','dot_xpos_pix','dot_ypos_pix', 'is_specialcore'};

% Store if requested
if params.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(params.stim.dot.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',params.stim.dot.stimfile,datestr(now,30))),'dot','mask','info','-v7.3');
    
    saveDir = fileparts(fullfile(params.stim.dot.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',params.stim.dot.infofile,datestr(now,30))))
end




%% Visualize if requested
if verbose
    makeprettyfigures;
    if store_imgs
        saveFigDir = fullfile(vcd_rootPath,'figs',params.disp.name,'dot','visual_checks');
        if ~exist(saveFigDir,'dir'), mkdir(saveFigDir); end
    end
    %% Make PNG image of dot alone:
    if store_imgs
        imwrite(dot, fullfile(vcd_rootPath,'figs',params.disp.name,'dot',sprintf('%04d_singledot.png', params.stim.dot.unique_im_nrs_core(1))));
    end
    
    %% Visualize effect of alpha transparency mask
    figure(99); clf;
    subplot(131); imagesc(dot); colormap gray; axis image; set(gca, 'CLim', [1 255]);
    title('single dot'); xlabel('pixels'); ylabel('pixels')
    subplot(132); imagesc(mask); colormap gray; axis image;  set(gca, 'CLim', [1 255]);
    title('alpha mask'); xlabel('pixels'); ylabel('pixels')
    subplot(133); imagesc(dot, 'AlphaData',mask); colormap gray; axis image;  set(gca, 'CLim', [1 255]);
    title('dot+alpha mask'); xlabel('pixels'); ylabel('pixels')

    if store_imgs 
        print(fullfile(saveFigDir,'singledot_alphamask'),'-dpng','-r150');
    end
    
    %% Visualize reference location vs offset locations for WM images
    cmap = [0,0,0; lines(4)];
    for ii = params.stim.dot.unique_im_nrs_core
        
        idx1 = find(info.unique_im==ii);
        idx2 = find(info.angle_i==info.angle_i(idx1) &  sum(info.delta_deg==params.stim.dot.delta_from_ref,2));
        
        figure(1); clf;
        pax = polaraxes;
        polarscatter(info.angle_rad([idx1;idx2])', info.eccen_deg([idx1;idx2])',params.stim.dot.radius_pix,cmap,'LineWidth',3);
        pax.LineWidth = 1;
        pax.ThetaZeroLocation = 'top';
        pax.ThetaDir = 'clockwise';
        pax.FontSize = 20;
        if store_imgs
            title(sprintf('Single dot location %d + ref: [%d,%d,0,%d,%d]',ii, params.stim.dot.delta_from_ref));
            print(fullfile(saveFigDir,sprintf('%04d_singledot', ii)),'-dpng','-r150');
        end
    end
    
    
    %% All dot locations per hemifield with "fake" fixation circle
    bckground = uint8(ones(params.disp.h_pix,params.disp.w_pix))*params.stim.bckgrnd_grayval;
    im1 = repmat(bckground,[1 1 3]);
    dot_halfsz = (size(dot,1)/2);

    for ai = 1:length(info.angle_rad(1:16))
        
        [x_shift,y_shift] = pol2cart(info.angle_rad(ai),params.stim.dot.iso_eccen);
    
        ys = params.disp.yc - round(y_shift*params.disp.ppd);
        xs = params.disp.xc - round(x_shift*params.disp.ppd);
        dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz -1);
        dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz -1);
    
        im1(dot_coords_y, dot_coords_x,:) = dot;
    end
    
    figure;
    imshow(im1,[1 255]);
    hold on;
    h0 = drawcircle('Center',[params.disp.xc,params.disp.yc],'Radius',3,'color', [1 1 1],'LineWidth',1);
    h0.InteractionsAllowed = 'none';
    title('single dot locations - both hemifields');
    axis image;
    if store_imgs
%         set(gcf, 'InvertHardCopy', 'off');
        print(fullfile(saveFigDir,'singledot_all_loc'),'-dpng','-r150');
    end
    
    %% Visualize individual dot locations
    im1 =double(ones(params.disp.h_pix,params.disp.w_pix))*0.5; % gray background
    dot_ref_locs = [0, params.stim.dot.delta_from_ref];
    figure(99); 
    for aa = 1:length(info.dot_xpos_pix)
        
        im_nr = info.unique_im(aa);
        xpos_dots_to_plot = info.dot_xpos_pix(aa);
        ypos_dots_to_plot = info.dot_ypos_pix(aa);
        dlta = find(info.delta_deg(aa)==[0,params.stim.dot.delta_from_ref]);
        
        clf;
        imshow(im1,[0 1]);
        hold all;
        h0 = drawcircle('Center',[params.disp.xc,params.disp.yc],'Radius',3,'color', [1 1 1],'LineWidth',1);
        h0.InteractionsAllowed = 'none';
        h1 = drawcircle('Center',[xpos_dots_to_plot,ypos_dots_to_plot],'Radius',params.stim.dot.radius_pix,'color',[1 1 1],'EdgeAlpha',0);
        h1.InteractionsAllowed = 'none';
        h1.FaceAlpha=1;
        
        title(sprintf('Dot stim nr #%03d; Angle: %2.2f; Delta: %2.2f',info.unique_im(aa), info.angle_deg(aa),dot_ref_locs(dlta)), 'FontSize',20);
        axis image tight;
%         set(gcf, 'InvertHardCopy', 'off');
        if store_imgs
            print(fullfile(saveFigDir,sprintf('%04d_singledot_delta%02d', im_nr,dlta-1)),'-dpng','-r150');
        end
    end
    
    
end

return





