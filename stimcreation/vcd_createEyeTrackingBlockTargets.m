function [sac_im,pupil_im_white,pupil_im_black] = vcd_createEyeTrackingBlockTargets(params, verbose, store_imgs)
% VCD function to create stimuli for eye tracking bloc
% 
%   [sac_im,pupil_im_white,pupil_im_black] =
%        vcd_createEyeTrackingBlockTargets(params, verbose, store_imgs)
% 
% There are 5 saccade targets (params.exp.block.nr_of_saccades): central, 
% left, right, up, down, and 1 pupil trial with a mid-gray central fixation 
% target on a black and white background.
%
% An eye tracking block has the following sequence;
%   1. a 1-second fixation period (central fixation target on a mid-gray luminance)
%   2. 5 x 2 second = 6 seconds of saccades (mimicing EL HV5 grid, ±4 deg 
%      (EIZOFLEXSCAN) or ±3 deg (BOLDSCREEN) in all directions)
%   3. a 2-seconds rest period (central fixation target on a mid-gray luminance)
%   4. a 4-seconds pupil trial: 3-s black adaptation, 1-s white screen to evoke max pupil response.
%
% The coordinates of each saccade image define the top left and bottom right 
% coordinates of our rectangle support:
% [top-left-x top-left-y bottom-right-x bottom-right-y].
%
% The distance between the center and 4 left/right/up/down targets are set
% as [xc,yc] ± 264 pixels (BOLDscreen) or ± 256 pixels (EIZOFLEXSCAN). 
% This results in dots at the following pixel coordinates for the 5 targets
% BOLDscreen target rect coordinates in pixels 
% [x1,y1,x1,y2] = [top-left-x, top-left-y, bottom-right-x bottom-right-y]:
%
%                    [948,264,972,288]
% [684,528,708,552]  [948,528,972,552]   [1212,528,1236,552]
%                    [948,792,972,816]
%
% EIZOFLEXSCAN target rect coordinates in pixels 
% [x1,y1,x1,y2] = [top-left-x, top-left-y, bottom-right-x bottom-right-y]:
%
%                    [951,335,969,353]
% [695,591,713,609]  [951,591,969,609]   [1207,591,1225,609]
%                    [951,847,969,865]
%

% EMPIRICAL target distance:
% * BOLDscreen: 264 pixels, which corresponds to 2.994 degrees.
% * PP room EIZOFLEX: 256 pixels, which corresponds to 4.0061 degrees.
%
% See vcd_setEyelinkParams.m for other parameters regarding Eyelink.
%
% INPUTS:
%  params           : (struct) parameter struct, which should contain the following fields:
%   * exp.block.nr_of_saccades         :  nr of saccade targets (positive integral number) (default is 5)
%   * stim.bckgrnd_grayval             :  mid-gray luminance level (default is 128)
%   * stim.el.point2point_distance_deg :  desired target distance in degrees from center of the screen (default = 4 deg for EIZOFLEXSCAN and 3 deg for BOLDSCREEN)
%   * stim.el.point2point_distance_pix :  desired target distance in pixels from center of the screen (default is 264 pixels for BOLD screen and 256 pixels for PP room eizoflex)
%   * stim.el.total_target_diam_pix    :  total diameter of target (inner circle + outer rim) in pixels (same as thick fixation circle, 22 pixels for BOLDscreen)
%   * stim.el.target_center_diam_pix   :  diameter of the inner circle of the target in pixels (same as fixation circle, 12 pixels for BOLDscreen)
%   * disp.h_pix                       :  height of the display in pixels
%   * disp.w_pix                       :  width of the display in pixels
% verbose           : (logical) show debug figures
% store_imgs        : (logical) store stimuli and debug figures as pngs 
% 
% OUTPUTS:
%  sac_im           : (uint8) saccade stimuli (disp.h_pix x disp.w_pix x 3 x nr of saccade targets)
%  pupil_im_white   : (uint8) white background pupil trial stimulus (disp.h_pix x disp.w_pix x 3)
%  pupil_im_black   : (uint8) black background pupil trial stimulus (disp.h_pix x disp.w_pix x 3)
%
% Written by E. Kupers @ UMN 2025/05

% Define params
target_locations = params.exp.block.nr_of_saccades;
bckground_gray   = uint8(ones(params.disp.h_pix,params.disp.w_pix))*params.stim.bckgrnd_grayval;

% Define folders to store stimulus matlab file and pngs.
if store_imgs
    saveStimMatFileDir = fullfile(vcd_rootPath,'workspaces','stimuli',params.disp.name);
    saveFigDir = fullfile(params.saveFigsFolder,'eye');
    if ~exist(saveFigDir,'file'), mkdir(saveFigDir); end
    if ~exist(saveStimMatFileDir,'file'), mkdir(saveStimMatFileDir); end
end

% Create image support for saccade target
support_x = 2*params.stim.fix.dotcenterdiam_pix;
support_y = support_x; 

% Where to insert luminance val? (divide diam by 2 to get radius, which is
% expected by makecircleimage)
fixationmask_rimthick  = find(makecircleimage(support_x, params.stim.fix.dotthickborderdiam_pix/2) - ...
                              makecircleimage(support_x, params.stim.fix.dotcenterdiam_pix/2));  
                       
% Everything is initially gray
fixation_rimthick0_black_on_gray = params.stim.bckgrnd_grayval*ones(support_x*support_y, 3);
fixation_rimthick0_gray_on_white = 255*ones(support_x*support_y, 3);
fixation_rimthick0_gray_on_black = zeros(support_x*support_y, 3);

% add thin black rim 
fixation_rimthick0_black_on_gray(fixationmask_rimthick,:)   = zeros(size(fixationmask_rimthick,1), 3);  % add black rim
fixation_rimthick0_gray_on_white(fixationmask_rimthick,:)   = params.stim.bckgrnd_grayval*ones(size(fixationmask_rimthick,1), 3);  % add gray rim
fixation_rimthick0_gray_on_black(fixationmask_rimthick,:)   = params.stim.bckgrnd_grayval*ones(size(fixationmask_rimthick,1), 3);  % add gray rim

% reshape
fixation_rimthick1_black_on_gray       = reshape(fixation_rimthick0_black_on_gray,[support_x, support_y, 3]);
fixation_rimthick1_gray_on_white       = reshape(fixation_rimthick0_gray_on_white,[support_x, support_y, 3]);
fixation_rimthick1_gray_on_black       = reshape(fixation_rimthick0_gray_on_black,[support_x, support_y, 3]);


%% Place saccade targets on grey background

sac_im      = [];
dot_sz      = support_x;
half_dot_sz = dot_sz/2;

% Saccade target coordinates are in the format:
% [top-left-x, top-left-y, bottom-right-x bottom-right-y]
% We follow PTB's convention where we start at at the outer corner of the 
% top left pixel (0,0). We subtract one pixel from the support size because 
% the outer side of the bottom right pixel will otherwise be counted as the 
% 19th pixel in image space, resulting in a mismatch between the size of 
% destination rectangle [19x19] allocated for the target image [18x18].
target_rect = [0 0 dot_sz-1 dot_sz-1]; 

for nr_loc = 1:target_locations
    
    % reset background
    im1 = repmat(bckground_gray,[1 1 3]);
    
    switch nr_loc
        
        case 1 % center
            x1 = params.disp.xc - half_dot_sz; % top-left-x
            y1 = params.disp.yc - half_dot_sz; % top-left-y
            x2 = params.disp.xc + half_dot_sz; % bottom-right-x
            y2 = params.disp.yc + half_dot_sz; % bottom-right-y
            
            % center the target rect on a particular location of the screen.
            coords = CenterRect(target_rect,[x1 y1 x2 y2]); % same as CenterRect(target_rect,[0 0 params.disp.w_pix,params.disp.h_pix]);
            im1(coords(2):coords(4), coords(1):coords(3),:) = fixation_rimthick1_black_on_gray; % subtract one because coords are counting the outer side of last pixel
            
        case 2 % left
            x1 = params.disp.xc - half_dot_sz - params.stim.el.point2point_distance_pix; 
            x2 = params.disp.xc + half_dot_sz - params.stim.el.point2point_distance_pix;
            y1 = params.disp.yc - half_dot_sz; 
            y2 = params.disp.yc + half_dot_sz;
            coords = CenterRect(target_rect,[x1 y1 x2 y2]);
            im1(coords(2):coords(4), coords(1):coords(3),:) = fixation_rimthick1_black_on_gray;
            
        case 3 % right
            x1 = params.disp.xc - half_dot_sz + params.stim.el.point2point_distance_pix; 
            x2 = params.disp.xc + half_dot_sz + params.stim.el.point2point_distance_pix;
            y1 = params.disp.yc - half_dot_sz; 
            y2 = params.disp.yc + half_dot_sz;
            coords = CenterRect(target_rect,[x1 y1 x2 y2]);
            im1(coords(2):coords(4), coords(1):coords(3),:) = fixation_rimthick1_black_on_gray;
            
        case 4 % up
            x1 = params.disp.xc - half_dot_sz; 
            x2 = params.disp.xc + half_dot_sz;
            y1 = params.disp.yc - half_dot_sz - params.stim.el.point2point_distance_pix; 
            y2 = params.disp.yc + half_dot_sz - params.stim.el.point2point_distance_pix;
            coords = CenterRect(target_rect,[x1 y1 x2 y2]);
            im1(coords(2):coords(4), coords(1):coords(3),:) = fixation_rimthick1_black_on_gray;
            
        case 5 % down
            x1 = params.disp.xc - half_dot_sz; 
            x2 = params.disp.xc + half_dot_sz;
            y1 = params.disp.yc - half_dot_sz + params.stim.el.point2point_distance_pix; 
            y2 = params.disp.yc + half_dot_sz + params.stim.el.point2point_distance_pix;
            coords = CenterRect(target_rect,[x1 y1 x2 y2]);
            im1(coords(2):coords(4), coords(1):coords(3),:) = fixation_rimthick1_black_on_gray;
            
    end
    
    sac_im = cat(4,sac_im,im1); 

    % Visualize image when requested
    if verbose
        figure(1); clf
        imshow(im1,[1 255]);
        axis image
    end
    
    % Store image when requested
    if store_imgs
        imwrite(im1, fullfile(saveFigDir,sprintf('vcd_eyetarget_%02d.png', nr_loc)));
    end
end

%% Pupil trial
bckground_white = uint8(ones(params.disp.h_pix,params.disp.w_pix))*255;
bckground_black = uint8(zeros(params.disp.h_pix,params.disp.w_pix));

pupil_im_white = repmat(bckground_white,[1 1 3]);
pupil_im_black = repmat(bckground_black,[1 1 3]);

x1 = params.disp.xc - half_dot_sz;
y1 = params.disp.yc - half_dot_sz;
x2 = params.disp.xc + half_dot_sz;
y2 = params.disp.yc + half_dot_sz;
coords = CenterRect(target_rect,[x1 y1 x2 y2]);
pupil_im_white(coords(2):coords(4), coords(1):coords(3),:) = fixation_rimthick1_gray_on_white;
pupil_im_black(coords(2):coords(4), coords(1):coords(3),:) = fixation_rimthick1_gray_on_black;

% Visualize image when requested
if verbose
    figure(2); clf
    imshow(pupil_im_white,[1 255]);
    axis image
    
    figure(3); clf
    imshow(pupil_im_black,[1 255]);
    axis image
end

% Store image when requested
if store_imgs
    save(fullfile(saveStimMatFileDir,sprintf('eye_%s%s.mat', params.disp.name, datestr(now,30))), 'sac_im','pupil_im_white','pupil_im_black');
    imwrite(pupil_im_white, fullfile(saveFigDir,'vcd_eyepupil_white.png'));
    imwrite(pupil_im_black, fullfile(saveFigDir,'vcd_eyepupil_black.png'));   
end

    

