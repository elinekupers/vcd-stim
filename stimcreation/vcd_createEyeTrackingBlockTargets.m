function [sac_im,pupil_im_white,pupil_im_black] = vcd_createEyeTrackingBlockTargets(params)
% VCD function to create stimuli for eye tracking block
% 
%   [sac_im,pupil_im_white,pupil_im_black] = vcd_createEyeTrackingBlockTargets(params)
% 
% There are 5 saccade targets (params.exp.block.nr_of_saccades): central, 
% left, right, up, down, and 1 pupil trial with a mid-gray central fixation 
% target on a black and white background.
%
% An eye tracking block has the following sequence;
%   1. a 1-second fixation period (central fixation target on a mid-gray luminance)
%   2. 5 x 1.2 second = 6 seconds of saccades (mimicing EL HV5 grid,±3 deg in all directions)
%   3. a 2-seconds rest period (central fixation target on a mid-gray luminance)
%   4. a 4-seconds pupil trial: 3-s black adaptation, 1-s white screen to evoke max pupil response.
%
% The distance between the center and 4 left/right/up/down
% points are set as [xc,yc] ± 265 pixels (BOLDscreen)
% or ± 194 pixels (EIZOFLEXSCAN). This results in dots at the following
% BOLDscreen coordinates in pixels:
%                  [x3,y3]=[0,375]
% [x1,y1]=[695,0]  [x0,y0]=[960,640]   [x2,y2]=[1225,0]
%                  [x4,y4]=[0,905]
%
% EIZOFLEXSCAN coordinates in pixels:
%                  [x3,y3]=[0,406]
% [x1,y1]=[766,0]  [x0,y0]=[960,600]   [x2,y2]=[1154,0]
%                  [x4,y4]=[0,794]
%
% EMPIRICAL target distance:
% * BOLDscreen: 265 pixels, which corresponds to 3.0059 degrees.
% * PP room EIZOFLEX: 194 pixels, which corresponds to 3.0139 degrees.
%
% See vcd_setEyelinkParams.m for other parameters regarding Eyelink.
%
% INPUTS:
%  params           : (struct) parameter struct, which should contain the following fields:
%   * exp.block.nr_of_saccades         :  nr of saccade targets (positive integral number) (default is 5)
%   * stim.bckgrnd_grayval             :  mid-gray luminance level (default is 128)
%   * stim.el.point2point_distance_deg :  desired target distance in degrees from center of the screen (default = 3 deg)
%   * stim.el.point2point_distance_pix :  desired target distance in pixels from center of the screen (default is  265 pixels for BOLD screen and 194 pixels for PP room eizoflex)
%   * stim.el.total_target_diam_pix    :  total diameter of target (inner circle + outer rim) in pixels (same as thick fixation circle, 22 pixels for BOLDscreen)
%   * stim.el.target_center_diam_pix   :  diameter of the inner circle of the target in pixels (same as fixation circle,10 pixels for BOLDscreen)
%   * disp.h_pix                       :  height of the display in pixels
%   * disp.w_pix                       :  width of the display in pixels
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
saveStimDir = fullfile(vcd_rootPath,'workspaces','stimuli',params.disp.name);
saveFigDir = fullfile(vcd_rootPath,'figs',params.disp.name,'eye');
if ~exist(saveFigDir,'file'), mkdir(saveFigDir); end
if ~exist(saveStimDir,'file'), mkdir(saveStimDir); end

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
dot_halfsz     = size(fixation_rimthick1_black_on_gray,1)/2;
sac_im         = [];

for nr_loc = 1:target_locations
    
    % reset background
    im1 = repmat(bckground_gray,[1 1 3]);
    
    switch nr_loc
        
        case 1 % center
            ys = params.disp.yc;
            xs = params.disp.xc;
            dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz -1);
            dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz -1);
            
            im1(dot_coords_y, dot_coords_x,:) = fixation_rimthick1_black_on_gray;
        case 2 % left
            ys = params.disp.yc;
            xs = params.disp.xc - params.stim.el.point2point_distance_pix;
            dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz -1);
            dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz -1);
            
            im1(dot_coords_y, dot_coords_x,:) = fixation_rimthick1_black_on_gray;
        case 3 % right
            ys = params.disp.yc;
            xs = params.disp.xc + params.stim.el.point2point_distance_pix;
            dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz -1);
            dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz -1);
            
            im1(dot_coords_y, dot_coords_x,:) = fixation_rimthick1_black_on_gray;
        case 4 % up
            ys = params.disp.yc - params.stim.el.point2point_distance_pix;
            xs = params.disp.xc;
            dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz -1);
            dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz -1);
            
            im1(dot_coords_y, dot_coords_x,:) = fixation_rimthick1_black_on_gray;
        case 5 % down
            ys = params.disp.yc + params.stim.el.point2point_distance_pix;
            xs = params.disp.xc;
            dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz -1);
            dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz -1);
            
            im1(dot_coords_y, dot_coords_x,:) = fixation_rimthick1_black_on_gray;
    end
    
    sac_im = cat(4,sac_im,im1); 

    % Visualize image when requested
    if params.verbose
        figure(1); clf
        imshow(im1,[1 255]);
        axis image
    end
    
    % Store image when requested
    if params.store_imgs
        imwrite(im1, fullfile(saveFigDir,sprintf('vcd_eyetarget_%02d.png', nr_loc)));
    end
    
    
end

%% Pupil trial
bckground_white = uint8(ones(params.disp.h_pix,params.disp.w_pix))*255;
bckground_black = uint8(zeros(params.disp.h_pix,params.disp.w_pix));

pupil_im_white = repmat(bckground_white,[1 1 3]);
pupil_im_black = repmat(bckground_black,[1 1 3]);

ys = params.disp.yc;
xs = params.disp.xc;
dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz -1);
dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz -1);

pupil_im_white(dot_coords_y, dot_coords_x,:) = fixation_rimthick1_gray_on_white;
pupil_im_black(dot_coords_y, dot_coords_x,:) = fixation_rimthick1_gray_on_black;

% Visualize image when requested
if params.verbose
    figure(2); clf
    imshow(pupil_im_white,[1 255]);
    axis image
    
    figure(3); clf
    imshow(pupil_im_black,[1 255]);
    axis image
end

% Store image when requested
if params.store_imgs
    save(fullfile(saveStimDir,sprintf('eye_%s%s.mat', params.disp.name, datestr(now,30))), 'sac_im','pupil_im_white','pupil_im_black');
    imwrite(pupil_im_white, fullfile(saveFigDir,'vcd_eyepupil_white.png'));
    imwrite(pupil_im_black, fullfile(saveFigDir,'vcd_eyepupil_black.png'));   
end

    

