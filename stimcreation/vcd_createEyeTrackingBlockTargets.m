function vcd_createEyeTrackingBlockTargets(params)

% exp.block.eye_gaze_fix_ID         = 990; % fixation target
% exp.block.eye_gaze_sac_target_ID  = 991:995; % central, left, right, up, down.
% exp.block.eye_gaze_pupil_ID       = 996; % white then black
% % eye gaze block
% exp.block.nr_of_saccades      = 5;
% exp.block.eye_gaze_fix0       = presentationrate_hz * 1.0; % start with 1 second fixation period
% exp.block.eye_gaze_sac_target = presentationrate_hz * 1.2; % then 5x1.2 = 6 seconds of saccades (mimicing EL HV5 grid,±3 deg in all directions)
% exp.block.eye_gaze_fix1       = presentationrate_hz * 2.0; % then a 2-seconds rest trial
% exp.block.eye_gaze_pupil      = presentationrate_hz .* [3.0,1.0]; % then a 4-seconds pupil trial: 3-s black adaptation, 1-s white screen to evoke max pupil response.
% exp.block.total_eyetracking_block_dur = sum([exp.block.eye_gaze_fix0, ...
%     exp.block.eye_gaze_sac_target*exp.block.nr_of_saccades, ...
%     exp.block.eye_gaze_fix1, ...
%     exp.block.eye_gaze_pupil]);
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
% EMPIRICAL target distance:
% * BOLDscreen: 265 pixels, which corresponds to 3.0059 degrees.
% * PP room EIZOFLEX: 194 pixels, which corresponds to 3.0139 degrees.
% See vcd_setEyelinkParams.m for other parameters regarding Eyelink.
% stim.el.point2point_distance_deg = 3.0;                                % desired target distance (in deg) from fixation
% stim.el.point2point_distance_pix = round((stim.el.point2point_distance_deg*disp_params.ppd/2))*2; % desired target distance in pixels
% stim.el.total_target_diam_pix    = stim.fix.dotthickborderdiam_pix;    % same as thick fixation circle (22 pixels for BOLDscreen)
% stim.el.target_center_diam_pix   = stim.fix.dotcenterdiam_pix;         % same as inner fixation circle (10 pixels for BOLDscreen)
% 

support_x = 2*params.stim.fix.dotcenterdiam_pix;
support_y = support_x; 

target_im  = uint8(zeros([support_x, support_y, 3])); 

x = [1:(2*params.stim.fix.dotcenterdiam_pix)]-params.stim.fix.dotcenterdiam_pix;
[XX,~] = meshgrid(x,x);

% Where to insert luminance val? (divide diam by 2 to get radius, which is
% expected by makecircleimage)
fixationmask_inner    = find(makecircleimage(support_x, params.stim.fix.dotcenterdiam_pix/2));
fixationmask_rimthick  = find(makecircleimage(support_x, params.stim.fix.dotthickborderdiam_pix/2) - ...
                           makecircleimage(support_x,   params.stim.fix.dotcenterdiam_pix/2));  
                       
% everything is initially gray
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

saveStimDir = fullfile(vcd_rootPath,'workspaces','stimuli',params.disp.name);
saveFigDir = fullfile(vcd_rootPath,'figs',params.disp.name,'eye');
if ~exist(saveFigDir), mkdir(saveFigDir); end
if ~exist(saveStimDir), mkdir(saveStimDir); end

target_locations = 5;

bckground_gray = uint8(ones(params.disp.h_pix,params.disp.w_pix))*params.stim.bckgrnd_grayval;
dot_halfsz = size(fixation_rimthick1_black_on_gray,1)/2;
sac_im = [];

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

    if params.verbose
        figure(1); clf
        imshow(im1,[1 255]);
        axis image
    end
    
    if params.store_imgs
        imwrite(im1, fullfile(saveFigDir,sprintf('vcd_eyetarget_%02d.png', nr_loc)));
    end
    
    
end

% Pupil trial
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

if params.verbose
    figure(2); clf
    imshow(pupil_im_white,[1 255]);
    axis image
    
    figure(3); clf
    imshow(pupil_im_black,[1 255]);
    axis image
end

if params.store_imgs
    imwrite(pupil_im_white, fullfile(saveFigDir,'vcd_eyepupil_white.png'));
    imwrite(pupil_im_black, fullfile(saveFigDir,'vcd_eyepupil_black.png'));   
end

if params.store_imgs
    save(fullfile(saveStimDir,sprintf('eye_%s%s.mat', params.disp.name, datestr(now,30))), 'sac_im','pupil_im_white','pupil_im_black');
end
    

