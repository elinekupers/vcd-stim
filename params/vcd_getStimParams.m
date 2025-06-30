function stim = vcd_getStimParams(varargin)
% VCD function to get stimulus parameters:
%
%   stim = vcd_getStimParams(['disp_name','7TAS_BOLDSCREEN32',] ...
%                            ['load_params',false,] ...
%                            ['store_params',true,], ...
%                            ['saveInfoDir', fullfile(vcd_rootPath,'workspaces','info')], ...
%                            ['verbose',true]);
%
% Stimulus params such as size and location will depend on display params.
%
% All parameters are deterministic, meaning there is NO randomization 
% component involved in creating these stimulus parameters.
%
% INPUTS:
%  disp_name      : Display name to load params (see vcd_getDisplayParams.m)
%                   Default: '7TAS_BOLDSCREEN32'
%  load_params    : Load prior stored parameters or not. Default: true
%  store_params   : Store generated parameters or not. Default: true
%  save_info_dir  : Where to store generated parameters?. 
%                   Default: fullfile(vcd_rootPath,'workspaces','info')
%  verbose        : Print info in the command window (true) or not (false).
%                   Default: true
%
% OUTPUT:
%  stim           : struct with stimulus params, including:
%                   * fix (fixation dot)
%                   * bckground (pink noise background).
%                   * cd (contrast decrement)
%                   * el (eyelink)
%                   * gabor
%                   * rdk (random dot motion kinetograms)
%                   * dot (single, simple dot)
%                   * obj (complex objects)
%                   % ns (natural scenes)
%
% Written by Eline Kupers November 2024 (kupers [at] umn [dot] edu)

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addParameter('disp_name'    , '7TAS_BOLDSCREEN32', @(x) any(strcmp(x,{'7TAS_BOLDSCREEN32', 'KKOFFICE_AOCQ3277', 'PPROOM_EIZOFLEXSCAN', 'EKHOME_ASUSVE247','CCNYU_VIEWPIXX3D'})));                   
p0.addParameter('load_params'  , true, @islogical);                    
p0.addParameter('store_params' , true, @islogical); 
p0.addParameter('save_info_dir', fullfile(vcd_rootPath,'workspaces','info'), @ischar);
p0.addParameter('verbose'      , true, @islogical); 

% Parse inputs
p0.parse(varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

% Define stimulus classes
stim_class = {'gabor','rdk','dot','obj','ns'};

%% Obtain display params
disp_params = vcd_getDisplayParams(disp_name);

%% Load params from stored file if requested
if load_params
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('stim_%s*.mat',disp_name)));
    if ~isempty(d)
        if verbose
            fprintf('[%s]: Found %d stim params .mat file(s)\n',mfilename,length(d));
            if length(d) > 1
                warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
            end
            fprintf('[%s]: Loading stim params .mat file: %s\n', mfilename, d(end).name);
        end
        load(fullfile(d(end).folder,d(end).name),'stim');
    else
        error('[%s]: Can''t find stim params file!\n', mfilename)
    end
else 
    %% We will load params if requested
    if verbose, fprintf('[%s]: Define stim params\n', mfilename); end
    
    % Setup struct
    stim = struct();
    
    %% GENERAL PARAMS
    stim.store_imgs          = true;                                       % Store images when creating them? (see s_createStimuli.m)
    
    % we want to present images at a rate of 60 Hz (every 16.6667 ms)
    stim.presentationrate_hz = 60;                                         % rate of stimulus presentation (Hz)
    stim.f2f                 = disp_params.refresh_hz*(1/stim.presentationrate_hz); % frame-to-frame duration: nr of monitor refreshes in between presentation frames
                                                                           % we want to wait to achieve 60 Hz presentation rate. 
                                                                           % so this number is 2 for a monitor refresh rate of 120 Hz and 1 for 60 Hz. 

    % Ensure we have a round number of monitor refreshes to achieve 30 Hz presentation rate  
    assert(isint(disp_params.refresh_hz/stim.presentationrate_hz));
    
    stim.framedur_s          = stim.f2f*(1/disp_params.refresh_hz);        % duration for each presentation frame (should be 16.67 ms)
    assert(isequal(stim.framedur_s, (1/stim.presentationrate_hz)));        % check if this is 16.67 ms.
    
    stim.bckgrnd_grayval     = ceil(255/2);                                % background color (middle of 1-255 pixel lum)
    assert(isequal(stim.bckgrnd_grayval, 128));                            % we want 128 to be mean luminance
    
    stim.scfactor            = 1;                                          % no global scaling (only for specific stimulus images)

    % **** TEMPORAL **** 
    stim.stimdur_frames      = stim.presentationrate_hz * 1.0;             % nr of presentation frames (each 16.67 ms), currently 2 s (7TAS has 120 frames = 2 sec, PP room has 60 frames = 2 sec)
    
    % Define default params to inherit (or overwrite) by STIM CLASSES BELOW

    % **** SPATIAL ****
    % EMPIRICAL parafoveal stimulus locations:
    %  * BOLDscreen: 352 pixels (4.0059 deg). 
    %  * PProom EIZOFLEX: 258 pixels (4.0082 deg). 
    x0_deg                   = [-4 4];                                     % Desired x-center location for left right stim apertures (degrees)
    y0_deg                   = [0 0];                                      % Desired y-center location for left right stim apertures (degrees)
    x0_pix                   = round((x0_deg.*disp_params.ppd)/2)*2;       % Desired x-center location for left & right stim apertures (pixels). 
    y0_pix                   = round((y0_deg.*disp_params.ppd)/2)*2;       % Desired y-center location for left & right stim apertures (pixels).
    
    % EMPIRICAL parafoveal circular aperture for gabors, rdks, and objects:
    %  * BOLDscreen: 354 pixels (4.0059 deg). 
    %  * PProom EIZOFLEX: 258 pixels (4.0082 deg).
    parafov_circle_diam_deg  = 4;                                          % Desired  parafoveal circular diameter aperture (degrees).
    parafov_circle_diam_pix  = round((parafov_circle_diam_deg*disp_params.ppd)/2)*2;  % Desired parafoveal circular diameter aperture (pixels).
    
    % EMPIRICAL center square image for scenes:
    %  * BOLDscreen: 742 pixels (8.3965 degrees). NB: slightly fewer pixels than original NSD stim size = 744 pixels with BOLDscreen Nova1x32
    %  * PProom EIZOFLEX: 541 pixels (8.3965 degrees).
    ctr_square_deg           = 8.4;                                         % desired center square side length in degrees. 
    ctr_square_pix           = round((ctr_square_deg*disp_params.ppd)/2)*2; % desired center square side length in pixels
    
    
    %% NOISE BACKGROUND Puzzle piece
    % GENERAL
    stim.bckground.stimfile         = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('bckgrnd_%s',disp_params.name)); % mat file
    stim.bckground.infofile         = fullfile(vcd_rootPath,'workspaces','info',sprintf('bckgrnd_info_%s',disp_params.name)); % csv file
    
    % SPATIAL
    stim.bckground.alpha            = 1;                                    % exponent to apply to the amplitude spectrum (i.e. 1/f^alpha).  default: 1.
    stim.bckground.num              = 1;                                    % number of unique background images desired.  default: 1.
    stim.bckground.mode             = 0;                                    % mode of pinknoise function, 0 means fixed amplitude spectrum + random phase
    stim.bckground.std_clip_range   = 3.5;                                  % when converting pixel values to 1-255 range, how many std's of image values do we allow before clipping the range
    
    if verbose
        fprintf('*** %s SCREEN SIZE (height x width): FoV: [%2.2f,%2.2f] degrees visual angle. Resolution = [%02d,%02d] pixels.***\n', disp_params.name, disp_params.h_deg,disp_params.w_deg ,disp_params.h_pix,disp_params.w_pix);
    end
    %% FIXATION CIRCLE
    
    % GENERAL
    stim.fix.stimfile               = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('fix_%s',disp_params.name)); % mat file
    stim.fix.infofile               = fullfile(vcd_rootPath,'workspaces','info',sprintf('fix_info_%s',disp_params.name)); % csv file
    
    % TEMPORAL
    stim.fix.dotmeanchange          = 1.4*stim.presentationrate_hz;         % (time frames) dot changes occur every 1.4 seconds
    stim.fix.dotchangeplusminus     = 0*stim.presentationrate_hz;           % (time frames) earliest and latests time that dot changes occur. 2 seconds means [-1:1] from meanchange
    stim.fix.dres                   = [];                                   % rescale factor as a fraction between 0-1  
    stim.fix.fixsoafun              = @() round(stim.fix.dotmeanchange);    % frequency of fixation circle luminance change (in time frames)
    
    % SPATIAL
    stim.fix.dotcenterdiam_deg      = 0.14;                                 % (deg) idealized diameter of inner fixation circle
    stim.fix.dotthinborderdiam_deg  = 0.20;                                 % (deg) idealized diameter of thick fixation circle rim
    stim.fix.dotthickborderdiam_deg = 0.25;                                 % (deg) idealized diameter of thin fixation circle rim
    
    % THIN RIM Empirical diameters are:
    % * BOLDscreen: 12 pixels for center + 6 pixels for rim (split across both sides) = 18 pixels ~ 0.2037 deg. 
    % * PP room EIZOFLEX screen = 9 pixels for center + 4 pixels for rim (split across both sides) = 13 pixels ~ 0.2020 deg.
    
    % THICK RIM Empirical diameters are:
    % * BOLDscreen: 12 pixels for center + 10 pixels for rim (split across both sides) = 22 pixels ~ 0.249 deg. 
    % * PP room EIZOFLEX screen = 9 pixels for center + 7 pixels for rim (split across both sides) = 16 pixels ~ 0.2486 deg
    stim.fix.dotcenterdiam_pix      = round(stim.fix.dotcenterdiam_deg * disp_params.ppd);      % inner fixaton cirle diameter in pixels (boldscreen: 12 pixels, pproom: 9 pixels)
    stim.fix.dotthinborderdiam_pix  = round(stim.fix.dotthinborderdiam_deg * disp_params.ppd);  % total fixation circle diameter with thin border in pixels (during ITI/IBI) 
    stim.fix.dotthickborderdiam_pix = round(stim.fix.dotthickborderdiam_deg * disp_params.ppd); % total fixation circle diameter with thick border in pixels (during trial)
    
    % FIX CIRCLE APPEARANCE: LUMINANCE, OPACITY & COLOR
    lumminmaxstep_norm              = [0.15,0.85,5];                        % min, max, and nr of dot luminance values normalized to range between [0-1],
    dotlum_lin                      = linspace(lumminmaxstep_norm(1),lumminmaxstep_norm(2),lumminmaxstep_norm(3)); % normalized dot gray luminance levels (0-1):
    dotlum_sq                       = round(255*(dotlum_lin.^2));           % dot luminance values adjusted for monitor with linearized gamma, ranging between [1-255],
    stim.fix.dotlum                 = [dotlum_sq, stim.bckgrnd_grayval];
    stim.fix.dotopacity             = 0.5;                                 % dot and border have 50% opacity
    stim.fix.color                  = [255, 255, 255; 255 0 0; 0 0 0];     % white and red (for spatial cue) and black (pre eyetracking run) 

    if verbose
        fprintf('*** FIXATION CIRCLE: thin rim diameter = %1.2f deg (inner center diameter = %d pixels, inner+thin rim = %d pixels) ***\n', ...
            stim.fix.dotthinborderdiam_deg, stim.fix.dotcenterdiam_pix,stim.fix.dotthinborderdiam_pix);
        fprintf('*** FIXATION CIRCLE: thick rim diameter = %1.2f deg (inner center diameter = %d pixels, inner+thick rim = %d pixels) ***\n', ...
            stim.fix.dotthickborderdiam_deg, stim.fix.dotcenterdiam_pix,stim.fix.dotthickborderdiam_pix);
    end
    
    %% CONTRAST DECREMENT -- TEMPORAL MODULATION FUNCTION
    % For the CD task crossings, we make the stimulus dip in its contrast 
    % using sharp step function from mean contrast 100% --> 70%.
    %
    % * UNITS: the temporal function is in units of time frames, where each
    % time point has a duration of time frames (see presentationrate_hz).
    % We treat timepoint t=0 as the onset time of a frame flip, and t=1 is
    % the offset of the first frame flip and the onset of the second frame
    % flip (so in between monitor refreshes). 
    % 
    % * DURATION: For the contrast decrement, we implement an temporal
    % contrast modulation function in the shape of step function. The
    % support of the modulation function is 6 frames [1 1 1 0.7 0.7 0.7].
    % The contrast modulation function is a "one-directional" step
    % function: once the contrast is changed, it will stay decremented
    % until the end of the stimulus duration. The onset time of contrast
    % decrement is variable, and is sampled uniformly between meanchange ±
    % changeplusminus (500 ± 300 ms) by the function handle "cdsoafun". 
    % The earliest onset time is the start of the 12th time frame after 
    % stimulus onset (or 200 ms, assuming a 60 Hz presentation rate). 
    % The latest onset time is the start of the 48th time frame after
    % stimulus onset (or 800 ms, assuming a 60 Hz presentation rate).
    % The output of feval(cdsoafun) refers to the onset of the time frame 
    % when the contrast is changed from baseline (1) to mn_cd (0.7). In the 
    % case of our step function, the onset time refers to the third frame 
    % of the contrast modulation function (t_cmodfun). This means the
    % stimulus presentation code will ignore the support values = 1;
    % 
    % * MODULATION: How quickly the contrast decrement change occurs is 
    % determined by the % modulation. Current desired modulation is 30% 
    % from the mean luminance. We will only modulate the cued stimulus
    % location for classic stimulus classes.
    stim.cd.t_support_N            = 6;                                     % support of temporal modulatin function in number of time frames (one time frame = 16.67 ms)
    stim.cd.meanchange             = round(stim.presentationrate_hz * 0.5); % mean onset of temporal modulation function in time frames (500 ms = 30 time frames)  
    stim.cd.changeplusminus        = round(stim.presentationrate_hz * 0.3); % range of onsets: 500 ± 300 ms (so min = 200 ms, max = 800 ms)  
    stim.cd.min_cd                 = 0.2;                                   % stimulus contrast is reduced by 20% of mean luminance at lowest point of temporal gaussian window (note: this corresponds to subtracting a contrast fraction of 10.^(log10(c)-0.1))
    stim.cd.cdsoafun               = @() round(stim.cd.meanchange + stim.cd.changeplusminus*(2*(rand-.5))); % onset of temporal contrast modulation function 
    
    % Create 1D step function downward
    t_step                         = ones(1,stim.cd.t_support_N);           % units of presentation 16.67 ms frames
    t_step(round(stim.cd.t_support_N/2):end) = 1-stim.cd.min_cd;            % set second half of step function to stim.cd.min_cd
    stim.cd.t_cmodfun              = t_step; 

    %% EYETRACKING BLOCK PARAMS
    % Each run starts with an eyetracking "block", which mimics the 
    % eyelink calibration/validation points on the display, and evokes a 
    % pupil constriction response by adapting pupil to black screen and 
    % flashing a white screen for 1 second.
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
    % EMPIRICAL target distance:
    % * BOLDscreen: XXX pixels, which corresponds to [4] degrees.
    % * PP room EIZOFLEX: XXX pixels, which corresponds to [4] degrees.
    % See vcd_setEyelinkParams.m for other parameters regarding Eyelink.
    stim.el.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('eye_%s',disp_params.name)); % mat file
    stim.el.point2point_distance_deg = 4.0;                                % desired target distance (in deg) from fixation 
    stim.el.point2point_distance_pix = round((stim.el.point2point_distance_deg*disp_params.ppd/2))*2; % desired target distance in pixels
    stim.el.total_target_diam_pix    = stim.fix.dotthickborderdiam_pix;    % same as thick fixation circle (22 pixels for BOLDscreen)
    stim.el.target_center_diam_pix   = stim.fix.dotcenterdiam_pix;         % same as inner fixation circle (10 pixels for BOLDscreen)

    
    %% STIM PARAMS
    for ii = 1:length(stim_class)
        
        p = [];
        
        switch stim_class{ii}
            
            case {'gabor',1}
                
                % Where to store stimulus images?
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('gabor_%s',disp_params.name)); % mat file
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('gabor_info_%s',disp_params.name)); % csv file
                
                % GENERAL
                p.unique_im_nrs_core   = [1:24];                                     % Unique image nrs associated with the CORE 24 Gabors

                % TEMPORAL -- Fixed params
                p.duration        = stim.stimdur_frames;                             % frames (nr of monitor refreshes)

                % SPATIAL -- Fixed params
                p.img_sz_deg      = parafov_circle_diam_deg;                    % desired height (and width) of stimulus support (deg)
                p.img_sz_pix      = parafov_circle_diam_pix;                    % desired height (and width) of square stimulus support (pix)
                p.og_res_stim     = p.img_sz_pix;                               % resolution of stored dot stimuli (in pixels)
                p.dres            = (( (p.img_sz_pix/disp_params.ppd) /disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                p.iscolor         = false;                                      % use color or not [[[IF WE USE COLOR: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]]
                p.square_pix_val  = false;
                
                % Empirical gauss window std is:
                % * BOLDscreen:   44 pixels (0.4979 deg)
                % * EIZOFLEXSCAN: 32 pixels (0.4971 deg)
                p.gauss_std_deg   = 0.5;                                        % desired standard deviation of gaussian window in degrees
                p.gauss_std_pix   = round((p.gauss_std_deg * disp_params.ppd)/2)*2; % standard deviation of gaussian window in pixels.
                p.sf_cpd          = 1;                                          % spatial frequency (cycles/deg)
                p.cycles_per_pix  = 1/disp_params.ppd;                          % nr of cycles per pixel (1 cpd = 0.0113 cyc/pix for BOLD screen)
                p.x0_deg          = x0_deg;                                     % x-center loc in deg (translation from 0,0)
                p.y0_deg          = y0_deg;                                     % y-center loc in deg (translation from 0,0)
                p.x0_pix          = x0_pix;                                     % x-center loc in pix (translation from 0,0)
                p.y0_pix          = y0_pix;                                     % y-center loc in pix (translation from 0,0)
                p.ph_deg          = [0:(180/2):359];                            % 4 quadrature Gabor phases (deg) (0,90,180,270)
                
                % SPATIAL -- Manipulated params
                p.contrast        = [0.05, 0.20, 0.8];                          % Michelson contrasts [0-1] (fraction)
                p.num_ori         = 8;                                          % gabor orientations (deg), 0 = 12 o'clock
                p.ori_deg         = [0:(180/p.num_ori):179]+(0.5*(180/p.num_ori)); % rotate half shift away from vertical, to avoid ill-defined response options                                                                 
                                                                                % will return [11.25 33.75 56.25 78.75 101.25 123.75 146.25 168.75]
                
                % ensure orientations have equal distance from cardinal meridians
                assert(isequal(sort(abs(0-p.ori_deg),2),sort(abs(180-p.ori_deg),2)))
                assert(isequal(sort(abs(90-p.ori_deg(1:(p.num_ori/2))),2),sort(abs(90-p.ori_deg(((p.num_ori/2)+1):p.num_ori)),2)))
                                                                                
                % WORKING MEMORY: gabor orientation deltas for test images
                p.delta_from_ref        = [-16, -8, 8, 16];                     % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                                                                                % the bigger the delta, the easier the judgement in a trial
                p.unique_im_nrs_wm_test = [111:142];                            % unique image nrs associated with the 32 WM test images

                
                % check if all test images for WM have unique orientations
                tmp = p.ori_deg + [0, p.delta_from_ref]';
                tmp(tmp<0)   = tmp(tmp<0)+180;
                tmp(tmp>180) = tmp(tmp>180)-180;
                assert(isequal(length(unique(tmp(:))), length(tmp(:))));
                clear tmp

                % LTM
                p.ltm_pairs = [];
                
                % IMAGERY
                p.unique_im_nrs_specialcore   = p.unique_im_nrs_core(17:end);               % SELECTED UNIQUE IMAGES (SUBSET of all 24) (only high contrast)
                p.imagery_sz_deg              = 5.658;                                      % QUIZ DOT PARAMS (STIM 2) desired diameter (deg) of the quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix              = (round(p.img_sz_deg* disp_params.ppd)/2)*2; % QUIZ DOT PARAMS (STIM 2) diameter of quiz dot image in pixels (ensure even nr of pixels)                
                p.unique_im_nrs_img_test      = [423:582];                                  % Unique image nrs associated with the 8*20=160 IMG gabor test dot images
                p.imagery_quiz_images         = [ones(1,10),2.*ones(1,10)];                 % quiz dots overlap (1) or not (2)

                % Add params to struct
                stim.gabor = p;
                
                
            case {'rdk',2}
                % Where to store stimulus images?
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('rdk_%s',disp_params.name)); % mat file
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('rdk_info_%s',disp_params.name)); % csv file
                
                % GENERAL
                p.unique_im_nrs_core = [25:48];                                  % Unique image nrs associated with the CORE 24 RDK stimuli
                p.iscolor            = false;                                    % Use color or not? 
                if p.iscolor
                    p.square_pix_val = true;                                     % [IF YES: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]
                else
                    p.square_pix_val = false;
                end

                % SPATIAL -- general stim params
                p.img_sz_deg      = parafov_circle_diam_deg;                      % stimulus aperture diameter (deg)
                p.img_sz_pix      = parafov_circle_diam_pix;                      % stimulus aperture diameter (pix)
                p.og_res_stim     = p.img_sz_pix;                                 % resolution of stored dot stimuli (in pixels)
                p.dres            = (( (p.img_sz_pix/disp_params.ppd) /disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                p.x0_deg          = x0_deg;                                       % Desired x-center loc of stimulus in deg (translation from 0,0)
                p.y0_deg          = x0_deg;                                       % Desired y-center loc of stimulus in deg (translation from 0,0)
                p.x0_pix          = x0_pix;                                       % x-center loc in pix (translation from 0,0)
                p.y0_pix          = y0_pix;                                       % y-center loc in pix (translation from 0,0)
                
                % SPATIAL -- RDK specific
                p.dots_size_deg   = 0.068/2;                                      % single dot radius in deg (we go by radius because that is what "drawcircle" expects
                p.dots_size_pix   = round(p.dots_size_deg*disp_params.ppd);       % single dot radius in pixels
                p.dots_color      = [255 255 255; 1 1 1];                    % 50:50 white:black, color in RGB and converted to [0-1] as expected by stimulus creation function
                p.dots_density    = 15.9;                                         % density of dots within circular aperture (dots/deg^2)
                p.max_dots_per_frame = round(p.dots_density*(pi*(p.img_sz_deg/2)^2)); % how many individual dots within a square support. Density is 15.9 dots / deg^2   (Number is similar to Kiani lab, rokers lab aims for 150) and roughly matches to nr of pixels in aperture
                p.dots_contrast   = 1;                                            % Michelson [0-1] (fraction)
                
                % TEMPORAL -- fixed params
                p.duration        = stim.stimdur_frames;                          % frames (nr of monitor refreshes)
                p.dots_speed      = (8*disp_params.ppd)/stim.presentationrate_hz; % speed in pixels per frame (currently set to 8 deg/s). For reference: Kiani lab uses usually 5 to 10 deg/s. Rokers lab uses 5 deg/s.
                p.dots_interval   = 1;                                            % update dots every frame   (For reference: Kiani's 75 hz refresh rate + interval = 3 -->  25 frames/sec)
                p.dots_lifetime   = 0.1 * stim.presentationrate_hz;               % 6 frames (0.1 seconds)
                
                % TEMPORAL -- manipulated params
                p.dots_coherence  = [0.3, 0.65, 1.0];                             % fraction of coherent moving dots. Kiani lab uses usually one of these [0 0.032 0.064 0.128 0.256 0.512]
                
                % SPATIAL -- manipulated params
                p.num_mot_dir      = 8;                                           % number of sampled motion directions
                p.dots_direction   = [33.75:(360/p.num_mot_dir):359];             % sample direction of coherent motion from [0-359] in deg (0 deg is aligned with 12 o'clock)
                                                                                  % turns out to be: [33.7500, 78.7500. 123.7500, 168.7500, 213.7500, 258.7500, 303.7500, 348.7500] deg
                % ensure equal distance from cardinal meridians
                assert(isequal(abs(0-p.dots_direction(1:(p.num_mot_dir/2))), abs(180-p.dots_direction(((p.num_mot_dir/2)+1):p.num_mot_dir))));
                assert(isequal(abs(90-p.dots_direction(1:(p.num_mot_dir/2))), abs(270-p.dots_direction(((p.num_mot_dir/2)+1):p.num_mot_dir))));

                
                % WORKING MEMORY: motion direction deltas for test images
                p.delta_from_ref        = [-20, -10, 10, 20];                       % How much should stim iso-eccen loc deviate from reference (WM: double epochs)
                p.unique_im_nrs_wm_test = [143:174];                                % Unique image nrs associated with the 32 WM RDK test stimuli
 
                % check if all test images for WM have unique orientations
                tmp          = p.dots_direction+[0, p.delta_from_ref]';
                tmp(tmp<0)   = tmp(tmp<0)+360;
                tmp(tmp>360) = tmp(tmp>360)-360;
                assert(isequal(length(unique(tmp(:))), length(tmp(:))));
                clear tmp 

                % LTM
                p.ltm_pairs = [];
                
                % IMAGERY: 
                p.unique_im_nrs_specialcore    = p.unique_im_nrs_core(17:end);       % SELECTED UNIQUE IMAGES (SUBSET of all 24) (only high coherence) 
                p.imagery_sz_deg               = 5.658;                              % desired diameter (degree) of the second, quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix               = (round(p.img_sz_deg* disp_params.ppd)/2)*2; % diameter of quiz dot image (pixels) (ensure even nr of pixels)
                p.unique_im_nrs_img_test  	   = [583:742];                          % Unique image nrs associated with the 8*20=160 IMG RDK test dot images
                p.imagery_quiz_images          = [ones(1,10),2.*ones(1,10)];         % quiz dots overlap (1) or not (2)
                
                % Add params to struct
                stim.rdk = p;
                
            case {'dot',3}
                % Where to store stimulus images?
                p.stimfile  = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('dot_%s',disp_params.name)); % mat file
                p.infofile  = fullfile(vcd_rootPath,'workspaces','info',sprintf('dot_info_%s',disp_params.name)); % csv file
                p.iscolor   = false;                                        % Use color or not? 
                if p.iscolor
                    p.square_pix_val     = true;                            % [IF YES: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]
                else
                    p.square_pix_val     = false;
                end
                
                % GENERAL
                p.unique_im_nrs_core   = [49:64];                                     % Unique image nrs associated with the 16 single dot stimuli
                
                % TEMPORAL
                p.duration        = stim.stimdur_frames;                              % frames (nr of monitor refreshes)
                
                % SPATIAL
                % Empirical single dot radius is:
                % * BOLDscreen:   44 pixels (0.4979 deg)
                % * EIZOFLEXSCAN: 32 pixels (0.4971 deg)
                p.img_sz_deg      = 1.0;                                     	 % desired spatial support of dot image in deg.
                p.img_sz_pix      = round((p.img_sz_deg * disp_params.ppd)/2)*2; % spatial support of dot image in pixels (BOLDscreen: 88 pixels. EIZOFLEXSCAN: 64 pixels)
                p.og_res_stim     = p.img_sz_pix;                                % resolution of stored dot stimuli (in pixels)
                p.dres            = (( (p.img_sz_pix/disp_params.ppd) /disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                p.radius_deg      = 0.25;                                        % desired dot radius in deg
                p.radius_pix      = round((p.radius_deg * disp_params.ppd)/2)*2; % dot radius in pix
                p.color           = [255 255 255];                               % white in uint8 RGB
                p.contrast        = 1;                                           % Michelson [0-1] (fraction)
                p.square_pix_val  = false;                                       % no need to square pixel luminance values for linear CLUT
                
                % Alpha mask
                p.alpha_mask_diam_pix = p.radius_pix+1;                          % add one pixel for alpha mask (EK: why 1?)

                % Savitzky-Golay sliding polynomial filter (anti-aliasing
                % circle edge) This filter is the product of a 2D low-pass
                % Butterworth filters: one that cuts out the high
                % frequencies in x direction, one that cuts high
                % frequencies in the y direction.  
                p.antialias.fcutoff_x         = 0.4;                             % Cutoff frequencies for x direction (dB)                                   
                p.antialias.fcutoff_y         = 0.4;                             % Cutoff frequencies for y direction (dB)  
                p.antialias.butterworth_order = 2;                               % Order of the Butterworth filter 
                p.antialias.tapered_pad_pix   = 6;                               % Width of tapered padding (pixels)
                
                % DOT LOCATIONS
                p.num_loc         = 16;                                          % number of equally spaced dot angles (deg), 0 = 12 o'clock
                % idealized angle of center dot loc in deg, 0 deg = 12 o'clock. 
                % rotate half a shift away from vertical, to avoid
                % ill-defined response options. Angles start with 8 left dots 
                % (from ~11 o'clock to ~7 o'clock_, then 8 right dots (from 
                % ~5 o'clock to ~1 o'clock.)
                p.ang_deg = fliplr([0:(360/p.num_loc):359]+(0.5*(360/p.num_loc))); % single dot location angle (deg) 
                % location angles turn out to be: [348.75, 326.25, 303.75,
                % 281.25, 258.75, 236.25, 213.75, 191.25, 168.75, 146.25,
                % 123.75, 101.25, 78.75, 56.25, 33.75, 11.25] deg
                
                p.iso_eccen = 4.0;                                         % desired iso-eccentric dot location (for all angles)
                p.eccen_deg = repmat(p.iso_eccen,1,length(p.ang_deg));     % desired eccen of center dot loc in deg (translation from center screen 0,0)
                
                % ensure equal distance from cardinal meridians
                assert(isequal(abs(180-p.ang_deg(1:(p.num_loc/2))), abs(0-p.ang_deg(((p.num_loc/2)+1):p.num_loc))));
                assert(isequal(abs(270-p.ang_deg(1:(p.num_loc/2))), abs(90-p.ang_deg(((p.num_loc/2)+1):p.num_loc))));
                
                % Convert dot polar angle coords to cartesian coords
                [x,y] = pol2cart(deg2rad(p.ang_deg-90),p.eccen_deg);               % subtract 90 deg to ensure 0 deg is now 12 o'clock in x,y-pixel space
                p.x0_pix        = disp_params.xc + round(x * disp_params.ppd);     % desired x-center loc in pix (translation from upper left corner [0,0])
                p.y0_pix        = disp_params.yc + round(y * disp_params.ppd);     % desired y-center loc in pix (translation from upper left corner [0,0])
                
                % WORKING MEMORY: dot angle position deltas for test images
                p.delta_from_ref             = [-12, -6, 6, 12]; %[-16, -8, 8, 16];                   % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                                                                                   % the bigger the delta, the easier the trial
                p.unique_im_nrs_wm_test      = [175:238];                          % Unique image nrs associated with the 64 WM DOT test stimuli
                p.min_ang_distance_test_stim = 20;                                 % minimum angular distance (deg) between two wm test images; to avoid that test dot stimuli will overlap.
                
                % check if all test images for WM have unique angles
                tmp          = p.ang_deg+[0, p.delta_from_ref]';
                tmp(tmp<0)   = tmp(tmp<0)+360;
                tmp(tmp>360) = tmp(tmp>360)-360;
                assert(isequal(length(unique(tmp(:))), length(tmp(:))));
                clear tmp
                
                % Convert delta offset from polar angle coords to cartesian coords
                for dd = 1:length(p.delta_from_ref)
                    p.ang_deg_delta(dd,:)   = p.ang_deg + p.delta_from_ref(dd);    % desired angle offset for dot loc in deg (translation from center screen 0,0)
                    p.eccen_deg_delta(dd,:) = p.eccen_deg;                         % desired eccen for dot loc in deg (translation from center screen 0,0)
                    
                    [x_d,y_d] = pol2cart(deg2rad(p.ang_deg_delta(dd,:)-90),p.eccen_deg_delta(dd,:)); % again, note the -90
                    p.x0_pix_delta(dd,:) = disp_params.xc + round(x_d * disp_params.ppd);  % x-center ref dot loc in pix (translation from upper left corner [0,0])
                    p.y0_pix_delta(dd,:) = disp_params.yc + round(y_d * disp_params.ppd);  % y-center ref dot loc in pix (translation from upper left corner [0,0])
                end

                % LTM PAIR
                p.ltm_pairs          = [];                                       %%

                % IMAGERY: 
                p.unique_im_nrs_specialcore    = p.unique_im_nrs_core(1:2:p.num_loc);             % SELECTED UNIQUE IMAGES (SUBSET of all 16 dots)
                
                % IMAGERY: QUIZ DOT PARAMS
                p.imagery_sz_deg   = [disp_params.w_deg/2, disp_params.h_deg];   % desired diameter (deg) of the second, quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix   = [disp_params.xc,disp_params.h_pix];         % diameter of quiz dot image (pixels) (we already ensured even nr of pixels in display params function)
                p.unique_im_nrs_img_test = [743:902];                            % Unique image nrs associated with the 8*20=160 IMG DOT test dot images
                p.imagery_quiz_images    = [ones(1,10),2.*ones(1,10)];           % quiz dots overlap (1) or not (2)
                
                % Add params to struct
                stim.dot = p;
                
            case {'obj',4}
                
                % GENERAL
                p.indivobjfile   = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name, 'vcd_objects_2degstep_lumcorrected');
                p.stimfile       = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name, sprintf('object_%s',disp_params.name)); % mat-file where to store stimulus images?
                p.infofile       = fullfile(vcd_rootPath,'workspaces','info',sprintf('object_info_%s',disp_params.name));  % csv-file Where to find stimulus info?

                % GENERAL
                p.unique_im_nrs_core     = [65:80];                             %#ok<*NBRAK> % Unique image nrs associated with the 16 single dot stimuli
                p.iscolor = false;                                          % Use color or not? 
                if p.iscolor
                    p.square_pix_val     = true;                             % [IF YES: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]
                else
                    p.square_pix_val     = false;
                end
                
                % TEMPORAL
                p.duration       = stim.stimdur_frames;                          % frames (nr of monitor refreshes)
                
                % SPATIAL
                p.contrast      = 1;                                        % Michelson [0-1] (fraction) 
                p.x0_deg        = x0_deg;                                   % desired x-center loc in deg (translation from center screen 0,0)
                p.y0_deg        = y0_deg;                                   % desired y-center loc in deg (translation from center screen 0,0)
                p.x0_pix        = x0_pix;                                   % x-center loc in pix (translation from center screen 0,0)
                p.y0_pix        = y0_pix;                                   % y-center loc in pix (translation from center screen 0,0)
                
                p.og_res_stim_total_sz     = 1024;                          % original resolution of object stimuli
                p.og_res_stim_target_sz    = parafov_circle_diam_pix;       % object stimuli have a target size of 4 dva = 354 pixels (7TAS_BOLDSCREEN32) or 258 pixels (PP room)
                p.og_res_stim_deg          = parafov_circle_diam_deg;       % corresponding to 4 deg

                % we calculate the img size and scale factor relative to 1024x1024 original image size, 
                % as the object preprocessing script has already scaled the raw image to get desired object size. 
                p.img_sz_pix               = (p.og_res_stim_total_sz/p.og_res_stim_target_sz)*parafov_circle_diam_pix; 
                p.img_sz_deg               = parafov_circle_diam_deg;            % height (or width) of square stimulus support (deg)
                p.dres                     = p.img_sz_pix / p.og_res_stim_total_sz; % height (or width) of square stimulus support (pix)

                p.crop_img_sz_pix          = p.og_res_stim_total_sz/2; % height (and width) of cropped square stimulus support (512 x 512 pixels).

                
                % OBJECT FACING ROTATION ANGLES
                p.num_unique_objects       = 16;                            % nr of rotations (deg), 0 = rightward facing, 180 = leftward facing, 90 = forward facing
                % Constraints: 
                % * Do not include 0-25 deg and 155-180 deg to avoid edge cases in WM
                % * Have equal nr of "sideways" and "forward" facing objects.
                % * Rotations need to be at least ±5 degrees away from decision bounds (45 and 135 degrees, as well as 90 degrees)
                % * if possible, object rotations have equal distance from rotation category border (45 and 135 degrees).
                % * NOT POSSIBLE: equal distance between object rotations.
                % * We only pick even numbers because that's the granularity
                %   of rotation steps we have for raw images.
                % * We ignore 45±5 degrees and 135±5 degrees, otherwise PC task will be too difficult.
                facing_dir_deg             = cat(2, [26, 32, 36, 40], ...    4 sideways between 26-40 degrees. 
                                                    [50, 54, 70, 84], ...    4 forward between 50-84 degrees.
                                                    [96, 110, 126, 130], ... 4 forward between 96-130 degrees 
                                                    [140, 144, 148, 154]); % 4 sideways between 140-154 degrees.
                facing_dir_deg             = ceil(facing_dir_deg/2)*2; % force integrals of 2 as original images come in steps of 2 degrees

                % ensure equal distance from cardinal meridians
                assert(isequal(abs(facing_dir_deg(1:(p.num_unique_objects/2))-90), fliplr(abs(90-facing_dir_deg(((p.num_unique_objects/2)+1):p.num_unique_objects)))));
                assert(isequal(abs(45-facing_dir_deg(1:(p.num_unique_objects/2))), fliplr(abs(135-facing_dir_deg(((p.num_unique_objects/2)+1):p.num_unique_objects)))));
                
                % we carefully assign a unique rotation to an object. Less
                % than 10 degrees from decision boundary is considered a
                % "hard" trial, the others are considered "easy"
                obj_idx                    = [7,3,14, ...   70  deg: damon (easy),  36  deg: lisa (hard),  144 deg: sophia (hard)
                                              2,13,9,5, ... 32  deg: parrot (easy), 140 deg: cat (hard),    96 deg: bear (easy),   50 deg: giraffe (hard)
                                              12, 1, ...    130 deg: drill (hard),  26  deg: brush (easy)
                                              8, 11, ...    84  deg: bus (easy),    126 deg: suv (hard)
                                              4, 10, ...    40  deg: pizza (hard),  110 deg: banana (easy)
                                              16, 6,15]; %  154 deg: church (easy), 54  deg: house (hard), 148 deg: watertower (easy)
                                          
                assert(isequal([1:p.num_unique_objects],sort(obj_idx)))
                p.facing_dir_deg           = facing_dir_deg(obj_idx);
                % facing_dir_deg =
                    % 1     70  (easy) damon
                    % 2     36  (hard) lisa
                    % 3    144  (hard) sophia
                    % 4     32  (easy) parrot
                    % 5    140  (hard) cat
                    % 6     96  (easy) bear
                    % 7     50  (hard) giraffe
                    % 8    130  (hard) drill
                    % 9     26  (easy) brush
                    % 10    84  (easy) bus
                    % 11   126  (hard) suv
                    % 12    40  (hard) pizza
                    % 13   110  (easy) banana
                    % 14   154  (easy) church 
                    % 15    54  (hard) house
                    % 16   148  (easy) watertower
                
                
                % For 20% of the PC-OBJ trials, we will use a catch facing
                % direction (in the condition_master this is logged as
                % "is_objcatch" = true or false). About 2/3 of those
                % is_objcatch should be from the opposite facing direction.
                
                % We have 91 possible rotations per object ((180/2)+1). One
                % rotations will be the base rotations for the core object,
                % which leaves 91-1 = 90 possible rotations for objcatch.
                % From these 90 possible rotations, we sample 18 catch
                % rotations. (as 18 is easily divisible by 3). We then
                % sample catch rotations such that 1/3 from same facing
                % direction as the core rotation, and 2/3 from other facing 
                % direction than the core rotation.
                nr_objcatch_rotations = 90/5; % 18 possible catch rotations for each of the 16 individual objects
                p.catch_rotation      = zeros(p.num_unique_objects,nr_objcatch_rotations);  
                forward_facing_rotations(1,:)  = 0:2:44;    % degrees (note: we exclude 45 degrees)
                forward_facing_rotations(2,:)  = 136:2:180; % degrees (note: we exclude 135 degrees)
                sideways_facing_rotations(1,:) = 46:2:88;   % degrees (note: we exclude 45 and 90 degrees)
                sideways_facing_rotations(2,:) = 92:2:134;  % degrees (note: we exclude 90 and 135 degrees)
                
                % how many forward/sideways catch rotations do we sample?
                p.distribution_objcatch = round(nr_objcatch_rotations.*[1/3, 2/3]); % 1/3 from same facing direction, and 2/3 from other facing direction

                % catch rotation is in (absolute) degrees using the same
                % convention as core object base rotation: 0 = rightward facing, 180 = leftward facing, 90 = forward facing
                p.catch_rotation = [12    22    26    30    36    96   102   104   110   112   134   138   142   146   160   170   176   178; ... damon
                                    62    76    94    96    98   100   104   110   112   120   124   134   150   156   162   164   172   178; ... lisa
                                     2    14    18    24    36    38    46    56    58    60    88   106   110   114   118   126   128   130; ... sophia
                                    46    48    52    68    72    78    84   106   110   120   128   132   146   148   154   162   178   180; ... parrot
                                     2     6    20    26    32    44    56    64    80    86    88    98   100   104   112   114   126   132; ... cat
                                     0     2     8    12    20    22    26    30    36    46    48    50    68    74    82   140   158   164; ... bear
                                     6    18    20    22    32    94   100   110   118   124   132   138   142   152   158   168   170   178; ... giraffe
                                     2     6    18    20    26    28    36    42    52    56    60    74    80    82   136   142   160   170; ... drill
                                    46    48    82    88    96   104   110   112   122   124   128   134   138   142   144   148   164   172; ... brush
                                     6    12    18    20    30    40   100   102   104   114   116   132   154   156   162   166   176   178; ... bus
                                     2    18    32    34    44    60    64    76    78    82    88   136   138   140   142   154   166   174; ... suv
                                    58    66    70    94    96   102   104   108   110   120   122   132   136   146   156   160   162   166; ... pizza
                                     4    12    32    34    50    52    66    72    86    88   138   142   150   154   158   168   174   178; ... banana
                                     2     4    16    32    34    36    48    52    54    60    62    72    76    86    92   106   108   110; ... church
                                     0     6    12    22    28    30    32    38    44   102   106   110   112   122   130   142   150   156; ... house
                                     0    10    14    34    38    42    50    56    62    64    76    88    94   116   118   124   132   134]; % watertower
                                         
               p.unique_im_nrs_objcatch = [1423:(1423+length(p.catch_rotation(:))-1)]; % unique image nrs for object catch rotations are 1423:1710;
                           
                 % Above sampling of catch rotations was created by the following code:
                
                % % %      for ff = 1:p.num_unique_objects
                % % %          % if the base rotation is forward facing
                % % %          if ismember(p.facing_dir_deg(ff),sideways_facing_rotations)
                % % %              % and in the first row of forward rotations
                % % %              if p.facing_dir_deg(ff) <= max(sideways_facing_rotations(1,:))
                % % %                  sampleSW_same = setdiff(sideways_facing_rotations(1,:),p.facing_dir_deg(ff));
                % % %                  sampleSW_diff = sideways_facing_rotations(2,:);
                % % %                  % or the second first row of forward rotations
                % % %              elseif p.facing_dir_deg(ff) >= max(sideways_facing_rotations(1,:)) && p.facing_dir_deg(ff) <= max(sideways_facing_rotations(2,:))
                % % %                  sampleSW_same = setdiff(sideways_facing_rotations(2,:),p.facing_dir_deg(ff));
                % % %                  sampleSW_diff = sideways_facing_rotations(1,:);
                % % %              end
                % % %              % we sample all forward facing rotations
                % % %              sampleFF = forward_facing_rotations;
                % % % 
                % % %              % select 1/3 sideways rotations that are 90 degrees
                % % %              % away and 2/3 of facing forward rotations for
                % % %              % objcatch
                % % %              p.catch_rotation(ff,:) = sort(cat(2, datasample(sampleSW_diff,p.distribution_objcatch(1),'Replace',false), datasample(sampleFF(:),p.distribution_objcatch(2),'Replace',false)'));
                % % % 
                % % %              % if the base rotation is sideways facing
                % % %          elseif ismember(p.facing_dir_deg(ff),forward_facing_rotations)
                % % %              % and in the first row of sideways rotations
                % % %              if p.facing_dir_deg(ff) <= max(forward_facing_rotations(1,:))
                % % %                  sampleFF_same = setdiff(forward_facing_rotations(1,:),p.facing_dir_deg(ff));
                % % %                  sampleFF_diff = forward_facing_rotations(2,:);
                % % %                  % or the second first row of sideways rotations
                % % %              elseif p.facing_dir_deg(ff) >= max(forward_facing_rotations(1,:)) && p.facing_dir_deg(ff) <= max(forward_facing_rotations(2,:))
                % % %                  sampleFF_same = setdiff(forward_facing_rotations(2,:),p.facing_dir_deg(ff));
                % % %                  sampleFF_diff = forward_facing_rotations(1,:);
                % % %              end
                % % %              % we sample all sideways facing rotations
                % % %              sampleSW = sideways_facing_rotations(:);
                % % % 
                % % %              % select 1/3 sideways rotations that are 90 degrees
                % % %              % away and 2/3 of facing forward rotations for
                % % %              % objcatch
                % % %              p.catch_rotation(ff,:) = sort(cat(2, datasample(sampleFF_diff,p.distribution_objcatch(1),'Replace',false), datasample(sampleSW(:),p.distribution_objcatch(2),'Replace',false)'));
                % % % 
                % % %          end
                % % % 
                % % %          % ensure the catch rotations are not the same as the core
                % % %          % rotations
                % % %          assert(all(p.catch_rotation(ff,:)~=p.facing_dir_deg(ff)))
                % % %      end

                % Define the 5 superordinate, 1-3 basic, and 16 subordinate categories
                p.super_cat          = {'human','animal','object','food','place'};     
                
                p.basic_cat{1}       = {'facemale','facefemale','facefemale'};
                p.basic_cat{2}       = {'small','small','large','large'};
                p.basic_cat{3}       = {'tool','tool','vehicle','vehicle'};
                p.basic_cat{4}       = {'manmade','produce'};
                p.basic_cat{5}       = {'building','building','building'};
                
                p.sub_cat{1}         = {'damon','lisa','sophia'};
                p.sub_cat{2}         = {'parrot','cat','bear','giraffe'};
                p.sub_cat{3}         = {'drill','brush','bus','suv'};
                p.sub_cat{4}         = {'pizza','banana'};
                p.sub_cat{5}         = {'church','house','watertower'};
                
                % Define the affordances for each object subcategory (for HOW task) 
                p.affordance{1}      = {'greet','greet','greet'};
                p.affordance{2}      = {'greet','greet','observe','observe'};
                p.affordance{3}      = {'grasp','grasp','enter','enter'};
                p.affordance{4}      = {'grasp','grasp'};
                p.affordance{5}      = {'enter','enter','observe'};

                % WORKING MEMORY: dot angle position deltas for test
                % images.
                % NOTE that WM have do not have "unique" rotations (i.e.,
                % no other core object has that rotation, given the limited 
                % nr of rotations we can use, and the use of objcatch
                % rotations).
                p.delta_from_ref         = [-24, -12, 12, 24];                  % Relative rotation from reference image for WM test image
                                                                                %  the bigger the delta, the easier the trial. 
                                                                                % for 1-89 deg rotations: Negative values are leftwards, positive values is rightwards
                                                                                % for 91-180 deg rotations: Negative values are rightward, positive values is leftward
                p.unique_im_nrs_wm_test  = [239:302];                           %  Unique image nrs associated with the 64 WM OBJ test stimuli
                                                                            
                                   
                % IMAGERY: 
                p.unique_im_nrs_specialcore = p.unique_im_nrs_core([1,3,5,7,8,10,12,14]);   % 8 SELECTED IMAGES USED (SUBSET of 16 IMAGES)--these are hand picked! damon, sophia, cat, giraffe, drill, bus, pizza, church
                
                % IMAGERY QUIZ DOT PARAMS
                p.imagery_sz_deg         = 5.658;                              % desired diameter (degree) of the second, quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix         = (round(p.img_sz_deg* disp_params.ppd)/2)*2; % pixel diameter of quiz dot image (ensure even nr of pixels)
                p.unique_im_nrs_img_test = [903:1062];                         % Unique image nrs associated with the 8*20=160 IMG OBJ test dot images
                p.imagery_quiz_images    = [ones(1,10),2.*ones(1,10)];         % quiz dots overlap (1) or not (2)
                
                % LTM PAIR
                p.ltm_pairs          = [];                                     %% TODO: list of core stim numbers of same/other stim classes in the order of obj core images
                
                % Add params to struct
                stim.obj = p;
                
            case {'ns',5}
                
                % GENERAL
                p.core_png_folder    = fullfile(vcd_rootPath,'workspaces','stimuli','RAW','vcd_natural_scenes'); % folder, where to find core scene pngs?
                p.wmtest_png_folder  = fullfile(p.core_png_folder, 'wm_test');                                   % folder, where to find individual wm test pngs?
                p.ltmlure_png_folder = fullfile(p.core_png_folder, 'ltm_novel_lures');                           % folder, where to find individual novel ltm lure pngs?
                p.imgtest_png_folder = fullfile(p.core_png_folder, 'img_test');                                  % folder, where to find individual imagery test pngs?
                
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('scene_%s',disp_params.name)); % prefix to mat file of preprocessed stimulus images
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('scene_info_%s',disp_params.name));                % csv-file with stimulus info
                
                p.unique_im_nrs_core      = [81:110];                                                         % Unique image nrs associated with the 30 core scene stimuli

                p.iscolor = true;                                                                             % Use color or not? 
                if p.iscolor, p.square_pix_val = true;                                                        % [IF YES: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]                                                                
                else, p.square_pix_val = false; end

                % TEMPORAL
                p.duration       = stim.stimdur_frames;                                                            % frames (nr of monitor refreshes)
                
                % SPATIAL
                % For reference: we resize scenes to 8.4 x 8.4 degrees (as this is what has been used in NSD) 
                % NSD scene resolution with the old 7T screen setup resulted in 714 x 714 pixels
                % For VCD, 8.4 x 8.4 degrees  results in: 
                % PP room Eizoflexscan: 537 x 537 pixels, or 8.4034 x 8.4034 degrees.
                % 7TAS BOLDSCREEN: XX x XX pixels, or XX x XX degrees.
                p.og_res_stim    = 425;                                                                       % Original stimulus reslution is 425 x 425 (pixels)
                p.img_sz_deg     = ctr_square_deg;                                                            % desired height (or width) of square stimulus support (deg)
                p.img_sz_pix     = round(p.img_sz_deg.*disp_params.ppd);                                      % height (and width) of square stimulus support (pixels). 
                p.dres           = p.img_sz_pix/p.og_res_stim;                                                % Scale factor for imresize
                p.x0_deg         = 0;                                                                         % desired x-center loc in degrees (translation from 0,0)
                p.y0_deg         = 0;                                                                         % desired y-center loc in degrees (translation from 0,0)
                p.x0_pix         = round((p.x0_deg.*disp_params.ppd)/2)*2;                                    % x-center loc in pixels (translation from 0,0)
                p.y0_pix         = round((p.y0_deg.*disp_params.ppd)/2)*2;                                    % y-center loc in pixels (translation from 0,0)
                
                % SEMANTIC CATEGORY INFORMATION
                
                % 5 superordinate semantic object categories
                p.super_cat     = {'human','animal','object','food','place'};      
                
                %  x 2 basic semantic object categories (scene location)
                p.basic_cat{1}  = {'indoor','outdoor'};
                p.basic_cat{2}  = {'indoor','outdoor'};
                p.basic_cat{3}  = {'indoor','outdoor'};
                p.basic_cat{4}  = {'indoor','outdoor'};
                p.basic_cat{5}  = {'indoor','outdoor'};
                
                %  x 3 subordinate categories (object spatial location)
                p.sub_cat{1,1}  = {'face1_left','face2_center','face3_right'};
                p.sub_cat{1,2}  = {'face1_left','face2_center','face3_right'};
                
                p.sub_cat{2,1}  = {'cat1_left','cat2_center','cat3_right'};
                p.sub_cat{2,2}  = {'giraffe1_left','giraffe2_center','giraffe3_right'};
               
                p.sub_cat{3,1}  = {'vase1_left','vase2_center','vase3_right'};
                p.sub_cat{3,2}  = {'bus1_left','bus2_center','bus3_right'};
                
                p.sub_cat{4,1}  = {'donut1_left','donut2_center','donut3_right'};
                p.sub_cat{4,2}  = {'banana1_left','banana2_center','banana3_right'};
                
                p.sub_cat{5,1}  = {'bathroom1_left','bathroom2_center','bathroom3_right'};
                p.sub_cat{5,2}  = {'building1_left','building2_center','building3_right'};
                
                % 4 affordance categories
                p.affordance{1,1} = {'observe','greet','greet'}; % people indoor
                p.affordance{1,2} = {'walk','greet','greet'};    % people outdoor
                p.affordance{2,1} = {'greet','greet','greet'};   % cats
                p.affordance{2,2} = {'observe','observe','observe'}; % giraffes
                p.affordance{3,1} = {'grasp','grasp','grasp'}; % object indoor
                p.affordance{3,2} = {'walk','observe','walk'}; % object outdoor
                p.affordance{4,1} = {'grasp','grasp','grasp'}; % food indoor
                p.affordance{4,2} = {'grasp','grasp','observe'}; % food outdoor
                p.affordance{5,1} = {'walk','walk','walk'}; % bathrooms indoor
                p.affordance{5,2} = {'observe','walk','walk'};  % streets outdoor

                % FOR WM task crossing, we have manipulated the original
                % NSD image by adding or removing something in the image.
                % These changes can be obvious (easy) or subtle (hard) to
                % detect:
                p.change_im                 = [-2,-1,1,2]; % where ±2 = easy, ±1 = hard, (-) = remove (+) = add.
                p.change_im_name            = {'easy_remove','hard_remove','hard_add','easy_add'};
                p.unique_im_nrs_wm_test     = [303:422];                               % Unique image nrs associated with the 120 WM NS changed stimuli
                
                % IMAGERY 
                p.unique_im_nrs_specialcore = p.unique_im_nrs_core([2,4,5,8,10,11,13,15,18,20,21,23,26,27,30]); % Half of the images will be used for  IMG/LTM pairing (these are carefully handpicked! see scene_info csv file)
                p.imagery_sz_deg            = p.img_sz_deg;                                                     % desired diameter (degree) of the second, quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix            = p.img_sz_pix;                                                     % diameter of quiz dot image (pixels) (we already ensured even nr of pixels)
                p.unique_im_nrs_img_test    = [1063:1362];                                                      % Unique image nrs associated with the 15*20=300 IMG NS test dot images
                p.imagery_quiz_images       = [ones(1,10),2.*ones(1,10)];                                       % Quiz dots overlap (1) or not (2)
                
                % FOR LTM incorrect trials, we have very similar looking images called "lures":
                p.lure_im                   = {'lure01', 'lure02', 'lure03', 'lure04'};
                p.unique_im_nrs_ltm_lures   = [1363:1422];                                                      % Unique image nrs associated with the 15*4=60 WM NS lure images
                
                % LTM PAIRED ASSOCIATES
                p.ltm_pairs                 = [];

                % Add params to struct
                stim.ns = p;
        end
        if verbose
            % Print out stimulus params
            fprintf('*** %s:\tstimulus size = %2.2f deg (%3.2f pixels) ***\n', upper(stim_class{ii}), p.img_sz_deg, p.img_sz_pix);
        end
    end
    
    stim.all_core_im_nrs         = sort(cat(2, stim.gabor.unique_im_nrs_core, stim.rdk.unique_im_nrs_core, stim.dot.unique_im_nrs_core, stim.obj.unique_im_nrs_core, stim.ns.unique_im_nrs_core));
    stim.all_specialcore_im_nrs  = sort(cat(2, stim.gabor.unique_im_nrs_specialcore, stim.rdk.unique_im_nrs_specialcore, stim.dot.unique_im_nrs_specialcore, stim.obj.unique_im_nrs_specialcore, stim.ns.unique_im_nrs_specialcore));
    stim.all_wm_test_im_nrs      = sort(cat(2, stim.gabor.unique_im_nrs_wm_test, stim.rdk.unique_im_nrs_wm_test, stim.dot.unique_im_nrs_wm_test, stim.obj.unique_im_nrs_wm_test, stim.ns.unique_im_nrs_wm_test));
    stim.all_ltm_pairs           = sort(cat(2, stim.gabor.ltm_pairs, stim.rdk.ltm_pairs, stim.dot.ltm_pairs, stim.obj.ltm_pairs, stim.ns.ltm_pairs));
    stim.all_img_test_im_nrs     = sort(cat(2, stim.gabor.unique_im_nrs_img_test, stim.rdk.unique_im_nrs_img_test, stim.dot.unique_im_nrs_img_test, stim.obj.unique_im_nrs_img_test, stim.ns.unique_im_nrs_img_test));
    stim.all_ltm_lure_im_nrs     = sort(stim.ns.unique_im_nrs_ltm_lures);
    stim.all_objectcatch_im_nrs  = sort(stim.obj.unique_im_nrs_objcatch);
    stim.all_test_im_nrs         = sort(cat(2, stim.all_wm_test_im_nrs, stim.all_img_test_im_nrs, stim.all_ltm_lure_im_nrs, stim.all_objectcatch_im_nrs));
    stim.all_im_nrs              = sort(cat(2, stim.all_core_im_nrs, stim.all_test_im_nrs));

    % Tell the user more info
    if verbose
        fprintf('*** ALL STIMULI: Using a stimulus duration of %d frames (%3.2f seconds), where one frame relates to %d monitor refreshes (%d Hz) ***\n', ...
            p.duration, p.duration*stim.framedur_s, stim.f2f, disp_params.refresh_hz);
        fprintf('*** Number of stimuli: \n')
        fprintf('\t %d total (%03d-%03d) \n',length(stim.all_im_nrs),min(stim.all_im_nrs), max(stim.all_im_nrs)) 
        fprintf('\t %d core (%03d-%03d) \n',length(stim.all_core_im_nrs),min(stim.all_core_im_nrs), max(stim.all_core_im_nrs)) 
        fprintf('\t %d WM test (%03d-%03d) \n',length(stim.all_wm_test_im_nrs),min(stim.all_wm_test_im_nrs), max(stim.all_wm_test_im_nrs))
        fprintf('\t %d IMG test (%03d-%03d) \n',length(stim.all_img_test_im_nrs),min(stim.all_img_test_im_nrs), max(stim.all_img_test_im_nrs))
        fprintf('\t %d LTM novel lures (%03d-%03d) \n',length(stim.all_ltm_lure_im_nrs),min(stim.all_ltm_lure_im_nrs), max(stim.all_ltm_lure_im_nrs))
        fprintf('\t %d OBJ catch images (%03d-%03d) \n', length(stim.all_objectcatch_im_nrs), min(stim.all_objectcatch_im_nrs), max(stim.all_objectcatch_im_nrs));
        fprintf('\t %d special core\n',length(stim.all_specialcore_im_nrs))
        fprintf('\t %d LTM pairs \n',length(stim.all_ltm_pairs)) 
        
    end
    
    
    % do some checks
    assert(isequal(1:length(stim.all_core_im_nrs),stim.all_core_im_nrs));
    assert(all(ismember(stim.all_specialcore_im_nrs,stim.all_core_im_nrs)));
%     assert(all(ismember(stim.all_ltm_pairs,stim.all_core_im_nrs)));
    assert(all(~ismember(stim.all_wm_test_im_nrs,stim.all_core_im_nrs)));
    assert(all(~ismember(stim.all_img_test_im_nrs,stim.all_core_im_nrs)));
    assert(all(~ismember(stim.all_ltm_lure_im_nrs,stim.all_core_im_nrs)));

    
    % Store params if requested
    if store_params
        if verbose, fprintf('[%s]: Storing params..\n',mfilename); end
        if ~exist(save_info_dir,'dir'); mkdir(save_info_dir); end
        save(fullfile(save_info_dir,sprintf('stim_%s_%s.mat',disp_params.name,datestr(now,30))),'stim')
    end
end


return