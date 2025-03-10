function stim = vcd_getStimParams(type, disp_params,load_params,store_params, overwrite_randomized_params)

if ~exist('type','var') || isempty(type)
    type = {'all'};
end

if ischar(type)
    type = {type};
end

if any(strcmp(type{:},'all'))
    type = {'gabor','rdk','dot','cobj','ns'};
end

if ~exist('disp_params','var') || isempty(store_params) || ~isfield(disp_params, 'refresh_hz')
    disp_name = '7TAS_BOLDSCREEN32'; % or 'KKOFFICEQ3277' or psychophys room??
    disp_params = vcd_getDisplayParams(disp_name);
end

if ~exist('store_params','var') || isempty(store_params)
    store_params = true;
end

if ~exist('load_params','var') || isempty(load_params)
    load_params = true;
end

if ~exist('overwrite_randomized_params','var') || isempty(overwrite_randomized_params)
    overwrite_randomized_params = false;
end

if load_params
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('stim_%s*.mat',disp_params.name)));
    if ~isempty(d)
        if length(d) > 1
            warning('[%s]: Multiple .mat files! Will pick the most recent one', mfilename);
        end
        load(fullfile(d(end).folder,d(end).name),'stim');
    else
        error('[%s]: Can''t find stim params file!', mfilename)
    end
else
    
    
    % Setup struct
    stim = struct();
    
    %% GENERAL PARAMS
    
    stim.store_imgs                 = true;                                     % Store images when creating them? (see s_createStimuli.m)
    
    % we want to present images at a rate of 30 Hz, so for a monitor
    % refresh rate of 120 Hz, this number is 4, for 60 Hz this is 2. 
    stim.presentationrate_hz        = 30;                                       % rate of frame presentation (Hz)
    
    % ensure we have a round number of monitor refreshes to achieve 30 Hz presentation rate  
    assert(isint(disp_params.refresh_hz/stim.presentationrate_hz));
    
    stim.framedur_s                 = (disp_params.refresh_hz/stim.presentationrate_hz) * (1/disp_params.refresh_hz); % frames per second (seconds)
    stim.bckgrnd_grayval            = uint8(1+255/2);                           % background color (middle of 1-255 pixel lum)
    stim.scfactor                   = 1;                                        % no scaling
    
    % Default params to inherit (or overwrite)
    stimdur_frames                  = 60;                                       % nr of (33 ms) presentation frames (60 frames = 2 sec)
    x0_deg                          = [-4 4];                                   % desired x-center location for left right stim apertures (degrees)
    y0_deg                          = [0 0];                                    % desired y-center location for left right stim apertures (degrees)
    
    x0_pix                          = round((x0_deg.*disp_params.ppd)/2)*2;       % 354 pixels (4.0059 deg) for BOLD screen, empirical x-center location for left right stim apertures (pixels)
    y0_pix                          = round((y0_deg.*disp_params.ppd)/2)*2;       % 354 pixels (4.0059 deg) for BOLD screen, empirical y-center location for left right stim apertures (pixels)
    
    % Parafoveal circular aperture (for gabors, rdk, and cobj)
    parafov_circle_diam_deg         = 4;                                        % desired  parafoveal circular diameter aperture (empirical is 4.0059 degrees)
    parafov_circle_diam_pix         = round((parafov_circle_diam_deg*disp_params.ppd)/2)*2;  % parafoveal circular diameter aperture (pixels), 354 pixels
    
    % Center square image (for ns)
    ctr_square_deg                  = 8.4;                                       % desired center square side length in degrees, same as NSD. Empirical is 8.3965 degrees.
    ctr_square_pix                  = round((ctr_square_deg*disp_params.ppd)/2)*2; % center square side length in pixels (742 pixels)
    
    
    %% NOISE BACKGROUND Puzzle piece
    % general
    stim.bckground.stimfile         = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('bckgrnd_%s',disp_params.name)); % mat file
    % stim.bckground.infofile         = fullfile(vcd_rootPath,'workspaces','info',sprintf('bckgrnd_info_%s',disp_params.name)); % csv file
    
    stim.bckground.alpha            = 1; %<alpha> (optional) is the exponent to apply to the amplitude spectrum (i.e. 1/f^alpha).  default: 1.
    stim.bckground.num              = 1; %<num> (optional) is the number of images desired.  default: 1.
    stim.bckground.mode             = 0; %<mode> (optional) is means fixed amplitude spectrum + random phase
    stim.bckground.std_clip_range   = 3.5; % std of image values, that defines the range we use to rescale image to 1 255.
    
    %% FIXATION DOT
    
    % general
    stim.fix.stimfile               = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('fix_%s',disp_params.name)); % mat file
    stim.fix.infofile               = fullfile(vcd_rootPath,'workspaces','info',sprintf('fix_info_%s',disp_params.name)); % csv file
    
    % TEMPORAL
    stim.fix.dotmeanchange          = 90;                                   % (nr 33ms frames), on average, dot changes occur every 3 seconds
    stim.fix.dotchangeplusminus     = 6;                                    % (nr 33ms frames), earliest and latests time that dot changes occur. 2 seconds means [-1:1] from meanchange
    stim.fix.dres                   = [];                                   % rescale factor (fraction)
        
    % SPATIAL
    stim.fix.dotcenterdiam_pix      = 18;                                   % dot diameter in pixels (18 pixels)
    stim.fix.dotcenterdiam_deg      = stim.fix.dotcenterdiam_pix/disp_params.ppd; % 0.2037 deg, dot diameter
    stim.fix.dotthinborderdiam_pix  = stim.fix.dotcenterdiam_pix+6;         % pixel-width for dot thin border (during ITI/IBI)
    stim.fix.dotthickborderdiam_pix  = stim.fix.dotcenterdiam_pix+10;       % pixel-width for dot border (during trial)
    
    stim.fix.lumminmaxstep          = [42,212,5];                       % min and max luminance of dot [1-255],
    stim.fix.dotlum                 = uint8(linspace(stim.fix.lumminmaxstep(1),stim.fix.lumminmaxstep(2),stim.fix.lumminmaxstep(3))); % dot gray levels
    stim.fix.dotopacity             = 0.5;                              % dot and border have 50% opacity
    
    fprintf('*** FIXATION MARK: fixation center diam = %d, thick rim = %d, thick rim = %d pixels ***\n', ...
        stim.fix.dotcenterdiam_pix,stim.fix.dotthinborderdiam_pix,stim.fix.dotthickborderdiam_pix);
    
    %% CONTRAST DECREMENT -- INVERTED GAUSSIAN TEMPORAL WINDOW
    % We treat timepoint t=0 as the onset time of a frame flip, and t=1 is 
    % the offset of the first frame flip and the onset of the second frame 
    % flip (so in between flips). Each time point has a duration of 33.33 ms 
    % (see framedur_s). 
    % For the contrast decrement, we implement an inverted gaussian with a
    % support of 15 time points (a 1/4 of the total stimulus duration),
    % thus t=[0,14]. How quickly the contrast decrement change occurs is 
    % determined by the std. Here, we pick std = 3 time points (3*33.33 ms). 
    % This means that peak contrast decrement occurs at t=7 (233.33 ms).
    stim.cd.t_gausswin_N            = 15;                                   % 15 number of timepoints for gaussian time window (contrast decrement)
    stim.cd.t_gausswin_std          = 3;                                    % standard devation of gaussian window in time (frames)
    stim.cd.meanchange              = 1/stim.framedur_s;                     % mean of gaussian window in time (30 frames = 1 sec)  
    stim.cd.changeplusminus         = (0.5/stim.framedur_s)-1;               % plus or minus this amount (14 frames = 0.46 sec)  
    stim.cd.max_cd                  = 0.5;                                  % we reduce contrast by 50% of mean luminance
    
    % Create 1D gaussian
    t_support           = linspace(-stim.cd.t_gausswin_N / 2, stim.cd.t_gausswin_N / 2, stim.cd.t_gausswin_N);
    t_gauss             = exp(-t_support .^ 2 / (2 * stim.cd.t_gausswin_std ^ 2));
    
    % invert and scale it
    t_gauss             = 1-(t_gauss*stim.cd.max_cd);
    stim.cd.t_gauss     = t_gauss;

    %% EYELINK PARAMS
    % we manually place the eyelink calibration/validation points on the
    % display. The distance between the center and 4 left/right/up/down
    % points are set as [xc,yc] ± 350 pixels. This results in dots at
    % [x1,y1]=[610,0],
    % [x2,y2]=[1310,0],
    % [x3,y3]=[0,190],
    % [x4,y4]=[0,890] pixels
    stim.el.point2point_distance = 354; % pixels (4.0059 deg for BOLDscreen)

    %% STIM PARAMS
    for ii = 1:length(type)
        
        p = [];
        
        switch type{ii}
            
            case {'gabor',1}
                
                % Where to store stimulus images?
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('gabor_%s',disp_params.name)); % mat file
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('gabor_info_%s',disp_params.name)); % csv file
                
                % TEMPORAL
                % Fixed params
                p.duration        = stimdur_frames;                                    % frames (nr of monitor refreshes)

                % SPATIAL
                % Fixed params
                p.img_sz_deg      = parafov_circle_diam_deg;                    % height (and width) of stimulus support (deg)
                p.img_sz_pix      = parafov_circle_diam_pix;                    % height (and width) of square stimulus support (pix)
                p.og_res_stim     = p.img_sz_pix;                               % resolution of stored dot stimuli (in pixels)
                p.dres            = (( (p.img_sz_pix/disp_params.ppd) /disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                p.iscolor         = false;                                      % use color or not [[[IF WE USE COLOR: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]]
                p.square_pix_val  = false;
                
                p.gauss_std_deg   = 0.5;                                        % desired standard deviation of gaussian window in degrees, empirical is 0.4979 deg for BOLD screen.
                p.gauss_std_pix   = round((p.gauss_std_deg * disp_params.ppd)/2)*2; % standard deviation of gaussian window in pixels.
                
                p.sf_cpd          = 4;                                          % spatial frequency (cycles/deg)
                p.cycles_per_pix  = (p.sf_cpd/p.img_sz_deg)/disp_params.ppd;           % nr of cycles per image (pix)
                
                p.x0_deg          = x0_deg;                                     % x-center loc in deg (translation from 0,0)
                p.y0_deg          = y0_deg;                                     % y-center loc in deg (translation from 0,0)
                p.x0_pix          = x0_pix;                                     % x-center loc in pix (translation from 0,0)
                p.y0_pix          = y0_pix;                                     % y-center loc in pix (translation from 0,0)
                
                p.ph_deg          = [0:(180/4):179];                            % 4 quadrature Gabor phases
                
                % Manipulated params
                p.contrast        = [0.05, 0.10, 1];                            % Michelson contrasts [0-1] (fraction)
                p.ori_jitter_sd   = 2;                                          % std of normal distribution to sample orientation jitter
                p.ori_jitter_mu   = 1;                                          % mean of normal distribution to sample orientation jitter
                p.num_ori         = 8;                                          % orientation "bins" from which we create final gabor orientations (deg), 0 = 12 o'clock
                p.ori_bins        = [10:(180/p.num_ori):170];                   % rotate 10 deg away from vertical, to avoid ill-defined response options
                                                                                % [10 32.5000 55 77.5000 100 122.5000 145 167.5000]
                if overwrite_randomized_params 
                    p.ori_jitter  = p.ori_jitter_mu + (p.ori_jitter_sd.*randn(1,p.num_ori)); % add a small amount of jitter around the orientation
                    p.ori_deg     = round(p.ori_bins + p.ori_jitter);           % final gabor orientations
                else
                    p.ori_jitter  = [2.0753 4.6678 -3.5177 2.7243 1.6375 -1.6154 0.1328 1.6852];
                    p.ori_deg     = [12 37 51 80 102 121 145 169];           % final gabor orientations
                end
                
                p.delta_from_ref  = [-15, -5, 5, 15];                           % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                % the bigger the delta, the easier the trial
                % Add params to struct
                stim.gabor = p;
                
                
            case {'rdk',2}
                % Where to store stimulus images?
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('rdk_%s',disp_params.name)); % mat file
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('rdk_info_%s',disp_params.name)); % csv file
                
                % TEMPORAL
                p.duration        = stimdur_frames;                                   % frames (nr of monitor refreshes)
                p.dots_coherence  = [0.064, 0.128, 0.512];                     % fraction of dots. Kiani lab uses usually one of these [0 0.032 0.064 0.128 0.256 0.512]
                p.dots_speed      = 5;                                         % pixels/frame? Kiani lab uses usually 5 to 10
                p.dots_interval   = 1;                                         % framedur_s interval by which dots update (so 30 disp.framedur_s / 1 interval = 30 frames/sec) 
                                                                               % currently set to 30 frames per second to approx Kiani's 75 hz refresh rate / 3 interval = 25 frames/sec
                
                % SPATIAL
                p.img_sz_deg      = parafov_circle_diam_deg;                    % stimulus aperture diameter (deg)
                p.img_sz_pix      = parafov_circle_diam_pix;                    % stimulus aperture diameter (pix)
                p.og_res_stim     = p.img_sz_pix;                               % resolution of stored dot stimuli (in pixels)
                p.dres            = (( (p.img_sz_pix/disp_params.ppd) /disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                p.iscolor         = false;                                      % use color or not [[[IF WE USE COLOR: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]]
                p.square_pix_val  = false;

                p.num_mot_dir      = 8;
                p.motdir_bins      = [0:(360/p.num_mot_dir):359]+10;             % sample direction of coherent motion from [0-359] in deg (0 deg is aligned with 12 o'clock)
                                                                                 % [10 55 100 145 190 235 280 325]
                p.motdir_jitter_sd = 2;                                          % std of normal distribution to sample orientation jitter
                p.motdir_jitter_mu = 1;                                          % mean of normal distribution to sample orientation jitter
                
                if overwrite_randomized_params 
                    p.motdir_jitter    = p.motdir_jitter_mu + (p.motdir_jitter_sd.*randn(1,p.num_mot_dir)); % add a small amount of jitter around the orientation
                    p.dots_direction   = round(p.motdir_bins + p.motdir_jitter);     % final gabor orientations
                else
                    p.motdir_jitter    = [8.1568 6.5389 -1.6998 7.0698 2.4508 0.8739 2.4295 0.5901]; 
                    p.dots_direction   = [18 62 98 152 192 236 282 326];
                end
                
                % RDK specific
                p.dots_density     = 16.7;                                      % dots/deg??
                p.dots_size        = 3;                                         % radius in pixels??
                p.dots_color       = [255 255 255; 1 1 1]./255;                  % 50:50 white:black, color in RGB [0-1]
                p.max_dots_per_frame = 200;                                     % from Kiani lab (roughly matches to nr of pixels in aperture)
                p.dots_contrast    = 1;                                         % Michelson [0-1] (fraction)
                
                p.x0_deg          = x0_deg;                                     % x-center loc in deg (translation from 0,0)
                p.y0_deg          = x0_deg;                                     % y-center loc in deg (translation from 0,0)
                p.x0_pix          = x0_pix;                                     % x-center loc in pix (translation from 0,0)
                p.y0_pix          = y0_pix;                                     % y-center loc in pix (translation from 0,0)
                
                p.delta_from_ref  = [-15, -5, 5, 15];                           % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                
                % Add params to struct
                stim.rdk = p;
                
            case {'dot',3}
                % Where to store stimulus images?
                p.stimfile        = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('dot_%s',disp_params.name)); % mat file
                p.infofile        = fullfile(vcd_rootPath,'workspaces','info',sprintf('dot_info_%s',disp_params.name)); % csv file
                p.iscolor         = false;                                      % use color or not [[[IF WE USE COLOR: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]]
                
                % TEMPORAL
                p.duration        = stimdur_frames;                                   % frames (nr of monitor refreshes)
                
                % SPATIAL
                p.img_sz_deg      = 1.1;                                     	 % desired spatial support of dot image in deg, empirical is 1.1090 deg for BOLD screen
                p.img_sz_pix      = round((p.img_sz_deg * disp_params.ppd)/2)*2; % spatial support of dot image in pixels (98 pixels)
                p.og_res_stim     = p.img_sz_pix;                                % resolution of stored dot stimuli (in pixels)
                p.dres            = (( (p.img_sz_pix/disp_params.ppd) /disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                p.radius_deg      = 0.5;                                         % desired dot radius in deg (empirical is 0.4979 for BOLDscreen)
                p.radius_pix      = round((p.radius_deg * disp_params.ppd)/2)*2; % dot radius in pix (44 pixels for BOLDscreen)
                p.color           = [255 255 255];                               % white in uint8 RGB
                p.contrast        = 1;                                           % Michelson [0-1] (fraction)
                p.square_pix_val  = false;

                p.num_loc         = 16;                                          % orientation "bins" from which we create final gabor orientations (deg), 0 = 12 o'clock
                p.loc_jitter_sd   = 1;                                           % std of normal distribution to sample orientation jitter
                p.loc_jitter_mu   = 1;                                           % mean of normal distribution to sample orientation jitter
                p.loc_bins        = [3:(180/p.num_loc):176];                     % rotate 10 deg away from vertical, to avoid ill-defined response options
                
                if overwrite_randomized_params 
                    p.loc_jitter      = p.loc_jitter_sd + (p.loc_jitter_mu.*randn(1,p.num_loc)); % add a small amount of jitter around the orientation
                    p.loc_deg         = round(p.loc_bins + p.loc_jitter) + 90;   % final gabor orientations to make North 0;
                else
                    p.loc_jitter      = [0.8759,2.4897,2.4090,2.4172,1.6715,-0.2075,1.7172,2.6302,1.4889, 2.0347,1.7269,0.6966,1.2939,0.2127,1.8884,-0.1471];
                    p.loc_deg         = [4 17 28 39 50 59 72 84 94 106 117 127 139 149 162 172]; 
                end
                
                p.iso_eccen       = 4.5;
                [x,y] = pol2cart(deg2rad(p.loc_deg),repmat(p.iso_eccen,1,length(p.loc_deg)));
                p.x0_deg          = x;                                        % x-center loc in deg (translation from 0,0)
                p.y0_deg          = y;                                        % y-center loc in deg (translation from 0,0)
                p.x0_pix          = round(p.x0_deg * disp_params.ppd);        % x-center loc in pix (translation from 0,0)
                p.y0_pix          = round(p.y0_deg * disp_params.ppd);        % y-center loc in pix (translation from 0,0)
                                
                p.delta_from_ref  = [-15, -5, 5, 15];                           % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                % the bigger the delta, the easier the trial
                % Add params to struct
                stim.dot = p;
                
            case {'cobj',4}
                % GENERAL
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('object_%s',disp_params.name)); % mat-file Where to store stimulus images?
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('object_info_%s',disp_params.name));  % csv-file Where to find stimulus info?
                p.iscolor         = false;                                      % use color or not [[[IF WE USE COLOR: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]]
                
                
                % TEMPORAL
                p.duration      = stimdur_frames;                                     % frames (nr of monitor refreshes)
                
                % SPATIAL
                p.contrast      = 1;                                            % Michelson [0-1] (fraction)
                p.x0_deg        = x0_deg;                                       % x-center loc in deg (translation from 0,0)
                p.y0_deg        = y0_deg;                                       % y-center loc in deg (translation from 0,0)
                p.x0_pix        = x0_pix;                                       % x-center loc in pix (translation from 0,0)
                p.y0_pix        = y0_pix;                                       % y-center loc in pix (translation from 0,0)
                
                
                p.og_res_stim_total_sz  = 1024;                                 % original resolution of object stimuli
                p.og_res_stim_target_sz = 354;                                  % these og stimuli have a target size of 354 pixels = 4 dva
                
                p.og_res_stim_deg   = parafov_circle_diam_deg;                  % corresponding to 4 deg
                
                p.img_sz_pix        = (p.og_res_stim_total_sz/p.og_res_stim_target_sz)*parafov_circle_diam_pix;                      % scale factor to apply
                p.img_sz_deg        = parafov_circle_diam_deg;                  % height (or width) of square stimulus support (deg)
                p.dres              = p.img_sz_pix / p.og_res_stim_total_sz;    % height (or width) of square stimulus support (pix)
                
                
                p.num_unique_objects = 16;                                      % orientation "bins" from which we create final gabor orientations (deg), 0 = 12 o'clock
                p.facing_dir_deg     = repmat(90 + [-34, 34],1,p.num_unique_objects/2); % rotate 10 deg away from canonical view
                
                p.super_cat          = {'human','animal','object','place'};     % 4 superordinate categories
                
                p.basic_cat{1}       = {'facemale','facefemale','facefemale'};
                p.basic_cat{2}       = {'small','small','big','big'};
                p.basic_cat{3}       = {'tool','tool','food','food','vehicle','vehicle'};
                p.basic_cat{4}       = {'building','building','building'};
                
                p.sub_cat{1}         = {'damon','lisa','sophia'};
                p.sub_cat{2}         = {'parrot','cat','bear','giraffe'};
                p.sub_cat{3}         = {'drill','brush','pizza','banana','bus','suv'};
                p.sub_cat{4}         = {'church','house','watertower'};
                
                p.square_pix_val     = true;
                p.delta_from_ref     = [-8, -4, 4, 8];                        % how much should stim pose rotate from reference (WM: for double epochs)
                % the bigger the delta, the easier the trial. Negative is counter-clockwise, positive is clockwise
                % Add params to struct
                stim.cobj = p;
                
            case {'ns',5}
                % GENERAL
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('scene_%s',disp_params.name)); % mat-file Where to store stimulus images?
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('scene_info_%s',disp_params.name)); % csv-file Where to find stimulus info?
                p.iscolor = true;                                               % Use color or not? [IF YES: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]
                
                % TEMPORAL
                p.duration              = stimdur_frames;                                       % frames (nr of monitor refreshes)
                
                % SPATIAL
                p.og_res_stim           = 425;                                            % original resolution of NSD stimuli
                p.rz_res_stim           = 714;                                            % resized resolution of NSD stimuli
                p.dres                  = ((ctr_square_deg/disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                %             assert(isequal(stim.ns.rz_res_stim, round(stim.ns.img_sz_pix * -stim.dres))); % check: slightly larger than OG NSD stim size?? = 744 pixels 7TAS BOLDscreen Nova1x32
                
                p.img_sz_deg            = ctr_square_deg;                                 % height (or width) of square stimulus support (deg)
                p.img_sz_pix            = ceil(p.og_res_stim.*p.dres);                    % height (or width) of square stimulus support (pix)
                p.square_pix_val        = true;

                
                
                p.x0_deg                = 0;                                       % desired x-center loc in deg (translation from 0,0)
                p.y0_deg                = 0;                                       % desired y-center loc in deg (translation from 0,0)
                p.x0_pix                = round((p.x0_deg.*disp_params.ppd)/2)*2;   % x-center loc in pix (translation from 0,0)
                p.y0_pix                = round((p.y0_deg.*disp_params.ppd)/2)*2;   % y-center loc in pix (translation from 0,0)
                
                p.super_cat     = {'human','animal','food','object','place'};  % 5 superordinate categories
                
                p.basic_cat{1}  = {'indoor','outdoor'};
                p.basic_cat{2}  = {'indoor','outdoor'};
                p.basic_cat{3}  = {'indoor','outdoor'};
                p.basic_cat{4}  = {'indoor','outdoor'};
                p.basic_cat{5}  = {'indoor','outdoor'};
                
                p.sub_cat{1,1}  = {'face1_left','face2_center','face3_right'};
                p.sub_cat{1,2}  = {'face1_left','face2_center','face3_right'};
                
                p.sub_cat{2,1}  = {'cat1_left','cat2_center','cat3_right'};
                p.sub_cat{2,2}  = {'giraffe1_left','giraffe2_center','giraffe3_right'};
                
                p.sub_cat{3,1}  = {'donut1_left','donut2_center','donut3_right'};
                p.sub_cat{3,2}  = {'banana1_left','banana2_center','banana3_right'};
                
                p.sub_cat{4,1}  = {'vase1_left','vase2_center','vase3_right'};
                p.sub_cat{4,2}  = {'bus1_left','bus2_center','bus3_right'};
                
                p.sub_cat{5,1}  = {'bathroom1_left','bathroom2_center','bathroom3_right'};
                p.sub_cat{5,2}  = {'building1_left','building2_center','building3_right'};
                
                p.change_im     = {'easy_added',  'hard_added','easy_removed'  'hard_removed'};
                p.lure_im       = {'lure1',  'lure2', 'lure3', 'lure4'};
                
                % Add params to struct
                stim.ns = p;
                
                fprintf('*** %s: dres (scale factor) = %.4f ***\n',upper(type{ii}),stim.ns.dres);
        end
        
        fprintf('*** %s: Using a stimulus size of %2.2f deg (%3.2f pixels). ***\n', upper(type{ii}), p.img_sz_deg, p.img_sz_pix);
        fprintf('*** %s: Using a stimulus duration of %d frames (%3.2f seconds), where 1 bin is %d frames of 1/%d Hz ***\n', ...
            upper(type{ii}), p.duration, p.duration*stim.framedur_s, stim.presentationrate_hz, disp_params.refresh_hz);
    end

    
    if store_params
        fprintf('[%s]:Storing params..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        save(fullfile(saveDir,sprintf('stim_%s_%s.mat',disp_params.name,datestr(now,30))),'stim')
    end
end


return