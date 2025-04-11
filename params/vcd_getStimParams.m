function stim = vcd_getStimParams(varargin)
% VCD function to get stimulus parameters:
%
%   stim = vcd_getStimParams(type, disp_params,load_params,store_params, overwrite_randomized_params)
%
% Stimulus params such as size and location will depend on display params.
%
% !!WARNING!! There is a randomization component involved in creating some
% stimuli (e.g., orientation of gabor stimuli or dot locations). If you
% don't want this, this leave the fifth argument
% (overwrite_randomized_params) empty (default is set to false) or set to
% false. If you do want regenerate probabilistic params, set the fifth
% argument to true and some stimulus values will change.
%
% INPUTS:
%  stim_class                  : Stimulus class to load params, choose from 
%                                 'gabor','rdk','dot','obj','ns','all' (default is 'all') 
%  disp_name                   : Display name to load params (see vcd_getDisplayParams.m)
%                                 Default: '7TAS_BOLDSCREEN32'
%  load_params                 : Load prior stored parameters or not. Default: true
%  store_params                : Store generated parameters or not. Default: true
%  overwrite_randomized_params : Overwrite stored parameters and regenerate 
%                                 params with probabilistic/randomized element.
%                                 Default: false 
%
% OUTPUT:
%  stim                        : struct with stimulus params, including:
%                                  * fix (fixation dot)
%                                  * bckground (pink noise background).
%                                  * cd (contrast decrement)
%                                  * el (eyelink)
%                                  * gabor
%                                  * rdk (random dot motion kinetograms)
%                                  * dot (single, simple dot)
%                                  * obj (complex objects)
%                                  % ns (natural scenes)
%
% Written by Eline Kupers November 2024 (kupers [at] umn [dot] edu)

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addParameter('stim_class'                 , {'all'}, @(x) ischar(x) || any(strcmp(x, {'all','gabor','rdk','dot','obj','ns'})));
p0.addParameter('disp_name'                  , '7TAS_BOLDSCREEN32', @(x) any(strcmp(x,{'7TAS_BOLDSCREEN32', 'KKOFFICE_AOCQ3277', 'PPROOM_EIZOFLEXSCAN', 'EKHOME_ASUSVE247'})));                   
p0.addParameter('load_params'                , true   , @islogical);                    
p0.addParameter('store_params'               , true   , @islogical); 
p0.addParameter('overwrite_randomized_params', false  , @islogical); 

% Parse inputs
p0.parse(varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

if ischar(stim_class)
    stim_class = {stim_class};
end

if any(strcmp(stim_class{:},'all'))
    stim_class = {'gabor','rdk','dot','obj','ns'};
end

% Retrieve display params
disp_params = vcd_getDisplayParams(disp_name);

%% Load params from stored file if requested
if load_params
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('stim_%s*.mat',disp_name)));
    if ~isempty(d)
        fprintf('[%s]: Found %d stim params .mat file(s)\n',mfilename,length(d));
        if length(d) > 1
            warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
        end
        fprintf('[%s]: Loading stim params .mat file: %s\n', mfilename, d(end).name);
        load(fullfile(d(end).folder,d(end).name),'stim');
    else
        error('[%s]: Can''t find stim params file!\n', mfilename)
    end
else
    fprintf('[%s]: Define stim params and %s overwrite randomized params\n', mfilename, choose(overwrite_randomized_params,'will','will NOT'));
    
    % Setup struct
    stim = struct();
    
    %% GENERAL PARAMS
    stim.store_imgs          = true;                                       % Store images when creating them? (see s_createStimuli.m)
    
    % we want to present images at a rate of 30 Hz
    stim.presentationrate_hz = 30;                                         % rate of stimulus presentation (Hz)
    stim.f2f                 = disp_params.refresh_hz*(1/stim.presentationrate_hz); % nr of monitor refreshes we want to wait to achieve 30 Hz presentation rate. 
                                                                           % so for a monitor refresh rate of 120 Hz, this number is 4, for 60 Hz this is 2. 

    % ensure we have a round number of monitor refreshes to achieve 30 Hz presentation rate  
    assert(isint(disp_params.refresh_hz/stim.presentationrate_hz));
    
    stim.framedur_s          = stim.f2f*(1/disp_params.refresh_hz);        % duration for each presentation frame (should be 33 ms)
    stim.bckgrnd_grayval     = ceil(255/2);                                % background color (middle of 1-255 pixel lum)
    stim.scfactor            = 1;                                          % no scaling
    
    % Default params to inherit (or overwrite)
    stimdur_frames           = stim.presentationrate_hz * 2.0;             % 2 s nr of (33 ms) presentation frames (60 frames = 2 sec)
    x0_deg                   = [-4 4];                                     % desired x-center location for left right stim apertures (degrees)
    y0_deg                   = [0 0];                                      % desired y-center location for left right stim apertures (degrees)
    
    x0_pix                   = round((x0_deg.*disp_params.ppd)/2)*2;       % Empirical x-center location for left right stim apertures (pixels): 354 pixels (4.0059 deg) for BOLD screen, 258 pixels for PP room EIZO screen. 
    y0_pix                   = round((y0_deg.*disp_params.ppd)/2)*2;       % Empirical y-center location for left right stim apertures (pixels): 354 pixels (4.0059 deg) for BOLD screen, 258 pixels for PP room EIZO screen.
    
    % Parafoveal circular aperture (for gabors, rdk, and obj)
    parafov_circle_diam_deg  = 4;                                          % desired  parafoveal circular diameter aperture (empirical is 4.0059 degrees)
    parafov_circle_diam_pix  = round((parafov_circle_diam_deg*disp_params.ppd)/2)*2;  % parafoveal circular diameter aperture (pixels), 354 pixels
    
    % Center square image (for ns)
    ctr_square_deg           = 8.4;                                         % desired center square side length in degrees, same as NSD. Empirical is 8.3965 degrees.
    ctr_square_pix           = round((ctr_square_deg*disp_params.ppd)/2)*2; % center square side length in pixels (742 pixels)
    
    
    %% NOISE BACKGROUND Puzzle piece
    % general
    stim.bckground.stimfile         = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('bckgrnd_%s',disp_params.name)); % mat file
    stim.bckground.infofile         = fullfile(vcd_rootPath,'workspaces','info',sprintf('bckgrnd_info_%s',disp_params.name)); % csv file
    
    stim.bckground.alpha            = 1;                                    % exponent to apply to the amplitude spectrum (i.e. 1/f^alpha).  default: 1.
    stim.bckground.num              = 1;                                    % number of unique background images desired.  default: 1.
    stim.bckground.mode             = 0;                                    % mode of pinknoise function, 0 means fixed amplitude spectrum + random phase
    stim.bckground.std_clip_range   = 3.5;                                  % when converting pixel values to 1-255 range, how many std's of image values do we allow before clipping the range
    
    fprintf('*** %s SCREEN SIZE (hxw): FoV: [%2.2f,%2.2f] degrees visual angle. Resolution = [%02d,%02d] pixels.***\n', disp_params.name, disp_params.h_deg,disp_params.w_deg ,disp_params.h_pix,disp_params.w_pix);
 
    %% FIXATION DOT
    
    % General
    stim.fix.stimfile               = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('fix_%s',disp_params.name)); % mat file
    stim.fix.infofile               = fullfile(vcd_rootPath,'workspaces','info',sprintf('fix_info_%s',disp_params.name)); % csv file
    
    % TEMPORAL
    stim.fix.dotmeanchange          = 1.4*stim.presentationrate_hz;         % (nr 33ms frames), on average, dot changes occur every 1.4 seconds
    stim.fix.dotchangeplusminus     = 0*stim.presentationrate_hz;           % (nr 33ms frames), earliest and latests time that dot changes occur. 2 seconds means [-1:1] from meanchange
    stim.fix.dres                   = [];                                   % rescale factor (fraction)
        
    % SPATIAL
    stim.fix.dotcenterdiam_deg      = 0.14;                                 % (deg) idealized diameter of inner fixation circle
    stim.fix.dotthinborderdiam_deg  = 0.20;                                 % (deg) idealized diameter of thick fixation circle rim
    stim.fix.dotthickborderdiam_deg = 0.25;                                 % (deg) idealized diameter of thin fixation circle rim
    % empirical diameters are boldscreen 12+6 pixels = 18 ~ 0.2037 deg. pproom eizo screen = 9+4 pixels = 13 ~ 0.2020 deg
    % empirical diameters are boldscreen 12+10 pixels = 22 ~ 0.249 deg. pproom eizo screen = 9+7 pixels = 16 pixels ~ 0.2486 deg
    stim.fix.dotcenterdiam_pix      = round(stim.fix.dotcenterdiam_deg * disp_params.ppd);      % inner fixaton cirle diameter in pixels (boldscreen: 12 pixels, pproom: 9 pixels)
    stim.fix.dotthinborderdiam_pix  = round(stim.fix.dotthinborderdiam_deg * disp_params.ppd);  % total fixation circle diameter with thin border in pixels (during ITI/IBI) 
    stim.fix.dotthickborderdiam_pix = round(stim.fix.dotthickborderdiam_deg * disp_params.ppd); % total fixation circle diameter with thick border in pixels (during trial)
    

    stim.fix.lumminmaxstep          = [41,213,5];                          % min, max, and nr of dot luminance values [1-255],
    stim.fix.dotlum                 = linspace(stim.fix.lumminmaxstep(1),stim.fix.lumminmaxstep(2),stim.fix.lumminmaxstep(3)); % dot gray luminance levels (1-255): 41 84 127 170 213
    stim.fix.dotopacity             = 0.5;                                 % dot and border have 50% opacity
    stim.fix.color                  = [255, 255, 255; 255 0 0];            % white and red (for spatial cue)
   
    fprintf('*** FIXATION CIRCLE: inner center diameter = %d, inner+thin rim = %d, inner+thick rim = %d pixels ***\n', ...
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
    stim.cd.t_gausswin_N            = 15;                                   % 15 number of presentation frames (33 ms) for gaussian time window (contrast decrement)
    stim.cd.t_gausswin_std          = 3;                                    % standard devation of gaussian window in time (presentation frames)
    stim.cd.meanchange              = stim.presentationrate_hz * 1.0;       % mean of gaussian window in time (30 frames = 1 sec)  
    stim.cd.changeplusminus         = (0.5/stim.framedur_s)-1;              % plus or minus this amount (14 frames = 0.46 sec)  
    stim.cd.max_cd                  = 0.2;                                  % stimulus contrast is reduced by 20% of mean luminance at lowest point of temporal gaussian window (note: this corresponds to subtracting a contrast fraction of 10.^(log10(c)-0.1))
    stim.cd.prob                    = 0.5;                                  % 50% probability that a trial will have a luminance change
    
    % Create 1D gaussian
    t_support           = linspace(-stim.cd.t_gausswin_N / 2, stim.cd.t_gausswin_N / 2, stim.cd.t_gausswin_N);
    t_gauss             = exp(-t_support .^ 2 / (2 * stim.cd.t_gausswin_std ^ 2));
    
    % invert and scale it
    t_gauss             = 1-(t_gauss*stim.cd.max_cd);
    stim.cd.t_gauss     = t_gauss;

    %% EYELINK PARAMS
    % we manually place the eyelink calibration/validation points on the
    % display. The distance between the center and 4 left/right/up/down
    % points are set as [xc,yc] ± 350 pixels. This results in dots at the
    % following screen coordinates in pixels:
    %                  [x3,y3]=[0,190]
    % [x1,y1]=[610,0]  [x0,y0]=[0,0  ]   [x2,y2]=[1310,0]
    %                  [x4,y4]=[0,890] 
    stim.el.point2point_distance_deg = 3.0;                                % desired target distance (in deg) from fixation (empirical is 3.0059 degrees)
    stim.el.point2point_distance_pix = round((stim.el.point2point_distance_deg*disp_params.ppd/2))*2; % pixels

    %% STIM PARAMS
    for ii = 1:length(stim_class)
        
        p = [];
        
        switch stim_class{ii}
            
            case {'gabor',1}
                
                % Where to store stimulus images?
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('gabor_%s',disp_params.name)); % mat file
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('gabor_info_%s',disp_params.name)); % csv file
                
                % TEMPORAL
                % Fixed params
                p.duration        = stimdur_frames;                             % frames (nr of monitor refreshes)

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
                p.cycles_per_pix  = (p.sf_cpd/p.img_sz_deg)/disp_params.ppd;    % nr of cycles per image (pix)
                
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
                    p.ori_deg     = [12 37 51 80 102 121 145 169];              % final gabor orientations
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
                p.duration        = stimdur_frames;                            % frames (nr of monitor refreshes)
                p.dots_coherence  = [0.064, 0.128, 0.512];                     % fraction of dots. Kiani lab uses usually one of these [0 0.032 0.064 0.128 0.256 0.512]
                p.dots_speed      = 5;                                         % pixels/s?   Kiani lab uses usually 5 to 10
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
                    p.dots_direction   = round(p.motdir_bins + p.motdir_jitter);     % final gabor orientations. Note this is only for one hemifield (left), will be mirrored for right.
                else
                    p.motdir_jitter    = [8.1568 6.5389 -1.6998 7.0698 2.4508 0.8739 2.4295 0.5901]; 
                    p.dots_direction   = [18 62 98 152 192 236 282 326];
                end
                
                p.dots_direction = p.dots_direction-90;                         % subtract 90 deg to ensure 0 deg is now 12 o'clock in x,y-pixel space  
                
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
                p.img_sz_deg      = 1.0;                                     	 % desired spatial support of dot image in deg, empirical is 0.9958 deg for BOLD screen
                p.img_sz_pix      = round((p.img_sz_deg * disp_params.ppd)/2)*2; % spatial support of dot image in pixels (88 pixels)
                p.og_res_stim     = p.img_sz_pix;                                % resolution of stored dot stimuli (in pixels)
                p.dres            = (( (p.img_sz_pix/disp_params.ppd) /disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                p.radius_deg      = 0.5;                                         % desired dot radius in deg (empirical is 0.4979 for BOLDscreen)
                p.radius_pix      = round((p.radius_deg * disp_params.ppd)/2)*2; % dot radius in pix (44 pixels for BOLDscreen)
                p.color           = [255 255 255];                               % white in uint8 RGB
                p.contrast        = 1;                                           % Michelson [0-1] (fraction)
                p.square_pix_val  = false;                                       % no need to square pixel luminance values for linear CLUT
                
                % Alpha mask
                p.alpha_mask_diam_pix = p.radius_pix+1;                          % add one pixel for alpha mask (EK: why 1?)

                p.num_loc         = 16;                                          % orientation "bins" from which we create final gabor orientations (deg), 0 = 12 o'clock
                p.loc_jitter_sd   = 1;                                           % std of normal distribution to sample orientation jitter
                p.loc_jitter_mu   = 1;                                           % mean of normal distribution to sample orientation jitter
                p.loc_bins        = 10:(360/p.num_loc):350;                      % rotate 2 deg away from vertical, to avoid ill-defined response options and have integer angle differences between locations
                
                if overwrite_randomized_params 
                    p.loc_jitter   = p.loc_jitter_mu + (p.loc_jitter_sd.*randn(1,p.num_loc)); % add a small amount of jitter around the orientation
                    p.loc_deg      = round(p.loc_bins + p.loc_jitter);        % final gabor orientations (add 90 deg to make 0 deg 12 o'clock / north / upper vertical     
                    assert(all(p.loc_jitter < 3));                  % make sure all jitter is small (within 3 deg), so we don't get crazy outliers due to random sampling 
                else
%                     p.loc_jitter_R    = [0.8759,2.4897,2.4090,2.4172,1.6715,-0.2075,1.7172,2.6302,1.4889, 2.0347,1.7269,0.6966,1.2939,0.2127,1.8884,-0.1471];
%                     p.loc_deg_R       = [4 17 28 39 50 59 72 84 94 106 117 127 139 149 162 172]; 
%                     p.loc_jitter_L    = [0.2685,1.5216,1.0164,1.4870,-0.9407,0.7928,1.3945,1.1130,0.3466,-0.3815,0.9270,1.6167,0.0299,1.5483,0.6760,-0.9645];
%                     p.loc_deg_L       = [183,196,207,218,227,240,252,263,273,284,296,308,318,331,341,351]; 
                    p.loc_jitter    = [0.2539,2.8806,-1.0658,1.4559,0.4893,1.7060,0.5951,1.8281,1.4718,1.8239,0.7439,1.5946,-0.3104,0.8743,0.8019,2.1788];
                    p.loc_deg       = [10,35,54,79,100,124,146,169,191,214,236,259,280,303,326,350]; 

                end
                
                p.iso_eccen           = 4.0;                                       % idealized iso-eccentric dot location (for all angles)
                p.ang_deg             = p.loc_deg;                                 % idealized angle of center dot loc in deg (translation from center screen 0,0)
                p.eccen_deg           = repmat(p.iso_eccen,1,length(p.loc_deg));   % idealized eccen of center dot loc in deg (translation from center screen 0,0)
                
                % convert polar angle coords to cartesian coords
                [x,y] = pol2cart(deg2rad(p.ang_deg+90),p.eccen_deg);               % add 90 deg to ensure 0 deg is now 12 o'clock in x,y-pixel space
                p.x0_pix        = disp_params.xc - round(x * disp_params.ppd);     % idealized x-center loc in pix (translation from upper left corner [0,0])
                p.y0_pix        = disp_params.yc - round(y * disp_params.ppd);     % idealized y-center loc in pix (translation from upper left corner [0,0])
                 
%                 figure(1); clf; scatter(p.x0_pix,p.y0_pix,'bo');
%                 hold on; plot([disp_params.xc disp_params.xc],[100 1000],'k')
%                 hold on; plot([400 1400],[disp_params.yc disp_params.yc],'k')
%                 
                p.delta_from_ref  = [-15, -5, 5, 15];                                % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                % do the same for delta offsets convert polar angle coords to cartesian coords
                for dd = 1:length(p.delta_from_ref)
                    p.ang_deg_delta(dd,:)   = p.ang_deg + p.delta_from_ref(dd);  % idealized angle offset for dot loc in deg (translation from center screen 0,0)
                    p.eccen_deg_delta(dd,:) = p.eccen_deg;                       % idealized eccen for dot loc in deg (translation from center screen 0,0)
                    
                    [x_d,y_d] = pol2cart(deg2rad(p.ang_deg_delta(dd,:)+90),p.eccen_deg_delta(dd,:)); % again, note the +90
                    p.x0_pix_delta(dd,:) = disp_params.xc - round(x_d * disp_params.ppd);      % idealized x-center ref dot loc in pix (translation from upper left corner [0,0])
                    p.y0_pix_delta(dd,:) = disp_params.yc - round(y_d * disp_params.ppd);      % idealized y-center ref dot loc in pix (translation from upper left corner [0,0])
                    
                end
%                 figure(1); 
%                 hold on; scatter(p.x0_pix_delta,p.y0_pix_delta,'rx');
                 
                
                
                % the bigger the delta, the easier the trial
                % Add params to struct
                stim.dot = p;
                
            case {'obj',4}
                % GENERAL
                p.indivobjfile = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,'vcd_objects_2degstep_lumcorrected');
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',sprintf('object_%s',disp_params.name)); % mat-file where to store stimulus images?
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('object_info_%s',disp_params.name));  % csv-file Where to find stimulus info?
                p.iscolor         = false;                                      % use color or not [[[IF WE USE COLOR: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]]
                
                
                % TEMPORAL
                p.duration      = stimdur_frames;                                     % frames (nr of monitor refreshes)
                
                % SPATIAL
                p.contrast      = 1;                                            % Michelson [0-1] (fraction) 
                p.x0_deg        = x0_deg;                                       % x-center loc in deg (translation from center screen 0,0)
                p.y0_deg        = y0_deg;                                       % y-center loc in deg (translation from center screen 0,0)
                p.x0_pix        = x0_pix;                                       % x-center loc in pix (translation from center screen 0,0)
                p.y0_pix        = y0_pix;                                       % y-center loc in pix (translation from center screen 0,0)
                p.og_res_stim_total_sz  = 1024;                                 % original resolution of object stimuli

                p.og_res_stim_target_sz = parafov_circle_diam_pix;              % object stimuli have a target size of 4 dva = 354 pixels (7TAS_BOLDSCREEN32) or 258 pixels (PP room)
                p.og_res_stim_deg   = parafov_circle_diam_deg;                  % corresponding to 4 deg
                
                p.img_sz_pix        = (p.og_res_stim_total_sz/p.og_res_stim_target_sz)*parafov_circle_diam_pix;                      % scale factor to apply
                p.img_sz_deg        = parafov_circle_diam_deg;                  % height (or width) of square stimulus support (deg)
                p.dres              = p.img_sz_pix / p.og_res_stim_total_sz;    % height (or width) of square stimulus support (pix)
                
                
                p.num_unique_objects = 16;                                      % orientation "bins" from which we create final gabor orientations (deg), 0 = 12 o'clock
                p.canonical_facing_dir_deg = repmat(90 + [-45, 45],1,p.num_unique_objects/2); % rotate 10 deg away from canonical view
                p.rot_bins      = linspace(10,170,p.num_unique_objects);
                p.rot_bins      = p.rot_bins(1:2:end);
                p.rot_jitter_sd = 1;                                           % std of normal distribution to sample rotation jitter
                p.rot_jitter_mu = 1;                                           % mean of normal distribution to sample rotation jitter
                
                if overwrite_randomized_params 
                    p.rot_jitter       = p.rot_jitter_mu + (p.rot_jitter_sd.*randn(1,p.num_unique_objects)); % add a small amount of jitter around the rotation
                    p.facing_dir_deg   = round(p.rot_bins + p.rot_jitter);     % final facing direction for all objects
                    p.facing_dir_deg   = ceil(p.facing_dir_deg/2)*2;           % ensure we only deal with even rotations
                    p.facing_dir_deg   = p.facing_dir_deg(shuffle_concat([1:p.num_unique_objects],1)); % shuffle order of facing directions
                else
                    p.rot_jitter       = [0.2280    0.7061    2.8998    0.5255    1.2055    0.8602   -0.1122  -0.6935    0.8379    2.1032    1.0862   -0.3847    1.5598    0.6225  0.2114  0.1853]; 
                    p.facing_dir_deg   = [64 44 150 34 96 128 84 22 10 170 108 74 118 160 54 140 64 44 150 34 96 128 84 22 10 170 108 74 118 160 54 140];
                end
                
                % Define the 4 superordinate, 8 basic, and 16 subordinate categories
                p.super_cat          = {'human','animal','object','place'};     
                
                p.basic_cat{1}       = {'facemale','facefemale','facefemale'};
                p.basic_cat{2}       = {'small','small','big','big'};
                p.basic_cat{3}       = {'tool','tool','food','food','vehicle','vehicle'};
                p.basic_cat{4}       = {'building','building','building'};
                
                p.sub_cat{1}         = {'damon','lisa','sophia'};
                p.sub_cat{2}         = {'parrot','cat','bear','giraffe'};
                p.sub_cat{3}         = {'drill','brush','pizza','banana','bus','suv'};
                p.sub_cat{4}         = {'church','house','watertower'};
                
                % Define the affordances for each object subcategory (for HOW task) 
                p.affordance{1}       = {'greet','greet','greet'};
                p.affordance{2}       = {'greet','greet','observe','observe'};
                p.affordance{3}       = {'grasp','grasp','grasp','grasp','enter','enter'};
                p.affordance{4}       = {'enter','enter','observe'};

                p.square_pix_val     = true;
                p.delta_from_ref     = [-8, -4, 4, 8];                        % Stim pose rotation from reference for WM test image
                                                                              % the bigger the delta, the easier the trial. 
                                                                              % for 1-89 deg rotations: Negative values are leftwards, positive values is rightwards
                                                                              % for 91-180 deg rotations: Negative values are rightward, positive values is leftward
                % Add params to struct
                stim.obj = p;
                
            case {'ns',5}
                % GENERAL
                p.indivscenefile = fullfile(vcd_rootPath,'workspaces','stimuli','vcd_natural_scenes'); % mat-file Where to load indiv pngs?
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
                
                p.affordance{1,1} = {'walk','greet','greet'};
                p.affordance{1,2} = {'greet','observe','observe'};
                p.affordance{2,1} = {'greet','greet','greet'};
                p.affordance{2,2} = {'observe','observe','observe'};
                p.affordance{3,1} = {'grasp','grasp','grasp'};
                p.affordance{3,2} = {'grasp','grasp','observe'};
                p.affordance{4,1} = {'grasp','grasp','grasp'};
                p.affordance{4,2} = {'walk','observe','walk'};
                p.affordance{5,1} = {'walk','walk','walk'};
                p.affordance{5,2} = {'observe','walk','walk'};
                
                p.ltm_im        = [2,4,5,8,10,11,13,15,18,20,21,23,26,27,30];
                
                p.change_im     = {'easy_add', 'hard_add','easy_remove', 'hard_remove'};
                p.lure_im       = {'lure01', 'lure02', 'lure03', 'lure04'};
                
                % Add params to struct
                stim.ns = p;
                
                fprintf('*** %s: dres (scale factor) = %.4f ***\n',upper(stim_class{ii}),stim.ns.dres);
        end
        
        fprintf('*** %s: Using a stimulus size of %2.2f deg (%3.2f pixels). ***\n', upper(stim_class{ii}), p.img_sz_deg, p.img_sz_pix);
        fprintf('*** %s: Using a stimulus duration of %d frames (%3.2f seconds), where one frame relates to %d monitor refreshes (%d Hz) ***\n', ...
            upper(stim_class{ii}), p.duration, p.duration*stim.framedur_s, stim.f2f, disp_params.refresh_hz);
    end

    
    if store_params
        fprintf('[%s]:Storing params..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        save(fullfile(saveDir,sprintf('stim_%s_%s.mat',disp_params.name,datestr(now,30))),'stim')
    end
end


return