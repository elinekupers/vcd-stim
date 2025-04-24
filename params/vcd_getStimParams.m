function stim = vcd_getStimParams(varargin)
% VCD function to get stimulus parameters:
%
%   stim = vcd_getStimParams(['stim_class',{'all'},] ...
%                            ['disp_name','7TAS_BOLDSCREEN32',] ...
%                            ['load_params',false,] ...
%                            ['store_params',true,] ...
%                            ['overwrite_randomized_params',false]);
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

%% Obtain display params
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
    %% We will load params if requested
    fprintf('[%s]: Define stim params and %s overwrite randomized params\n', mfilename, choose(overwrite_randomized_params,'will','will NOT'));
    
    % Setup struct
    stim = struct();
    
    %% GENERAL PARAMS
    stim.store_imgs          = true;                                       % Store images when creating them? (see s_createStimuli.m)
    
    % we want to present images at a rate of 30 Hz
    stim.presentationrate_hz = 30;                                         % rate of stimulus presentation (Hz)
    stim.f2f                 = disp_params.refresh_hz*(1/stim.presentationrate_hz); % frame-to-frame duration: nr of monitor refreshes 
                                                                           % we want to wait to achieve 30 Hz presentation rate. 
                                                                           % so this number is 4 for a monitor refresh rate of 120 Hz and 2 for 60 Hz. 

    % Ensure we have a round number of monitor refreshes to achieve 30 Hz presentation rate  
    assert(isint(disp_params.refresh_hz/stim.presentationrate_hz));
    
    stim.framedur_s          = stim.f2f*(1/disp_params.refresh_hz);        % duration for each presentation frame (should be 33 ms)
    assert(isequal(stim.framedur_s, (1/stim.presentationrate_hz)));        % check if this is 33 ms.
    
    stim.bckgrnd_grayval     = ceil(255/2);                                % background color (middle of 1-255 pixel lum)
    assert(isequal(stim.bckgrnd_grayval, 128));                            % we want 128 to be mean luminance
    
    stim.scfactor            = 1;                                          % no global scaling (only for specific stimulus images)
    
    % Default params to inherit (or overwrite)
    
    % **** TEMPORAL **** 
    stimdur_frames           = stim.presentationrate_hz * 2.0;             % nr of (33 ms) presentation frames that results in 2 s (60 frames = 2 sec)
    
    % **** SPATIAL ****
    
    % EMPIRICAL parafoveal stimulus locations:
    %  * BOLDscreen: 354 pixels (4.0059 deg). 
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
    % general
    stim.bckground.stimfile         = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('bckgrnd_%s',disp_params.name)); % mat file
    stim.bckground.infofile         = fullfile(vcd_rootPath,'workspaces','info',sprintf('bckgrnd_info_%s',disp_params.name)); % csv file
    
    stim.bckground.alpha            = 1;                                    % exponent to apply to the amplitude spectrum (i.e. 1/f^alpha).  default: 1.
    stim.bckground.num              = 1;                                    % number of unique background images desired.  default: 1.
    stim.bckground.mode             = 0;                                    % mode of pinknoise function, 0 means fixed amplitude spectrum + random phase
    stim.bckground.std_clip_range   = 3.5;                                  % when converting pixel values to 1-255 range, how many std's of image values do we allow before clipping the range
    
    fprintf('*** %s SCREEN SIZE (hxw): FoV: [%2.2f,%2.2f] degrees visual angle. Resolution = [%02d,%02d] pixels.***\n', disp_params.name, disp_params.h_deg,disp_params.w_deg ,disp_params.h_pix,disp_params.w_pix);
 
    %% FIXATION CIRCLE
    
    % General
    stim.fix.stimfile               = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('fix_%s',disp_params.name)); % mat file
    stim.fix.infofile               = fullfile(vcd_rootPath,'workspaces','info',sprintf('fix_info_%s',disp_params.name)); % csv file
    
    % TEMPORAL
    stim.fix.dotmeanchange          = 1.4*stim.presentationrate_hz;         % nr 33ms frames, dot changes occur every 1.4 seconds
    stim.fix.dotchangeplusminus     = 0*stim.presentationrate_hz;           % nr 33ms frames, earliest and latests time that dot changes occur. 2 seconds means [-1:1] from meanchange
    stim.fix.dres                   = [];                                   % rescale factor as a fraction between 0-1  
        
    % SPATIAL
    stim.fix.dotcenterdiam_deg      = 0.14;                                 % (deg) idealized diameter of inner fixation circle
    stim.fix.dotthinborderdiam_deg  = 0.20;                                 % (deg) idealized diameter of thick fixation circle rim
    stim.fix.dotthickborderdiam_deg = 0.25;                                 % (deg) idealized diameter of thin fixation circle rim
    
    % THIN RIM Empirical diameters are:
    % * BOLDscreen: 12 pixels for center + 6 pixels for rim (split across both sides) = 18 pixels ~ 0.2037 deg. 
    % * PP room EIZOFLEX screen = 9 pixels for center + 4 pixels for rim (split across both sides) = 13 pixels ~ 0.2020 deg.
    
    % THICK RIM Empirical diameters are:
    % * BOLDscreen: 12 pixels for center + 10 pixels for rim (split across both sides) = 22 ~ 0.249 deg. 
    % * PP room EIZOFLEX screen = 9 pixels for center + 7 pixels for rim (split across both sides) = 16 pixels ~ 0.2486 deg
    stim.fix.dotcenterdiam_pix      = round(stim.fix.dotcenterdiam_deg * disp_params.ppd);      % inner fixaton cirle diameter in pixels (boldscreen: 12 pixels, pproom: 9 pixels)
    stim.fix.dotthinborderdiam_pix  = round(stim.fix.dotthinborderdiam_deg * disp_params.ppd);  % total fixation circle diameter with thin border in pixels (during ITI/IBI) 
    stim.fix.dotthickborderdiam_pix = round(stim.fix.dotthickborderdiam_deg * disp_params.ppd); % total fixation circle diameter with thick border in pixels (during trial)
    
    % FIX CIRCLE LUMINANCE & COLOR
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
    t_support           = linspace(-stim.cd.t_gausswin_N / 2, stim.cd.t_gausswin_N / 2, stim.cd.t_gausswin_N);  % [-7:1:7] in units of presentation 33 ms frames
    t_gauss             = exp(-t_support .^ 2 / (2 * stim.cd.t_gausswin_std ^ 2)); % 1D gaussian with a sd of 3 presentation frames
    
    % invert and scale it
    t_gauss             = 1-(t_gauss*stim.cd.max_cd); % invert gaussian, we start from 1, then dip to 0.5 and go back to 1. 
    stim.cd.t_gauss     = t_gauss;

    %% EYELINK PARAMS
    % we manually place the eyelink calibration/validation points on the
    % display. The distance between the center and 4 left/right/up/down
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
    stim.el.point2point_distance_deg = 3.0;                                % desired target distance (in deg) from fixation 
    stim.el.point2point_distance_pix = round((stim.el.point2point_distance_deg*disp_params.ppd/2))*2; % desired target distance in pixels
    
    
    %% STIM PARAMS
    for ii = 1:length(stim_class)
        
        p = [];
        
        switch stim_class{ii}
            
            case {'gabor',1}
                
                % Where to store stimulus images?
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('gabor_%s',disp_params.name)); % mat file
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('gabor_info_%s',disp_params.name)); % csv file
                
                % GENERAL
                p.unique_im_nrs   = [1:24];                                     % Unique image nrs associated with the CORE 24 Gabors

                % TEMPORAL
                % Fixed params
                p.duration        = stimdur_frames;                             % frames (nr of monitor refreshes)

                % SPATIAL
                % Fixed params
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
                
                p.sf_cpd          = 4;                                          % spatial frequency (cycles/deg)
                p.cycles_per_pix  = (p.sf_cpd/p.img_sz_deg)/disp_params.ppd;    % nr of cycles per image (pix)
                
                p.x0_deg          = x0_deg;                                     % x-center loc in deg (translation from 0,0)
                p.y0_deg          = y0_deg;                                     % y-center loc in deg (translation from 0,0)
                p.x0_pix          = x0_pix;                                     % x-center loc in pix (translation from 0,0)
                p.y0_pix          = y0_pix;                                     % y-center loc in pix (translation from 0,0)
                
                p.ph_deg          = [0:(180/4):179];                            % 4 quadrature Gabor phases
                
                % Manipulated params
                p.contrast        = [0.05, 0.10, 1];                            % Michelson contrasts [0-1] (fraction)
                p.num_ori         = 8;                                          % gabor orientations (deg), 0 = 12 o'clock
                p.ori_deg         = [0:(180/p.num_ori):179]+(0.5*(180/p.num_ori)); % rotate half shift away from vertical, to avoid ill-defined response options                                                                 
                                                                                % will return [11.25 33.75 56.25 78.75 101.25 123.75 146.25 168.75]
                
                % ensure orientations have equal distance from cardinal meridians
                assert(isequal(sort(abs(0-p.ori_deg),2),sort(abs(180-p.ori_deg),2)))
                assert(isequal(sort(abs(90-p.ori_deg(1:(p.num_ori/2))),2),sort(abs(90-p.ori_deg(((p.num_ori/2)+1):p.num_ori)),2)))
                                                                                
                % WORKING MEMORY: gabor orientation deltas for test images
                p.delta_from_ref  = [-15, -5, 5, 15];                           % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                                                                                % the bigger the delta, the easier the judgement in a trial
                p.unique_im_nrs_WM = [111:206];                                 % Unique image nrs associated with the 96 WM test images

                % check if all test images for WM have unique orientations
                tmp = p.ori_deg + [0, p.delta_from_ref]';
                tmp(tmp<0)   = tmp(tmp<0)+180;
                tmp(tmp>180) = tmp(tmp>180)-180;
                assert(isequal(length(unique(tmp(:))), length(tmp(:))));
                clear tmp

                % IMAGERY
                p.imagery_im_nr   = p.unique_im_nrs(17:end);                    % SELECTED UNIQUE IMAGES (SUBSET of all 24) (only high contrast)
                p.imagery_sz_deg  = 5.658;                                      % QUIZ DOT PARAMS (STIM 2) desired diameter (deg) of the quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix  = (round(p.img_sz_deg* disp_params.ppd)/2)*2; % QUIZ DOT PARAMS (STIM 2) diameter of quiz dot image in pixels (ensure even nr of pixels)
                
                p.unique_im_nrs_IMG_test  = [610:669];                          % Unique image nrs associated with the 8*20=60 IMG gabor test dot images

                % Add params to struct
                stim.gabor = p;
                
                
            case {'rdk',2}
                % Where to store stimulus images?
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('rdk_%s',disp_params.name)); % mat file
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('rdk_info_%s',disp_params.name)); % csv file
                
                % GENERAL
                p.unique_im_nrs    = [25:48];                                  % Unique image nrs associated with the CORE 24 RDK stimuli

                % TEMPORAL
                p.duration        = stimdur_frames;                            % frames (nr of monitor refreshes)
                p.dots_coherence  = [0.064, 0.128, 0.512];                     % fraction of coherent moving dots. Kiani lab uses usually one of these [0 0.032 0.064 0.128 0.256 0.512]
                p.dots_speed      = (5*disp_params.ppd)/stim.presentationrate_hz; % pixels/frames (or 5 deg/s). For reference: Kiani lab uses usually 5 to 10 deg/s. Rokers lab uses 5 deg/s.
                p.dots_interval   = 1;                                         % update dots every frame   (for references: Kiani's 75 hz refresh rate / 3 interval = 25 frames/sec)
                p.dots_lifetime   = 0.1 * stim.presentationrate_hz;            % 3 frames / 0.1 seconds   
                
                % SPATIAL
                p.img_sz_deg      = parafov_circle_diam_deg;                    % stimulus aperture diameter (deg)
                p.img_sz_pix      = parafov_circle_diam_pix;                    % stimulus aperture diameter (pix)
                p.og_res_stim     = p.img_sz_pix;                               % resolution of stored dot stimuli (in pixels)
                p.dres            = (( (p.img_sz_pix/disp_params.ppd) /disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                p.iscolor         = false;                                      % use color or not [[[IF WE USE COLOR: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]]
                p.square_pix_val  = false;

                p.num_mot_dir      = 8;                                         % number of sampled motion directions
                p.dots_direction   = [0:(360/p.num_mot_dir):359]+(0.5*(360/p.num_mot_dir)); % sample direction of coherent motion from [0-359] in deg (0 deg is aligned with 12 o'clock)
                                                                                % turns out to be: [22.5   67.5  112.5  157.5  202.5  247.5  292.5 337.5]
                % ensure equal distance from cardinal meridians
                assert(isequal(abs(0-p.dots_direction(1:(p.num_mot_dir/2))), abs(180-p.dots_direction(((p.num_mot_dir/2)+1):p.num_mot_dir))));
                assert(isequal(abs(90-p.dots_direction(1:(p.num_mot_dir/2))), abs(270-p.dots_direction(((p.num_mot_dir/2)+1):p.num_mot_dir))));
                                
                % RDK specific
                p.dots_size        = 3;                                         % single dot radius in pixels??
                p.dots_color       = [255 255 255; 1 1 1]./255;                 % 50:50 white:black, color in RGB and converted to [0-1] as expected by stimulus creation function
                p.max_dots_per_frame = 200;                                     % how many dots within a square support from (number is similar to Kiani lab, rokers lab aims for 150) and roughly matches to nr of pixels in aperture
                p.dots_contrast    = 1;                                         % Michelson [0-1] (fraction)
                p.x0_deg           = x0_deg;                                    % Desired x-center loc of stimulus in deg (translation from 0,0)
                p.y0_deg           = x0_deg;                                    % Desired y-center loc of stimulus in deg (translation from 0,0)
                p.x0_pix           = x0_pix;                                    % x-center loc in pix (translation from 0,0)
                p.y0_pix           = y0_pix;                                    % y-center loc in pix (translation from 0,0)
                
                % WORKING MEMORY: motion direction deltas for test images
                p.delta_from_ref   = [-15, -5, 5, 15];                          % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                p.unique_im_nrs_WM = [207:302];                                 % Unique image nrs associated with the 96 WM RDK test stimuli
 
                % check if all test images for WM have unique orientations
                tmp          = p.dots_direction+[0, p.delta_from_ref]';
                tmp(tmp<0)   = tmp(tmp<0)+360;
                tmp(tmp>360) = tmp(tmp>360)-360;
                assert(isequal(length(unique(tmp(:))), length(tmp(:))));
                clear tmp 

                % IMAGERY: 
                p.imagery_im_nr    = p.unique_im_nrs(17:end);                   % SELECTED UNIQUE IMAGES (SUBSET of all 24) (only high coherence) 
                p.imagery_sz_deg      = 5.658;                                  % desired diameter (degree) of the second, quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix      = (round(p.img_sz_deg* disp_params.ppd)/2)*2; % diameter of quiz dot image (pixels) (ensure even nr of pixels)
                p.unique_im_nrs_IMG_test  = [670:729];                          % Unique image nrs associated with the 8*20=60 IMG RDK test dot images

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
                p.unique_im_nrs   = [49:64];                                     % Unique image nrs associated with the 16 single dot stimuli
                
                % TEMPORAL
                p.duration        = stimdur_frames;                              % frames (nr of monitor refreshes)
                
                % SPATIAL
                % Empirical single dot radius is:
                % * BOLDscreen:   44 pixels (0.4979 deg)
                % * EIZOFLEXSCAN: 32 pixels (0.4971 deg)
                p.img_sz_deg      = 1.0;                                     	 % desired spatial support of dot image in deg.
                p.img_sz_pix      = round((p.img_sz_deg * disp_params.ppd)/2)*2; % spatial support of dot image in pixels (BOLDscreen: 88 pixels. EIZOFLEXSCAN: 64 pixels)
                p.og_res_stim     = p.img_sz_pix;                                % resolution of stored dot stimuli (in pixels)
                p.dres            = (( (p.img_sz_pix/disp_params.ppd) /disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                p.radius_deg      = 0.5;                                         % desired dot radius in deg
                p.radius_pix      = round((p.radius_deg * disp_params.ppd)/2)*2; % dot radius in pix
                p.color           = [255 255 255];                               % white in uint8 RGB
                p.contrast        = 1;                                           % Michelson [0-1] (fraction)
                p.square_pix_val  = false;                                       % no need to square pixel luminance values for linear CLUT
                
                % Alpha mask
                p.alpha_mask_diam_pix = p.radius_pix+1;                          % add one pixel for alpha mask (EK: why 1?)

                % DOT LOCATIONS
                p.num_loc         = 16;                                          % number of equally spaced dot angles (deg), 0 = 12 o'clock
                p.ang_deg         = [0:(360/p.num_loc):359]+(0.5*(360/p.num_loc)); % idealized angle of center dot loc in deg, 0 deg = 12 o'clock. rotate half a shift away from vertical, to avoid ill-defined response options
                                                                                   % location angles turn out to be: 11.25 33.75 56.25 78.75 101.25 123.75 146.25 11.25 33.75 56.25 78.75 101.25 123.75 146.25 326.25 348.75
                p.iso_eccen       = 4.0;                                         % desired iso-eccentric dot location (for all angles)
                p.eccen_deg       = repmat(p.iso_eccen,1,length(p.ang_deg));     % desired eccen of center dot loc in deg (translation from center screen 0,0)
                
                % ensure equal distance from cardinal meridians
                assert(isequal(abs(0-p.ang_deg(1:(p.num_loc/2))), abs(180-p.ang_deg(((p.num_loc/2)+1):p.num_loc))));
                assert(isequal(abs(90-p.ang_deg(1:(p.num_loc/2))), abs(270-p.ang_deg(((p.num_loc/2)+1):p.num_loc))));
                
                % Convert dot polar angle coords to cartesian coords
                [x,y] = pol2cart(deg2rad(p.ang_deg+90),p.eccen_deg);               % add 90 deg to ensure 0 deg is now 12 o'clock in x,y-pixel space
                p.x0_pix        = disp_params.xc - round(x * disp_params.ppd);     % desired x-center loc in pix (translation from upper left corner [0,0])
                p.y0_pix        = disp_params.yc - round(y * disp_params.ppd);     % desired y-center loc in pix (translation from upper left corner [0,0])
                
                % WORKING MEMORY: dot angle position deltas for test images
                p.delta_from_ref  = [-15, -5, 5, 15];                              % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                                                                                   % the bigger the delta, the easier the trial
                p.unique_im_nrs_WM = [303:366];                                    % Unique image nrs associated with the 64 WM DOT test stimuli
                                
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
                    
                    [x_d,y_d] = pol2cart(deg2rad(p.ang_deg_delta(dd,:)+90),p.eccen_deg_delta(dd,:)); % again, note the +90
                    p.x0_pix_delta(dd,:) = disp_params.xc - round(x_d * disp_params.ppd);  % x-center ref dot loc in pix (translation from upper left corner [0,0])
                    p.y0_pix_delta(dd,:) = disp_params.yc - round(y_d * disp_params.ppd);  % y-center ref dot loc in pix (translation from upper left corner [0,0])
                end

                % LTM PAIR
                p.ltm_pairs          = [];                                       %%

                % IMAGERY: 
                p.imagery_im_nr    = p.unique_im_nrs(1:2:p.num_loc);             % SELECTED UNIQUE IMAGES (SUBSET of all 16 dots)
                
                % IMAGERY: QUIZ DOT PARAMS
                p.imagery_sz_deg   = [disp_params.w_deg/2, disp_params.h_deg];   % desired diameter (deg) of the second, quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix   = [disp_params.xc,disp_params.h_pix];         % diameter of quiz dot image (pixels) (we already ensured even nr of pixels in display params function)
                p.unique_im_nrs_IMG_test = [730:789];                            % Unique image nrs associated with the 8*20=60 IMG DOT test dot images
                
                % Add params to struct
                stim.dot = p;
                
            case {'obj',4}
                
                % GENERAL
                p.indivobjfile   = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name, 'vcd_objects_2degstep_lumcorrected');
                p.stimfile       = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name, sprintf('object_%s',disp_params.name)); % mat-file where to store stimulus images?
                p.infofile       = fullfile(vcd_rootPath,'workspaces','info',sprintf('object_info_%s',disp_params.name));  % csv-file Where to find stimulus info?

                % GENERAL
                p.unique_im_nrs      = [65:80];                             % Unique image nrs associated with the 16 single dot stimuli
                p.iscolor = false;                                          % Use color or not? 
                if p.iscolor
                    p.square_pix_val     = true;                             % [IF YES: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]
                else
                    p.square_pix_val     = false;
                end
                
                % TEMPORAL
                p.duration       = stimdur_frames;                          % frames (nr of monitor refreshes)
                
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

                % OBJECT FACING ROTATION ANGLES
                p.num_unique_objects       = 16;                            % nr of rotations (deg), 0 = rightward facing, 180 = leftward facing, 90 = forward facing
                p.facing_dir_deg           = linspace(10,170,p.num_unique_objects+1); % do not include 0-10, 170-180 deg to avoid edge cases in WM
                p.facing_dir_deg(p.facing_dir_deg==90) = [];                % remove 90 to avoid ill-defined rotation.
                p.facing_dir_deg           = ceil(p.facing_dir_deg/2)*2;    % force integrals of 2 as original images come in steps of 2 degrees
                
                % ensure equal distance from cardinal meridians
                assert(isequal(abs(p.facing_dir_deg(1:(p.num_unique_objects/2))-90), fliplr(abs(90-p.facing_dir_deg(((p.num_unique_objects/2)+1):p.num_unique_objects)))));
                assert(isequal(abs(45-p.facing_dir_deg(1:(p.num_unique_objects/2))), fliplr(abs(135-p.facing_dir_deg(((p.num_unique_objects/2)+1):p.num_unique_objects)))));
                
                % we carefully assign a unique rotation to an object
                obj_idx                    = [6,4,14,3,9,12,8,2,1,16,10,7,11,15,5,13]; 
                p.facing_dir_deg           = p.facing_dir_deg(obj_idx);
                
                % Define the 4 superordinate, 8 basic, and 16 subordinate categories
                p.super_cat          = {'human','animal','object','place'};     
                
                p.basic_cat{1}       = {'facemale','facefemale','facefemale'};
                p.basic_cat{2}       = {'small','small','large','large'};
                p.basic_cat{3}       = {'tool','tool','food','food','vehicle','vehicle'};
                p.basic_cat{4}       = {'building','building','tower'};
                
                p.sub_cat{1}         = {'damon','lisa','sophia'};
                p.sub_cat{2}         = {'parrot','cat','bear','giraffe'};
                p.sub_cat{3}         = {'drill','brush','pizza','banana','bus','suv'};
                p.sub_cat{4}         = {'church','house','watertower'};
                
                % Define the affordances for each object subcategory (for HOW task) 
                p.affordance{1}      = {'greet','greet','greet'};
                p.affordance{2}      = {'greet','greet','observe','observe'};
                p.affordance{3}      = {'grasp','grasp','grasp','grasp','enter','enter'};
                p.affordance{4}      = {'enter','enter','observe'};

                % WORKING MEMORY: dot angle position deltas for test images
                p.delta_from_ref     = [-8, -4, 4, 8];                      % Relative rotation from reference image for WM test image
                                                                            %  the bigger the delta, the easier the trial. 
                                                                            % for 1-89 deg rotations: Negative values are leftwards, positive values is rightwards
                                                                            % for 91-180 deg rotations: Negative values are rightward, positive values is leftward
                p.unique_im_nrs_WM   = [367:430];                           % Unique image nrs associated with the 64 WM OBJ test stimuli
                                                                            
                % check if all test images for WM have unique rotations
                tmp          = p.facing_dir_deg+[0, p.delta_from_ref]';
                assert(isequal(length(unique(tmp(:))), length(tmp(:))));
                clear tmp
                                                       
                % IMAGERY: 
                p.imagery_im_nr = p.unique_im_nrs([1,3,5,7,8,10,12,14]);   % SELECTED IMAGES USED (SUBSET of 16 IMAGES)--these are hand picked!
                
                % IMAGERY QUIZ DOT PARAMS
                p.imagery_sz_deg      = 5.658;                              % desired diameter (degree) of the second, quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix      = (round(p.img_sz_deg* disp_params.ppd)/2)*2; % pixel diameter of quiz dot image (ensure even nr of pixels)
                p.unique_im_nrs_IMG_test = [790:849];                      % Unique image nrs associated with the 8*20=60 IMG OBJ test dot images

                % LTM PAIR
                p.ltm_pairs          = [];                                       %%
                
                % Add params to struct
                stim.obj = p;
                
            case {'ns',5}
                
                % GENERAL
                p.indivscenefile = fullfile(vcd_rootPath,'workspaces','stimuli','RAW','vcd_natural_scenes'); % mat-file Where to load indiv pngs?
                p.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli',disp_params.name,sprintf('scene_%s',disp_params.name)); % mat-file Where to store stimulus images?
                p.infofile = fullfile(vcd_rootPath,'workspaces','info',sprintf('scene_info_%s',disp_params.name)); % csv-file Where to find stimulus info?
                p.unique_im_nrs      = [81:110];                             % Unique image nrs associated with the 16 single dot stimuli

                p.iscolor = true;                                            % Use color or not? 
                if p.iscolor
                    p.square_pix_val     = true;                             % [IF YES: MAKE SURE TO SQUARE IMAGE VALS FOR CLUT]
                else
                    p.square_pix_val     = false;
                end
                
                % TEMPORAL
                p.duration              = stimdur_frames;                          % frames (nr of monitor refreshes)
                
                % SPATIAL
                p.og_res_stim           = 425;                                     % original resolution of NSD stimuli
                p.rz_res_stim           = 714;                                     % resized  resolution of NSD stimuli
                p.dres                  = ((ctr_square_deg/disp_params.h_deg * disp_params.h_pix) / p.og_res_stim);  % scale factor to apply
                
                p.img_sz_deg            = ctr_square_deg;                          % desired height (or width) of square stimulus support (deg)
                p.img_sz_pix            = ceil(p.og_res_stim.*p.dres);             % height (or width) of square stimulus support (pix)
                p.square_pix_val        = true;                                    % MAKE SURE TO SQUARE IMAGE VALS FOR LINEAR CLUT

                p.x0_deg                = 0;                                       % desired x-center loc in deg (translation from 0,0)
                p.y0_deg                = 0;                                       % desired y-center loc in deg (translation from 0,0)
                p.x0_pix                = round((p.x0_deg.*disp_params.ppd)/2)*2;  % x-center loc in pix (translation from 0,0)
                p.y0_pix                = round((p.y0_deg.*disp_params.ppd)/2)*2;  % y-center loc in pix (translation from 0,0)
                
                % 5 superordinate semantic object categories
                p.super_cat     = {'human','animal','food','object','place'};      
                
                % 2 basic semantic object categories (scene location)
                p.basic_cat{1}  = {'indoor','outdoor'};
                p.basic_cat{2}  = {'indoor','outdoor'};
                p.basic_cat{3}  = {'indoor','outdoor'};
                p.basic_cat{4}  = {'indoor','outdoor'};
                p.basic_cat{5}  = {'indoor','outdoor'};
                
                % 3 subordinate categories (object spatial location)
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
                
                % 4 affordance categories
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
                

                % FOR WM task crossing, we have manipulated the original
                % NSD image by adding or removing something in the image.
                % These changes can be obvious (easy) or subtle (hard) to
                % detect:
                p.change_im         = {'easy_add', 'hard_add','easy_remove', 'hard_remove'};
                p.unique_im_nrs_WM  = [431:550];                            % Unique image nrs associated with the 120 WM NS changed stimuli

                % FOR LTM incorrect trials, we have very similar looking images called "lures":
                p.lure_im       = {'lure01', 'lure02', 'lure03', 'lure04'};
                p.unique_im_nrs_LTM_lures  = [551:610];                    % Unique image nrs associated with the 15*4=60 WM NS lure images

                % LTM PAIR
                p.ltm_pairs     = [];                                       %%
                
                % IMAGERY 
                p.imagery_im_nr          = p.unique_im_nrs([2,4,5,8,10,11,13,15,18,20,21,23,26,27,30]); % Half of the images will be used for IMG/LTM pairing (these are carefully handpicked! see scene_info csv file)
                p.imagery_sz_deg         = p.img_sz_deg;                    % desired diameter (degree) of the second, quiz dots image in an imagery trial to encourage subjects to create a vidid mental image.
                p.imagery_sz_pix         = p.img_sz_pix;                    % diameter of quiz dot image (pixels) (we already ensured even nr of pixels)
                p.unique_im_nrs_IMG_test = [850:1149];                     % Unique image nrs associated with the 15*20=300 IMG NS test dot images
                
                % Add params to struct
                stim.ns = p;
                
                fprintf('*** %s: dres (scale factor) = %.4f ***\n',upper(stim_class{ii}),stim.ns.dres);
        end
        % Print out stimulus params
        fprintf('*** %s: Using a stimulus size of %2.2f deg (%3.2f pixels). ***\n', upper(stim_class{ii}), p.img_sz_deg, p.img_sz_pix);
        fprintf('*** %s: Using a stimulus duration of %d frames (%3.2f seconds), where one frame relates to %d monitor refreshes (%d Hz) ***\n', ...
            upper(stim_class{ii}), p.duration, p.duration*stim.framedur_s, stim.f2f, disp_params.refresh_hz);
    end

    % Store params if requested
    if store_params
        fprintf('[%s]: Storing params..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
        save(fullfile(saveDir,sprintf('stim_%s_%s.mat',disp_params.name,datestr(now,30))),'stim')
    end
end


return