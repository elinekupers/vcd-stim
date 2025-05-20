%% s_createStim.m
%
% Stand-alone script to create and store the stimuli presented in the VCD
% core experiment, as well as the background and small fixation circle at
% the center of the screen. This script will do the following: 
%
% % 1. Load or define display parameters.
%      Users can pick from a list of 4 environments:
%      * '7TAS_BOLDSCREEN32'  : BOLDscreen at the CMRR's 7 Tesla Actively Shielded MRI.
%      * 'PPROOM_EIZOFLEXSCAN': Eizo Flexscan monitor at the CMRR's psychophysics lab.
%      * 'KKOFFICE_AOCQ3277'  : AOC monitor in Kay office at CMRR (only used for testing purposes).
%      * 'EKHOME_ASUSVE247'   : ASUS monitor in Eline's home (only used for testing purposes).
%   This adjust stimulus parameters such that stimuli are the intended size.
% 2. Create background images (BCKGRND)
% 3. Create small, central, fixation circle images (FIX)
% 4. Create gabor images (GBR)
% 5. Create rdk movies (RDK)
% 6. Create single dot images (DOT)
% 7. Create object images (OBJ)
% 8. Create scene images (NS)

%% %%%%%%%%%%%%%%%%%%%
%%%%%% PARAMETERS %%%% 
%%%%%%%%%%%%%%%%%%%%%%

params = struct();
params.verbose      = true; % visualize stimuli (true) or not (false)
params.store_imgs   = true; % store visualization figures (true) or not (false)

% Get display params
% Users can pick from a list of 4 environments (or add their own):
%  * '7TAS_BOLDSCREEN32'  : BOLDscreen at the CMRR's 7 Tesla Actively Shielded MRI.
%  * 'PPROOM_EIZOFLEXSCAN': Eizo Flexscan monitor at the CMRR's psychophysics lab.
%  * 'KKOFFICE_AOCQ3277'  : AOC monitor in Kay office at CMRR (only used for testing purposes).
%  * 'EKHOME_ASUSVE247'   : ASUS monitor in Eline's home (only used for testing purposes).
% This function will use field of view and native refresh rate of monitor
% to adjust stimulus pixel size and frame duration.
dispname    = 'PPROOM_EIZOFLEXSCAN';
params.disp = vcd_getDisplayParams(dispname);

% Where to store visualization of stimuli (debug figures and PNGs)?
saveFigsFolder = fullfile(vcd_rootPath,'figs',dispname); 
if ~exist(saveFigsFolder,'dir'); mkdir(saveFigsFolder); end

% Get stimulus parameters
% When loading stored parameters from file, it will search for a file called
% "stim_<dispname>_*.mat" in fullfile(vcd_rootPath,'workspaces','info').
% When storing generated parameters, we save the parameters as a .mat file
% in "stim_<dispname>_YYYYMMDDTHHMMSS.mat" in fullfile(vcd_rootPath,'workspaces','info').
params.load_params  = false;  % if true, we load from file. if false, define params.
params.store_params = true;  % if false, we don't store params. if true, we store mat file in fullfile(vcd_rootPath,'workspaces','info')
                                       
% Reset random number generator with arbitrary number (based on system
% clock). Early versions of Matlab use different generators for rand and
% randn, hence we do it separately for each.
rand('seed', sum(100*clock));
params.rng.rand_seed  = rng;   % store rand seed
randn('seed', sum(100*clock)); 
params.rng.randn_seed  = rng;  % store randn seed

%% Define/Load stimulus params 
% Load stimulus params, all parameters are deterministic (i.e., there is no
% randomization involved in setting the stimulus parameters). This function
% well get appropriate display params by calling vcd_getDisplayParams.m.

params.stim   = vcd_getStimParams('disp_name', params.disp.name, ...       
                                  'load_params', params.load_params, ...  
                                  'store_params', params.store_params);  
                              
%% Define/Load experiment session params
params.exp    = vcd_getSessionParams('disp_name', params.disp.name, ...
                                'presentationrate_hz',params.stim.presentationrate_hz, ...
                                'load_params', params.load_params, ...
                                'store_params', params.store_params);
%% %%%%%%%%%%%%%%%%
%%%%%% STIMULI %%%% 
%%%%%%%%%%%%%%%%%%%

%% Background
% We create a pinknoise (1/f) image for the entire monitor size, and
% superimpose a mean gray luminance cutout. The center of this cutout can
% be adjusted with pixoffset parameter.
% If params.verbose = true, this function will make a PNG of each
% background image, and a figure plotting image with image nr and axes.
% If params.store_imgs = true, this function will store these figures in
% fullfile(vcd_rootPath,'figs',dispname,'background').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
% * cutout_type : what stimulus apertures to use to create cutout shape, choose from:
%   - 'puzzle'    (union of center square and 2 parafoveal stimulus apertures in 2-stimulus array)
%   - 'dotring'   (iso-eccentric ring that the single dot lives on) 
%   - 'comb'      (union of puzzle and dotring)
% * rim width   : what rim do you want, choose from: 
%   - 'skinny'    (no additional "buffer zone" between  and pink noise background
%   - 'fat'       (with additional "buffer zone": 1 degree added on each side of the cutout 
% * num         : number of unique noise images (should be an integer, default is 1)
% * pixoffset   : relative offset of [x,y] center in pixels from the native 
%                 center of the monitor (BOLDscreen width 540 x height 960
%                 pixels). Default is [0,0] pixels.
%
% OUTPUT: 
% * bckgrnd_im  : (uint8) background images, height in pixels x width in pixels x 3 (rbg) x number of images (int) 
%                 BOLDscreen dimensions are: height (1080 pixels) x width (1920 pixels)

gaptype     = 'circle';
borderwidth = 'fat';
if strcmp(dispname,'PPROOM_EIZOFLEXSCAN')
    % 14 runs * 1 session x 2 session type ("BEHAVIOR001")
    num = sum(params.exp.session.behavior.n_runs_per_session); 
elseif strcmp(dispname,'7TAS_BOLDSCREEN32')
    % 258 MRI runs: 10 runs x 2 (WIDE01A + WIDE01B) + 10 runs x 25 sessions (DEEP001-0025) + 4 runs x 2 sessions (DEEP26A/26Bp)
    num = sum(params.exp.session.mri.wide.n_runs_per_session) + sum(params.exp.session.mri.deep.n_runs_per_session);
else
    num = 1;
end
bckgrnd_im  = vcd_pinknoisebackground(params, ...
                                     'gaptype', gaptype, ...
                                     'borderwidth', borderwidth,...
                                     'num', num, ...
                                     'pixoffset', [0,0]); 

%% Eyetracking targets

vcd_createEyeTrackingBlockTargets(params)
                                 
%% Fixation circle
% vcd_fixationDot function creates 25 types of fixation circles, the full
% crossing between 5 inner circle luminance levels and 5 rim types. 
% If params.verbose = true, this function will make a PNG of each
% fixation dot, and a figure with all dots in a 5x5 array.
% If params.store_imgs = true, this function will store these figures in
% fullfile(vcd_rootPath,'figs',dispname,'fix').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
%
% OUTPUTS:
% * fix_im      : (uint8) unique fixation circle images, 5D array: 
%                   height (24 pixels) x width (24 pixels) x 3 (rgb) 
%                   x 5 luminance levels x 5 dot rims types 
%                   Rim types are 1: thin white, 2: thick white, 
%                   3: thick-red left, 4: thick-red right, 5: thick-red both
% * fix_mask    : (uint8) alpha transparency masks, 4D array: 
%                   height (24 pixels) x width (24 pixels) x 2 dot rims types (thin, thick)

[fix_im, fix_mask, fix_info] = vcd_fixationDot(params);

%% Gabors
% vcd_gabor function creates 120 Gabor stimuli: 24 core and 96 working
% memory test images.
% If params.verbose = true, this function will make a PNG of each
% gabor image, and a figure plotting gabors with image nr and axes as well 
% as histograms of the pixel luminance.
% If params.store_imgs = true, this function will store these figures in
% fullfile(vcd_rootPath,'figs',dispname,'gabor').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
%
% OUTPUTS:
% * gabors      : (uint8) 120 unique gabor images, 5D array: 
%                   height (354 pixels) x width (354 pixels) x 3 (rgb) 
%                   x 24 unique images (8 ori x 3 contrasts) 
%                   x 5 orientation tilt offsets (0 + -15, -5, +5, +15 deg)
% * masks       : (uint8) alpha transparency masks used by Psychtoolbox to
%                   crop the edges of the square support, 4D array:
%                   height (354 pixels) x width (354 pixels) 
%                   x 24 unique images x 5 tilt offsets
% * info        : table with stimulus information matching the gabor array
[gabors, masks, info] = vcd_gabor(params);

%% RDKs (Random Dot motion Kinetograms)
% vcd_rdk function creates 120 RDK movies: 24 core and 96 working memory
% test images.
%
% If params.verbose = true, this function will make a PNG of each RDK movie
% frame, a MP4 movie per RDK, and figure plotting dot motion vectors for
% each movie frame. 
% If params.store_imgs = true, this function will store
% these figures in fullfile(vcd_rootPath,'figs',dispname,'rdk').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
%
% OUTPUTS:
% * rdks        : (cell) 3D cell array with RDK stimuli: 
%                  8 motion directions x 3 coherence levels x 5 motion direction offsets (0, -15, -5, +5, +15 deg). 
%                  Each cell contains movie frames (uint8):
%                  For BOLDscreen:
%                   height (548 pixels) x width (548 pixels for BOLDscreen, 400 pixels for EIZOflexscan) 
%                   x 3 (rgb) x 30 frames (total of 2 s, 33 ms per frame).
% * masks       : (cell) 3D cell with alpha transparency masks:
%                   8 motion directions x 3 coherence levels x 5 motion direction offsets (0, -15, -5, +5, +15 deg). 
%                   Each cell contains one uint8 mask image: 
%                   For BOLDscreen:
%                     height (548 pixels) x width (548 pixels)
% * info        : table with stimulus information matching the rdk array 

[rdks, masks, info] = vcd_rdk(params);

%% Single dot
% vcd_singledot function creates 1 dot that will be used for 16 core and 64
% working memory test images.
%
% If params.verbose = true, this function will make a PNG for each dot 
% location (core and WM test), and figures to check dot WM test locations 
% relative to core image and core dot locations relative to one another.
% If params.store_imgs = true, this function will store
% these figures in fullfile(vcd_rootPath,'figs',dispname,'dot').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
%
% OUTPUTS:
% * single_dot  : (uint8) matrix with single dot image:
%                 For BOLDscreen:
%                   height (94 pixels) x width (94 pixels) x 3 (rgb).
% * masks       : (uint8) matrix with single alpha transparency mask:
%                 For BOLDscreen:
%                   height (94 pixels) x width (94 pixels)
% * info        : table with stimulus information about the dot locations

[single_dot, masks, info] = vcd_singledot(params);

%% Objects
% vcd_objects function creates 80 objects will be used for 16 core and 64
% working memory test images.
%
% If params.verbose = true, this function will make a PNG for each object 
% (core and WM test) and figures with axes/titles to check image nrs,
% rotation.
% If params.store_imgs = true, this function will store
% these figures in fullfile(vcd_rootPath,'figs',dispname,'obj').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
%
% OUTPUTS:
% * objects     : (uint8) is a 5D array:
%                 For BOLDscreen:
%                   height (1024 pixels) x width (1024 pixels) x 3 (rgb) 
%                   x 16 object categories (subordinate level) 
%                   x 4 rotation offsets (0, -8, -4, +4, +8 deg). 
% * masks       : (uint8) is a 4D array containing alpha transparency mask:
%                 For BOLDscreen:
%                   height (1024 pixels) x width (1024 pixels) 
%                   x 16 object categories (subordinate level) 
%                   x 5 rotation offsets (0, -8, -4, +4, +8 deg). 
%  info         :   Loaded csv from workspaces/stimuli/ with png file names
%                   and object category information. 

[objects, masks, info] = vcd_objects(params);

%% Natural scenes
% vcd_naturalscenes function loads, resizes (and squares pixel values if
% using linearized display) for 30 core scenes, 120 working memory test
% images, 80 long-term memory lure images.  Note: There are no alpha masks
% for this stimulus class.
%
% If params.verbose = true, this function will make a PNG for each scene 
% (core, WM test, and LTM lute) and figures with axes/titles to check image 
% nrs, pixel values and scene size.
% If params.store_imgs = true, this function will store
% these figures in fullfile(vcd_rootPath,'figs',dispname,'ns').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
%
% OUTPUTS:
% * scenes      : (uint8) core scene images, 6D array:
%                 For BOLDscreen:
%                   height (743 pixels) x width ( 743 pixels) x 3 (rgb)
%                   x 5 superordinate semantic object categories (human, animal, object, food, place)
%                   x 2 scene locations (indoor/outdoor)
%                   x 3 obj locations (left/middle/right)
% * ltm_lures   : (uint8) long term memory lure scenes, 7D array:
%                 For BOLDscreen:
%                   height (743 pixels) x width (743 pixels) x 3 (rgb)
%                   x 5 superordinate semantic object categories 
%                   x 2 scene locations (indoor/outdoor)
%                   x 3 obj locations (left/middle/right) 
%                   x 4 lure types:
%                       1: very similar/difficult
%                       2: somewhat similar
%                       3: somewhat different
%                       4: very different/easy
% * wm_im       : (uint8) working memory test scenes, 7D array:
%                 For BOLDscreen:
%                   height (743 pixels) x width (743 pixels) x 3 (rgb) 
%                   x 5 superordinate semantic object categories 
%                   x 2 scene locations (indoor/outdoor) 
%                   x 3 obj locations (left/middle/right)
%                   x 4 change types:
%                       1:easy add -- scene is altered by adding something big/obvious 
%                       2:hard add -- scene is altered by adding something small/subtle
%                       3:easy remove -- scene is altered by removing something big/obvious
%                       4:hard remove -- scene is altered by removing something small/subtle

[scenes, ltm_lures, wm_im, info] = vcd_naturalscenes(params);

