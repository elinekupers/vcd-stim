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

% If params.verbose = true and params.store_imgs =  true, every stimulus
% creation function will store pngs. You can define here the path to store 
% the pngs of stimuli in "saveFigsFolder" (both debug figures and PNGs).
saveFigsFolder = fullfile(vcd_rootPath,'figs',dispname); 
if ~exist(saveFigsFolder,'dir'); mkdir(saveFigsFolder); end
params.saveFigsFolder = saveFigsFolder;

% Get stimulus parameters
% When loading stored parameters from file (params.load_params = true), 
% vcd_getStimParams will search for a file called
% "stim_<dispname>_*.mat" in fullfile(vcd_rootPath,'workspaces','info').
% When storing generated parameters, we save the parameters as a .mat file
% in "stim_<dispname>_YYYYMMDDTHHMMSS.mat" in fullfile(vcd_rootPath,'workspaces','info').
params.load_params  = false;  % if true, we load from file. if false, define params.
params.store_params = true;  % if false, we don't store params. if true, we store mat file in fullfile(vcd_rootPath,'workspaces','info')
                                       
% Reset random number generator with arbitrary number (based on system
% clock). Early versions of Matlab use different generators for rand and
% randn, hence we do it separately for each. (Needed for RDKs).
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


%% Eyetracking targets
% There are 5 saccade targets (params.exp.block.nr_of_saccades): central, 
% left, right, up, down, and 1 pupil trial with a mid-gray central fixation 
% target on a black and white background.
%
% The distance between the center and 4 left/right/up/down targets are set
% as [xc,yc] ± 264 pixels (BOLDscreen) or ± 192 pixels (EIZOFLEXSCAN). 
% This results in dots at the following pixel coordinates for the 5 targets
% BOLDscreen coordinates in pixels:
%                    [x3,y3]=[960,376]
% [x1,y1]=[696,640]  [x0,y0]=[960,640]   [x2,y2]=[1224,640]
%                    [x4,y4]=[960,904]
%
% EIZOFLEXSCAN coordinates in pixels:
%                    [x3,y3]=[960,408]
% [x1,y1]=[768,600]  [x0,y0]=[960,600]   [x2,y2]=[1156,600]
%                    [x4,y4]=[960,792]
%
% EMPIRICAL target distance:
% * BOLDscreen: 264 pixels, which corresponds to 2.9936 degrees.
% * PP room EIZOFLEX: 192 pixels, which corresponds to 3.0046 degrees.
%
% INPUTS:
%  params           : (struct) parameter struct, which should contain the following fields:
%
% OUTPUTS:
%  sac_im           : (uint8) saccade stimuli (height in pixels x width in pixels x 3 x 5)
%  pupil_im_white   : (uint8) white background pupil trial stimulus (height in pixels x width in pixels x 3)
%  pupil_im_black   : (uint8) black background pupil trial stimulus (height in pixels x width in pixels x 3)

[sac_im,pupil_im_white,pupil_im_black] = vcd_createEyeTrackingBlockTargets(params);
                                 
%% Fixation circle
% vcd_fixationDot function creates 30 different fixation circle images, 
% which is the full crossing between 6 inner circle luminance levels and 5 
% rim types. If params.verbose = true, this function will make a PNG of 
% each fixation dot, and a figure with all dots in a 6x5 array.
% If params.store_imgs = true, this function will store these figures in
% fullfile(vcd_rootPath,'figs',dispname,'fix').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
%
% OUTPUTS:
% * fix_im      : (uint8) unique fixation circle images, 5D array: 
%                   height (24 pixels) x width (24 pixels) x 3 (rgb) 
%                   x 6 luminance levels x 5 dot rims types 
%                   Rim types are 1: thin white, 2: thick white, 
%                   3: thick-red left, 4: thick-red right, 5: thick-red both
% * fix_mask    : (uint8) alpha transparency masks, 4D array: 
%                   height (24 pixels) x width (24 pixels) x 2 dot rims types (thin, thick)

[fix_im, fix_mask, fix_info] = vcd_fixationDot(params);

%% Gabors
% vcd_gabor function creates 56 Gabor stimuli: 24 core and 32 working
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
% * gabors      : (uint8) 56 unique gabor images, 6D array: 
%                   height (352/256 pixels for MRI/PProom) 
%                   x width (352/256 pixels for MRI/PProom 
%                   x 3 (rgb) 
%                   x 8 orientations
%                   x 3 contrasts 
%                   x 5 orientation tilt offsets (0, -16, -8, +8, +16 deg)
%                   Note that dims gabors(:,:,:,:,[1,2],:) are empty
%                   as all wm test stimuli use the highest contrast level.
% * masks       : (uint8) 56 alpha transparency masks used by Psychtoolbox to
%                   crop the edges of the square support, 5D array:
%                   height (352/256 pixels for MRI/PProom) 
%                   x width (352/256 pixels for MRI/PProom) 
%                   x 8 orientations
%                   x 3 contrasts 
%                   x 5 orientation tilt offsets (0, -16, -8, +8, +16 deg).
%                   Note that dims masks(:,:,:,[1,2],:) are empty
%                   as all wm test stimuli use the highest contrast level.
% * info        : table with stimulus information matching the gabor array
[gabors, ~, ~] = vcd_gabor(params);

%% RDKs (Random Dot motion Kinetograms)
% vcd_rdk function creates 56 RDK movies: 24 core and 32 working memory
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
%                  8 motion directions x 3 coherence levels x 5 motion direction offsets (0, -20, -10, +10, +20 deg). 
%                  Each cell contains movie frames (uint8):
%                  For BOLDscreen:
%                   height (BOLDscreen: 546 pixels, EIZOflexscan: 396 pixels)
%                   x width (BOLDscreen: 546 pixels, EIZOflexscan: 396 pixels) 
%                   x 3 (rgb) x 30 frames (total of 2 s, 33 ms per frame).
% * masks       : (cell) 3D cell with alpha transparency masks:
%                   8 motion directions x 3 coherence levels x 5 motion direction offsets (0, -20, -10, +10, +20 deg). 
%                   Each cell contains one uint8 mask image: 
%                   For BOLDscreen:
%                     height (544 or 396 pixels) x width (544 or 396 pixels)
% * info        : table with stimulus information matching the rdk array 

[rdks, ~, ~] = vcd_rdk(params);

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
%                   height (BOLDscreen: 94 pixels, Eizoflexscan: 70 pixels)
%                   x width (BOLDscreen: 94 pixels, Eizoflexscan: 70 pixels) 
%                   x 3 (rgb).
% * masks       : (uint8) matrix with single alpha transparency mask:
%                 For BOLDscreen:
%                   height (BOLDscreen: 94 pixels, Eizoflexscan: 70 pixels) 
%                   x width (BOLDscreen: 94 pixels, Eizoflexscan: 70 pixels)
% * info        : table with stimulus information about the dot locations

[single_dot, ~, ~] = vcd_singledot(params);

%% Objects
% vcd_objects function creates 80 objects will be used for 16 core and 64
% working memory test images.
%
% If params.verbose = true, this function will make a PNG for each object 
% (core and WM test) and figures with axes/titles to check image nrs,
% rotation.
% If params.store_imgs = true, this function will store
% these figures in fullfile(vcd_rootPath,'figs',<dispname>,'obj').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
%
% OUTPUTS:
% * objects     : (uint8) is a 5D array:
%                 For BOLDscreen:
%                   height (1024 pixels) x width (1024 pixels) x 3 (rgb) 
%                   x 16 object categories (subordinate level) 
%                   x 4 rotation offsets (0, -24, -12, +12, +24 deg). 
% * masks       : (uint8) is a 4D array containing alpha transparency mask:
%                 For BOLDscreen:
%                   height (1024 pixels) x width (1024 pixels) 
%                   x 16 object categories (subordinate level) 
%                   x 5 rotation offsets (0, -24, -12, +12, +24 deg). 
% * info        : (table) information about object png filenames,
%                  category information, and object rotation. Also stored
%                  as csv info file in params.stim.obj.infofile.

[objects, ~, ~] = vcd_objects(params);

%% Natural scenes
% vcd_naturalscenes function loads, resizes (and squares pixel values if
% using linearized display) for 30 core scenes, 120 working memory test
% images, 80 long-term memory lure images.  Note: There are no alpha masks
% for this stimulus class.
%
% If params.verbose = true, this function will make a PNG for each scene 
% (core, WM test, and LTM novel lure scenes) and figures with axes/titles 
% to check image nrs, pixel values and scene size.
% If params.store_imgs = true, this function will store
% these figures in fullfile(vcd_rootPath,'figs',dispname,'ns').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
%
% OUTPUTS:
% * scenes      : (uint8) core scene images, 6D array:
%                 For :
%                   height (BOLDscreen: 743 pixels, Eizoflexscan: 541 pixels)
%                   x width (BOLDscreen: 743 pixels, Eizoflexscan: 541 pixels)
%                   x 3 (rgb)
%                   x 5 superordinate semantic object categories (human, animal, object, food, place)
%                   x 2 scene locations (indoor/outdoor)
%                   x 3 obj locations (left/middle/right)
% * ltm_lures   : (uint8) long term memory lure scenes, 7D array:
%                   height (BOLDscreen: 743 pixels, Eizoflexscan: 541 pixels)
%                   x width (BOLDscreen: 743 pixels, Eizoflexscan: 541 pixels)
%                   x 3 (rgb)
%                   x 5 superordinate semantic object categories 
%                   x 2 scene locations (indoor/outdoor)
%                   x 3 obj locations (left/middle/right) 
%                   x 4 lure types:
%                       1: very similar/difficult
%                       2: somewhat similar
%                       3: somewhat different
%                       4: very different/easy
% * wm_im       : (uint8) working memory test scenes, 7D array:
%                   height (BOLDscreen: 743 pixels, Eizoflexscan: 541 pixels)
%                   x width (BOLDscreen: 743 pixels, Eizoflexscan: 541 pixels)
%                   x 3 (rgb)
%                   x 5 superordinate semantic object categories 
%                   x 2 scene locations (indoor/outdoor) 
%                   x 3 obj locations (left/middle/right)
%                   x 4 change types:
%                       1:easy add -- scene is altered by adding something big/obvious 
%                       2:hard add -- scene is altered by adding something small/subtle
%                       3:easy remove -- scene is altered by removing something big/obvious
%                       4:hard remove -- scene is altered by removing something small/subtle
% * info        : (table) information about scene png filenames and
%                  category information. Also stored
%                  as csv info file in params.stim.ns.infofile.

[scenes, ltm_lures, wm_im, ~] = vcd_naturalscenes(params);




%% -- Background -- OBSOLETE!!!
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
%
% * gaptype     : what stimulus apertures to use to create cutout shape, choose from:
%   - 'none'        : gray mean luminance background, no pink noise
%   - 'circle'      : circular shape 
%   - 'puzzle'      : union of center square and 2 parafoveal stimulus apertures in 2-stimulus array)
%   - 'dotring'     : iso-eccentric ring that the single dot lives on.
%   - 'comb'        : union of puzzle and dotring)
%
% * borderwidth : how much buffer zone do you want between stimulus edge and pink noise, choose from: 
%   - 'skinny'    (no additional "buffer zone" between and pink noise background
%   - 'fat'       (with additional "buffer zone": 1 degree added on each side of the cutout 
%   if gaptype = 'none', borderwidth will be ignored.
%
% * num         : number of unique noise images (should be an integer, default is 1)
% * pixoffset   : relative offset of [x,y] center in pixels from the native 
%                 center of the monitor (BOLDscreen width 540 x height 960
%                 pixels). Default is [0,0] pixels.
%
% OUTPUT: 
% * bckgrnd_im  : (uint8) background images, height in pixels x width in pixels x 3 (rbg) x number of images (int) 
%                 BOLDscreen dimensions are: height (1080 pixels) x width (1920 pixels)

% Define inputs
% gaptype     = 'none';
% borderwidth = 'none';
% 
% if strcmp(dispname,'PPROOM_EIZOFLEXSCAN')
%     % 14 runs * 1 session x 2 session type ("BEHAVIOR001")
%     num = 1; %sum(params.exp.session.behavior.n_runs_per_session); 
% elseif strcmp(dispname,'7TAS_BOLDSCREEN32')
%     % 258 MRI runs: 10 runs x 2 (WIDE01A + WIDE01B) + 10 runs x 25 sessions (DEEP001-0025) + 4 runs x 2 sessions (DEEP26A/26Bp)
%     num = 1; %sum(params.exp.session.mri.wide.n_runs_per_session) + sum(params.exp.session.mri.deep.n_runs_per_session);
% else
%     num = 1;
% end
% bckgrnd_im  = vcd_pinknoisebackground(params, ...
%                                      'gaptype', gaptype, ...
%                                      'borderwidth', borderwidth,...
%                                      'num', num, ...
%                                      'pixoffset', [0,0]); 