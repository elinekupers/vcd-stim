%% s_createStim.m
%
% Stand-alone script to create and store the stimuli presented in the VCD
% core experiment, as well as the background and small fixation circle at
% the center of the screen. This script will do the following: 
%
% % 1. Define display to load correct parameters.
%      Users can pick from a list of 5 displays (or add their own in vcd_getDisplayParams):
%      * '7TAS_BOLDSCREEN32'  : BOLDscreen at the CMRR's 7 Tesla Actively Shielded MRI.
%      * 'PPROOM_EIZOFLEXSCAN': Eizo Flexscan monitor at the CMRR's psychophysics lab.
%      * 'KKOFFICE_AOCQ3277'  : AOC monitor in Kay office at CMRR (only used for testing purposes).
%      * 'EKHOME_ASUSVE247'   : ASUS monitor in Eline's home (only used for testing purposes).
%      * 'CCNYU_VIEWPIXX3D'   : ViewPixx monitor in Curtis Lab psychophysics room.
%      The vcd_getDisplayParams function will use field of view and native 
%      refresh rate of monitor to adjust stimulus pixel size and time frame duration.
% 2. Create eyetracking targets (ET)
% 3. Create small, central, fixation circle images (FIX)
% 4. Create gabor images (GBR)
% 5. Create rdk movies (RDK)
% 6. Create single dot images (DOT)
% 7. Create object images (OBJ)
% 8. Create scene images (NS)


%% %%%%%%%%%%%%%%%%%%%
%%%%%% PARAMETERS %%%% 
%%%%%%%%%%%%%%%%%%%%%%

% Get general parameters
% * When loading stored parameters from file (params.load_params = true), 
% vcd_getStimParams will search for a file called "stim_<dispname>_*.mat" 
% in fullfile(vcd_rootPath,'workspaces','info').
% * When storing generated parameters, we save the parameters as a .mat file
% in "stim_<dispname>_YYYYMMDDTHHMMSS.mat" in fullfile(vcd_rootPath,'workspaces','info').
verbose      = true;  % visualize stimuli (true) or not (false)
store_imgs   = true;  % store visualization figures (true) or not (false)
load_params  = false; % if true, we load from file. if false, define params.
store_params = false; % if false, we don't store params. if true, we store mat file in fullfile(vcd_rootPath,'workspaces','info')

% Create params struct
params = struct();

% Get display params
dispname    = '7TAS_BOLDSCREEN32';
params.disp = vcd_getDisplayParams(dispname);

% If params. and params.store_imgs =  true, every stimulus
% creation function will store pngs. You can define here the path to store 
% the pngs of stimuli in "saveFigsFolder" (both debug figures and PNGs).
if verbose && store_imgs
    params.saveFigsFolder = fullfile(vcd_rootPath,'figs',dispname);
    if ~exist(params.saveFigsFolder,'dir'); mkdir(params.saveFigsFolder); end
else
    params.saveFigsFolder = [];
end
                                       
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
                                  'load_params', load_params, ...  
                                  'store_params', store_params);  
                              
%% Define/Load experiment session params
params.exp    = vcd_getSessionParams('disp_name', params.disp.name, ...
                                'presentationrate_hz',params.stim.presentationrate_hz, ...
                                'load_params', load_params, ...
                                'store_params', store_params);
%% %%%%%%%%%%%%%%%%
%%%%%% STIMULI %%%% 
%%%%%%%%%%%%%%%%%%%


%% Eyetracking targets
% There are 5 saccade targets (params.exp.block.nr_of_saccades): central, 
% left, right, up, down, and 1 pupil trial with a mid-gray central fixation 
% target on a black and white background.
%
% DESIRED target distance:
% Between the center and 4 left/right/up/down targets is 4 degrees visual
% angle. This results in [xc,yc] ± 352 pixels (7TAS BOLDscreen) or [xc,yc] 
% ± 265 pixels (PPROOM EIZOFLEXSCAN).
%
% EMPIRICAL target distance:
% * 7TAS BOLDscreen: 352 pixels, which corresponds to 3.9914 degrees.
% * PPROOM EIZOFLEXSCAN: 265 pixels, which corresponds to 4.0061 degrees.
%
% This results in dots at the following pixel coordinates for the 5 targets
% BOLDscreen target rect coordinates in pixels 
% [x1,y1,x1,y2] = [top-left-x, top-left-y, bottom-right-x bottom-right-y]:
%
%                    [949,177,972,200]
% [597,529,620,552]  [949,529,972,552]   [1301,529,1324, 552]
%                    [949,881,972,904]
%
% EIZOFLEXSCAN target rect coordinates in pixels 
% [x1,y1,x1,y2] = [top-left-x, top-left-y, bottom-right-x bottom-right-y]:
%
%                    [951,335,969,353]
% [695,591,713,609]  [951,591,969,609]   [1207,591,1225,609]
%                    [951,847,969,865]
%
% INPUTS:
% * params         : (struct) parameter struct.
% * verbose        : (logical) show debug figures
% * store_imgs     : (logical) store stimuli and debug figures as pngs 
%
% OUTPUTS:
% * sac_im         : (uint8) saccade stimuli
%                     height (BOLDscreen 1080 pixels; PProom: 1200 pixels)
%                     x width (BOLDscreen 1920 pixels; PProom: 1920 pixels)
%                     x 3 (RGB) 
%                     x 5 targets (1: center, 2: left, 3: right, 4: up, 5:
%                     down).
% * pupil_im_white : (uint8) white background pupil trial stimulus 
%                     height (BOLDscreen 1080 pixels; PProom: 1200 pixels)
%                     x width (BOLDscreen 1080 pixels; PProom: 1200 pixels)
%                     x 3 (RGB) 
% * pupil_im_black : (uint8) black background pupil trial stimulus
%                     height (BOLDscreen 1080 pixels; PProom: 1200 pixels)
%                     x width (BOLDscreen 1080 pixels; PProom: 1200 pixels)
%                     x 3 (RGB)

[sac_im, pupil_im_white, pupil_im_black] = vcd_createEyeTrackingBlockTargets(params, verbose, store_imgs);
                                 
%% Fixation circle
% vcd_fixationDot function creates 36 different fixation circle images, 
% which is the full crossing between 6 inner circle luminance levels and 6 
% rim types. If params.verbose = true, this function will make a PNG of 
% each fixation dot, and two figures with all dots in a 6x6 array (one 
% mimicking the 50% transparency, the other one not).
% If params.store_imgs = true, this function will store these figures in
% fullfile(vcd_rootPath,'figs',<dispname>,'fix').
%
% INPUTS:
% * params      : parameter struct (requires display and stimulus parameters)
% * verbose     : (logical) show debug figures
% * store_imgs  : (logical) store stimuli and debug figures as pngs 
%
% OUTPUTS:
% * fix         : (uint8) unique fixation circle images, 5D array: 
%                   height (BOLDscreen: 24 pixels) 
%                   x width (BOLDscreen: 24 pixels) 
%                   x 3 (rgb) 
%                   x 6 luminance levels: corresponding to uint8 values [6,
%                   27, 64, 116, 184, 128] (in that particular order).
%                   x 6 dot rims types, which are ordered as followed:
%                     1: thin rim white both, 
%                     2: thick rim white both, 
%                     3: thick rim red left, 
%                     4: thick rim red right, 
%                     5: thick rim red both
%                     6: thick rim black both
% * masks       : (uint8) alpha transparency masks, 4D array: 
%                   height (BOLDscreen: 24 pixels) 
%                   x width (BOLDscreen: 24 pixels) 
%                   x 6 dot rims types (same order as above)
% * info        : (table) 36x4 table with specs about the different types 
%                   of fixation circles. Stored as csv info file.
%                   Columns include:
%                     1: luminance  - uint8 luminance values
%                     2: rim_width  - 'thin' or 'thick'
%                     3: color_side - 'both', 'left', or 'right
%                     4: rim_color  - cell with double: [255, 255, 255]
%                     (white), [255 0 0] (red), or [0 0 0] (black).
[fix, ~, ~] = vcd_fixationDot(params, verbose, store_imgs); % outputs are [fix, masks, info]

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
% * verbose     : (logical) show debug figures
% * store_imgs  : (logical) store stimuli and debug figures as pngs 
%
% OUTPUTS:
% * gabors      : (uint8) 56 unique gabor images, 6D array: 
%                   height (BOLDscreen: 352 pixels; PProom: 256 pixels) 
%                   x width (BOLDscreen: 352 pixels; PProom: 256 pixels)
%                   x 3 (rgb) 
%                   x 8 orientations (11.25, 33.75, 56.25, 78.75, 101.25,
%                   123.75, 146.25, 168.75 deg)
%                   x 3 contrasts (low: 5%, medium: 20%, high: 80%)
%                   x 5 orientation tilt offsets (0, -16, -8, +8, +16 deg)
%                   Note that dims gabors(:,:,:,:,[1,2],[2:5]) are all zeros
%                   as wm test stimuli use the highest contrast level.
% * masks       : (uint8) 56 alpha transparency masks used by Psychtoolbox 
%                   to crop the edges of the square support, 5D array:
%                   height (BOLDscreen: 352 pixels; PProom: 256 pixels) 
%                   x width (BOLDscreen: 352 pixels; PProom: 256 pixels) 
%                   x 8 orientations
%                   x 3 contrasts 
%                   x 5 orientation tilt offsets (0, -16, -8, +8, +16 deg).
%                   Note that dims masks(:,:,:,[1,2],[2:5]) are all zeros
%                   as wm test stimuli use the highest contrast level.
% * info        : 56x12 table with stimulus information matching the gabor 
%                   array. Also stored as csv info file in 
%                   params.stim.gabor.infofile.  Columns include:
%                   1: unique_im  - (double) unique image nr for each Gabor: 
%                                   range 1-24, 111-142 for wm test images.
%                   2: stim_pos   - (double) stimulus position index. 
%                                   1=left, 2=right
%                   3: stim_pos_i - (cell) stimulus position, same as 
%                                   stim_pos_i but human-readable 
%                                   ({'left'} or {'right'})
%                   4: orient_deg - (double) tilt in degrees (0=12 o'clock).
%                   5: orient_i   - (double) same as orient_deg but indexed 
%                                   1: 11.25 deg through 8: 168.75 deg
%                   6: contrast   - (double) stimulus contrast level (fraction of 1),   
%                   7: contrast_i - (double) same as 'contrast' but indexed 1 
%                                   (0.05), 2 (0.2), or 3 (0.8).
%                   8: phase_deg  - (double) gabor phase (deg), one of four 
%                                   quadrature phases (0, 90, 180, 270).
%                   9: phase_i    - (double) same as phase_deg but indexed  
%                                   1 (0 deg) to 4 (270 deg).
%                   10: delta_deg - (double) Gabor orientation relative from
%                                   corresponding core Gabor: 0, -16, -8, +8, +16 (deg)   
%                   11: delta_i   - (double) same as delta_deg but indexed
%                                   0 (0 deg), 1 (-16 deg), 2 (-8 deg), 
%                                   3 (+8 deg), or 4 (+16 deg).
%                   12: is_specialcore - (logical) whether the Gabor 
%                                   stimulus is part of the subselected 
%                                   stimuli used in imagery and long-term 
%                                   memory task.
[gabors, ~, ~] = vcd_gabor(params, verbose, store_imgs); % outputs are [gabors, masks, info]

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
% * verbose     : (logical) show debug figures
% * store_imgs  : (logical) store stimuli and debug figures as pngs 
%
% OUTPUTS:
% * rdks        : (cell) 3D cell array with RDK stimuli: 
%                  8 motion directions (33.75, 78.75, 123.75, 168.75,
%                  213.75, 258.75, 303.75, 348.75 deg)
%                  x 3 coherence levels (0.3, 0.65, 1.0)
%                  x 5 motion direction offsets (0, -20, -10, +10, +20 deg). 
%                  Each cell contains movie frames (uint8), with dimensions:
%                  height (BOLDscreen: 546 pixels, EIZOflexscan: 396 pixels)
%                  x width (BOLDscreen: 546 pixels, EIZOflexscan: 396 pixels) 
%                  x 3 (rgb) 
%                  x 60 frames (total of 1 s, 16.67 ms per frame).
%                  Note that rdks{:,[1,2],[2:5]} are empty as all wm test
%                  stimuli use the highest coherence level.
% * masks       : (cell) 3D cell with alpha transparency masks:
%                  8 motion directions (33.75, 78.75, 123.75, 168.75,
%                  213.75, 258.75, 303.75, 348.75 deg)
%                  x 3 coherence levels (0.30, 0.65, 1.0)
%                  x 5 motion direction offsets (0, -20, -10, +10, +20 deg). 
%                  Each cell contains one uint8 mask image, with
%                  dimensions:
%                  height (BOLDscreen: 546 pixels, EIZOflexscan: 396 pixels)
%                  x width (BOLDscreen: 546 pixels, EIZOflexscan: 396 pixels)
%                  Note that masks{:,[1,2],[2:5]} are empty as all wm test 
%                  stimuli use the highest coherence level.
% * info        : 56x11 table with information matching the rdk array.
%                 Also stored as csv info file in params.stim.rdk.infofile.
%                 Columns include:
%                   1: unique_im  - (double) unique image nr for each RDK 
%                                   movie
%                   2: stim_pos_i - (double) stimulus position index. 1=left, 2=right
%                   3: stim_pos   - (cell) stimulus position, same as 
%                                   stim_loc_i but human-readable
%                   4: dot_motdir_deg - (double) motion direction in degrees 
%                                       (0 = 12 o'clock, 90 = 3 o'clock).
%                   5: dot_motdir_deg_i - (double) same as dot_motdir_deg but 
%                                   indexed 1 (33.75 deg) through 8 (348.75 deg)
%                   6: dot_coh    - (double) coherence level (fraction of 1),
%                   7: dot_coh_i  - (double) same as dot_coh but indexed 1 
%                                   0.30), 2 (0.65) or 3 (1.0).
%                   8: rel_motdir_deg - (double) motion direction relative 
%                                   from corresponding core RDK 0, -20, 
%                                   -10, +10, +20 (deg)
%                   9: rel_motdir_deg_i - (double) same as rel_motdir_deg 
%                                   but indexed 0 (0 deg), 1 (-20 deg), 
%                                   2 (-10 deg), 3 (+10 deg), or 4 (+20 deg).
%                   10: dot_pos   - {1×1 cell} [x,y] dot positions on each
%                                   frame relative to the center of the 
%                                   aperture (pixels). Dimensions for each
%                                   cell (dot by [x,y] by time): 
%                                   BOLDscreen: 478×2×60.
%                                   PProom: 480x2x60. 
%                   11: is_specialcore - (logical) whether the RDK movie is 
%                                   part of the subselected stimuli used in 
%                                   imagery and long-term memory task.
[rdks, ~, ~] = vcd_rdk(params, verbose, store_imgs); % outputs are [rdks, masks, info]

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
% * verbose     : (logical) show debug figures
% * store_imgs  : (logical) store stimuli and debug figures as pngs 
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
% * info        : 80x12 table with stimulus information about single dots.
%                 Also stored as csv info file in params.stim.dot.infofile.
%                 Columns include:
%                 1: unique_im  - (double) unique image nr for each dot 
%                                 location:  core: 49-64, wm test: 303-366.
%                 2: stim_pos   - (cell) stimulus position, same as stim_loc_i 
%                                 but human-readable ({'left'} or {'right'})
%                 3: stim_pos_i - (double) stimulus position index. 1=left, 2=right
%                 4: angle_deg  - (double) dot angle in degrees (0=12 o'clock).
%                 5: angle_i    - (double) same as angle_deg but indexed 1
%                                  (11.25 deg) through 8 (348.75 deg)
%                 6: eccen_deg  - (double) eccentricity of dot (deg) 
%                 7: angle_rad  - (double) dot angle in radians.
%                 8: delta_deg  - (double) dot angle relative from core dot 
%                                  angle: 0, -10, -5, +5, +10 (deg)
%                 9: delta_deg_i - (double) same as delta_deg but indexed 
%                                  0 (0 deg), 1 (-10 deg), 2 (-5 deg), 
%                                  3 (+5 deg), 4 (+10 deg).            
%                 10: dot_xpos_pix - (double) x-position of dot angle in 
%                                  full screen image space (pixels),
%                                  where [0,0] is the upper left pixel
%                 11: dot_ypos_pix - (double) y-position of dot angle in 
%                                  full screen image space (pixels),
%                                  where [0,0] is the upper left pixel
%                 12: is_specialcore - (logical) if dot location is part of
%                                  subselected stimuli used in imagery and 
%                                  long-term memory task crossing.
[single_dot, ~, ~] = vcd_singledot(params, verbose, store_imgs); % outputs are [single_dot, masks, info]

%% Objects
% vcd_objects function creates a total of 80 objects:
% - 16 core images 
% - 64 working memory test images.
% - 288 catch object images.
%
% Note that while BOLDscreen and PProom stimuli have the same number of
% pixels for width and height, the object stimuli themselves are scaled to
% different sizes within this 512x512 pixel support to accomplish 4 deg
% diameter stimuli.
% If params.verbose = true, this function will make a PNG for each object 
% (core and WM test) and figures with axes/titles to check image nrs,
% rotation.
% If params.store_imgs = true, this function will store
% these figures in fullfile(vcd_rootPath,'figs',<dispname>,'obj').
%
% INPUTS:
% * params        : parameter struct (requires display and stimulus parameters)
% * verbose       : (logical) show debug figures
% * store_imgs    : (logical) store stimuli and debug figures as pngs 
%
% OUTPUTS:
% * objects       : (uint8) is a 5D array:
%                   height (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x width (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x 3 (rgb) 
%                   x 16 object categories (subordinate level) 
%                   x 5 rotation offsets (0, -24, -12, +12, +24 deg). 
% * masks         : (uint8) is a 4D array containing alpha transparency mask: 
%                   height (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x width (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x 16 object categories (subordinate level) 
%                   x 5 rotation offsets (0, -24, -12, +12, +24 deg). 
% * objects_catch : (uint8) is a 5D array:
%                   height (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x width (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x 3 (rgb) 
%                   x 16 object categories (subordinate level) 
%                   x 18 catch rotations sampled from all rotation objects
%                   (0-180 degrees), except the core rotation of the
%                   object.
% * masks_catch   : (uint8) is a 4D array containing alpha transparency mask:
%                   height (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x width (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x 16 object categories (subordinate level) 
%                   x 18 catch rotations sampled from all rotation objects
%                   (0-180 degrees), except the core rotation of the
%                   object.
% * info          : (table) 386x12 info table about object stimuli. 
%                   Also stored as csv info file in params.stim.obj.infofile.
%                 Columns include:
%                 1: filename  - (cell) filename of original png,
%                 2: unique_im - (double) unique image nr for each object
%                 3: stim_pos_name - (cell) stimulus position ('left','right'),
%                 4: stim_pos  - (double) stimulus position index. 1=left, 2=right
%                 5: super_cat_name - (cell) superordinate semantic category label:
%                                {'human','animal','object','food','place'}
%                 6: super_cat - (double) same as superordinate, but indexed 
%                                 by nr 1:'human' through 5: 'place'
%                 7: basic_cat_name - (cell) basic semantic category label 
%                                for each superordinate semantic category category:
%                                {'facefemale','facefemale', 'small','large', 'tool',
%                                'vehicle', 'man-made','produce', 'building'};
%                 8: basic_cat - (double) same as basic, but indexed by nr
%                                1: facefemale, 2: facefemale; 1: small, 2: large; 
%                                1: tool, 2: vehicle; 1: man-made, 2: produce; 1: building.
%                 9: sub_cat_name - (cell) subordinate semantic category 
%                                label, name for each individual core object: {'damon',
%                                'lisa','sophia','parrot','cat','bear','giraffe',
%                                'drill','brush','pizza','banana','bus','suv''church',
%                                'house','watertower'};
%                 10: sub_cat - (double) same as subordinate, but indexed 
%                                 by nr 1:damon, 1:lisa, 2:sophia, 1:parrot, 2:cat, 
%                                 1:bear, 2: giraffe, 1:drill, 2: brush, 1:bus, 2: suv,
%                                 1:pizza, 1: banana, 1:church, 2: house, 3: watertower
%                 11: affordance_name - (cell) affordance label for each 
%                                object, names are 1-4: 'greet','grasp','enter','observe'
%                 12: affordance_cat - (double) same as affordance_name, but
%                                 computer-readable: 1: greet, 2: grasp, 3: enter, and 4: observe.
%                 13: is_specialcore - (logical) whether the object is part of the
%                                 subselected stimuli used in imagery and long-term memory task.
%                 14: base_rot - (double) rotation (in deg) of object shown in
%                                 the first stimulus array (prior to applying relative 
%                                 rotation ("rel_rot). "base_rot" is the same for core 
%                                 images, and only relevant for wm test images only. 
%                                 Values range from 0-180 deg (excluding 90 deg).
%                 15: rot_abs - (double) absolute rotation (deg) of object.
%                                 ranges from 0-180, (excluding 90 degrees).
%                                 0 deg = profile view facing to the right.
%                                 90 deg = forward facing view.
%                                 180 deg = profile view facing to the left.
%                 16: rot_rel - (double) relative rotation (deg) of object 
%                                 from core object rotation, ranges from 0, -24, -12, 
%                                 +12,+24 deg. Only relevant for wm test images.
%                 17: facing_dir - (cell) facing direction, 'forward' 
%                                 directions are between 46-134 degree. 'sideways' 
%                                 directions are between 0-44 and 136-180 degrees. 
%                                 Rotations of 45, 90, 135 deg are ill-defined and don't occur.
%                 18: rel_rot_name - (cell) relative rotation name for WM 
%                                 test images whether the relative rotation results in the
%                                 object facing more left or rightwards from the 
%                                 perspective of the observer. 'rightwards' is when 
%                                 absolute rotation is between 0-90 deg and relative 
%                                 rotation is > 0 deg, OR when absolute rotation is 
%                                 between 91-180 deg and relative rotation < 0 deg.
%                                 'leftwards' is when absolute rotation is between 0-90 
%                                 deg and relative rotation is < 0 deg, OR when
%                                 absolute rotation is between 91-180 deg and relative
%                                 rotation is > 0 deg.
%                 19: is_objectcatch - (logical) whether the object is part 
%                                 of the catch object rotations used in PC-OBJ task.
[objects, ~, ~, ~, ~] = vcd_objects(params, verbose, store_imgs); % outputs are [objects, masks, objects_catch, masks_catch, info]

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
% * verbose     : (logical) show debug figures
% * store_imgs  : (logical) store stimuli and debug figures as pngs 
%
% OUTPUTS:
% * scenes      : (uint8) core scene images, 6D array:
%                   height (BOLDscreen: 741 pixels, Eizoflexscan: 537 pixels)
%                   x width (BOLDscreen: 741 pixels, Eizoflexscan: 537 pixels)
%                   x 3 (rgb)
%                   x 5 superordinate semantic object categories (human, animal, object, food, place)
%                   x 2 scene locations (indoor/outdoor)
%                   x 3 obj locations (left/middle/right)
% * ltm_lures   : (uint8) long term memory lure scenes, 7D array:
%                   height (BOLDscreen: 741 pixels, Eizoflexscan: 541 pixels)
%                   x width (BOLDscreen: 741 pixels, Eizoflexscan: 541 pixels)
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
%                   height (BOLDscreen: 741 pixels, Eizoflexscan: 541 pixels)
%                   x width (BOLDscreen: 741 pixels, Eizoflexscan: 541 pixels)
%                   x 3 (rgb)
%                   x 5 superordinate semantic object categories 
%                   x 2 scene locations (indoor/outdoor) 
%                   x 3 obj locations (left/middle/right)
%                   x 4 change types:
%                       1:easy add -- scene is altered by adding something big/obvious 
%                       2:hard add -- scene is altered by adding something small/subtle
%                       3:easy remove -- scene is altered by removing something big/obvious
%                       4:hard remove -- scene is altered by removing something small/subtle
% * info        : (table) 210x14 info table about scenes. Also stored as 
%                       csv info file in params.stim.ns.infofile. Columns
%                       include:
%                  1: filename  - (cell) filename of original png.
%                  2: unique_im - (double) unique image nr for each object.
%                  3: stim_pos_name - (cell) name of stimulus position on 
%                                 the screen. All scenes are {'center'}.
%                  4: stim_pos  - (double) index of "stim_pos". All scenes 
%                                 are 3. (machine readable)
%                  5: super_cat_name - (cell) superordinate semantic 
%                                 category label of the dominant object in the scene:
%                                 {'human','animal','object','food','place'}
%                  6: super_cat - (double) same as super_cat_name, but 
%                                 indexed by nr 1:'human' through 5: 'place'
%                  7: basic_cat_name - (double) basic semantic category label, 
%                                 whether the scene is indoor or outdoor:
%                                 {'indoor', 'outdoor'}
%                  8: basic_cat - (double) same as basic_cat_name, but indexed 
%                                 by nr 1: indoor, 2: outdoor
%                  9: sub_cat_name - (cell) subordinate semantic category 
%                                 label, indicating where the majority of the dominant 
%                                 object in the scene is spatially located  relative
%                                 to central fixation circle: {'left','center','right'};
%                  10: sub_cat - (double) same as subordinate, but indexed 
%                                 by nr 1:left, 2: center, 3: right.
%                  11: affordance_name - (cell) scene affordance inferred 
%                                 by the experimenter {'greet','grasp','walk','observe'}
%                  12: affordance_cat - (double) same as affordance_name, but indexed 
%                                 by nr 1:greet, 2: grasp, 3: walk, 4:observe
%                  13: is_specialcore - (logical) whether the scene is part 
%                                 of the subselected stimuli used in imagery and
%                                 long-term memory task.
%                  14: is_lure - (logical) whether the scenes are lures for 
%                                 long-term memory task: ordered 1-4: 1 is most similar
%                                 to core scene and 4 is least similar.
[scenes, ~, ~, ~] = vcd_naturalscenes(params, verbose, store_imgs); % outputs are [scenes, ltm_lures, wm_im, info]



