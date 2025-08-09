function [data,all_images] = runme_vcdcore(subj_nr,ses_nr,ses_type,run_nr, dispName, varargin)
% Main wrapper function to run the core VCD experiment. This function can 
% run both behavioral and functional MRI versions of the VCD experiment.
%
%     [data,all_images] = runme_vcdcore(subj_nr,ses_nr,ses_type,run_nr,dispName, varargin)
%
% PURPOSE: 
% Execute a single 6.5-min run (behavioral or fMRI) of the core
% Visual Cognition Dataset (VCD) experiment. The behavioral experiment is
% meant to run on the CMRR psychophysics lab Eizoflex Scan monitor
% ('PPROOM_EIZOFLEXSCAN') and the fMRI experiment is meant to run on the
% CMRR 7 Tesla Actively Shielded (7TAS) MRI ('7TAS_BOLDSCREEN32'). This
% function will load in the images according to the time_table_master in
% ./workspaces/info/time_table_master_complete_*.mat, the requested session
% nr, session type, and run number.
% 
% The behavioral experiment has only 1 session (1), one version 
% (use session_type = 1), and 15 runs (use run_nr = 1-15).
% All subjects will see the same runs in the behavioral experiment.
%
% The fMRI experiment has 27 sessions (1-27), each session has 10 runs (use
% run_nr = 1-10). The first fMRI session is a wide session, and has two
% versions. Half of the subjects will see version 1A (use session_type = 1),
% the other half of the subjects will see version 1B (use session_type = 2).
% fMRI session 2-26 are deep sessions, where every subject will see the
% same version (session_type = 1). The last fMRI session (session_nr = 27)
% has two versions split across subjects Half of the subjetcs will see
% version 27A (use session_type = 1), the other half of the subjects will
% see version 27B (use session_type = 2). MRI runs range from 1-10.
% 
% The wrapper will store 2 files in a subject specific data folder for the 
% session that was run:  ./data/<environment>/subj###_ses##
%
% * behavioral .mat file containing triggers, subject's key presses and 
%   VBL monitor refresh timing: 
%   'behavior_<datetime>_vcd_subj<subject_nr>_ses<session_nr>_<session_type>_run<run_nr>.mat'
%   
% * eyetracking .edf file containing raw time series of gaze xy-position, 
%   pupil size, eye velocity, and sync messages:
%   'eye_<datetime>_vcd_subj<subject_nr>_ses<session_nr>_<session_type>_run<run_nr>.edf'
%
% <environment> is either:
%  * 'MRI'      (when using dispName = '7TAS_BOLDSCREEN32'), 
%  * 'BEHAVIOR' (when using dispName = 'PPROOM_EIZOFLEXSCAN'), 
% <subject_nr> is a zero-padded integral number ranging from 1 to 999
% <session_nr> is a zero-padded integral number ranging from 1 to 27
% <session_type> is either  "A" for session_type=1 or "B" for session_type=2
%   Note: most sessions only have session_type=1. The only sessions with 
%   two session types are MRI session 1 (the one and only wide session) 
%   and MRI session 27 (the last deep session). 
% <run_nr> is a zero-padded integral number ranging from 1 to 15.
%
% Setup display resolution:
% 7TAS_BOLDSCREEN32: (BOLDscreen + Nova1x32 headcoil + large mirror)
%  * viewing distance 183.5 cm; 
%  * monitor size (h x w): 39.29 x 69.84 cm; 
%  * monitor resolution (h x w): 1080 x 1920 pixels; 
%  * monitor refresh rate: 120 Hz. 
% PPROOM_EIZOFLEXSCAN: (EizoFlexScan + Bits#)
%  * viewing distance 99.0 cm; 
%  * monitor size (h x w): 32.5 x 52.0 cm; 
%  * monitor resolution (h x w): 1200 x 1920 pixels; 
%  * monitor refresh rate: 60 Hz. 
% 
% From resolution to degrees visual angle (DVA):  
% BOLDscreen height   : atan( (39.29 / 2) / 183.5) / pi*180*2 = 12.22 deg 
% BOLDscreen width    : atan( (69.84 / 2) / 183.5) / pi*180*2 = 21.55 deg 
% BOLDscreen pixels per deg for 4 deg eccentricity:  88.178261586694418 [precise]
% EIZOFlexscan height : atan( (32.5 / 2) / 99.0) / pi*180*2 = 18.64 deg 
% EIZOFlexscan width  : atan( (52.0 / 2) / 99.0) / pi*180*2 = 29.43 deg 
% EIZOFlexscan pixels per deg for 4 deg eccentricity: 63.902348145300280 [precise]
%
% For 7TAS BOLDscreen, a conservative estimate of what subject can see is 
% height: 1080 px, width: 1152 pix.
% For the width/horizontal dimension, this results in a total (jointly across two eyes):
%   atan(1152/1920 * (69.84 / 2) /183.5) / pi*180*2 = 13.03 deg
% For the height/vertical dimension, this results in a total (jointly across two eyes):
%   atan(1080/1080 * (39.29 / 2) /183.5) / pi*180*2 = 12.2 deg
% For reference: nsdheightdeg = 12.70 in degrees (BOLDscreen vertical, Nova1x32)
% 
% INPUTS
%  Input parser requires first 5 inputs (subj_nr, ses_nr, ses_type, run_nr, dispName).
%  The other input arguments are optional and will be set to default if 
%  undefined. To define optional input arguments, use: 'vararg', <val>.
%
% MANDATORY INPUTS
%   subj_nr              : subject number 
%   ses_nr               : session number  (numeric integer between 1-27)
%   ses_type             : session type number (numeric integer: 1=A and 2=B)
%   run_nr               : run number      (numeric integer between 1-10)
%   dispName             : name of display (string), this affects what stimuli and time_table_master are loaded
%                          Choose: '7TAS_BOLDSCREEN32'   - BOLD screen at the 7TAS MRI scanner
%                                  'KKOFFICE_AOSQ3277'   - external monitor in kendrick's CMRR office
%                                  'EKHOME_ASUSVE247'    - external monitor at Eline's home
%                                  'PPROOM_EIZOFLEXSCAN' - EizoFlexscan monitor @ CMRR's psychophysics room
%                                  'CCNYU_VIEWPIXX3D'    - ViewPixx monitor in Clay Curtis' psychophysics room     
% OPTIONAL INPUTS
%   [wantsynctest]       : if true, we do the PTB monitor sync test. 
%                          Default: true
%   [loadparams]         : if true, load stored parameter values. 
%                          Default: true
%   [storeparams]        : if true, store created "stim" and "exp" parameter values. 
%                          Default: false
%   [savestim]           : if true, store stimuli in matlab file prior to ptb flipping. File is ~25-50 MB!
%                          Default: false
%   [loadstimfromrunfile]: if true, load subject single run file with stimuli (to save time a minute and rerun the same run) 
%                          Default: false
%   [verbose]            : if true, print all the text in the command window, if false, we do the minimum stuff
%                          Default: false
%   [storeimgs]          : if true, we save debug figures as .png (set to false when running the actual experiment!!)
%                          Default: false
%   [offsetpix]          : offset of center [x,y]-coordinate in pixels. 
%                          Default: No offset [0,0]
%   [movieflip]          : flip presented stimulus display up-down (first input, first dimension of pixel image) or left-right (second input, second dimension of pixel image). 
%                          Default: no flipping: [0,0].
%   [savedatafolder]     : folder where subject's button presses and stim timing are stored. 
%                          Default: [], which will turn into: sprintf('vcd_subj%03d_ses%02d',subj_nr,sesID);
%   [subjfilename]       : pre-fix for behavioral and eyelink mat-file. 
%                          Default: [], which will turn into: sprintf('%s_vcd_subj%03d_ses%02d_%s_run%02d.mat',ts0,subj_nr,ses_nr,ses_type,run_nr);
%   [wanteyetracking]    : if true, we will initialize eyetracking. 
%                          Default: false;
%   [wantdatabypass]     : if true, we will skip actually running the experiment and just generate dummy behavioral .mat files. 
%                          Default: false;
%   [ptbMaxVBLstd]       : allowable error (std in seconds) in the 50 VBL time stamp measurements taken by PTB during the synctest. 
%                          If the measurement error is larger than the desired monitor refresh rate, psychtoolbox throws an error during the synctest.
%                          Default: 0.0009 s. Larger standard deviations are more forgiving and less likely for psychtoolbox to throw an error.
%   [all_images]         : struct with all VCD images, to avoid extra load time.
%   [timetable_file]     : mat-file specifying the time_table_master based on the subject's shuffled condition_master table.
%                          Default: []. No file means that the runme_vcdcore code will ask user if they want to select an existing file or generate a new time_table_master file on the fly.
%   [stimfolder]         : folder where pre-made stimuli are stored.
%                          Default: fullfile(vcd_rootPath,'workspaces','stimuli')
%   [instrfolder]        : folder where the instruction ong and txt files are stored.
%                          Default: fullfile(vcd_rootPath,'workspaces','instructions')
%   [infofolder]         : folder where pre-made stimuli are stored.
%                          Default: fullfile(vcd_rootPath,'workspaces','info')
%   [exp_env]            : argument to vcd_startup.m. Default: []. Choose between 1:'7tas', 2:'cmrr pproom', 2:'nyu pproom', 4:'other'
%   [is_demo]            : if true, we will treat this as a demo run. if false, we treat this as a normal run.
%                          Default: false;
%
% OUTPUTS:
%   data                 : struct with behavioral button presses and monitor
%                          refresh rate timing, as well as other parameters.
%   all_images           : struct with all VCD stimuli and instruction images 
%                          in uint8 pixels, to avoid additional loading time 
%                          when executing multiple runs.
%
% EXAMPLES:
%  [data,all_images] = runme_vcdcore(1, 1, 1, 1, '7TAS_BOLDSCREEN32'); 
%  [data,all_images] = runme_vcdcore(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN'); % default is with sync test and no eyetracking
%  [data,all_images] = runme_vcdcore(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wantsynctest', true, 'ptbMaxVBLstd', 0.0006); 
%  [data,all_images] = runme_vcdcore(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wantsynctest', true, 'wanteyetracking', true); % sync test + eyetracking
%  [data,all_images] = runme_vcdcore(1, 1, 1, 1, '7TAS_BOLDSCREEN32'  , 'wantsynctest', false);
%  [data,all_images] = runme_vcdcore(1, 1, 1, 1, 'KKOFFICE_AOCQ3277'  , 'wantsynctest', true);
%  [data,all_images] = runme_vcdcore(1, 1, 1, 1, 'EKHOME_ASUSVE247'   , 'wantsynctest', false);
%  [data,all_images] = runme_vcdcore(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'timetable_file',[])
%  [data,all_images] = runme_vcdcore(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN','timetable_file', fullfile(vcd_rootPath,'data','BEHAVIOR','vcd_subj000','vcd_subj000_time_table_master_PPROOM_EIZOFLEXSCAN_20250610T110319.mat'))
%  [data,all_images] = runme_vcdcore(1, 1, 1, 2, 'PPROOM_EIZOFLEXSCAN','timetable_file', fullfile(vcd_rootPath,'data','BEHAVIOR','vcd_subj001','vcd_subj001_time_table_master_PPROOM_EIZOFLEXSCAN_20250611T170541.mat'), 'all_images', all_images)
% 
%
% DEPENDENCIES:
%  * Psychtoolbox-3 (https://github.com/Psychtoolbox/Psychtoolbox-3) 
%       7TAS   : v. 3.0.14 commit ef093cbf296115badddb995fa06452e34c8c7d02 (origin/master). build date Nov 17 19:46:08 2020 +0100
%       PP room: v. 3.0.14 build date 0ct 3, 2017
%
%  * knkutils (https://github.com/cvnlab/knkutils)
%     commit 27dd66770edcf0ef3adbd73e1892678a275e2383 (origin/master)
% 
% AUTHOR:
%  Written by Eline Kupers @ UMN (kupers@umn.edu) (2024,2025)
%
% HISTORY:
%  - 2018/12/25 - initial runnsd version 
%  - 2024/11/06 - initial vcd version with some stuff adopted from runnsd
%  - 2024/12/01 - first version committed to github 

close all; clc;

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p = inputParser;
p.addRequired ('subj_nr'            , @isnumeric); % subject number 
p.addRequired ('ses_nr'             , @isnumeric); % session number 
p.addRequired ('ses_type'           , @isnumeric); % session type (1=A or 2=B) 
p.addRequired ('run_nr'             , @isnumeric); % run number
p.addRequired('dispName'            , @(x) any(strcmp(x, {'7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','PPROOM_EIZOFLEXSCAN','EKHOME_ASUSVE247', 'CCNYU_VIEWPIXX3D'}))); % display name
p.addParameter('wantsynctest'       , true    , @islogical);  % whether we want to run the PTB sync test or not
p.addParameter('loadparams'         , false   , @islogical);  % whether we want to load stim or exp params from the mat file in the info folder
p.addParameter('storeparams'        , false   , @islogical);  % whether we want to load stim or exp params as a mat file in the info folder
p.addParameter('savestim'           , false   , @islogical);  % whether we want to store matlab file with stimuli and timing
p.addParameter('loadstimfromrunfile', false   , @islogical);  % whether we want to load stim from run file
p.addParameter('verbose'            , false   , @islogical);  % whether to print out text in command window.   
p.addParameter('storeimgs'          , false   , @islogical);  % whether we want to store debug figures
p.addParameter('offsetpix'          , [0 0]   , @isnumeric);  % [x,y]-offset of screen in pixels relative to center. [0 0] means no offset, centered stimuli. [10 20] means move 10-px right, 20-px down
p.addParameter('movieflip'          , [0 0]   , @isnumeric);  % whether to flip up-down, whether to flip left-right: % up/down, left/right
p.addParameter('savedatafolder'     , ''      , @ischar);     % place to store data with today's date
p.addParameter('subjfilename'       , ''      , @ischar);     % subject specific preamble given to behavioral and eyelink file names
p.addParameter('wanteyetracking'    , false   , @islogical);  % whether to try to hook up to the eyetracker
p.addParameter('wantdatabypass'     , false   , @islogical);  % whether to skip the experiment and just save dummy .mat file
p.addParameter('ptbMaxVBLstd'       , 0.0009  , @isnumeric);  % what standard deviation for screen flip duration do we allow?
p.addParameter('all_images'         , struct(), @isstruct);   % preloaded all_images in a single struct (to save time)
p.addParameter('timetable_file'     , ''      , @ischar);     % what subject time_table_master file are we loading? file should exist in savedatafolder  
p.addParameter('stimfolder'         , fullfile(vcd_rootPath,'workspaces','stimuli')     , @ischar); % Where do the stimulus mat files live?
p.addParameter('instrfolder'        , fullfile(vcd_rootPath,'workspaces','instructions'), @ischar); % Where do the task instruction txt and png files live?
p.addParameter('infofolder'         , fullfile(vcd_rootPath,'workspaces','info')        , @ischar); % where do the *_info.csv file(s) live?
p.addParameter('exp_env'            , []        , @isnumeric); % input argument given to vcd_startup to add paths to toolbox dependencies.  Default: [] will ask experimenter for environment. Other options are 1:'7tas', 2:'cmrr pproom', 2:'nyu pproom', 4:'other'
p.addParameter('deviceNr'           , -3        , @isnumeric); % kbWait/kbCheck input device number to listen to. Default = -3 will listen to all devices. Previously: vcd_checkDevices(params.deviceNr, params.device_check);
p.addParameter('device_check'       , 'both'    , @char);      % what type of devices do we want to check for button presses: 'external','internal', or 'both'
p.addParameter('triggerkey'         , {'5','t'}, @(x) iscell(x) || isstring(x)) % key(s) that starts the experiment
p.addParameter('triggerkeyname'     , '''5'' or ''t''', @isstring) % for display only
p.addParameter('userkeys'           , {'1','2','3','4'}, @(x) iscell(x) || isstring(x)) % key(s) that participants are expected to push.
p.addParameter('store_imgs'         , false     , @islogical)  % whether to save figures locally (true) or not (false).                      
p.addParameter('is_demo'            , false     , @islogical)  % whether this is a demo run (true) or not (false).
% Parse inputs
p.parse(subj_nr, ses_nr, ses_type, run_nr, dispName, varargin{:});

% Rename variables into struct
params = p.Results;

%% Do startup

vcd_startup(params.exp_env); % function will navigate to root of vcd-stim code folder

%% Check environment and input parameters:
       
% Are we running the behavioral or MRI experiment? Or are we testing on an office monitor?
if strcmp(params.dispName,'7TAS_BOLDSCREEN32')
    params.env_type = 'MRI';
elseif ismember(params.dispName,{'PPROOM_EIZOFLEXSCAN','CCNYU_VIEWPIXX3D','KKOFFICE_AOCQ3277','EKHOME_ASUSVE247'})
    params.env_type = 'BEHAVIOR';
end

      
% Can we actually run this experiment?
assert(params.subj_nr>=0 && (params.subj_nr<=999));
if strcmp(params.env_type,'BEHAVIOR')
    params.is_wide  = false;
    if params.is_demo
        if params.run_nr~=1, error('[%s]: Run number can only be 1 for demo runs',mfilename); end
        if ~ismember(params.ses_nr,[1:3]), error('[%s]: Session number can only be 1, 3, or 3 for demo runs',mfilename); end
    else
        assert(params.run_nr>=1 && params.run_nr<=15);
        assert(isequal(params.ses_nr,1));
    end
    assert(isequal(params.ses_type,1));
    
elseif strcmp(params.env_type,'MRI')
    assert(params.ses_nr>=1 && params.ses_nr<=31);
    if params.ses_nr == 1, params.is_wide = true; end
    if ismember(params.ses_nr, [1,31]), assert(ismember(params.ses_type, [1,2]));
    else, assert(isequal(params.ses_type,1)); end
    assert(params.run_nr>=1 && params.run_nr<=10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%% Deal with folders and filenames %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(params.savedatafolder)
    params.savedatafolder = fullfile(vcd_rootPath,'data',params.env_type,sprintf('vcd_subj%03d_ses%02d',params.subj_nr, params.ses_nr));
end
if ~exist(params.savedatafolder, 'dir'); mkdir(params.savedatafolder); end


if isempty(params.timetable_file)
   load_existing_timetable = input('No timetable_file was specified, do you want to regenerate the file?  1: YES   2: NO \n');

   if load_existing_timetable == 1
       fprintf('OK, will regenerate new time table for this subject''s session and run.\n')
   elseif load_existing_timetable == 2
       tmp_timetable_dir = fullfile(vcd_rootPath,'data',params.env_type, sprintf('vcd_subj%03d',params.subj_nr));
       if strcmp(params.dispName,'CCNYU_VIEWPIXX3D')
           tempfiles = matchfiles(fullfile(tmp_timetable_dir,sprintf('vcd_subj%03d_time_table_master_%s%s*.mat',params.subj_nr, choose(params.is_demo,'demo_',''),'PPROOM_EIZOFLEXSCAN')));
       else
           tempfiles = matchfiles(fullfile(tmp_timetable_dir,sprintf('vcd_subj%03d_time_table_master_%s%s%s*.mat',params.subj_nr, choose(params.is_wide,'wide_',''), choose(params.is_demo,'demo_',''),params.dispName)));
       end
       % If user doesn't want to generate a time table and we can't find any, end the run gracefully
       if isempty(tempfiles)
           fprintf('\nIt appears you have NO time tables for this specific subject and session. Will abort the run.\n');
           data = []; all_images = struct;
           return
       else
           tempfiles_short = strrep(tempfiles, tmp_timetable_dir, '.');
           fprintf('\nIt appears you have the following time tables already:\n');
           
           for mm = 1:length(tempfiles_short)
               fprintf('  %d: %s\n', mm, tempfiles_short{mm})
           end
           
           timetable_idx = input('Type nr of time table file you want to use: \n');
           if ismember(timetable_idx, 1:length(tempfiles)) % check if we can load the user's requested file               
               params.timetable_file = tempfiles{timetable_idx};
               fprintf('\nWill use timetable_file = %s\n',params.timetable_file);
           else
               % If user selects a file number that doesn't exists, we also end the run gracefully
               fprintf('\nRequested timetable_file doesn''t exist. Please rerun and pick an existing time table file.\n');
               data = []; all_images = struct;
               return
           end
           clear timetable_idx tempfiles_short tmp_timetable_dir load_existing_timetable
       end
   end

end

% Where do we store behavioral and eyetracking data?
ts0 = gettimestring;
if isempty(params.subjfilename)
    params.behavioralfilename  = sprintf('behavior_%s_vcd_subj%03d_ses%02d_%s_run%02d.mat',...
        ts0,params.subj_nr,params.ses_nr,choose(params.ses_type==1,'A','B'),params.run_nr);
    if params.wanteyetracking
        params.eyelinkfilename = sprintf('eye_%s_vcd_subj%03d_ses%02d_%s_run%02d.edf',...
            ts0,params.subj_nr,params.ses_nr,choose(params.ses_type==1,'A','B'),params.run_nr);
    else
        params.eyelinkfilename = '';
    end
else
    params.behavioralfilename  = fullfile(params.savedatafolder,sprintf('%s_%s.mat',ts0,params.subjfilename));
    if params.wanteyetracking
        params.eyelinkfilename = fullfile(params.savedatafolder,sprintf('%s_%s.edf',ts0,params.subjfilename));
    else
        params.eyelinkfilename = '';
    end
end
    

%% %%%%%%%%%%%%%%%%%%%%%%%% DON'T EDIT BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tell operator what experiment are running
fprintf('[%s]: Running VCD core %s experiment: subj_nr %03d - session %02d %s - run %02d \n', ...
    mfilename, params.env_type, params.subj_nr,params.ses_nr,choose(params.ses_type==1,'A','B'),params.run_nr)
fprintf('[%s]: %s eyetracking \n', mfilename, choose(params.wanteyetracking,'YES','NO'))
fprintf('[%s]: Running experiment with images optimized for %s\n', mfilename, params.dispName)
if isempty(params.timetable_file)
    fprintf('[%s]: Subject''s time table file was NOT specified. Will create one on the fly. \n', mfilename)
else
    fprintf('[%s]: Will load subject''s time table file:\n', mfilename)
    fprintf('[%s]: \t %s\n',mfilename, params.timetable_file)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Run experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for optional inputs use: params.('var') = <val>
[data, all_images] = vcd_singleRun(params);

return
    
