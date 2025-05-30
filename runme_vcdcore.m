function data = runme_vcdcore(subj_nr,ses_nr,ses_type,run_nr, dispName, varargin)
% Main wrapper function to run the core VCD experiment. This function can 
% run both behavioral and functional MRI versions of the VCD experiment.
%
%     data = runme_vcdcore(subj_nr,ses_nr,ses_type,run_nr,dispName, varargin)
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
% versions. Half of the subjetcs will see version 1A (use session_type = 1),
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
%  * 'TEST'     (when using dispName = 'KKOFFICE_AOCQ3277' or 'EKHOME_ASUSVE247')
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
% BOLDscreen pixels per deg: 88.3702666667131 [precise]
% EIZOFlexscan height : atan( (32.5 / 2) / 99.0) / pi*180*2 = 18.64 deg 
% EIZOFlexscan width  : atan( (52.0 / 2) / 99.0) / pi*180*2 = 29.43 deg 
% EIZOFlexscan pixels per deg: 64.3673991642835 [precise]
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
%  Input parser requires first three inputs (subj_nr, sesID, run_nr).
%  The other input arguments are optional and will be set to default if 
%  undefined. To define optional input arguments, use: 'vararg', <val>.
%
%   subj_nr           : subject number 
%   ses_nr            : session number  (numeric integer between 1-27)
%   ses_type          : session type number (numeric integer: 1=A and 2=B)
%   run_nr            : run number      (numeric integer between 1-10)
%   dispName          : name of display (string), this affects what stimuli and time_table_master are loaded
%                       Choose: '7TAS_BOLDSCREEN32'   - BOLD screen at the 7TAS MRI scanner
%                               'KKOFFICE_AOSQ3277'   - external monitor in kendrick's CMRR office
%                               'EKHOME_ASUSVE247'    - external monitor at Eline's home
%                               'PPROOM_EIZOFLEXSCAN' - CMRR's Psychophysics room monitor                      
%   [wantsynctest]    : if true, we do the PTB monitor sync test. 
%                        Default: false
%   [loadparams]      : if true, load stored parameter values. 
%                        Default: true
%   [storeparams]     : if true, store created parameter values. 
%                        Default: true
%   [savestim]        : if true, store stimuli in matlab file prior to ptb flipping. File is ~25-50 MB!
%                        Default: false
%   [loadstimfromrunfile]  : if true, load subject single run file with stimuli (to save time a minute and rerun the same run) 
%                        Default: false
%   [offsetpix]       : offset of center [x,y]-coordinate in pixels. 
%                        Default: No offset [0,0]
%   [movieflip]       : flip presented stimulus display left-right (first argument) or up-down (second argument). 
%                        Default: no flipping: [0,0].
%   [stimdir]         : folder where pre-made stimuli are stored.
%                        Default: fullfile(vcd_rootPath,'workspaces','info')
%   [savedatadir]     : folder where subject's button presses and stim timing are stored. 
%                        Default: [], which will turn into: sprintf('vcd_subj%03d_ses%02d',subj_nr,sesID);
%   [subjfilename]    : pre-fix for behavioral mat-file. 
%                        Default: [], which will turn into: sprintf('%s_vcd_subj%03d_ses%02d_run%02d.mat',ts0,subj_nr,sesID,run_nr);
%   [wanteyetracking] : if true, we will initialize eyetracking. 
%                        Default: false;
%   [ptbMaxVBLstd]    : allowable error (std in seconds) in the 50 VBL time stamp measurements taken by PTB during the synctest. 
%                       If the measurement error is larger than the desired monitor refresh rate, psychtoolbox throws an error during the synctest.
%                        Default: 0.0004 s. Larger standard deviations are more forgiving and less likely for psychtoolbox to throw an error.
%
% OUTPUTS:
%   data              : struct with behavioral button presses and monitor
%                        refresh rate timing, as well as other parameters.
%
% EXAMPLES:
%  runme_vcdcore(1, 1, 1, 1, '7TAS_BOLDSCREEN32'); 
%  runme_vcdcore(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN'); % default is with sync test and no eyetracking
%  runme_vcdcore(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wantsynctest', true, 'ptbMaxVBLstd', 0.0006); 
%  runme_vcdcore(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wantsynctest', true, 'wanteyetracking', true); % sync test + eyetracking
%  runme_vcdcore(1, 1, 1, 1, '7TAS_BOLDSCREEN32'  , 'wantsynctest', false);
%  runme_vcdcore(1, 1, 1, 1, 'KKOFFICE_AOCQ3277'  , 'wantsynctest', true);
%  runme_vcdcore(1, 1, 1, 1, 'EKHOME_ASUSVE247'   , 'wantsynctest', false);
%
% DEPENDENCIES:
%  * Psychtoolbox-3 (https://github.com/Psychtoolbox/Psychtoolbox-3) 
%  7TAS: v. 3.0.14
%      commit ef093cbf296115badddb995fa06452e34c8c7d02 (origin/master)
%      Date:   Tue Nov 17 19:46:08 2020 +0100
%  PP room: v. 3.0.14 — build date 0ct 3, 2017
%
%  * knkutils (https://github.com/cvnlab/knkutils)
%     commit 27dd66770edcf0ef3adbd73e1892678a275e2383 (origin/master)
% 
% AUTHOR:
%  Written by Eline Kupers @ UMN (kupers@umn.edu)
%
% HISTORY:
%  - 2018/12/25 - initial runnsd version 
%  - 2024/11/06 - initial vcd version with some stuff adopted from runnsd
%  - 2024/12/01 - first version committed to github 

close all; clc;
vcd_startup;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DON'T EDIT BELOW

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p = inputParser;
p.addRequired ('subj_nr'        , @isnumeric); % subject number 
p.addRequired ('ses_nr'         , @isnumeric); % session number 
p.addRequired ('ses_type'       , @isnumeric); % session type (1=A or 2=B) 
p.addRequired ('run_nr'         , @isnumeric); % nun number
p.addRequired('dispName'        , @(x) any(strcmp(x, {'7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','PPROOM_EIZOFLEXSCAN','EKHOME_ASUSVE247'})))
p.addParameter('wantsynctest'       , true, @islogical);
p.addParameter('loadparams'         , true, @islogical);
p.addParameter('storeparams'        , true, @islogical);
p.addParameter('savestim'           , false, @islogical);
p.addParameter('loadstimfromrunfile', false, @islogical);
p.addParameter('offsetpix'          , [0 0], @isnumeric); % [x,y]
p.addParameter('movieflip'          , [0 0], @isnumeric); % up/down, left/right
p.addParameter('stimDir'            , fullfile(vcd_rootPath,'workspaces','info'), @ischar);
p.addParameter('savedatadir'        , [], @ischar);
p.addParameter('subjfilename'       , [], @ischar);
p.addParameter('wanteyetracking'    , false, @islogical);
p.addParameter('ptbMaxVBLstd'       , [], @isnumeric);

% Parse inputs
p.parse(subj_nr, ses_nr, ses_type, run_nr, dispName, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p.Results.%s;', rename_me{ff},rename_me{ff})]); %#ok<NBRAK>
end
clear rename_me ff p

cd(vcd_rootPath);
addpath(genpath(pwd));

% Subject anon function
subjfun = @(x) find(ismember([1:999],x)); %#ok<NBRAK> % EK Question: how do we map subject ID to subject numbers 1-999?

%%%%%%%%%%%%%%%% need to get info from user
% ask the user what to run
if ~exist('subj_nr','var') || isempty(subj_nr)
  subj_nr = input('What is the subj number? ');
end

if ~exist('ses_nr','var') || isempty(ses_nr)
    % ask the user more stuff
    ses_nr  = input('What session nr (1-27)? ');
end

if ismember(ses_nr,[1,27])
    if ~exist('ses_type','var') || isempty(ses_type)
        % ask the user more stuff
        ses_type  = input('What session type (1 = A and 2 = B)? ');
    end
end

if ~exist('run_nr','var') || isempty(run_nr)
    run_nr = input('What run number (1-15)? ');
end

if ~exist('dispName','var')
    dispName = input('What monitor are you using? (press 1 for 7TAS, 2 for Psychophysics room, 3 for Office) ');
end

% Check environment: are we running the behavioral or MRI experiment? 
% Or are we testing on a different/office monitor?
if isnumeric(dispName)
    if isequal(dispName,1)
        env_type = 'MRI';
    elseif isequal(dispName,2)
        env_type = 'BEHAVIOR';
    else
        env_type = 'TEST';
    end
elseif ischar(dispName)
    if strcmp(dispName,'7TAS_BOLDSCREEN32')
        env_type = 'MRI';
    elseif strcmp(dispName,'PPROOM_EIZOFLEXSCAN')
        env_type = 'BEHAVIOR';
    else
        env_type = 'TEST';
    end
end
      

% Give some hints to the operator
if isempty(savedatadir) %#ok<NODEF>
    savedatadir = fullfile(vcd_rootPath,'data',env_type,sprintf('vcd_subj%03d_ses%02d',subj_nr, ses_nr));
end
if ~exist(savedatadir, 'dir'); mkdir(savedatadir); end

tempfiles = matchfiles(fullfile(savedatadir,sprintf('*_vcd_subj%03d_ses%02d_%s_run%02d.mat',subj_nr, ses_nr, choose(ses_type==1,'A','B'), run_nr)));
fprintf('\nIt appears you have completed the following runs already:\n');
tempfiles %#ok<NOPRT>


if ~exist('wanteyetracking','var')
    wanteyetracking  = input('Do you want to use eyetracking? (press 1 for yes, 0 for no) ');
end

% Tell operator what experiment are running
fprintf('[%s]: Running VCD core %s experiment: subj_nr %03d - session %02d %s - run %02d \n', mfilename, env_type, subj_nr,ses_nr,choose(ses_type==1,'A','B'),run_nr)
fprintf('[%s]: %s eyetracking \n', mfilename, choose(wanteyetracking,'YES','NO'))
fprintf('[%s]: Running experiment with images optimized for %s\n', mfilename, dispName)


% Deal with folders and filenames
instructionsDir = fullfile(vcd_rootPath, 'workspaces','instructions');

ts0 = gettimestring;
if isempty(subjfilename)
    behavioralfilename  = sprintf('behavior_%s_vcd_subj%03d_ses%02d_%s_run%02d.mat',ts0,subj_nr,ses_nr,choose(ses_type==1,'A','B'),run_nr);
    if wanteyetracking
        eyelinkfilename = sprintf('eye_%s_vcd_subj%03d_ses%02d_%s_run%02d.edf',ts0,subj_nr,ses_nr,choose(ses_type==1,'A','B'),run_nr);
    else
        eyelinkfilename = '';
    end
else
    behavioralfilename  = fullfile(savedatadir,sprintf('%s_%s.mat',ts0,subjfilename));
    if wanteyetracking
        eyelinkfilename = fullfile(savedatadir,sprintf('%s_%s.edf',ts0,subjfilename));
    else
        eyelinkfilename = '';
    end
end
    

% Do some checks, can we actually run this experiment?
assert(~isempty(subjfun(subj_nr)));
assert(subjfun(subj_nr)>=1 && subjfun(subj_nr)<=999);
if strcmp(env_type,'BEHAVIOR')
    assert(run_nr>=1 && run_nr<=15);
    assert(isequal(ses_type,1));
    assert(isequal(ses_nr,1));
elseif strcmp(env_type,'MRI')
    assert(ses_nr>=1 && ses_nr<=27);
    if ismember(ses_nr, [1,27]), assert(ismember(ses_type, [1,2]));
    else, assert(isequal(ses_type,1)); end
    if ses_nr == 27, assert(run_nr>=1 && run_nr<=4);
    else, assert(run_nr>=1 && run_nr<=10); end
end


%% run experiment
% for optional inputs use: 'var',<val>
data = vcd_singleRun(subj_nr, ses_nr, ses_type, run_nr, dispName, ... % mandatory inputs
    'env_type',        env_type, ... % optional inputs
    'wantsynctest',    wantsynctest, ...
    'behaviorfile',    behavioralfilename, ...
    'eyelinkfile',     eyelinkfilename, ...
    'savedatadir',     savedatadir, ...
    'wanteyetracking', wanteyetracking, ...
    'loadparams',      loadparams, ...
    'storeparams',     storeparams, ...
    'offsetpix',       offsetpix, ...
    'movieflip',       movieflip, ...
    'instrtextdir',    instructionsDir, ...
    'savestim',        savestim, ...
    'loadstimfromrunfile', loadstimfromrunfile, ...
    'ptbMaxVBLstd',    ptbMaxVBLstd); 

