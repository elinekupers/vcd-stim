%% s_createDesignMatrix.m
%
% Script to create and store a giant matrix (or struct??) that defines the 
% the order in which VCD-core stimuli and tasks are shown across subject 
% at 4 levels of the experiment:
%  1: subject session (1.5 hrs), contains multiple runs (highest level) 
%  2: runs (6 mins), contains multiple blocks, interspersed with rest
%     periods 
%  3: block (1-1.5 mins), contains multiple trials of a single 
%     stimulus-task crossing. blocks are interspersed with variable
%     (5-9 sec) Inter-Block-Intervals (IBIs), which are "rest" periods where
%     subject is presented with background+fixation dot (and continues to 
%     fixate). Number of trials per block depends on the trial type (8 
%     trials for single epoch trials, 4 trials for double epoch trials). 
%     Trials within a block are interspersed with variable (0.2-1.6 sec)
%     Inter-Trial-Intervals (ITIs).
%  4: trial (single: 4.2 sec or double: 14.2 sec, without ITI). Structure
%     sequence of the following events: 
%     SINGLE EPOCH:
%       * trial start (fixation dot rim thickening) (0.4 sec)
%       * covert spatial attention cue (0.8 sec)
%       * stimulus array (epoch 1) (2.0 sec)
%       * response window (1.0 sec)
%       * trial end + ITI (fixation dot rim thinning (0.2-1.6 sec)
%     DOUBLE EPOCH:
%       * trial start (fixation dot rim thickening) (0.4 sec)
%       * covert spatial attention cue (0.8 sec)
%       * stimulus array (epoch 1) (2.0 sec)
%       * inter-stimulus interval (ISI) (8.0 sec)
%       * stimulus array (epoch 2) (2.0 sec)
%       * response window (1.0 sec)
%       * trial end + ITI (fixation dot rim thinning (0.2-1.6 sec)
%
% Requirements:
%   * vcd-stim code repo (github.com/elinekupers/vcd-stim)
%
%
% Written by Eline Kupers @ UMN
% History:
% v0.0: 11/2024 - first attempt
% v0.1: 03/2025 - refactoring code

%% %%%%%%%%%%%%%%%%%%%
%%%%%% PARAMETERS %%%% 
%%%%%%%%%%%%%%%%%%%%%%
params = struct();
params.verbose        = false; % visualize stimuli or not
params.store_imgs     = true; % store visualization figures
params.saveFigsFolder = fullfile(vcd_rootPath,'figs'); % where to store visualization figures

% Get display params
dispname    = 'PPROOM_EIZOFLEXSCAN'; % Choose from: '7TAS_BOLDSCREEN32' or 'KKOFFICE_AOCQ3277' or 'EKHOME_ASUSVE247' or 'PPROOM_EIZOFLEXSCAN'
params.disp = vcd_getDisplayParams(dispname);

% Infer session type
if strcmp(dispname,'7TAS_BOLDSCREEN32')
    session_type = 'MRI';
elseif strcmp(dispname,'PPROOM_EIZOFLEXSCAN')
    session_type = 'BEHAVIOR';
else
    session_type = 'MRI';
end

% Get stimulus parameters
params.load_params                 = true; % load stored params
params.store_params                = true;

% SETUP RNG
params.rng.rand_seed = sum(100*clock);
rand('seed', params.rng.rand_seed);
params.rng.randn_seed = sum(100*clock);
randn('seed', params.rng.randn_seed);

%% Define/Load stimulus params 
% Input 1: Display name to load disp params struct (see vcd_getDisplayParams.m) 
% Input 2: Load prior stored parameters or not? (logical)
% Input 3: Store generated parameters or not? (logical)
params.stim   = vcd_getStimParams('disp_name', params.disp.name, ...
                            'load_params', params.load_params, ...
                            'store_params', params.store_params); 

%% Define/Load experiment session params
params.exp    = vcd_getSessionParams('disp_name', params.disp.name, ...
                                'presentationrate_hz',params.stim.presentationrate_hz, ...
                                'load_params', params.load_params, ...
                                'store_params', params.store_params);

%% Make/Load blocks with trials that sample unique stimuli from each class
% !!WARNING!! There is a randomization component involved in creating the
% trial sequence (e.g., order of unique images within a block). If you
% don't want this, set second input (load_params) to true and load an
% existing file
%
% This function contains the following important steps/functions:
% * vcd_defineUniqueImages
% * vcd_createConditionMaster
% * vcd_shuffleStimForTaskClass
params.verbose = false;
[params, condition_master, all_unique_im, all_cond] = ...
            vcd_createBlocksAndTrials(params,'load_params', params.load_params, ...
                                          'store_params', params.store_params, ...
                                          'session_type', session_type);
params.trials =  condition_master;

%% Create/Load miniblocks into runs and sessions, shuffle blocks within a run for each subject's run
% !!WARNING!! There is a randomization component involved in creating the
% block order within a run. If you don't want this, set second input
% (load_params) to true and load an existing file
%
% This function contains the following important steps/functions:
% * vcd_allocateBlocksToRuns
% * vcd_createRunTimeTables
% * vcd_addFIXandCDtoTimeTableMaster
[params,time_table_master] = vcd_createSessions(params,'load_params',  false, ...params.load_params, ...
                                                       'store_params', params.store_params, ...
                                                       'session_type', session_type);



