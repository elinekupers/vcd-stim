%% s_createDesignMatrix.m
%
% Script to create and store the entire VCD experiment as a MATLAB table
% including the order in which VCD-core stimuli and tasks are shown to  
% subject at the granularity of trial events.
%
% At a bird's eye level: there are sessions (~1.5-2 hrs), which contain
% multiple runs (~6.5 mins), which contain multiple blocks interspersed 
% with rest periods of a gray (mean luminance) screen + fixation circle
% called Inter-Block-Intervals (IBIs, 5-9 sec) or pre/post-run rest
% periods.
%
% Each block is 44.5 seconds (single-stimulus presentation trials) or 60.5 
% seconds (double-stimulus presentation trials) and contains multiple 
% trials of a single stimulus-task crossing. For a given trial, the subject 
% is presented with gray background + fixation circle that the subject is 
% instructed to continuously fixate at throughout the run.
%
% A single block contains starts with a task cue (4 seconds, followed by a 
% 1-s blank screen), followed by either 8 trials (for single-stimulus 
% presentation blocks, also known as trial_type 1) or 4 trials (double-
% stimulus presentation blocks, also known as trial_type 2). Trials within
% a block are interspersed with variable (0-1.5 sec) Inter-Trial-Intervals
% (ITIs).
% 
% A single trial (single: 4.5 sec or double: 13.5 sec, without ITI)
% contains the following events: 
%     SINGLE stimulus-presentation trial:
%       * covert spatial attention cue (0.5 sec)
%       * gap (gray screen             (0.5 sec)
%       * stimulus array               (1.0 sec)
%       * response window              (2.5 sec)
%       * ITI                          (0-1.5 sec)
%     DOUBLE stimulus-presentation trial:
%       * covert spatial attention cue (0.5 sec)
%       * gap (gray screen             (0.5 sec)
%       * stimulus array 1             (1.0 sec)
%       * delay                        (8.0 sec)
%       * stimulus array 2             (1.0 sec)
%       * response window              (2.5 sec)
%       * ITI                          (0-1.5 sec)
%
% Code dependencies:
%   * vcd-stim code repo (github.com/elinekupers/vcd-stim)
%   * knkutils code repo (github.com/cvnlab/knkutils)
%
%
% Written by Eline Kupers @ UMN
% History:
% v0.0: 11/2024 - first attempt
% v0.1: 03/2025 - refactoring code
% v0.2: 06/2025 - minor tweaks

%% %%%%%%%%%%%%%%%%%%%
%%%%%% PARAMETERS %%%% 
%%%%%%%%%%%%%%%%%%%%%%

verbose      = true;  % print text and visualize stimuli (true) or not (false)?
store_imgs   = true;  % store debug figures (true) or not (false)?
load_params  = false; % load stored params (true) or recreate them (false)?
store_params = true;  % save stored params as mat file (true) or not (false)?

params = struct();
params.saveFigsFolder = fullfile(vcd_rootPath,'figs'); % where to store visualization figures

% Get display params
% Choose from: '7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','EKHOME_ASUSVE247','PPROOM_EIZOFLEXSCAN','CCNYU_VIEWPIXX3D'
dispname    = '7TAS_BOLDSCREEN32'; 
params.disp = vcd_getDisplayParams(dispname);

% Infer environment type
if strcmp(dispname,'7TAS_BOLDSCREEN32')
    env_type = 'MRI';
elseif ismember(dispname,{'PPROOM_EIZOFLEXSCAN','KKOFFICE_AOCQ3277','EKHOME_ASUSVE247','CCNYU_VIEWPIXX3D'})
    env_type = 'BEHAVIOR';
end

% Flag if this is a demo run (true) or not (false)..
params.is_demo = false;

% Flag if this is a wide MRI experiment (true) or deep MRI experiment (false)..
params.is_wide = false;

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
                                  'load_params', load_params, ...
                                  'store_params', store_params); 

%% Define/Load experiment session params
params.exp    = vcd_getSessionParams('disp_name', params.disp.name, ...
                                'presentationrate_hz',params.stim.presentationrate_hz, ...
                                'load_params', load_params, ...
                                'store_params', store_params);

%% Create unique conditions
% !!WARNING!! There is a randomization component involved in creating the
% conditions (i.e., order of trials within a session). If you don't want
% this, set second input (load_params) to true and load an existing
% condition_master.
%
% This function will create the unique stimuli and conditions for each
% task-stim crossing occuring across all sessions (MRI or BEHAVIOR).
%
% This function contains the following important steps/subfunctions:
% * vcd_defineUniqueImages
% * vcd_createConditionMaster
% * vcd_shuffleStimForTaskClass
% * vcd_getCorrectButtonResponse
% * (no subfunction): Select object catch trials for OBJ-PC crossing
% * vcd_balanceButtonCorrectPresses
% * vcd_allocateBlocksToRuns
% * vcd_determineContrastDecrementChangeTrials

[params, condition_master, all_unique_im, all_cond] = ...
            vcd_createConditions(params,  'load_params',  load_params, ...
                                          'store_params', store_params, ...
                                          'env_type', env_type);

%% Create/Load miniblocks into runs and sessions, shuffle blocks within a run for each subject's run
% !!WARNING!! There is a randomization component involved in creating the
% block order within a run. If you don't want this, set second input
% (load_params) to true and load an existing file
%
% This function contains the following important steps/functions:
% * vcd_createRunTimeTables
% * vcd_createRunFrames
[params,condition_master_shuffled,time_table_master_shuffled, all_subj_run_frames] = ...
    vcd_createSessions(params,'condition_master',condition_master,...
                               'load_params',  load_params, ...
                               'store_params', store_params, ...
                               'env_type', env_type);



