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
params.verbose        = true; % visualize stimuli or not
params.store_imgs     = true; % store visualization figures
params.saveFigsFolder = fullfile(vcd_rootPath,'figs'); % where to store visualization figures

% Get display params
dispname = '7TAS_BOLDSCREEN32'; % Choose from: '7TAS_BOLDSCREEN32' or 'KKOFFICE_AOCQ3277' or 'EKHOME_ASUSVE247' or 'PPROOM_EIZOFLEXSCAN'
params.disp   = vcd_getDisplayParams(dispname);

% Get stimulus parameters
params.load_params                 = true; % load stored params
params.store_params                = true;
params.overwrite_randomized_params = false; % DO NOT OVERWRITE PARAMS

% SETUP RNG (EK Q: why not set rng itself? will that work for older matlab versions?)
params.rng.rand_seed = sum(100*clock);
rand('seed', params.rng.rand_seed);
params.rng.randn_seed = sum(100*clock);
randn('seed', params.rng.randn_seed);

%% Define/Load stimulus params 
% !!WARNING!! There is a randomization component involved in creating some
% stimuli (e.g., orientation of gabor stimuli or dot locations). If you
% don't want this, this leave the fifth argument
% (overwrite_randomized_params) empty (default is set to false) or set to
% false.
%
% If you do want regenerate probabilistic params, set the fifth argument to
% true and some stimulus values will change.
%
% Input 1: Stimulus class, choose from
% 'gabor','rdk','dot','obj','ns','all' (default is 'all') 
% Input 2: Display name to load disp params struct (see vcd_getDisplayParams.m) 
% Input 3: Load prior stored parameters or not? (logical)
% Input 4: Store generated parameters or not? (logical)
% Input 5: Overwrite stored parameters and regenerate probabilistic params?
% (logical)
params.stim   = vcd_getStimParams('stim_class','all',...
                            'disp_name', params.disp.name, ...
                            'load_params', params.load_params, ...
                            'store_params', params.store_params, ...
                            'overwrite_randomized_params', params.overwrite_randomized_params); 

%% Define/Load experiment session params
params.exp    = vcd_getSessionParams('disp_name', params.disp.name, ...
                                'presentationrate_hz',params.stim.presentationrate_hz, ...
                                'load_params', params.load_params, ...
                                'store_params', params.store_params);

%% Make/Load blocks with trials that sample unique stimuli from each class
% !!WARNING!! There is a randomization component involved in creating the 
% trial sequence (e.g., order of unique images within a block). If you
% don't want this, set second input (load_params) to true.
[params, condition_master, all_unique_im, all_cond] = ...
            vcd_createBlocksAndTrials(params,'load_params', params.load_params, ...
                                          'store_params', params.store_params);
params.trials =  condition_master;

%% Create/Load miniblocks into runs and sessions, shuffle blocks within a run for each subject's run
% !!WARNING!! There is a randomization component involved in creating the
% block order within a run. If you don't want this, set second input
% (load_params) to true.
[params,time_table_master] = vcd_createSessions(params,'load_params',false, ...
                                          'store_params', params.store_params);
                                      
                                      
%% Select unique image for a given subject run

% test on subject 1, run 1 for now
for sj = 1:params.exp.total_subjects
    for ses = 1:params.exp.n_sessions
        for rr = 1:params.exp.n_runs_per_session
        
            subj_run = time_table_master((time_table_master.subj_nr==sj & time_table_master.run_nr==rr & time_table_master.session_nr==ses),:);
            
            % run_images is a cell matrix with trials x locations (1:l, 2:r)
            [run_images, run_alpha_masks] = vcd_getImageOrderSingleRun(params, subj_run, sj, ses, rr, ...
                'load_params',false,'store_params', true);
            
            % run_images is a cell matrix with trials x locations (1:l, 2:r)
            [run_images, run_alpha_masks] = vcd_expandImageOrderSingleRun_30Hz(params, subj_run, sj, ses, rr, ...
                'load_params',false,'store_params', true);
            
        end
    end
end
    




