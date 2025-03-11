%% s_createDesignMatrix.m

% Stand-alone script to create and store the order in which stimuli are
% shown for each trial, miniblock, run, and session in VCD core experiment.

%% %%%%%%%%%%%%%%%%%%%
%%%%%% PARAMETERS %%%% 
%%%%%%%%%%%%%%%%%%%%%%

verbose        = true; % visualize stimuli or not
p.store_imgs   = true; % store visualization figures
saveFigsFolder = fullfile(vcd_rootPath,'figs'); % where to store visualization figures

% Get display params
dispname = 'KKOFFICE_AOCQ3277'; %'7TAS_BOLDSCREEN32'; % or 'KKOFFICE_AOCQ3277' or 'EKHOME_ASUSVE247' or 'PPROOM_EIZOFLEXSCAN'
p.disp   = vcd_getDisplayParams(dispname);

% Get stimulus parameters
p.load_params                 = true; % load stored params
p.store_params                = true;
p.overwrite_randomized_params = false; % DO NOT OVERWRITE PARAMS

% SETUP RNG 
rand('seed', sum(100*clock));
randn('seed', sum(100*clock));
p.rng.rand = rand;
p.rng.randn = randn;

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
% This function will Input 1: Stimulus type, choose from
% 'gabor','rdk','dot','cobj','ns','all' (default is 'all') Input 2: Display
% params struct (see vcd_getDisplayParams.m) Input 3: Load prior stored
% parameters or not. Input 4: Store generated parameters or not. Input 5:
% Overwrite stored parameters and regenerate probabilistic params
p.stim   = vcd_getStimParams('all',p.disp,p.load_params,p.store_params, p.overwrite_randomized_params); 

%% Define/Load experiment session params
p.exp    = vcd_getSessionParams(p,p.load_params,p.store_params);

%% Make/Load miniblocks with trials that sample unique stimuli from each class
% !!WARNING!! There is a randomization component involved in creating the 
% trial sequence (e.g., order of unique images within a miniblock). If you
% don't want this, set second input (load_params) to true.
p.trials = vcd_makeTrials(p,p.load_params,p.store_params);

%% Create/Load miniblocks into runs and sessions, shuffle blocks within a run for each subject's run
% !!WARNING!! There is a randomization component involved in creating the
% miniblock order within a run. If you don't want this, set second input
% (load_params) to true.
subject_sessions = vcd_createSessions(p,false,p.store_params);
