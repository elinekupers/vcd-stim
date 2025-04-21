function exp_session = vcd_getSessionParams(varargin)
% VCD function to get experimental session parameters related to block and
% run timing, order, and condition nr's.
%
%   exp_session = vcd_getSessionParams(disp_name,load_params,store_params)
%
% Stimulus params such as size and location will depend on display params.
%
% INPUTS:
%  [disp_name]           : (optional) Display name to load params (see vcd_getDisplayParams.m)
%                           Default: '7TAS_BOLDSCREEN32'
%  [presentationrate_hz] : (optional) Nr of stimulus frames presented per
%                           second,as we use a frame-locked stimulus presentation protocol, 
%                           Default: 30 frames per second
%  [load_params]         : (optional) Load prior stored parameters or not. 
%                           Default: true
%  [store_params]        : (optional) Store generated parameters or not.
%                           Will store params in fullfile(vcd_rootPath,'workspaces','info');
%                           Default: true
%
% OUTPUT:
%  exp_session           : struct with stimulus params, including:
%                           * priority_stim_manip
%                           * tktktkt
%
% Written by Eline Kupers November 2024 (kupers [at] umn [dot] edu)

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addParameter('disp_name'             , '7TAS_BOLDSCREEN32', @(x) any(strcmp(x,{'7TAS_BOLDSCREEN32', 'KKOFFICE_AOCQ3277', 'PPROOM_EIZOFLEXSCAN', 'EKHOME_ASUSVE247'})));                   
p0.addParameter('presentationrate_hz'   , 30     , @isnumeric);
p0.addParameter('load_params'           , true   , @islogical);                    
p0.addParameter('store_params'          , true   , @islogical); 

% Parse inputs
p0.parse(varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

%% Load params if requested
if load_params
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('exp_session_%s*.mat',disp_name)));
    if ~isempty(d)
        fprintf('[%s]: Found %d exp params .mat file(s)\n',mfilename,length(d));
        if length(d) > 1
            warning('[%s]: Multiple .mat files! Will pick the most recent one', mfilename);
        end
        fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
        load(fullfile(d(end).folder,d(end).name),'exp_session');
    else
        error('[%s]: Can''t find experiment session params file!\n', mfilename)
    end
else
    fprintf('[%s]: Define exp params\n', mfilename);
    
    % We will create a big struct where we separate params based on the
    % experimental hierarchy, i.e.: if it is at the session level, run
    % level, block, or trial level.
    exp_session = struct('session',[],'run',[],'block',[],'trial', []);
    
    %% %%%% STIMULUS - TASK CROSSINGS %%%%
    
    % Define big stim-task crossing table
    exp_session.stimClassLabels = {'gabor','rdk','dot','obj','ns'};
    exp_session.taskClassLabels = {'fix','cd','scc','pc','wm','ltm','img','what','where','how'};
    
    % Define stim-task crossings (master table)
    exp_session.crossings = false(length(exp_session.stimClassLabels),length(exp_session.taskClassLabels));
    
    exp_session.stimTaskLabels = cell(length(exp_session.taskClassLabels),length(exp_session.stimClassLabels));
    for row = 1:size(exp_session.crossings,1)
        for col = 1:size(exp_session.crossings,2)
            if strcmp(exp_session.taskClassLabels{col},'scc') % SCC will mix 4 stim classes (Gabors, RDKs, dot, obj)
                exp_session.stimTaskLabels{col,row} = sprintf('%s-all',lower(exp_session.taskClassLabels{col}));
            else
                exp_session.stimTaskLabels{col,row} = sprintf('%s-%s',lower(exp_session.taskClassLabels{col}),lower(exp_session.stimClassLabels{row}));
            end
        end
    end
    
    % Remove any empty labels
    exp_session.stimTaskLabels(cellfun(@isempty, exp_session.stimTaskLabels)) = [];
    
    % Set Classic block
    exp_session.crossings(1:5,1:7) = true; % we keep scc x gabor/rdk/dot/obj despite that they are technically the same crossing. 
    exp_session.crossings(5,3) = false; % we remove scc form NS.

    % Set Naturalistic tail
    exp_session.crossings(4:5,8:10) = true;
    exp_session.crossings(4,9) = false;
    
    % Organize stim-task labels according to the crossings
    exp_session.stimTaskLabels = exp_session.stimTaskLabels(exp_session.crossings');
    
    %% %%%% SESSION PARAMS %%%%
    exp_session.n_unique_trial_repeats = 8;   % 8 allows for allocation of all trials across blocks and runs.
    exp_session.n_sessions             = 5;   % 40; 
    exp_session.n_runs_per_session     = 10;  % right now we have 11x 5.5 min runs per session
    exp_session.TR                     = 1.6; % seconds
    exp_session.total_subjects         = 3;   % 3 subjects for now.. EK: we probably want to separate wide and deep subjects
    
    % %%%% SESSION %%%%
    exp_session.session.wideSessions     = [1,2]; % Sessions dedicated to wide subject sampling
    exp_session.session.baselineSessions = [1:6]; % Sessions dedicated to establish baseline
    
    % Response for LTM/IMAG (= 2 repeats of unique image conditions)
    exp_session.session.task_start       = [1,1,1,1,1,6,6,1,1,1]; % When do we start sampling the tasks (LTM/IMG have later starts)
    
     
    %% %%%% RUN %%%%
    % general
    exp_session.run.run_type1 = [7; 0]; % single-stim, double-stim blocks within a run, total run time is 227 TRs or 363 s
    exp_session.run.run_type2 = [4; 2];
    exp_session.run.run_type3 = [1; 4];
    exp_session.run.run_type4 = [0; 5]; 
   
    % timing
    exp_session.run.pre_blank_dur     = presentationrate_hz * 4.0; % pre-run blank period: 4 seconds in number of presentation frames
    exp_session.run.post_blank_dur    = presentationrate_hz * 12.2; % 12 seconds in number of presentation frames
    exp_session.run.total_run_dur     = presentationrate_hz * 363.2; % 363.2 s or 227 volumes (1.6 s TR)
    
    assert(isint(exp_session.run.total_run_dur/exp_session.TR)); % ensure this is an integer nr of TRs

    %% %%%% BLOCK PARAMS %%%%
    
    % general
    exp_session.block.n_trials_single_epoch = 8; % number of trials per block when we only have a single stimulus epoch 
    exp_session.block.n_trials_double_epoch = 4; % number of trials per block when we only have a two stimulus epochs (less trials because each trial is longer)
    
    % eye gaze block
    exp_session.block.nr_of_saccades      = 5;
    exp_session.block.eye_gaze_fix0       = presentationrate_hz * 1.0; % start with 1 second fixation period
    exp_session.block.eye_gaze_sac_target = presentationrate_hz * 1.2; % then 5x1.2 = 6 seconds of saccades (mimicing EL HV5 grid,±3 deg in all directions)
    exp_session.block.eye_gaze_fix1       = presentationrate_hz * 2.0; % then a 2-seconds rest trial
    exp_session.block.eye_gaze_pupil      = presentationrate_hz .* [3.0,1.0]; % then a 4-seconds pupil trial: 3-s black adaptation, 1-s white screen to evoke max pupil response.
    exp_session.block.total_eyetracking_block_dur = sum([exp_session.block.eye_gaze_fix0, ...
                                                        exp_session.block.eye_gaze_sac_target*exp_session.block.nr_of_saccades, ...
                                                        exp_session.block.eye_gaze_fix1, ...
                                                        exp_session.block.eye_gaze_pupil]);
    
    exp_session.block.eye_gaze_fix_ID         = 990;
    exp_session.block.eye_gaze_sac_target_ID  = [991:995]; % central, left, right, up, down.
    exp_session.block.eye_gaze_pupil_ID       = 996;
    
    exp_session.run.actual_task_dur = exp_session.run.total_run_dur - exp_session.block.total_eyetracking_block_dur - exp_session.run.pre_blank_dur - exp_session.run.post_blank_dur; % nr of presentation frames we actually spend doing the experiment
    
    % event IDs
    exp_session.block.stim_epoch1_ID        = 91; % generic stim ID
    exp_session.block.stim_epoch2_ID        = 92; % generic stim ID
    exp_session.block.response_ID           = 93; % Time for subject to respond
    exp_session.block.trial_start_ID        = 94; % Fixation dot thickening
    exp_session.block.spatial_cue_ID        = 95; % Fixation dot turning black on L/R/both sides
    exp_session.block.delay_ID              = 96; % Delay period between two stimulus epochs
    exp_session.block.task_cue_ID           = 97; % Text on display to instruct subject
    exp_session.block.ITI_ID                = 98; % Inter-trial interval
    exp_session.block.IBI_ID                = 99; % Inter-block interval
    
    % Check if these IDs do not already exist in stim-task labels
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.block.response_ID)));
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.block.trial_start_ID)));
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.block.spatial_cue_ID)));
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.block.task_cue_ID)));
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.block.delay_ID)));
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.block.task_cue_ID)));
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.block.ITI_ID)));
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.block.IBI_ID)));
    
    
    % Timing
    exp_session.block.task_cue_dur        = presentationrate_hz * 2.0; % 2.0 seconds in number of presentation frames
    exp_session.block.IBI                 = presentationrate_hz * linspace(5,9,5); % [5:1:9] seconds Inter-block interval -- uniformly sample between [min,max]
    
    exp_session.block.total_single_epoch_dur =  presentationrate_hz * 42.0;  % 42.0 seconds in number of presentation frames (excl. IBI)
    exp_session.block.total_double_epoch_dur =  presentationrate_hz * 62.0;  % 62.0 seconds in number of presentation frames (excl. IBI)
    
    % Make we have integer number of frames
    assert(isint(exp_session.block.task_cue_dur));
    assert(all(isint(exp_session.block.IBI)));

    % In each run, we have manipulations that we prioritize to fully sample,
    % otherwise it is difficult to compare conditions (e.g., we want to sample
    % all contrast levels within the run).
    exp_session.priority_stim_manip = struct('name',{},'priority',{},'other',{});
    
    exp_session.priority_stim_manip(1).name     = {'gabor'};
    exp_session.priority_stim_manip(1).priority = {'contrast','delta_ref'}; % Priority manipulation
    exp_session.priority_stim_manip(1).other    = {'ori_bin'};                 % Other manipulations
    exp_session.priority_stim_manip(2).name     = {'rdk'};
    exp_session.priority_stim_manip(2).priority = {'coherence','delta_ref'};
    exp_session.priority_stim_manip(2).other    = {'ori_bin'};
    exp_session.priority_stim_manip(3).name     = {'dot'};
    exp_session.priority_stim_manip(3).priority = {'ori_bin'};
    exp_session.priority_stim_manip(3).other    = {};
    exp_session.priority_stim_manip(4).name     = {'obj'};
    exp_session.priority_stim_manip(4).priority = {'super_cat','basic_cat'};
    exp_session.priority_stim_manip(4).other    = {'sub_cat'};
    exp_session.priority_stim_manip(5).name     = {'ns'};
    exp_session.priority_stim_manip(5).priority = {'super_cat','basic_cat'};
    exp_session.priority_stim_manip(5).other    = {'sub_cat'}; % ERK: do we sample/index these or does it not matter since stim are selected to be balanced?
    
    
    
    %% %%%% TRIAL %%%%
    % general
    exp_session.trial.single_epoch_tasks = logical([1 1 1 1 0 0 0 1 1 1]);
    exp_session.trial.double_epoch_tasks = ~exp_session.trial.single_epoch_tasks;
    exp_session.trial.stim_LR_loc_cue    = logical([1 1 1 1 0]);
    
    % timing
    exp_session.trial.start_cue_dur       = presentationrate_hz * 0.4; % 12 x 33 ms frames = 0.4 seconds (thickening of dot rim)
    exp_session.trial.spatial_cue_dur     = presentationrate_hz * 0.8; % 24 x 33 ms frames = 0.8 seconds
    exp_session.trial.stim_array_dur      = presentationrate_hz * 2.0; % 60 x 33 ms frames = 2.0 seconds
    exp_session.trial.response_win_dur    = presentationrate_hz * 1.0; % 30 x 33 ms frames = 1.0 seconds
   
    exp_session.trial.totalITI            = presentationrate_hz .* [6.4, 3.2];
    exp_session.trial.ITI                 = presentationrate_hz .* [0.2:0.1:1.6]; % [6:2:48] frames corresponds to 0.2:0.1:1.6 seconds (thinning of dot rim)
    exp_session.trial.delay_dur           = presentationrate_hz * 8.0 ; % 240 x 33 ms frames = 8.0 seconds
    
    exp_session.trial.single_epoch_dur   = ...
        sum([exp_session.trial.start_cue_dur,... % frames
        exp_session.trial.spatial_cue_dur, ...
        exp_session.trial.stim_array_dur, ...
        exp_session.trial.response_win_dur]);
    
    exp_session.trial.double_epoch_dur   = ...
        sum([exp_session.trial.start_cue_dur,... % frames
        exp_session.trial.spatial_cue_dur, ...
        exp_session.trial.stim_array_dur, ...
        exp_session.trial.delay_dur, ...
        exp_session.trial.stim_array_dur, ...
        exp_session.trial.response_win_dur]);
    
    assert( nearZero(mod(exp_session.trial.single_epoch_dur / presentationrate_hz,1)))
    
    %% TASK SPECIFIC PROBABILITY
    
    % CD
    exp_session.trial.cd.prob_change                         = 0.5;  % chance of a contrast change
    
    % LTM
    exp_session.trial.ltm.prob_correct_pair                  = 0.5;  % chance of a given stim-stim pair in a trial is correct
    exp_session.trial.ltm.prob_incorrect_pair_same_stimclass = 0.25; % these are lures
    exp_session.trial.ltm.prob_incorrect_pair_diff_stimclass = 0.25; % these are non-lures 
    assert(sum([exp_session.trial.ltm.prob_correct_pair, exp_session.trial.ltm.prob_incorrect_pair_same_stimclass, exp_session.trial.ltm.prob_incorrect_pair_diff_stimclass])==1)

    exp_session.session.ltm.prob_new_pairing                = 0.2;   % chance that LTM stim A will be match to stim C (instead of stim B), in a given session
    exp_session.session.ltm.prob_pair_order_flip            = 0.2;   % chance that LTM stim A -> B will flip to B -> A in a given session
    
    % IMG
    exp_session.trial.img.test_task                         = 0.5;  % chance of a contrast change
    exp_session.trial.img.quiz_images                       = [ones(1,20),2.*ones(1,20)];  % quiz dots overlap (1) or not (2)


    
    %% Nr of blocks per sessions 
    
    % 2 wide sessions with 23 runs
    exp_session.session.nr_of_type1_runs([1,2])  = [6,6]; % 7 single-stim blocks / 0 double-stim blocks 
    exp_session.session.nr_of_type2_runs([1,2])  = [0,0]; % 4 single-stim blocks / 2 double-stim blocks
    exp_session.session.nr_of_type3_runs([1,2])  = [4,4]; % 1 single-stim blocks / 4 double-stim blocks
    exp_session.session.nr_of_type4_runs([1,2])  = [0,0]; % 0 single-stim blocks / 5 double-stim blocks
    
    exp_session.session.nr_of_type1_runs([3:5])  = [5,6,5]; % 7 single-stim blocks / 0 double-stim blocks 
    exp_session.session.nr_of_type2_runs([3:5])  = [0,0,0]; % 4 single-stim blocks / 2 double-stim blocks
    exp_session.session.nr_of_type3_runs([3:5])  = [5,4,5]; % 1 single-stim blocks / 4 double-stim blocks
    exp_session.session.nr_of_type4_runs([3:5])  = [0,0,0]; % 0 single-stim blocks / 5 double-stim blocks

    exp_session.session.nr_of_type1_runs([6:22])  = 3;
    exp_session.session.nr_of_type2_runs([6:22])  = 0;
    exp_session.session.nr_of_type3_runs([6:22])  = 0;
    exp_session.session.nr_of_type4_runs([6:22])  = 7;
    
    exp_session.session.nr_of_type1_runs([23:27])  = 4;
    exp_session.session.nr_of_type2_runs([23:27])  = 0;
    exp_session.session.nr_of_type3_runs([23:27])  = 0;
    exp_session.session.nr_of_type4_runs([23:27])  = 6;
    
    ses_blocks = zeros(size(exp_session.crossings,1),size(exp_session.crossings,2),exp_session.n_sessions);
    
    % sessions 1-2 are WIDE
    %                   fix cd scc  pc  wm ltm img what where  how  
    ses_blocks(:,:,1) = [1	 2	 2	 3	 3	 0	 0	 0	 0	   0; % Gabor:
                         1	 2	 2	 3	 4	 0	 0	 0	 0	   0; % RDK:
                         1	 1	 1	 2	 3	 0	 0	 0	 0	   0; % Dot:
                         1	 1	 1	 2	 2	 0	 0	 1	 0	   2; % Obj:
                         2	 3	 0	 3	 4	 0	 0	 3	 3	   3];% NS: 
                     
    %                   fix  cd   scc   pc    wm    ltm   img   what where  how                 
    ses_blocks(:,:,2) = [1    2    2	3	4	0	0	0	0	0;% Gabor:
                         1    2	   2	3	3	0	0	0	0	0;% RDK:
                         1    1	   1	2	2	0	0	0	0	0;% Dot:
                         1    1	   1	2	3	0	0	2	0	1; % Obj:
                         2    3    0    3	4	0	0	3	3	3];% NS: 

    % sessions 3-6 have no LTM and no IMG
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how                 
    ses_blocks(:,:,3) = [1      2       2	2	4	0	0	0	0	0;  % Gabor:
                         1      2       2	1	4	0	0	0	0	0;  % RDK:
                         1      1       1	1	3	0	0	0	0	0;  % Dot:
                         1      2       2	1	3	0	0	1	0	1;  % Obj:
                         3      3       0	2	6	0	0	3	2	2]; % NS:

    ses_blocks(:,:,4) = [1      2       2	2	4	0	0	0	0	0;  % Gabor:
                         1      2       2	2	4	0	0	0	0	0;  % RDK:
                         1      2       2	1	2	0	0	0	0	0;  % Dot:
                         1      2       2	1	2	0	0	2	0	1;  % Obj:
                         3      3       0	3	4	0	0	3	2	3]; % NS
                     
    ses_blocks(:,:,5) = [1      2       2	1	4	0	0	0	0	0;  % Gabor:
                         1      2       2	2	4	0	0	0	0	0;  % RDK:
                         1      2       2	1	3	0	0	0	0	0;  % Dot:
                         1      1       1	1	3	0	0	1	0	2;  % Obj:
                         2      2       0	3	6	0	0	2	3	2]; % NS: 
                     

                     
    % sessions 6- have all stim-task crossings
    % we alternate 3 vs 2 LTM & IMG blocks for Gabor & RDK
    % Dot: starting session 7, we have 2 IMG blocks every 3 sessions. 
    % Dot: starting session 10, we have 2 LTM blocks every 3 sessions. 
    ses_blocks(:,:,6) = [1     2     1     2     3     0     0     0     0     0; % Gabor:
                         1     2     1     2     3     0     0     0     0     0; % RDK:
                         1     1     1     1     2     0     0     0     0     0; % Dot:
                         1     1     1     2     2     0     0     2     0     2; % Obj:
                         1     2     0     2     2     0     0     2     1     1]; % NS:
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how                 
    ses_blocks(:,:,7) = [1     1     1     1     2     2     2     0     0     0; % Gabor: 
                         1     0     1     0     1     2     2     0     0     0; % RDK: 
                         1     1     1     1     1     2     2     0     0     0; % Dot: 
                         1     1     1     1     1     2     2     1     0     1; % Obj: -
                         1     1     0     0     1     2     2     1     1     1]; % NS: 

    %                   fix    cd   scc   pc    wm    ltm   img   what where  how                     
    ses_blocks(:,:,8) = [1     0     1     1     1     2     2     0     0     0; % Gabor: 
                         1     1     1     1     2     2     2     0     0     0; % RDK: 
                         1     1     1     1     1     2     1     0     0     0; % Dot: 
                         1     1     1     1     1     2     1     1     0     1; % Obj: 
                         1     1     0     1     1     2     2     1     0     1]; % NS: 

    %                   fix    cd   scc   pc    wm    ltm   img   what where  how                     
    ses_blocks(:,:,9) = [1     1     1     1     2     2     2     0     0     0; % Gabor: 
                         1     0     1     1     1     2     2     0     0     0; % RDK: 
                         1     1     1     1     1     2     2     0     0     0; % Dot: -
                         1     1     1     1     1     2     1     1     0     1; % Obj: -
                         1     0     0     0     2     2     2     0     1     0]; % NS: 

    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,10) =[1     0     1     1     1     2     2     0     0     0; % Gabor: 3 LTM blocks
                         1     1     1     1     2     2     2     0     0     0; % RDK: 3 IMG blocks 
                         1     1     1     0     1     2     1     0     0     0; % Dot: 
                         1     1     1     0     1     2     1     1     0     1; % Obj: skip PC block
                         1     1     0     1     1     2     2     1     0     1]; % NS: 2 IMG & 2 LTM blocks

    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,11) =[1     1     1     1     1     2     2     0     0     0; % Gabor: 3 IMG blocks
                         1     0     1     1     1     2     2     0     0     0; % RDK: 3 LTM blocks 
                         1     1     1     1     1     2     2     0     0     0; % Dot: 2 IMG blocks
                         1     1     1     1     1     2     1     1     0     1;
                         1     0     0     0     1     2     2     0     1     0];

    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,12) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];

    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,13) =[1     1     1     1     1     2     2     0     0     0;
                         1     0     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];

    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,14) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     1     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];

    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,15) =[1     1     1     1     1     2     2     0     0     0;
                         1     0     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];

    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,16) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];
                     
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,17) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];
                     
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,18) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];                     
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,19) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];
                     
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,20) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];
                     
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,21) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];                     
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,22) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];
                     
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,23) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];
                     
    %                   fix    cd   scc   pc    wm    ltm   img   what where  how
    ses_blocks(:,:,24) =[1     0     1     1     1     2     2     0     0     0;
                         1     1     1     1     1     2     2     0     0     0;
                         1     0     1     0     1     2     1     0     0     0; % Dot: skip PC block
                         1     0     1     0     1     2     1     1     0     1;
                         1     1     0     1     2     2     1     1     0     1];                     

    exp_session.ses_blocks = ses_blocks;
    
    % Remove scc duplicate stim task label
    scc_duplicates = find(~cellfun(@isempty, strfind(exp_session.stimTaskLabels,'scc-all')));
    exp_session.stimTaskLabels(scc_duplicates(2:end)) = [];

    % Remove scc duplicate stim task label
    ltm_duplicates = find(~cellfun(@isempty, strfind(exp_session.stimTaskLabels,'ltm-all')));
    exp_session.stimTaskLabels(ltm_duplicates(2:end)) = [];
    
    if store_params
        fprintf('[%s]:Storing session data..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir, sprintf('exp_session_%s_%s.mat',disp_name, datestr(now,30))),'exp_session','-v7.3');
    end
end


return
    
    


