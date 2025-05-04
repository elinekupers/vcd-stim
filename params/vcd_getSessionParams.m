function exp = vcd_getSessionParams(varargin)
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
%  [verbose]             : (optional) Print text in command window or not
%                           Default: true
%
% OUTPUT:
%  exp                   : struct with experimental session params, including:
%                           * session, run, block, trial
%                           * tktktkt
%
% Written by Eline Kupers November 2024 (kupers [at] umn [dot] edu)

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addParameter('disp_name'             , '7TAS_BOLDSCREEN32', @(x) any(strcmp(x,{'7TAS_BOLDSCREEN32', 'KKOFFICE_AOCQ3277', 'PPROOM_EIZOFLEXSCAN', 'EKHOME_ASUSVE247'})));                   
p0.addParameter('presentationrate_hz'   , 30     , @isnumeric);
p0.addParameter('load_params'           , true   , @islogical);                    
p0.addParameter('store_params'          , true   , @islogical);
p0.addParameter('verbose'               , true   , @islogical); 

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
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('exp_%s*.mat',disp_name)));
    if ~isempty(d)
        if verbose
            fprintf('[%s]: Found %d exp params .mat file(s)\n',mfilename,length(d));
            if length(d) > 1
                warning('[%s]: Multiple .mat files! Will pick the most recent one', mfilename);
            end
            fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
        end
        load(fullfile(d(end).folder,d(end).name),'exp');
    else
        error('[%s]: Can''t find experiment session params file!\n', mfilename)
    end
else
    if verbose, fprintf('[%s]: Define exp params\n', mfilename); end
    
    % We will create a big struct where we separate params based on the
    % experimental hierarchy, i.e.: if it is at the session level, run
    % level, block, or trial level.
    exp = struct('session',[],'run',[],'block',[],'trial', []);
    
    %% %%%% STIMULUS - TASK CROSSINGS %%%%
    
    % Define big stim-task crossing table
    exp.stimclassnames = {'gabor','rdk','dot','obj','ns'};
    exp.taskclassnames = {'fix','cd','scc','pc','wm','ltm','img','what','where','how'};
    
    % Define stim-task crossings (master table)
    exp.crossings = false(length(exp.stimclassnames),length(exp.taskclassnames));
    
    exp.crossingnames = cell(length(exp.taskclassnames),length(exp.stimclassnames));
    for row = 1:size(exp.crossings,1)
        for col = 1:size(exp.crossings,2)
            % SCC and LTM will mix 4 stim classes (GBR, RDK, DOT, OBJ), NS
            % doesn't cross with SCC and will be treated separately in LTM
            % crossing (so only scenes in a LTM-NS block).
            if ismember(exp.taskclassnames{col},{'scc','ltm'})
                exp.crossingnames{col,row} = sprintf('%s-all',lower(exp.taskclassnames{col}));
            else
                exp.crossingnames{col,row} = sprintf('%s-%s',lower(exp.taskclassnames{col}),lower(exp.stimclassnames{row}));
            end
        end
    end
    
    % Remove any empty labels
    exp.crossingnames(cellfun(@isempty, exp.crossingnames)) = [];
    
    % Set Classic block 
    % NOTE: We keep scc x gabor/rdk/dot/obj and ltm x gabor/rdk/dot/obj/ns
    % despite the fact they are technically the same crossing.
    % We do this because we still want to ensure we sample equal
    % nr of stimuli from each stimulus class for this particular crossing.
    exp.crossings(1:5,1:7) = true; 
    exp.crossings(5,3)     = false; % we remove SCC-NS crossing

    % Set Naturalistic tail
    exp.crossings(4:5,8:10) = true;
    exp.crossings(4,9)      = false; % we remove WHERE-OBJ crossing
    
    % Organize stim-task labels according to the crossings
    exp.crossingnames = exp.crossingnames(exp.crossings');
                                                     
    % Remove scc duplicate stim task label
    scc_duplicates = find(~cellfun(@isempty, strfind(exp.crossingnames,'scc-all')));
    exp.crossingnames(scc_duplicates(2:end)) = [];

    % Remove scc duplicate stim task label
    ltm_duplicates = find(~cellfun(@isempty, strfind(exp.crossingnames,'ltm-all')));
    exp.crossingnames(ltm_duplicates(2:end)) = [];

    %% %%%% SESSION PARAMS %%%%
    
    % To achieve least 10 repetitions for each unique condition, we need to create 11+ trial repeats.
    % These unique trial repeats will then be divided across all blocks and runs, until we run out...
    % Some crossings have more repetitions, such as FIX, as we want at
    % least one FIX block per session.
    
                                   % fix cd  scc  pc  wm ltm img what where how
    exp.n_unique_trial_repeats_mri = [19  11  11  11  11  23  23   0   0   0;  % gabor
                                      19  11  11  11  11  23  23   0   0   0;  % rdk
                                      27  11  11  11  11  23  23   0   0   0;  % dot
                                      27  11  11  11  11  23  23   11  0   11; % obj
                                      11  11  0   11  11  11  11   11  11  11]; % ns
                                  
                                        % fix  cd  scc  pc  wm ltm img what where how
    exp.n_unique_trial_repeats_behavior = [3	3	3	3	3	0	0	0	0	0; % Gabor
                                           3	3	3	3	6	0	0	0	0	0; % RDK
                                           2	2	2	2	2	0	0	0	0	0; % Dot
                                           2	2	2	2	2	0	0	2	0	2; % Obj
                                           2	4	0	4	8	0	0	4	4	4]; % NS                          
                              

    exp.TR                             = 1.6;                               % repetiton time of the MRI pulse sequence in seconds
    exp.total_subjects                 = 3;                                 % 3 subjects for now.. EK: we probably want to separate wide and deep subjects
    
    % %%%% SESSION %%%%
    
    % BEHAVIORAL
    exp.session.behavior.session_nrs        = 1;                            % Sessions dedicated to behavior
    exp.session.behavior.n_runs_per_session = 15;                           % there are 15x ~5 min runs per session

    % WIDE
    exp.session.wide.session_nrs            = 1:2;                          % Sessions dedicated to wide subject sampling (WIDE 1A and WIDE 1B)
    exp.session.wide.n_runs_per_session     = [10,10];                      % WIDE session 1A and 1B both have 10x 6.05 min runs per session

    % DEEP
    exp.session.deep.n_runs_per_session = [repmat(10,1,25),4,4];            % there are 10x 6.05 min runs for DEEP sessions 1-25, 26A and 26B have 4 runs each
    exp.session.deep.session_nrs        = 1:26;                             % Sessions dedicated to deep subject sampling
    exp.session.deep.baseline_sessions  = 1:4;                              % Deep sessions dedicated to establish baseline, prior to introducing LTM/IMG

    % General
    exp.session.n_behavioral_sessions   = length(exp.session.behavior.session_nrs);
    exp.session.n_deep_sessions         = length(exp.session.deep.session_nrs); 
    exp.session.n_wide_sessions         = length(exp.session.wide.session_nrs); 
    exp.session.n_total_sessions        = exp.session.n_deep_sessions + exp.session.n_wide_sessions;
    
    % What session do we introduce/start sampling the tasks?
    exp.session.behavior.task_start    = [1,1,1,1,1,99,99,1,1,1];           % Everything but LTM/IMG in BEHAVIOR 01
    exp.session.wide.task_start        = [1,1,1,1,1,99,99,1,1,1];           % Everything but LTM/IMG in both WIDE01A and WIDE01B
    exp.session.deep.task_start        = [1,1,1,1,1,5,5,1,1,1];             % LTM/IMG start at DEEP06 (the 6th session), everything else starts at DEEP01
    
     
    %% %%%% RUN %%%%
    % total run time is 227 TRs or 363 s
    
    % general 
    exp.run.run_type1 = [7; 0]; % single-stim, double-stim blocks within a run
    exp.run.run_type2 = [4; 2]; % single-stim, double-stim blocks within a run
    exp.run.run_type3 = [1; 4]; % single-stim, double-stim blocks within a run
    exp.run.run_type4 = [0; 5]; % single-stim, double-stim blocks within a run
   
    % timing MRI
    exp.run.pre_blank_dur_MRI     = presentationrate_hz * 4.0;   % pre-run blank period: 4 seconds in number of presentation frames
    exp.run.post_blank_dur_MRI    = presentationrate_hz * 12.2;  % 12 seconds in number of presentation frames
    exp.run.total_run_dur_MRI     = presentationrate_hz * 363.2; % 363.2 s or 227 volumes (1.6 s TR)    
    assert(isint(exp.run.total_run_dur_MRI/exp.TR)); % ensure duration results in an integer nr of TRs

    % timing BEHAVIOR
    exp.run.pre_blank_dur_BEHAVIOR     = presentationrate_hz * 4.0;   % pre-run blank period: 4 seconds in number of presentation frames
    exp.run.post_blank_dur_BEHAVIOR    = presentationrate_hz * 4;     % post-blank period: 4 seconds in number of presentation frames
    exp.run.total_run_dur_BEHAVIOR     = presentationrate_hz * 343.2; % total run duration 343.2 s
    
    %% %%%% BLOCK PARAMS %%%%
    
    % general
    exp.block.n_trials_single_epoch = 8; % number of trials per block when we only have a single stimulus epoch 
    exp.block.n_trials_double_epoch = 4; % number of trials per block when we only have a two stimulus epochs (less trials because each trial is longer)
    
    % eye gaze block
    exp.block.nr_of_saccades      = 5;
    exp.block.eye_gaze_fix0       = presentationrate_hz * 1.0; % start with 1 second fixation period
    exp.block.eye_gaze_sac_target = presentationrate_hz * 1.2; % then 5x1.2 = 6 seconds of saccades (mimicing EL HV5 grid,±3 deg in all directions)
    exp.block.eye_gaze_fix1       = presentationrate_hz * 2.0; % then a 2-seconds rest trial
    exp.block.eye_gaze_pupil      = presentationrate_hz .* [3.0,1.0]; % then a 4-seconds pupil trial: 3-s black adaptation, 1-s white screen to evoke max pupil response.
    exp.block.total_eyetracking_block_dur = sum([exp.block.eye_gaze_fix0, ...
                                                        exp.block.eye_gaze_sac_target*exp.block.nr_of_saccades, ...
                                                        exp.block.eye_gaze_fix1, ...
                                                        exp.block.eye_gaze_pupil]);
    
    exp.block.eye_gaze_fix_ID         = 990;
    exp.block.eye_gaze_sac_target_ID  = 991:995; % central, left, right, up, down.
    exp.block.eye_gaze_pupil_ID       = 996;
    
    exp.run.actual_task_dur_MRI = exp.run.total_run_dur_MRI - exp.block.total_eyetracking_block_dur - exp.run.pre_blank_dur_MRI - exp.run.post_blank_dur_MRI; % nr of presentation frames we actually spend doing the experiment
    exp.run.actual_task_dur_BEHAVIOR = exp.run.total_run_dur_BEHAVIOR - exp.block.total_eyetracking_block_dur - exp.run.pre_blank_dur_BEHAVIOR - exp.run.post_blank_dur_BEHAVIOR; % nr of presentation frames we actually spend doing the experiment

    % event IDs
    exp.block.stim_epoch1_ID        = 91; % generic stim ID
    exp.block.stim_epoch2_ID        = 92; % generic stim ID
    exp.block.response_ID           = 93; % Time for subject to respond
    exp.block.trial_start_ID        = 94; % Fixation dot thickening
    exp.block.spatial_cue_ID        = 95; % Fixation dot turning black on L/R/both sides
    exp.block.delay_ID              = 96; % Delay period between two stimulus epochs
    exp.block.task_cue_ID           = 97; % Text on display to instruct subject
    exp.block.ITI_ID                = 98; % Inter-trial interval
    exp.block.IBI_ID                = 99; % Inter-block interval
    
    % Check if these IDs do not already exist in stim-task labels
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.response_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.trial_start_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.spatial_cue_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.task_cue_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.delay_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.task_cue_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.ITI_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.IBI_ID)));
    
    
    % Timing
    exp.block.task_cue_dur           = presentationrate_hz * 2.0;               % 2.0 seconds in number of presentation frames
    exp.block.IBI_MRI                = presentationrate_hz * linspace(5,9,5);   % [5:1:9] seconds Inter-block interval -- uniformly sample between [min,max]
    exp.block.IBI_BEHAVIOR           = presentationrate_hz * 4;                 % 4 seconds inter-block interval
    exp.block.total_single_epoch_dur =  presentationrate_hz * 42.0;             % 42.0 seconds in number of presentation frames (excl. IBI)
    exp.block.total_double_epoch_dur =  presentationrate_hz * 62.0;             % 62.0 seconds in number of presentation frames (excl. IBI)
    
    % Make we have integer number of frames
    assert(isint(exp.block.task_cue_dur));
    assert(all(isint(exp.block.IBI_MRI))); assert(all(isint(exp.block.IBI_BEHAVIOR)));

    % In each run, we have manipulations that we prioritize to fully sample,
    % otherwise it is difficult to compare conditions (e.g., we want to sample
    % all contrast levels within the run).
    exp.priority_stim_manip = struct('name',{},'priority',{},'other',{});
    exp.priority_stim_manip(1).name     = {'gabor'};
    exp.priority_stim_manip(1).priority = {'contrast','delta_from_ref'};   % First Priority manipulation
    exp.priority_stim_manip(1).other    = {'orient_deg'};                  % Other manipulations
    exp.priority_stim_manip(2).name     = {'rdk'};
    exp.priority_stim_manip(2).priority = {'coherence'};                   % First Priority manipulation
    exp.priority_stim_manip(2).other    = {'ori_bin','delta_from_ref'};
    exp.priority_stim_manip(3).name     = {'dot'};
    exp.priority_stim_manip(3).priority = {'ang_deg'};                     % First Priority manipulation
    exp.priority_stim_manip(3).other    = {''};
    exp.priority_stim_manip(4).name     = {'obj'};
    exp.priority_stim_manip(4).priority = {'super_cat'};                   % First Priority manipulation
    exp.priority_stim_manip(4).other    = {'basic_cat','sub_cat'};
    exp.priority_stim_manip(5).name     = {'ns'};
    exp.priority_stim_manip(5).priority = {'super_cat'};                   % First Priority manipulation
    exp.priority_stim_manip(5).other    = {'basic_cat','sub_cat'}; 
    
    
    
    %% %%%% TRIAL %%%%
    % general
    exp.trial.single_epoch_tasks = logical([1 1 1 1 0 0 0 1 1 1]);
    exp.trial.double_epoch_tasks = ~exp.trial.single_epoch_tasks;
    exp.trial.stim_LR_loc_cue    = logical([1 1 1 1 0]);
    
    % timing
    exp.trial.start_cue_dur       = presentationrate_hz * 0.4; % 12 x 33 ms frames = 0.4 seconds (thickening of dot rim)
    exp.trial.spatial_cue_dur     = presentationrate_hz * 0.8; % 24 x 33 ms frames = 0.8 seconds
    exp.trial.stim_array_dur      = presentationrate_hz * 2.0; % 60 x 33 ms frames = 2.0 seconds
    exp.trial.response_win_dur    = presentationrate_hz * 1.0; % 30 x 33 ms frames = 1.0 seconds
   
    exp.trial.totalITI            = presentationrate_hz .* [6.4, 3.2];
    exp.trial.ITI                 = presentationrate_hz .* [0.2:0.1:1.6]; % [6:2:48] frames corresponds to 0.2:0.1:1.6 seconds (thinning of dot rim)
    exp.trial.delay_dur           = presentationrate_hz * 8.0 ; % 240 x 33 ms frames = 8.0 seconds
    
    exp.trial.single_epoch_dur   = ...
        sum([exp.trial.start_cue_dur,... % frames
        exp.trial.spatial_cue_dur, ...
        exp.trial.stim_array_dur, ...
        exp.trial.response_win_dur]);
    
    exp.trial.double_epoch_dur   = ...
        sum([exp.trial.start_cue_dur,... % frames
        exp.trial.spatial_cue_dur, ...
        exp.trial.stim_array_dur, ...
        exp.trial.delay_dur, ...
        exp.trial.stim_array_dur, ...
        exp.trial.response_win_dur]);
    
    assert( nearZero(mod(exp.trial.single_epoch_dur / presentationrate_hz,1)))
    
    % TASK SPECIFIC PROBABILITY for trials
    
    % CD
    exp.trial.cd.prob_change                         = 0.5;  % chance of a contrast change
    
    % LTM
    exp.trial.ltm.prob_correct_pair                  = 0.5;  % chance of a given stim-stim pair in a trial is correct
    exp.trial.ltm.prob_incorrect_pair_same_stimclass = 0.25; % these are lures
    exp.trial.ltm.prob_incorrect_pair_diff_stimclass = 0.25; % these are non-lures 
    assert(sum([exp.trial.ltm.prob_correct_pair, ...
                exp.trial.ltm.prob_incorrect_pair_same_stimclass, ...
                exp.trial.ltm.prob_incorrect_pair_diff_stimclass])==1)

    exp.trial.ltm.prob_new_pairing                   = 0;    % chance that LTM stim A will be match to stim C (instead of stim B), in a given session
    exp.trial.ltm.prob_pair_order_flip               = 0;    % chance that LTM stim A -> B will flip to B -> A in a given session
    
    % IMG
    exp.trial.img.test_task                          = 0.5;  % chance of a contrast change

    
    %% Nr of blocks per sessions 
    
    % 1 behavioral session with 1 run only (all subjects do the same)
    exp.session.behavior.nr_of_type1_runs(1)  = 9; % 7 single-stim blocks / 0 double-stim blocks 
    exp.session.behavior.nr_of_type2_runs(1)  = 0; % 4 single-stim blocks / 2 double-stim blocks
    exp.session.behavior.nr_of_type3_runs(1)  = 6; % 1 single-stim blocks / 4 double-stim blocks
    exp.session.behavior.nr_of_type4_runs(1)  = 0; % 0 single-stim blocks / 5 double-stim blocks
    
    % 1 wide session with 2 run options (50% of subjects will see A or 50% of subjects will see B)
    exp.session.wide.nr_of_type1_runs([1,2])  = [6,6]; % 7 single-stim blocks / 0 double-stim blocks 
    exp.session.wide.nr_of_type2_runs([1,2])  = [0,0]; % 4 single-stim blocks / 2 double-stim blocks
    exp.session.wide.nr_of_type3_runs([1,2])  = [4,4]; % 1 single-stim blocks / 4 double-stim blocks
    exp.session.wide.nr_of_type4_runs([1,2])  = [0,0]; % 0 single-stim blocks / 5 double-stim blocks

    % 1-4 Deep sessions with 3 run options, only WM, no LTM/IMG
    exp.session.deep.nr_of_type1_runs([1:4])  = [6,5,5,5]; % 7 single-stim blocks / 0 double-stim blocks 
    exp.session.deep.nr_of_type2_runs([1:4])  = [0,0,1,0]; % 4 single-stim blocks / 2 double-stim blocks
    exp.session.deep.nr_of_type3_runs([1:4])  = [4,5,4,5]; % 1 single-stim blocks / 4 double-stim blocks
    exp.session.deep.nr_of_type4_runs([1:4])  = [0,0,0,0]; % 0 single-stim blocks / 5 double-stim blocks

    % 5-21 Deep sessions with 2 run options, all crossings
    exp.session.deep.nr_of_type1_runs([5:21])  = 3;
    exp.session.deep.nr_of_type2_runs([5:21])  = 0;
    exp.session.deep.nr_of_type3_runs([5:21])  = 0;
    exp.session.deep.nr_of_type4_runs([5:21])  = 7;

    % 22-24 Deep sessions with 2 run options, all crossings, winding down on single-stim blocks  
    exp.session.deep.nr_of_type1_runs([22:24])  = 4;
    exp.session.deep.nr_of_type2_runs([22:24])  = 0;
    exp.session.deep.nr_of_type3_runs([22:24])  = 0;
    exp.session.deep.nr_of_type4_runs([22:24])  = 6;
    
    % 22-24 Deep sessions with 2 run options, all crossings, winding down on single-stim blocks      
    exp.session.deep.nr_of_type1_runs([25:27])  = 5;
    exp.session.deep.nr_of_type2_runs([25:27])  = 0;
    exp.session.deep.nr_of_type3_runs([25:27])  = 0;
    exp.session.deep.nr_of_type4_runs([25:27])  = 5;
    
    exp.session.behavior.ses_blocks = zeros(size(exp.crossings,1),size(exp.crossings,2),exp.session.n_behavioral_sessions);
    exp.session.wide.ses_blocks = zeros(size(exp.crossings,1),size(exp.crossings,2),exp.session.n_wide_sessions);
    exp.session.deep.ses_blocks = zeros(size(exp.crossings,1),size(exp.crossings,2),exp.session.n_deep_sessions);

    
    % sessions for behavior 0                      fix cd scc  pc  wm ltm img what where  how  
    exp.session.behavior.ses_blocks(:,:,1) =    [3	3	3	3	3	0	0	0	0	0; % Gabor
                                                    3	3	3	3	6	0	0	0	0	0; % RDK
                                                    2	2	2	2	2	0	0	0	0	0; % Dot
                                                    2	2	2	2	2	0	0	2	0	2; % Obj
                                                    2	4	0	4	8	0	0	4	4	4]; % NS

    % sessions WIDE 1A                           fix cd scc  pc  wm ltm img what where  how  
    exp.session.wide.ses_blocks(:,:,1) =   [2	2	2	2	4	0	0	0	0	0;
                                            2	2	2	2	4	0	0	0	0	0;
                                            1	1	1	2	2	0	0	0	0	0;
                                            1	1	2	2	2	0	0	2	0	2;
                                            2	3	0	2	4	0	0	3	2	3]; % NS

    % sessions WIDE 1A                           fix cd scc  pc  wm ltm img what where  how  
    exp.session.wide.ses_blocks(:,:,2) =   [2	2	2	2	4	0	0	0	0	0;
                                            2	2	2	2	4	0	0	0	0	0;
                                            1	1	1	2	2	0	0	0	0	0;
                                            1	1	2	2	2	0	0	2	0	2;
                                            2	3	0	3	4	0	0	2	3	2];% NS

    % Deep sessions 1-4 have no LTM/IMG   fix   cd   scc pc  wm  ltm img what where  how                 
    exp.session.deep.ses_blocks(:,:,1) = [  1	2	2	3	4	0	0	0	0	0;
                                            1	2	2	3	3	0	0	0	0	0;
                                            1	1	1	2	2	0	0	0	0	0;
                                            1	1	1	2	3	0	0	2	0	1;
                                            2	3	0	3	4	0	0	3	3	3];% NS

    %                                    fix   cd   scc pc  wm  ltm img what where  how                      
    exp.session.deep.ses_blocks(:,:,2) =   [1	2	2	2	4	0	0	0	0	0;
                                            1	2	2	1	4	0	0	0	0	0;
                                            1	1	1	1	4	0	0	0	0	0;
                                            1	2	2	1	3	0	0	1	0	1;
                                            3	2	0	2	5	0	0	3	2	2];% NS
                                        
    %                                   fix   cd   scc pc  wm  ltm img what where  how          
    exp.session.deep.ses_blocks(:,:,3) = [ 1	2	2	2	4	0	0	0	0	0;
                                            1	2	2	2	4	0	0	0	0	0;
                                            1	2	2	1	2	0	0	0	0	0;
                                            1	2	2	1	3	0	0	1	0	1;
                                            3	2	0	3	5	0	0	2	2	3];% NS

    %                                   fix   cd   scc pc  wm  ltm img what where  how          
    exp.session.deep.ses_blocks(:,:,4) = [  1	2	2	1	4	0	0	0	0	0
                                            1	2	2	2	4	0	0	0	0	0
                                            1	2	2	1	3	0	0	0	0	0
                                            1	1	1	1	3	0	0	1	0	2
                                            2	2	0	3	6	0	0	2	3	2];% NS

                                                             
    % Deep sessions 5- have all stim-task crossings
    %                                    fix   cd   scc pc  wm  ltm img what where  how          
    exp.session.deep.ses_blocks(:,:,5) = [  1	1	1	1	2	2	3	0	0	0;
                                            1	1	1	1	2	2	3	0	0	0;
                                            1	1	0	0	1	2	2	0	0	0;
                                            1	0	1	1	1	2	2	1	0	1;
                                            1	1	0	1	3	4	4	1	1	1];% NS

    %                                    fix   cd   scc pc  wm  ltm img what where  how           
    exp.session.deep.ses_blocks(:,:,6) = [  1	1	1	1	2	3	2	0	0	0;
                                            1	1	1	1	2	3	2	0	0	0;
                                            1	1	1	0	1	2	2	0	0	0;
                                            1	1	0	1	1	2	2	1	0	1;
                                            1	1	0	1	3	4	4	1	0	1];% NS

    %                                      fix   cd   scc pc  wm  ltm img what where  how          
    exp.session.deep.ses_blocks(:,:,7) = [  1	1	1	1	3	2	3	0	0	0;
                                            1	1	0	1	2	2	2	0	0	0;
                                            1	1	1	0	1	2	2	0	0	0;
                                            1	1	0	1	1	2	2	1	0	1;
                                            1	1	0	1	3	4	4	1	1	1];% NS

                                        
                                        
    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how
    exp.session.deep.ses_blocks(:,:,8) = [  1		1		1		1		2		3		2		0		0		0; % Gabor
                                            1		1		1		1		2		2		2		0		0		0; % RDK
                                            1		1		0		1		1		3		2		0		0		0; % Dot
                                            1		0		1		1		1		2		3		0		0		1; % Obj
                                            1		1		0		1		3		3		4		1		1		1];% NS

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how
    exp.session.deep.ses_blocks(:,:,9) = [  1		1		1		1		3		2		3		0		0		0; % Gabor
                                            1		1		1		1		2		2		2		0		0		0; % RDK
                                            1		1		0		1		2		2		2		0		0		0; % Dot
                                            1		0		1		1		1		2		2		1		0		0; % Obj
                                            1		1		0		1		3		4		3		1		1		1];% NS

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how
    exp.session.deep.ses_blocks(:,:,10) = [  1		1		1		1		2		3		2		0		0		0; % Gabor
                                            1		1		0		1		2		2		2		0		0		0; % RDK
                                            1		1		0		1		2		2		2		0		0		0; % Dot
                                            1		1		0		1		2		2		2		1		0		1; % Obj
                                            1		1		0		1		4		3		3		1		1		1];% NS

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how
    exp.session.deep.ses_blocks(:,:,11) = [  1		1		1		1		3		2		2		0		0		0; % Gabor
                                            1		1		1		1		2		2		2		0		0		0; % RDK
                                            1		0		1		1		2		2		2		0		0		0; % Dot
                                            1		1		0		1		2		2		2		1		0		0; % Obj
                                            1		1		0		1		3		3		4		1		1		1];% NS

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how
    exp.session.deep.ses_blocks(:,:,12) = [  1		1		1		1		2		2		2		0		0		0; % Gabor
                                            1		1		1		1		2		2		2		0		0		0; % RDK
                                            1		1		0		1		1		2		2		0		0		0; % Dot
                                            1		0		0		1		2		3		2		1		0		1; % Obj
                                            1		1		0		1		3		4		4		1		1		1];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how
    exp.session.deep.ses_blocks(:,:,13) = [  1		1		1		1		3		2		2		0		0		0; % Gabor
                                            1		1		1		1		2		2		2		0		0		0; % RDK
                                            1		0		0		1		2		2		2		0		0		0; % Dot
                                            1		1		0		1		2		2		2		1		0		0; % Obj
                                            1		1		0		1		3		4		3		2		1		1];% NS

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how
    exp.session.deep.ses_blocks(:,:,14) = [  1		1		1		1		2		2		2		0		0		0; % Gabor
                                            1		1		0		1		2		2		2		0		0		0; % RDK
                                            1		0		0		1		2		2		2		0		0		0;
                                            1		1		1		0		2		2		2		1		0		1;
                                            1		1		0		1		3		4		4		1		2		1];% NS

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how
    exp.session.deep.ses_blocks(:,:,15) = [  1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		0		1		0		2		2		2		0		0		0;
                                            1		1		1		0		2		2		2		1		0		1;
                                            1		1		0		1		3		4		4		1		1		1];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how
    exp.session.deep.ses_blocks(:,:,16) = [  1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		0		1		0		2		2		2		0		0		1;
                                            1		1		0		1		3		4		4		1		1		1];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how        
    exp.session.deep.ses_blocks(:,:,17) = [  1		1		1		1		3		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		0		1		0		1		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		1;
                                            1		1		0		1		3		4		4		1		1		1];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how       
    exp.session.deep.ses_blocks(:,:,18) = [  1		1		1		1		3		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		0		1		1		2		2		2		0		0		0;
                                            1		1		0		1		2		4		4		1		1		1];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how            
    exp.session.deep.ses_blocks(:,:,19) = [  1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		2		2		3		0		0		0;
                                            1		0		1		0		2		2		3		0		0		0;
                                            1		1		1		0		1		2		2		1		0		1;
                                            1		1		0		1		2		4		4		1		1		1];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how               
    exp.session.deep.ses_blocks(:,:,20) = [  1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		1		2		3		0		0		0;
                                            1		1		1		0		1		3		3		0		0		0;
                                            1		1		0		1		2		4		4		1		1		1];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how           
    exp.session.deep.ses_blocks(:,:,21) = [  1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		2		3		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		1		0		0		1		2		2		1		0		0;
                                            1		1		0		1		3		4		4		1		1		1];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how           
    exp.session.deep.ses_blocks(:,:,22) = [  1		1		1		1		1		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		1		2		2		0		0		0;
                                            1		1		1		1		1		2		2		1		0		1;
                                            2		2		0		2		2		4		3		1		2		1];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how           
    exp.session.deep.ses_blocks(:,:,23) = [  1		1		1		1		1		2		2		0		0		0;
                                            1		1		2		1		2		2		2		0		0		0;
                                            1		1		1		1		1		3		2		0		0		0;
                                            1		1		1		1		1		2		2		1		0		1;
                                            2		2		0		1		1		4		3		1		1		2];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how           
    exp.session.deep.ses_blocks(:,:,24) = [  1		1		1		1		1		2		2		0		0		0;
                                            1		1		1		1		2		2		2		0		0		0;
                                            1		1		1		1		1		2		2		0		0		0;
                                            1		1		1		1		1		2		2		1		0		1;
                                            2		2		0		2		1		4		4		1		1		2];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how           
    exp.session.deep.ses_blocks(:,:,25) = [  1		1		2		1		0		2		2		0		0		0;
                                            1		1		3		1		1		3		2		0		0		0;
                                            1		2		2		1		1		2		2		0		0		0;
                                            1		1		2		1		0		2		2		0		0		1;
                                            2		2		0		2		0		2		4		2		2		2];

    %                                      fix      cd      scc     pc      wm      ltm     img     what    where  how           
    exp.session.deep.ses_blocks(:,:,26) = [  1		1		2		2		3		1		1		0		0		0;
                                            1		2		2		2		5		1		2		0		0		0;
                                            1		0		1		0		0		2		2		0		0		0;
                                            1		1		0		0		0		2		2		0		0		0;
                                            2		2		0		2		2		3		3		2		2		2];
    

    
    if store_params
        if verbose, fprintf('[%s]: Storing session data..\n',mfilename); end
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir, sprintf('exp_%s_%s.mat',disp_name, datestr(now,30))),'exp','-v7.3');
    end
end


return
    
    


