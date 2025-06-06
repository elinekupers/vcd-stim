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
%                           Default: 60 frames per second
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
p0.addParameter('presentationrate_hz'   , 60     , @isnumeric);
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
    
    % GENERAL PARAMS
    exp.TR                                              = 1.6;              % repetiton time of the MRI pulse sequence in seconds
    exp.total_subjects                                  = 1;                % 1 subjects for now.. EVERYONE SHARES THE SAME EXPERIMENTAL RUNS
    
    % What session do we introduce/start sampling the tasks?
    exp.session.behavior.task_start                     = [1,1,1,1,1,99,99,1,1,1];  % Everything but LTM/IMG in BEHAVIOR 01
    exp.session.mri.task_start                          = [1,1,1,1,1,6,6,1,1,1];    % LTM/IMG start at DEEP05 (the 6th session), everything else starts at DEEP01
    
    % BEHAVIORAL
    exp.session.behavior.session_nrs                    = 1;                % Sessions dedicated to behavioral experiment
    exp.session.behavior.session_types                  = 1;                % Only 1 BEHAVIOR session type (no A/B)
    exp.session.behavior.n_runs_per_session             = 12;               % There are 12 runs for the BEHAVIOR session

    % WIDE: Includes WIDE01A/WIDE01B
    % WIDE session 1A and 1B both have 10x 6.05 min runs per session
    exp.session.mri.wide.session_nrs                    = [1,1];            % [Session nr,session type]
    exp.session.mri.wide.session_types                  = [1,2];            % Refers to WIDE 1A and WIDE 1B
    exp.session.mri.wide.n_runs_per_session             = [10,10];          

    % DEEP: Includes DEEP001 to DEEP025, DEEP26A and DEEP26B
    % there are 10x 6.05 min runs for DEEP sessions 1-25, 26A and 26B have 4 runs each
    exp.session.mri.deep.session_nrs                    = [1:26];           %#ok<*NBRAK> % Session nr, deep subject sampling
    exp.session.mri.deep.session_types                  = NaN(26,2);        % [Session nr,session type] 
    exp.session.mri.deep.session_types(:,1)             = 1;                % insert session type = 1 for column one
    exp.session.mri.deep.session_types(end,2)           = 2;                % insert session type = 2 for column two (only last session).
    
    exp.session.mri.deep.n_runs_per_session             = NaN(26,2);       
    exp.session.mri.deep.n_runs_per_session(1:end-1,1)  = 10; 
    exp.session.mri.deep.n_runs_per_session(end,:)      = [5,5];
    exp.session.mri.deep.baseline_sessions              = 1:4;              % Deep sessions dedicated to establish baseline, prior to introducing LTM/IMG

    % MORE GENERAL PARAMS
    exp.session.n_behavioral_sessions                   = unique(exp.session.behavior.session_nrs);
    exp.session.n_deep_sessions                         = unique(exp.session.mri.deep.session_nrs);
    exp.session.n_wide_sessions                         = unique(exp.session.mri.wide.session_nrs);  
    exp.session.n_mri_sessions                          = [exp.session.n_wide_sessions,exp.session.n_deep_sessions];
   

    %% %%%% RUN %%%%
    % total run time is 242 TRs or 363 s
    
    % general 
%     exp.run.run_type1 = [7; 0]; % single-stim, double-stim blocks within a run
%     exp.run.run_type2 = [4; 2]; % single-stim, double-stim blocks within a run
%     exp.run.run_type3 = [1; 4]; % single-stim, double-stim blocks within a run
%     exp.run.run_type4 = [0; 5]; % single-stim, double-stim blocks within a run
   

    % timing MRI
    exp.run.pre_blank_dur_MRI          = presentationrate_hz * 4.0;         % pre-run blank period: 4 seconds in number of presentation frames
    exp.run.post_blank_dur_MRI         = presentationrate_hz * 12;          % 12 seconds in number of presentation frames
    exp.run.total_run_dur_BEHAVIOR     = presentationrate_hz * 387.2;       % TOTAL RUN DUR = 387.2 s or 242 EPI volumes (1.6 s TR) or 23232 time frames
    
    % timing BEHAVIOR
    exp.run.pre_blank_dur_BEHAVIOR     = presentationrate_hz * 4.0;         % pre-run blank period: 4 seconds in number of presentation frames
    exp.run.post_blank_dur_BEHAVIOR    = presentationrate_hz * 12.0;        % post-blank period: 4 seconds in number of presentation frames
    exp.run.total_run_dur_MRI          = presentationrate_hz * 387.2;       % TOTAL RUN DUR = 387.2 s or 242 EPI volumes (1.6 s TR) or 23232 time frames
    assert(isintnearzero(exp.run.total_run_dur_MRI/exp.TR));                % ensure MRI run duration results in an integer nr of TRs

   

    
    
    %% %%%% BLOCK PARAMS %%%%
    
    
    % General params
    exp.block.n_trials_single_epoch = 8;                                    % number of trials per block when we only have a single stimulus epoch 
    exp.block.n_trials_double_epoch = 4;                                    % number of trials per block when we only have a two stimulus epochs (less trials because each trial is longer)
    
    % Eye gaze block durations
    exp.block.nr_of_saccades       = 5;
    
    % ----- timing ---
    exp.block.eye_gaze_fix0        = presentationrate_hz * 1.0;             % start with 1 second fixation period
    exp.block.eye_gaze_sac_target  = presentationrate_hz * 1.5;             % then 5x1.2 = 6 seconds of saccades (mimicing EL HV5 grid,±3 deg in all directions)
    exp.block.eye_gaze_fix1        = presentationrate_hz * 2.0;             % then a 2-seconds rest trial
    exp.block.eye_gaze_pupil_black = presentationrate_hz * 3.0;             % then a 4-seconds pupil trial: 3-s black adaptation, 1-s white screen to evoke max pupil response.
    exp.block.eye_gaze_pupil_white = presentationrate_hz * 1.0;             % followed by a 1-s white screen to evoke max pupil response.
    exp.block.eye_gaze_fix2        = presentationrate_hz * 4.0;             % we end with 4 second fixation period (note this is separate from the pre-stimulus blank period).
    exp.block.total_eyetracking_block_dur = sum([exp.block.eye_gaze_fix0, ...
                                                        exp.block.eye_gaze_sac_target*exp.block.nr_of_saccades, ...
                                                        exp.block.eye_gaze_fix1, ...
                                                        exp.block.eye_gaze_pupil_black,...
                                                        exp.block.eye_gaze_pupil_white, ...
                                                        exp.block.eye_gaze_fix2]);
    % eye gaze block event IDs
    exp.block.eye_gaze_fix_ID         = 990;
    exp.block.eye_gaze_sac_target_ID  = 991:995;                            % central, left, right, up, down (last 4 are relative to central fixation.)
    exp.block.eye_gaze_pupil_black_ID = 996;
    exp.block.eye_gaze_pupil_white_ID = 997;

    % event IDs
    exp.block.task_cue_ID           = 90; % Text on display to instruct subject
    exp.block.post_task_cue_ITI_ID  = 91; % ITI between task cue and first spatial cue. Subjects will see a fixation circle with a white, THICK rim.
    exp.block.spatial_cue_ID        = 92; % Fixation circle rim turns red on left/right/both sides  
    exp.block.pre_stim_blank_ID     = 93; % Blank period in between spatial cue offset and stimulus onset
    exp.block.stim_epoch1_ID        = 94; % Stim onset (1st interval)
    exp.block.stim_epoch2_ID        = 95; % Stim onset (2nd interval after delay)
    exp.block.delay_ID              = 96; % Delay period between two stimulus epochs
    exp.block.response_ID           = 97; % Time for subject to respond
    exp.block.ITI_ID                = 98; % Inter-trial interval
    exp.block.IBI_ID                = 99; % Inter-block interval
    
    % Check if these IDs do not already exist in stim-task labels
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.pre_stim_blank_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.response_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.spatial_cue_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.task_cue_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.delay_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.task_cue_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.ITI_ID)));
    assert(isempty(intersect([1:length(exp.crossingnames)],exp.block.IBI_ID)));
    
    % Calculate min and max run durations based on duration of events.
    exp.run.min_run_dur_MRI        = exp.block.total_eyetracking_block_dur + exp.run.pre_blank_dur_MRI + exp.run.post_blank_dur_MRI;
    exp.run.min_run_dur_BEHAVIOR   = exp.block.total_eyetracking_block_dur + exp.run.pre_blank_dur_BEHAVIOR + exp.run.post_blank_dur_BEHAVIOR;
    
    % nr of presentation frames we actually spend doing the experiment
    exp.run.actual_task_dur_MRI      = exp.run.total_run_dur_MRI - exp.run.min_run_dur_MRI;
    exp.run.actual_task_dur_BEHAVIOR = exp.run.total_run_dur_BEHAVIOR - exp.run.min_run_dur_BEHAVIOR;
    
    
    %% %%%% TRIAL PARAMS %%%%
    
    % General params
    
    % Define what task classes have single vs double stim presentations
    exp.trial.single_epoch_tasks = logical([1 1 1 1 0 0 0 1 1 1]);
    exp.trial.double_epoch_tasks = ~exp.trial.single_epoch_tasks;
    
    % Define what stim classes are classic (i.e., have stimuli on the left
    % and right of central fixation).
    exp.trial.stim_LR_loc_cue    = logical([1 1 1 1 0]);
    
    %%%%% Timing %%%%
    
    % RUN/BLOCK LEVEL
    exp.block.IBI_MRI                = presentationrate_hz * [5:0.5:9];     % seconds inter-block interval -- uniformly sample between [min,max]
    exp.block.IBI_BEHAVIOR           = presentationrate_hz * [5:0.5:9];     % seconds inter-block interval -- currently the same as MRI exp
    exp.block.total_single_epoch_dur = presentationrate_hz * 44.5;          % 44.5 seconds in number of presentation frames (excl. IBI)
    exp.block.total_double_epoch_dur = presentationrate_hz * 60.5;          % 60.5 seconds in number of presentation frames (excl. IBI)
    
    % Make we have integer number of frames
    assert(all(isint(exp.block.IBI_MRI))); 
    assert(all(isint(exp.block.IBI_BEHAVIOR)));

    % BLOCK/TRIAL LEVEL
    exp.trial.task_cue_dur            = presentationrate_hz * 4.0;          % 4.0 seconds in number of presentation time frames
    exp.trial.post_task_cue_ITI_dur   = presentationrate_hz * 1.0;          % 1.0 second x 16.67 ms = 60 time frames (thick dot rim) 
    exp.trial.pre_stim_blank_dur      = presentationrate_hz * 0.5;          % 0.5 second x 16.67 ms = 30 time frames (thick dot rim, used in between spatial cue and stim onset) 
    exp.trial.spatial_cue_dur         = presentationrate_hz * 0.5;          % 0.5 second x 16.67 ms = 30 time frames
    exp.trial.stim_array_dur          = presentationrate_hz * 1.0;          % 1.0 second x 16.67 ms = 60 time frames 
    exp.trial.response_win_dur        = presentationrate_hz * 2.5;          % 2.5 seconds x 16.67 ms = 150 time frames
   
    exp.trial.ITI_single_block        = presentationrate_hz .* [0, 0, 0.5, 0.5, 0.5, 1, 1]; % in units of [0,0,30,30,30,60,60,60] time frames (each time frame is 16.67 ms). Presentation code will shuffle these ITIs prior to allocation within a block
    exp.trial.ITI_double_block        = presentationrate_hz .* [0, 0.5, 1]; % in units of [0,30,60] time frames (each time frame is 16.67 ms). Presentation code will shuffle these ITIs prior to allocation within a block
    exp.trial.delay_dur               = presentationrate_hz * 8.0;          % 8.0 seconds x 16.67 ms = 240 frames
    
    exp.trial.single_epoch_dur   = ...  % time frames
        sum([exp.trial.spatial_cue_dur, ... 
        exp.trial.pre_stim_blank_dur, ...
        exp.trial.stim_array_dur, ...
        exp.trial.response_win_dur]);
    
    exp.trial.double_epoch_dur   = ...  % time frames
        sum([exp.trial.spatial_cue_dur, ...
        exp.trial.pre_stim_blank_dur, ...
        exp.trial.stim_array_dur, ...
        exp.trial.delay_dur, ...
        exp.trial.stim_array_dur, ...
        exp.trial.response_win_dur]);
    
    assert( nearZero(mod(exp.trial.single_epoch_dur / 0.5*presentationrate_hz,1))); % ensure our durations are in 0.5 second time chunks

    %% In each run, we have manipulations that we prioritize to fully sample,
    % To make it easier to compare conditions (e.g., we want to sample
    % all contrast levels within the run). We mention these manipulations
    % here for archival purposes. Code does not actually use this param.
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
    
    
    
    %% TASK SPECIFIC PARAMS 
    
    % CD
    exp.trial.cd.prob_change                         = 0.2;  % probability of a contrast change across all cued trials. Note that for the uncued side, we do not change contrast.
    
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
    exp.trial.img.test_task                          = 0.5;  % chance of a yes/no overlapping test dots
 
    
    %% Nr of blocks per run and session
    % BEHAVIOR:
    % * There is 1 behavioral session with 11 runs and only one session type
    % hence second column is NaN. 
    % * Behavioral session 1 contains FIX/CD/SCC/PC/WM/WHAT/WHERE/HOW, no LTM/IMG.
    %
    % MRI:
    % * First MRI session is a wide session (only one wide session total).
    % * Wide session has with 2 session types (50% of subjects will see A or 50% of subjects will see B)
    % * Wide session contains FIX/CD/SCC/PC/WM/WHAT/WHERE/HOW, no LTM/IMG.
    % * Deep sessions 1-4 contain FIX/CD/SCC/PC/WM/WHAT/WHERE/HOW, no LTM/IMG.
    % * Only one session type (all subjetcs will do the same experimental
    % runs), hence second column is NaN.
    % * Deep sessions 5-20 contain all crossings.
    % * Deep sessions 21-24 contain all crossings, winding down on
    % single-stim blocks.
    % * Deep session 25 contains all crossings, winding down on single-stim
    % blocks.
    % * Last deep session (26) has two session types (A and B), with majority NS-stimulus crossings.

    exp.session.behavior.nr_blocks_per_run           = [7,NaN]; % 7 blocks/run, one session types
    exp.session.mri.wide.nr_blocks_per_run           = [7,7];   % 7 blocks/run, two session types
    exp.session.mri.deep.nr_blocks_per_run           = NaN(26,2);
    exp.session.mri.deep.nr_blocks_per_run(1:25,:)   = repmat([7,NaN],25,1); % 7 blocks/run, one session types
    exp.session.mri.deep.nr_blocks_per_run(26,:)     = [7,7];   % 7 blocks/run, two session types
    
    %% SPECIFIC CROSSINGS SHOWN PER SESSION
    % The goal is to achieve least 10 repetitions for each unique
    % condition. In practice, this means that we need to create 11+ trial 
    % repeats (e.g., due to increased repetition of special core stimuli, 
    % or the fact that the nr of unique core stimuli is not perfectly 
    % divisible by the nr of blocks per trial.
    % These unique trial repeats will then be divided across all blocks and
    % runs, and then we cut whatever we don't need. Some crossings have
    % more repetitions, such as FIX, as we want at least one FIX block per
    % session.
    exp.session.behavior.ses_blocks = zeros(size(exp.crossings,1),size(exp.crossings,2),length(exp.session.n_behavioral_sessions)); % last dim represents session type
    exp.session.wide.ses_blocks     = zeros(size(exp.crossings,1),size(exp.crossings,2),length(exp.session.n_wide_sessions),2);
    exp.session.deep.ses_blocks     = zeros(size(exp.crossings,1),size(exp.crossings,2),length(exp.session.n_deep_sessions),2);

    % sessions for behavior 1                  fix cd  scc  pc  wm ltm img what where how  
    exp.session.behavior.ses_blocks(:,:,1,1) = [1	1  0.5	 3	 6	 0	 0	 0	 0	  0; % gabor
                                                1	1  0.5	 3	 6	 0	 0	 0	 0	  0; % rdk
                                                1	1  0.5	 2	 4	 0	 0	 0	 0	  0; % dot
                                                1	1  0.5	 2	 4	 0	 0	 2	 0	  2; % obj
                                                1	1  0	 4   8	 0	 0	 4	 4	  4]; % ns


    % sessions WIDE 1A                    fix cd scc  pc  wm ltm img what where  how  
    exp.session.mri.wide.ses_blocks(:,:,1,1) = [2	2	2	2	4	0	0	0	0	0; % gabor
                                                2	2	2	2	4	0	0	0	0	0;
                                                1	1	1	2	2	0	0	0	0	0;
                                                1	1	2	2	2	0	0	2	0	2;
                                                2	3	0	2	4	0	0	3	2	3]; % NS

    % sessions WIDE 1B                    fix cd scc  pc  wm ltm img what where  how  
    exp.session.mri.wide.ses_blocks(:,:,1,2) = [2	2	2	2	4	0	0	0	0	0; % gabor
                                                2	2	2	2	4	0	0	0	0	0;
                                                1	1	1	2	2	0	0	0	0	0;
                                                1	1	2	2	2	0	0	2	0	2;
                                                3	2	0	3	4	0	0	2	3	2];% NS

    % Deep sessions 1-4 have no LTM/IMG   fix   cd   scc pc  wm  ltm img what where  how                 
    exp.session.mri.deep.ses_blocks(:,:,1,1) = [1	2	2	3	4	0	0	0	0	0; % gabor
                                                1	2	2	3	3	0	0	0	0	0;
                                                1	1	1	2	2	0	0	0	0	0;
                                                1	1	1	2	3	0	0	2	0	1;
                                                2	3	0	3	4	0	0	3	3	3];% NS

    %                                    fix   cd   scc pc  wm  ltm img what where  how                      
    exp.session.mri.deep.ses_blocks(:,:,2,1) = [1	2	2	2	4	0	0	0	0	0; % gabor
                                                1	2	2	1	4	0	0	0	0	0;
                                                1	1	1	1	4	0	0	0	0	0;
                                                1	2	2	1	3	0	0	1	0	1;
                                                3	3	0	2	5	0	0	3	2	2];% NS

    %                                   fix   cd   scc pc  wm  ltm img what where  how          
    exp.session.mri.deep.ses_blocks(:,:,3,1) = [1	2	2	2	4	0	0	0	0	0;
                                                1	2	2	2	4	0	0	0	0	0;
                                                1	2	2	1	2	0	0	0	0	0;
                                                1	2	2	1	3	0	0	1	0	1;
                                                3	2	0	3	5	0	0	2	2	3];% NS

    %                                   fix   cd   scc pc  wm  ltm img what where  how          
    exp.session.mri.deep.ses_blocks(:,:,4,1) = [1	2	2	1	4	0	0	0	0	0;
                                                1	2	2	2	4	0	0	0	0	0;
                                                1	2	2	1	3	0	0	0	0	0;
                                                1	1	1	1	3	0	0	1	0	2;
                                                2	2	0	3	6	0	0	2	3	2];% NS

                                                             
    % Deep sessions 5-end have all stim-task crossings
    %                                    fix   cd   scc pc  wm  ltm img what where  how          
    exp.session.mri.deep.ses_blocks(:,:,5,1) = [1	1	1	1	2	2	3	0	0	0;
                                                1	1	1	1	2	2	3	0	0	0;
                                                1	1	0	0	1	2	2	0	0	0;
                                                1	0	1	1	1	2	2	1	0	1;
                                                1	1	0	1	3	4	4	1	1	1];% NS

    %                                      fix cd  scc pc  wm ltm img what where how           
    exp.session.mri.deep.ses_blocks(:,:,6,1) = [1	1	1	1	2	3	2	0	0	0;
                                                1	1	1	1	2	3	2	0	0	0;
                                                1	1	1	0	1	2	2	0	0	0;
                                                1	1	0	1	1	2	2	1	0	1;
                                                1	1	0	1	3	4	4	1	0	1];% NS

    %                                      fix cd  scc pc  wm ltm img what where how     
    exp.session.mri.deep.ses_blocks(:,:,7,1) = [1	1	1	1	3	2	3	0	0	0;
                                                1	1	0	1	2	2	2	0	0	0;
                                                1	1	1	0	1	2	2	0	0	0;
                                                1	1	0	1	1	2	2	1	0	1;
                                                1	1	0	1	3	4	4	1	1	1];% NS                     
                                        
    %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,8,1) = [1	1	1	1	2	3	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	1	0	1	1	3	2	0	0	0;
                                                1	0	1	1	1	2	3	1	0	1;
                                                1	1	0	1	2	4	4	1	1	0];% NS

    %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,9,1) = [1	1	1	1	2	2	3	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	1	0	1	2	2	2	0	0	0;
                                                1	0	1	1	1	2	2	1	0	0;
                                                1	1	0	1	3	4	4	1	1	1];% NS

    %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,10,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	0	1	2	2	2	0	0	0;
                                                1	1	0	1	2	2	2	0	0	0;
                                                1	1	0	1	2	2	2	1	0	1;
                                                1	1	0	1	3	4	4	1	1	1];% NS

    %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,11,1) =[1	1	1	1	2	2	3	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	0	1	1	1	3	2	0	0	0;
                                                1	1	0	1	1	2	2	1	0	0;
                                                1	1	0	1	3	4	4	1	1	1];% NS

    %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,12,1) =[1	1	1	1	2	3	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	1	0	1	1	2	2	0	0	0;
                                                1	0	0	1	2	2	2	1	0	1;
                                                1	1	0	1	3	4	4	1	1	1];

    %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,13,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	0	0	1	2	2	2	0	0	0;
                                                1	1	0	1	2	2	2	1	0	0;
                                                1	1	0	1	3	4	4	2	1	1];% NS

    %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,14,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	0	1	2	2	2	0	0	0;
                                                1	0	0	1	2	2	2	0	0	0;
                                                1	1	1	0	2	2	2	1	0	1;
                                                1	1	0	1	3	4	4	1	2	1];% NS

    %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,15,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	0	1	0	2	2	2	0	0	0;
                                                1	1	1	0	2	2	2	1	0	1;
                                                1	1	0	1	3	4	4	1	1	1];

   %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,16,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	0	1	0	2	2	2	0	0	1;
                                                1	1	0	1	3	4	4	1	1	1];

   %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,17,1) =[1	1	1	1	3	2	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	0	1	0	1	2	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	1;
                                                1	1	0	1	3	4	4	1	1	1];

    %                                      fix cd  scc pc  wm ltm img what where how
    exp.session.mri.deep.ses_blocks(:,:,18,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	2	3	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	0	1	1	2	2	2	0	0	0;
                                                1	1	0	1	2	4	4	1	1	1];

    %                                      fix cd  scc pc  wm ltm img what where how       
    exp.session.mri.deep.ses_blocks(:,:,19,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	2	3	3	0	0	0
                                                1	0	1	0	1	2	3	0	0	0
                                                1	1	1	0	1	2	2	1	0	1
                                                1	1	0	1	2	4	4	1	1	1];

    %                                      fix cd  scc pc  wm ltm img what where how          
    exp.session.mri.deep.ses_blocks(:,:,20,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	3	2	2	0	0	0;
                                                1	1	1	1	2	2	3	0	0	0;
                                                1	1	1	0	2	2	2	0	0	0;
                                                1	1	0	1	3	3	3	1	1	1];

    %                                      fix cd  scc pc  wm ltm img what where how  
    exp.session.mri.deep.ses_blocks(:,:,21,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	1	2	2	0	0	0;
                                                1	1	1	1	1	2	2	1	0	1;
                                                2	2	0	2	2	3	3	1	2	1];

    %                                      fix cd  scc pc  wm ltm img what where how    
    exp.session.mri.deep.ses_blocks(:,:,22,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	1	2	2	0	0	0;
                                                1	1	1	1	1	2	2	1	0	1;
                                                2	2	0	2	2	3	3	1	2	1];

    %                                      fix cd  scc pc  wm ltm img what where how          
    exp.session.mri.deep.ses_blocks(:,:,23,1) =[1	1	1	1	2	2	2	0	0	0;
                                                1	1	2	1	2	2	2	0	0	0;
                                                1	1	1	1	2	2	2	0	0	0;
                                                1	1	1	1	1	2	2	1	0	1;
                                                2	2	0	1	1	3	3	1	1	2];

    %                                     % fix cd  scc  pc  wm ltm img what where how 
    exp.session.mri.deep.ses_blocks(:,:,24,1) =[1	1	1	1	2	1	1	0	0	0;
                                                1	1	2	1	3	1	2	0	0	0;
                                                1	1	2	1	2	1	1	0	0	0;
                                                1	1	1	1	2	3	2	0	0	1;
                                                2	1	0	1	2	4	3	2	1	2];

    %                                     % fix cd  scc  pc  wm ltm img what where how       
    exp.session.mri.deep.ses_blocks(:,:,25,1) =[1	2	2	2	2	1	0	0	0	0;
                                                1	2	3	2	3	1	1	0	0	0;
                                                1	1	1	1	1	2	2	0	0	0;
                                                1	1	0	0	0	2	2	0	0	0;
                                                2	2	0	2	1	4	3	2	3	3];

        %                                     % fix cd  scc  pc  wm ltm img what where how           
    exp.session.mri.deep.ses_blocks(:,:,26,1) =[0	0	0	0	0	1	1	0	0	0;
                                                0	0	0	0	0	1	1	0	0	0;
                                                0	0	0	0	0	1	1	0	0	0;
                                                0	0	0	0	0	1	1	0	0	0;
                                                2	2	0	3	3	0	2	3	3	3];
                                        
                                             % fix cd  scc  pc  wm ltm img what where how         
    exp.session.mri.deep.ses_blocks(:,:,26,2) =[0	0	0	0	0	1	1	0	0	0;
                                                0	0	0	0	0	1	1	0	0	0;
                                                0	0	0	0	0	1	1	0	0	0;
                                                0	0	0	0	0	1	1	0	0	0;
                                                1	3	0	2	3	0	2	4	2	4];

    % Get a summary of how many trials we expect/allocate across sessions.
    exp.nr_unique_trials_per_crossing          = exp.crossings.* [24,24,16,16,30]'; % 24 trials per stim class, etc 
    exp.nr_unique_trials_per_crossing(1:4,1)   = exp.nr_unique_trials_per_crossing(1:4,1)*0.5; % fixation crossings have half the nr of unique trials.
    exp.nr_unique_trials_per_crossing(:,[6:7]) = (exp.nr_unique_trials_per_crossing(:,[6:7]).* [1/3, 1/3, 1/2, 1/2, 1/2]'); % ltm/img crossings have special core images only (1/3 or 1/2 nr of core), but x2.
    exp.nr_trials_per_block = exp.crossings.* cat(2,repmat(exp.block.n_trials_single_epoch,1,4),repmat(exp.block.n_trials_double_epoch,1,3),repmat(exp.block.n_trials_single_epoch,1,3));
    
    ses_type_1_trials = ((sum(exp.session.deep.ses_blocks(:,:,:,1),3) + sum(exp.session.wide.ses_blocks(:,:,:,1),3)) .* (exp.nr_trials_per_block));
    ses_type_2_trials = ((sum(exp.session.deep.ses_blocks(:,:,:,2),3) + sum(exp.session.wide.ses_blocks(:,:,:,2),3)) .* (exp.nr_trials_per_block));
    exp.n_unique_trial_repeats_mri = ceil((ses_type_1_trials + ses_type_2_trials) ./  exp.nr_unique_trials_per_crossing) ;
   
    exp.n_unique_trial_repeats_behavior = ceil((exp.session.behavior.ses_blocks .*  exp.nr_trials_per_block)  ./  exp.nr_unique_trials_per_crossing);                                 
    
    % Store params struct if requested
    if store_params
        if verbose, fprintf('[%s]: Storing session data..\n',mfilename); end
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir, sprintf('exp_%s_%s.mat',disp_name, datestr(now,30))),'exp','-v7.3');
    end
end


return
    
    


