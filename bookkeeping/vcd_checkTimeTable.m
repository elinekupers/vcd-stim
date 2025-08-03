function results = vcd_checkTimeTable(time_table_master,varargin)
% VCD bookkeeping function to check columns in the time table for any
% inconsistencies and deviations from expectations.
%
%    results = vcd_checkTimeTable(time_table_master,data_results)
%
% INPUTS:
%
%  time_table_master    : (table) subject-specific time table master
%                         (N events x 48 columns) with complete information
%                         about all the sessions at the granularity of sub
%                         trial events. Table can be one used for a
%                         behavioral or MRI version of the VCD experiment.
%                         Rows represent the individual trial events,
%                         columns code different aspects of the experiment
%                         (e.g., condition_nr, block_nr, etc.)
% [data_results]        : (optional) (struct) outputs of vcdbehavioralanalysis,
%                         concatenated at the base, across all runs of a given session.
% [is_wide]             : (optional) (logical) flag to tell function this
%                          is the MRI wide session. If not defined, we
%                          assume this is either all deep MRI sessions, or 
%                          the behavioral experiment.
% OUTPUTS:
%  results              : (struct) summary table with a several statistics
%                         regarding the time_table_master. Note that this
%                         summary table is only a fraction of the checks we
%                         perform on the time_table_master.
%   struct fields:
%        total_nr_of_sessions: (integer)    number of sessions in the time table file
%            total_nr_of_runs: (integer)    number of runs across all sessions
%          total_nr_of_blocks: (integer)    number of blocks across all sessions
%          total_nr_of_trials: (integer)    number of trials across all sessions
%                crossing_nrs: [N×1 double] crossing numbers across all sessions
%              crossing_names: {N×1 cell}   crossing names across all sessions
%        trials_per_taskclass: [tasksclasses x 1 double] number of trials for each of the 10 task classes
%              trials_per_run: [sessions x runs x 1 double] number of trials per run, for each session
%              cues_per_block: [sessions × runs × blocks × 3 double] count of is_cued for each block.
%                              For the 4th dimension, order is cue=1 (left), cue=2 (right), cue=3 (both),
%     nr_of_unique_conditions: [sessions x runs x 1 double] number of unique condition numbers for each run, for each session
%              blocks_per_run: [sessions x runs x 1 double] number of blocks per run, for each session
%               condition_nrs: {[L×1 double]  [R×1 double]} condition
%                               numbers for left/central (Lx1) and right
%                               (Rx1) condition_nr columns in the time table master
%             condition_names: {{L×1 cell}  {R×1 cell}} condition
%                               names for left/central (Lx1) and right
%                               (Rx1) condition_names columns in the time table master
%                  cued_count: [1x3 double] total number of counted is_cue
%                               values across all sessions. Order is cue=1 (left), cue=2 (right), cue=3 (both)
%
% We want to know checks like:
% -What is the total number of trials per run?
% -What is the total number of trials and blocks per stimulus class
% -What is the total number of trials and blocks per task class?
% -What is the total number of unique crossing across a session?
% -What is the total number of unique stimuli shown across a session?
% -How many distinct conditions do we have? (e.g., GBR-0001-L-CUED-PC, etc.),
% -How many repetitions do we have per distinct condition?
% -Does each stimulus class have at least one fixation block?
% -How many left/right/neutral spatial cues do we have within a block?
% -Are the spatial cue L/R balanced across blocks/runs/entire session for each crossing involving classic stimuli?
% -How many contrast levels, rdk coherence levels, dot angles, object & ns super/basic/sub categories do we sample across a session/run/block?
% -Do we sample all three coherence levels in an RDK block?
% -Do we sample all three contrast levels within a GBR block?
% -Do we sample all superordinate categories across OBJ and/or NS blocks?
% -Are button presses balanced for each crossing across the session?
%
% Is it true that:
% - we do not repeat any stimuli within a block?
% - 20% of CD trials are CD+?
% - 20% of OBJ-PC trials are PC+?
% - the middle letter in the condition_name matches is_cued?
% - the number in the condition_name matches stimulus_nr_left/right?
% - the left column of condition_names only contains 'L' and 'C' (or 'X') stimulus positions?
% - the right column of condition_names only contains 'R' (or 'X') stimulus positions?
%
% Does the table show the following yoked columns?
% * Session_level (No NaN’s)
%     * session_nr; session_type; run_nr; global_run_nr;
%     * event_start; event_dur; event_end; event_id; event_name;
% * Block level (NaN’s for IBIs/pre-/post-blank periods)
%     * block_nr; global_block_nr
% * Semiblock level (NaN’s for ITIs/IBIs/pre-/post-blank periods)
%     * crossing_nr; crossing_name; nr_fix_changes
% * Trial level (NaN’s for ITIs)
%     * trial_nr; global_trial_nr
%     * stim_class; stim_class_name
%     * task_class; task_class_name
%     * is_cued; trial_type
% * Condition_level:
%     * condition_nr; condition_name
%     * stim_nr_left; stim_nr_right
%     * correct_response; cd_start
%     * orient_dir; contrast; gbr_phase; rdk_coherence
%     * super_cat; super_cat_name
%     * basic_cat; basic_cat_name
%     * sub_cat; sub_cat_name
%     * affordance_cat; affordance_name
%     * stim2_im_nr; stim2_delta; stim2_orient_dir
%     * repeat_nr
%     * is_objectcatch; is_catch; is_special_core; is_lure
%
% Example:
% {
% subj_nr     = 999;
% ses_nr      = 1;
% data_dir    = fullfile(vcd_rootPath, 'data','BEHAVIOR',sprintf('vcd_subj%03d',subj_nr));
% subj_folder = sprintf('vcd_subj%03d_ses%02d',subj_nr,ses_nr);
% tt_file     = dir(fullfile(data_dir, sprintf('vcd_subj%03d_time_table_master_PPROOM_EIZOFLEXSCAN_*.mat',subj_nr)));
% for rr = 1:12
%     mat_file  = sprintf('behavior_*_vcd_subj%03d_ses%02d_A_run%02d.mat',subj_nr,ses_nr,rr);
%     dd        = dir(fullfile(vcd_rootPath,'data','BEHAVIOR',subj_folder,mat_file));
%     load(fullfile(dd(end).folder,dd(end).name));
%     performance = vcdbehavioralanalysis(fullfile(params.savedatafolder,params.behaviorfile));
%     behresults(rr) = performance;
% end
% load(fullfile(tt_file(end).folder, tt_file(end).name))
% vcd_checkTimeTable(time_table_master,behresults)
% }
%
% Written by E Kupers @ UMN 2025/06

% Check inputs
if nargin==1 && ~exist('behresults','var')
    behresults = [];
elseif nargin==2 && ~isempty(varargin{1})
    behresults = varargin{1};
elseif nargin==3 && ~isempty(varargin{2})
    is_wide = varargin{2};
end

if ~exist('is_wide','var')
    is_wide = false;
end

% Get params
params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);
params.exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);

% Get all condition names
all_cond_names = vcd_getConditionNames;

% %%%% Bird's eye view stats %%%%
results = struct();

% record results
results.total_nr_of_sessions     = unique(time_table_master.session_nr);
results.total_nr_of_sessiontypes = unique(time_table_master.session_type);
results.total_nr_of_runs         = max(time_table_master.global_run_nr);
results.total_nr_of_blocks       = max(time_table_master.global_block_nr);
results.total_nr_of_trials       = max(time_table_master.global_trial_nr);

% check if nr of sessions is as we expect
assert(isequal(results.total_nr_of_sessions,[1:max(results.total_nr_of_sessions)]'))

% if we only have one session nr, tell the user
if length(results.total_nr_of_sessions)==1
    if is_wide
        fprintf('[%s]: Only found 1 session number, assumes this is the wide MRI experiment.\n',mfilename)
        assert(isequal(length(results.total_nr_of_sessions), params.exp.session.n_wide_sessions));
        env_type = 'MRI';
    else
        fprintf('[%s]: Only found 1 session number, assumes this is the behavioral experiment.\n',mfilename)
        assert(isequal(length(results.total_nr_of_sessions), params.exp.session.n_behavioral_sessions));
        env_type = 'BEHAVIOR';
    end
elseif length(results.total_nr_of_sessions)>1
    fprintf('[%s]: Only found 1 session number, assuming this is the MRI experiment.\n',mfilename)
    assert(isequal(length(results.total_nr_of_sessions), params.exp.session.n_deep_sessions));
    env_type = 'MRI';
else
    error('[%s]: No session number(s) found?!',mfilename)
end

% check if global_run_nr ascends as we expect (not skipping any run nr)
assert(isequal(unique(time_table_master.global_run_nr),[1:results.total_nr_of_runs]'))

% check if global_block_nr ascends as we expect (not skipping any block nr)
assert(isequal(unique(time_table_master.global_block_nr(~isnan(time_table_master.global_block_nr))),[1:results.total_nr_of_blocks]'))

% check if global_trial_nr ascends as we expect (not skipping any trial nr)
assert(isequal(unique(time_table_master.global_trial_nr(~isnan(time_table_master.global_trial_nr))),[1:results.total_nr_of_trials]'))

% check if task classes are as we expect
assert(all(ismember(unique(time_table_master.task_class(~isnan(time_table_master.task_class))),1:length(params.exp.taskclassnames))));
assert(all(ismember(unique(time_table_master.task_class_name(~cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.task_class_name))),params.exp.taskclassnames)));

% check if stim classes are as we expect
assert(all(ismember(unique(time_table_master.stim_class(~isnan(time_table_master.stim_class))),[1:length(params.exp.stimclassnames),99])));
assert(all(ismember(unique(time_table_master.stim_class_name(~cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.stim_class_name))),params.exp.stimclassnames)));

% Check if crossing nrs match with parameter struct
if strcmp(env_type,'BEHAVIOR')
    ses_blocks = params.exp.session.behavior.ses_blocks;
elseif strcmp(env_type,'MRI')
    if is_wide
        ses_blocks = params.exp.session.mri.wide.ses_blocks;
    else
        ses_blocks = cat(2,params.exp.session.mri.wide.ses_blocks,params.exp.session.mri.deep.ses_blocks); % dim 1:stim, 2:task, 3:session, 4:session type (a/b)
    end
end

% Compare crossing names
expected_crossing_names = {};
for kk = 1:size(ses_blocks,3) % session nr
    for ll = 1:size(ses_blocks,4) % session type
        for mm = 1:size(ses_blocks,1) % stimulus classes
            tc_names = params.exp.taskclassnames(ses_blocks(mm,:,kk,ll)>0);
            for nn = 1:length(tc_names) %  task classes crossed with given stimulus class
                if strcmp(tc_names{nn},'scc') || strcmp(tc_names{nn},'ltm')
                    expected_crossing_names = cat(1,expected_crossing_names, sprintf('%s-%s',tc_names{nn},'all'));
                else
                    expected_crossing_names = cat(1,expected_crossing_names, sprintf('%s-%s',tc_names{nn},params.exp.stimclassnames{mm}));
                end
            end
        end
    end
end

expected_crossing_names = unique(expected_crossing_names);
expected_crossing_nrs = find(ismember(params.exp.crossingnames,expected_crossing_names'));

assert(isequal(expected_crossing_names, unique(time_table_master.crossing_name(~cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.crossing_name)))))
assert(isequal(expected_crossing_nrs, setdiff(unique(time_table_master.crossing_nr(~isnan(time_table_master.crossing_nr))),999)));

% record
results.crossing_nrs   = expected_crossing_nrs;
results.crossing_names = expected_crossing_names;

% find eyetracking block
eye_tracking_block = find(time_table_master.block_nr==999);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check where the NaNs are to see if columns are yoked as expected

% Session level columns should have no NaNs:
% session nr, session_type, run_nr, global_run_nr, event_start, event_dur,
% event_end, event_id, event_name
assert(isempty(find(isnan(time_table_master.session_nr)))); %#ok<*EFIND>
assert(isempty(find(isnan(time_table_master.session_type))));
assert(isempty(find(isnan(time_table_master.run_nr))));
assert(isempty(find(isnan(time_table_master.global_run_nr))));
assert(isempty(find(isnan(time_table_master.event_start))));
assert(isempty(find(isnan(time_table_master.event_dur))));
assert(isempty(find(isnan(time_table_master.event_end))));
assert(isempty(find(isnan(time_table_master.event_id))));
assert(isempty(find(cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.event_name)))); %#ok<*DISEQN>

% Block level columns should have NaNs during IBI, pre-/post-blank.
% Crossing nr/name also have NaNs during ITI and eyetracking block.
block_nan           = find(isnan(time_table_master.block_nr));
global_block_nr_nan = find(isnan(time_table_master.global_block_nr));
crossing_nr_nan     = find(isnan(time_table_master.crossing_nr)); % Crossing numbers are NaN during ITIs
crossing_name_nan   = find(cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.crossing_name));
iti_nan             = find(ismember(time_table_master.event_id,98));

assert(isequal(block_nan,setdiff(global_block_nr_nan,eye_tracking_block)))
assert(isequal(block_nan,setdiff(crossing_nr_nan,[iti_nan;eye_tracking_block]))); % exclude ITIs and eyetracking block from counted nans in crossing_nr to match block nans
assert(isequal(block_nan,setdiff(crossing_name_nan,[iti_nan;eye_tracking_block])))% exclude ITIs and eyetracking block from counted nans crossing_name to match block nans
assert(all(ismember(time_table_master.event_id(block_nan),[0,99]))); % 0 = prepost run, 99 = ibi
assert(all(ismember(time_table_master.event_name(block_nan),{'pre-blank','post-blank','IBI'}))); % 0 = prepost run, 99 = ibi

% The row before a fixation block starts or ends should be either labeled
% as IBI (99) or pre-/post-blank (0)
nr_fix_changes_nan  = find(~isnan(time_table_master.nr_fix_changes));
assert(all(ismember(time_table_master.event_id([nr_fix_changes_nan(1)-1, nr_fix_changes_nan(end)+1]),[0,99]))); % 0 = prepost run, 99 = ibi

% Trial level columns should have NaNs during ITI, IBI, pre-/post-blank, and
% eyetracking block.
trial_nan           = find(isnan(time_table_master.trial_nr));
global_trial_nan    = find(isnan(time_table_master.global_trial_nr));
stim_class_nan      = find(isnan(time_table_master.stim_class));
stim_class_name_nan = find(cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.stim_class_name(:,1)));
task_class_nan      = find(isnan(time_table_master.task_class));
task_class_name_nan = find(cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.task_class_name));
trial_type_nan      = find(isnan(time_table_master.trial_type));

assert(isequal(trial_nan,global_trial_nan))
assert(isequal(trial_nan,stim_class_nan))
assert(isequal(trial_nan,task_class_nan))
assert(isequal(trial_nan,trial_type_nan))
assert(isequal(trial_nan,stim_class_name_nan))
assert(isequal(trial_nan,task_class_name_nan))

assert(all(ismember(time_table_master.event_id(trial_nan),[0,99,98, 990:998]))) % 0 = prepost run, 99 = ibi, 98=iti, 990:998=eyetracking block
assert(all(ismember(time_table_master.event_name(trial_nan),{'pre-blank','post-blank','IBI','ITI','et_sac','et_fix','et_pupil_black','et_pupil_white'})));

% eyetracking block is followed by pre-blank
pre_post_run = find(diff(eye_tracking_block)>1);
assert(all(ismember(time_table_master.event_name(eye_tracking_block(pre_post_run)+1),{'pre-blank'})))

% post-blank is followed by 10 eyetracking trials
assert(all(ismember(time_table_master.event_name(eye_tracking_block(pre_post_run(2:end))-10),{'post-blank'})))

% last entry is post-blank
assert(all(ismember(time_table_master.event_name(end),{'post-blank'})))

% Condition level columns should have NaNs during all events except for
% stim1/stim2
stim_nr_left_nan        = find(isnan(time_table_master.stim_nr_left));
stim_nr_right_nan       = find(isnan(time_table_master.stim_nr_right));
cd_start_not_nan        = find(~isnan(time_table_master.cd_start));
is_catch_nan            = find(isnan(time_table_master.is_catch));

% Catch check
assert(isequal(is_catch_nan,stim_nr_left_nan));

% Object catch checks
is_objectcatch_nan      = find(isnan(time_table_master.is_objectcatch));
obj_catch_condition_nrs = (time_table_master.condition_nr(time_table_master.crossing_nr==find(ismember(params.exp.crossingnames,'pc-obj'))));
assert(isequal(length(obj_catch_condition_nrs(~isnan(obj_catch_condition_nrs))), sum(~isnan(time_table_master.is_objectcatch))))
assert(isequal( sum(isnan(obj_catch_condition_nrs)) + sum(isnan(time_table_master.is_objectcatch(time_table_master.crossing_nr==find(ismember(params.exp.crossingnames,'pc-obj'))==0))), ...
    length(is_objectcatch_nan)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trials_per_run       = zeros(results.total_nr_of_sessions,length(results.total_nr_of_sessiontypes),max(time_table_master.run_nr));
trials_per_taskclass = zeros(results.total_nr_of_sessions,length(results.total_nr_of_sessiontypes),max(time_table_master.run_nr));
blocks_per_run       = zeros(results.total_nr_of_sessions,length(results.total_nr_of_sessiontypes),max(time_table_master.run_nr));
trials_per_block     = zeros(results.total_nr_of_sessions,length(results.total_nr_of_sessiontypes),max(time_table_master.run_nr), 7);
nr_of_unique_conditions = zeros(results.total_nr_of_sessions,length(results.total_nr_of_sessiontypes),max(time_table_master.run_nr));
cues_per_block       = zeros(results.total_nr_of_sessions,length(results.total_nr_of_sessiontypes),max(time_table_master.run_nr), 7, 3); % assume max 7 blocks and 3 cuing options (left/right/neutral)

for ses = 1:results.total_nr_of_sessions
    for st = 1:results.total_nr_of_sessiontypes % session type
        curr_ses_table = time_table_master(time_table_master.session_nr == ses & time_table_master.session_type == st,:);
        
        % Stats per run:
        [~,run_start] = unique(curr_ses_table.run_nr);
        run_start = cat(1,run_start,size(curr_ses_table,1));
        
        for rr = 1:length(run_start)-1
            tmp = curr_ses_table.trial_nr(run_start(rr):(run_start(rr+1)-1));
            trials_per_run(ses,st,rr) = length(tmp(~isnan(tmp)));
            nr_of_unique_conditions(ses,st,rr) = length(unique(curr_ses_table.condition_nr(curr_ses_table.event_id==94 & ~isnan(curr_ses_table.condition_nr) & curr_ses_table.run_nr==rr)));
            
            block_nrs = unique(curr_ses_table.block_nr(curr_ses_table.run_nr==rr & curr_ses_table.event_id==94));
            blocks_per_run(ses,rr) = length(block_nrs);
            for bb = 1:length(block_nrs)
                n0 = histcounts(curr_ses_table.is_cued(curr_ses_table.run_nr == rr & curr_ses_table.block_nr==block_nrs(bb) & curr_ses_table.event_id==94),[1:4]);
                cues_per_block(ses,rr,bb,1:length(n0)) = n0;
                trial_nrs = curr_ses_table.trial_nr(curr_ses_table.run_nr == rr & curr_ses_table.block_nr==block_nrs(bb) & curr_ses_table.event_id==94);
                trials_per_block(ses,st,rr,bb) = length(trial_nrs);
                assert(isequal([1:length(trial_nrs)]',trial_nrs)); % increment count
                assert(any(ismember(length(trial_nrs),[4,8]))); % assert only 4 and 8 trials per block
                stim_left_nrs = curr_ses_table.stim_nr_left(curr_ses_table.run_nr == rr & curr_ses_table.block_nr==block_nrs(bb) & curr_ses_table.event_id==94);
                stim_right_nrs = curr_ses_table.stim_nr_right(curr_ses_table.run_nr == rr & curr_ses_table.block_nr==block_nrs(bb) & curr_ses_table.event_id==94);
                assert(isempty(find(diff(stim_left_nrs)==0))) % no stim nr repeat
                assert(isempty(find(diff(stim_right_nrs)==0))) % no stim nr repeat
            end
        end
        
        for rr = 1:length(run_start)-1
            tmp = curr_ses_table.task_class(run_start(rr):(run_start(rr+1)-1));
            trials_per_taskclass(ses,st,rr) = length(tmp(~isnan(tmp)));
        end
        
        assert(isequal(unique(trials_per_block(ses,st,rr,:)),[0,4,8]')); % we expect either no, 4 or 8 blocks per trial.
        assert(all(mod(squeeze(trials_per_taskclass(ses,st,:)),2)==0)); % we expect even nr of trials per taskclass
        assert(all(all(mod(squeeze(trials_per_block(ses,st,rr,:)),2)==0))); % we expect even nr of trials per block
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Condition specific stats
        
        all_repeated_im_nrs_ns = cell(length(results.total_nr_of_sessiontypes),10); % 10 runs
        
        
        
        %%% Loop over left/right columns
        for side = [1,2]
            
            % check conditon nrs and names
            condition_nan                   = find(isnan(curr_ses_table.condition_nr(:,side)));
            results.condition_nrs{ses,st,side}   = curr_ses_table.condition_nr(curr_ses_table.event_id==94 & ~isnan(curr_ses_table.condition_nr(:,side)),side);
            results.condition_names{ses,st,side} = curr_ses_table.condition_name(curr_ses_table.event_id==94 & ~isnan(curr_ses_table.condition_nr(:,side)),side);
            
            % we expect conditon nrs and names columns to be yoked in terms of NaNs
            assert(isequal(isnan(curr_ses_table.condition_nr(:,side)), cellfun(@(x) isequalwithequalnans(x,NaN), curr_ses_table.condition_name(:,side))))
            
            % we expect conditon nrs to be yoked in terms of NaNs for left/right
            % stimulus numbers
            if side==1
                assert(isequal(condition_nan,find(isnan(curr_ses_table.stim_nr_left))))
                filter_me = condition_nan;
                
            elseif side==2
                assert(isequal(condition_nan,find(isnan(curr_ses_table.stim_nr_right))))
                filter_me = ismember(curr_ses_table.stim_class,[1:5,99]);
            end
            
            % we should expect the same nr of NaNs as for contrast, repeat
            % nr and special core column
            assert(isequal(condition_nan, find(isnan(curr_ses_table.contrast(:,side)))));
            assert(isequal(condition_nan, find(isnan(curr_ses_table.repeat_nr(:,side)))));
            assert(isequal(condition_nan, find(isnan(curr_ses_table.is_special_core(:,side)))));
            
            % For WM cases, we check if this is true for those task-crossings only
            assert(isequal(isnan(curr_ses_table.condition_nr(curr_ses_table.task_class==5,side)), isnan(curr_ses_table.stim2_im_nr(curr_ses_table.task_class==5,side))))
            assert(isequal(isnan(curr_ses_table.condition_nr(curr_ses_table.task_class==5,side)), isnan(curr_ses_table.stim2_delta(curr_ses_table.task_class==5,side))))
            assert(isequal(isnan(curr_ses_table.condition_nr(curr_ses_table.task_class==5 & curr_ses_table.stim_class~=5,side)), ...
                isnan(curr_ses_table.stim2_orient_dir(curr_ses_table.task_class==5 & curr_ses_table.stim_class~=5,side)))); % exclude NS stim since they have no orientation/motion direction/angle change
            
            for sc = [1:length(params.exp.stimclassnames), 99] % 99 is 'all', for scc/ltm
                
                switch sc
                    case 1 % Gabor
                        gbr_phase_nan  = find(isnan(curr_ses_table.gbr_phase(strcmp(curr_ses_table.stim_class_name(:,side),'gabor'),side)));
                        assert(isequal(gbr_phase_nan, find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'gabor'),side)))));
                        
                        gbr_ori       = curr_ses_table.orient_dir(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'gabor'),side);
                        gbr_phase     = curr_ses_table.gbr_phase(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'gabor'),side);
                        gbr_contrast  = curr_ses_table.contrast(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'gabor'),side);
                        wm_test_im    = curr_ses_table.stim2_orient_dir(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'gabor'),side);
                        ref_im        = curr_ses_table.orient_dir(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'gabor'),side);
                        delta_im      = curr_ses_table.stim2_delta(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'gabor'),side);
                        
                        if side == 1
                            assert(all(ismember(gbr_ori,params.stim.gabor.ori_deg(1:2:end))))
                            assert(all(ismember(gbr_phase,params.stim.gabor.ph_deg(1:2:end))))
                            assert(isequal(unique(gbr_ori)',params.stim.gabor.ori_deg(1:2:end)))
                            assert(isequal(unique(gbr_phase)',params.stim.gabor.ph_deg(1:2:end)))
                        elseif side == 2
                            assert(all(ismember(gbr_ori,params.stim.gabor.ori_deg(2:2:end))))
                            assert(all(ismember(gbr_phase,params.stim.gabor.ph_deg(2:2:end))))
                            assert(isequal(unique(gbr_ori)',params.stim.gabor.ori_deg(2:2:end)))
                            assert(isequal(unique(gbr_phase)',params.stim.gabor.ph_deg(2:2:end)))
                        end
                        
                        assert(all(ismember(gbr_contrast, params.stim.gabor.contrast)))
                        assert(isequal(unique(gbr_contrast)',params.stim.gabor.contrast))
                        assert(all(ref_im~=wm_test_im))
                        assert(isequal(circulardiff(wm_test_im,ref_im, 360),delta_im))
                        
                    case 2 % RDK
                        rdk_coherence_nan   = find(isnan(curr_ses_table.rdk_coherence(strcmp(curr_ses_table.stim_class_name(:,side),'rdk'),side)));
                        assert(isequal(rdk_coherence_nan, find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'rdk'),side)))));
                        
                        rdk_motdir  = curr_ses_table.orient_dir(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'rdk'),side);
                        rdk_coh     = curr_ses_table.rdk_coherence(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'rdk'),side);
                        
                        wm_test_im  = curr_ses_table.stim2_orient_dir(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'rdk'),side);
                        ref_im      = curr_ses_table.orient_dir(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'rdk'),side);
                        delta_im    = curr_ses_table.stim2_delta(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'rdk'),side);
                        
                        if side == 1
                            assert(all(ismember(rdk_motdir,params.stim.rdk.dots_direction(1:2:end))));
                            assert(isequal(unique(rdk_motdir)',params.stim.rdk.dots_direction(1:2:end)));
                        elseif side == 2
                            assert(all(ismember(rdk_motdir,params.stim.rdk.dots_direction(2:2:end))));
                            assert(isequal(unique(rdk_motdir)',params.stim.rdk.dots_direction(2:2:end)));
                        end
                        assert(all(ismember(rdk_coh,params.stim.rdk.dots_coherence)))
                        assert(isequal(unique(rdk_coh)',params.stim.rdk.dots_coherence))
                        assert(all(ref_im~=wm_test_im))
                        assert(isequal(circulardiff(wm_test_im,ref_im, 360),delta_im))
                        
                    case 3 % dot
                        % there is no special dot column
                        
                        dot_ang    = curr_ses_table.orient_dir(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'dot'),side);
                        wm_test_im = curr_ses_table.stim2_orient_dir(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'dot'),side);
                        ref_im     = curr_ses_table.orient_dir(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'dot'),side);
                        delta_im   = curr_ses_table.stim2_delta(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'dot'),side);
                        
                        if side == 1
                            assert(all(ismember(dot_ang,params.stim.dot.ang_deg(1:8))))
                            assert(isequal(sort(unique(dot_ang),'descend')',params.stim.dot.ang_deg(1:8)))
                            
                            assert(all(ismember(delta_im,params.stim.dot.delta_from_ref)))
                            assert(isequal(unique(delta_im)',params.stim.dot.delta_from_ref))
                        elseif side == 2
                            assert(all(ismember(dot_ang,params.stim.dot.ang_deg(9:16))))
                            assert(isequal(sort(unique(dot_ang),'descend')',params.stim.dot.ang_deg(9:16)))
                            
                            assert(all(ismember(delta_im,params.stim.dot.delta_from_ref)))
                            assert(isequal(unique(delta_im)',params.stim.dot.delta_from_ref))
                        end
                        assert(all(ref_im~=wm_test_im))
                        diff_delta = circulardiff(wm_test_im,ref_im,360);
                        assert(isequal(diff_delta,delta_im))
                        
                        
                    case 4 % obj
                        super_cat_nan       = find(isnan(curr_ses_table.super_cat(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)));
                        basic_cat_nan       = find(isnan(curr_ses_table.basic_cat(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)));
                        sub_cat_nan         = find(isnan(curr_ses_table.sub_cat(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)));
                        affordance_cat_nan  = find(isnan(curr_ses_table.affordance_cat(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)));
                        super_cat_name_nan  = find(cellfun(@(x) isequalwithequalnans(x,NaN), curr_ses_table.super_cat_name(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)));
                        basic_cat_name_nan  = find(cellfun(@(x) isequalwithequalnans(x,NaN), curr_ses_table.basic_cat_name(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)));
                        sub_cat_name_nan    = find(cellfun(@(x) isequalwithequalnans(x,NaN), curr_ses_table.sub_cat_name(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)));
                        affordance_name_nan = find(cellfun(@(x) isequalwithequalnans(x,NaN), curr_ses_table.affordance_name(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)));
                        
                        assert(isequal(super_cat_nan,       find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)))));
                        assert(isequal(basic_cat_nan,       find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)))));
                        assert(isequal(sub_cat_nan,         find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)))));
                        assert(isequal(affordance_cat_nan,  find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)))));
                        assert(isequal(super_cat_name_nan,  find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)))));
                        assert(isequal(basic_cat_name_nan,  find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)))));
                        assert(isequal(sub_cat_name_nan,    find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)))));
                        assert(isequal(affordance_name_nan, find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side)))));
                        
                        obj_ori        = curr_ses_table.orient_dir(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side);
                        is_objectcatch = curr_ses_table.is_objectcatch(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'obj'));
                        is_cued        = curr_ses_table.is_cued(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'obj'));
                        is_objpc       = curr_ses_table.is_objectcatch(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'obj') & curr_ses_table.task_class==4);
                        
                        results.nr_of_objectcatch_trials(ses,st) = sum(is_objectcatch==1);
                        results.nr_of_objpc_trials(ses,st) = length(is_objpc);
                        
                        if side == 1
                            objc = (is_objectcatch==1 & is_cued==1);
                            if ~isempty(objc)
                                catch_ori = params.stim.obj.catch_rotation(1:2:end,:);
                                assert(all(ismember(obj_ori(objc),catch_ori(:)')));
                            end
                            assert(all(ismember(obj_ori(~objc),params.stim.obj.facing_dir_deg(1:2:end))))
                            assert(isequal(unique(obj_ori(~objc))',sort(params.stim.obj.facing_dir_deg(1:2:end))))
                        elseif side == 2
                            objc = (is_objectcatch==1 & is_cued==2);
                            if ~isempty(objc)
                                catch_ori = params.stim.obj.catch_rotation(2:2:end,:);
                                assert(all(ismember(obj_ori(objc),catch_ori(:)')));
                            end
                            assert(all(ismember(obj_ori(~objc)',sort(params.stim.obj.facing_dir_deg(2:2:end)))))
                            assert(isequal(unique(obj_ori(~objc))',sort(params.stim.obj.facing_dir_deg(2:2:end))))
                        end
                        % check that wm test images have an updated orientation
                        wm_test_im     = curr_ses_table.stim2_orient_dir(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side);
                        ref_im         = curr_ses_table.orient_dir(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side);
                        delta_im       = curr_ses_table.stim2_delta(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'obj'),side);
                        assert(all(ref_im~=wm_test_im))
                        assert(isequal(wm_test_im-ref_im,delta_im))
                        
                        
                    case 5 % ns
                        super_cat_nan       = find(isnan(curr_ses_table.super_cat(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)));
                        basic_cat_nan       = find(isnan(curr_ses_table.basic_cat(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)));
                        sub_cat_nan         = find(isnan(curr_ses_table.sub_cat(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)));
                        affordance_cat_nan  = find(isnan(curr_ses_table.affordance_cat(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)));
                        
                        super_cat_name_nan  = find(cellfun(@(x) isequalwithequalnans(x,NaN), curr_ses_table.super_cat_name(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)));
                        basic_cat_name_nan  = find(cellfun(@(x) isequalwithequalnans(x,NaN), curr_ses_table.basic_cat_name(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)));
                        sub_cat_name_nan    = find(cellfun(@(x) isequalwithequalnans(x,NaN), curr_ses_table.sub_cat_name(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)));
                        affordance_name_nan = find(cellfun(@(x) isequalwithequalnans(x,NaN), curr_ses_table.affordance_name(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)));
                        
                        assert(isequal(super_cat_nan,       find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)))));
                        assert(isequal(basic_cat_nan,       find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)))));
                        assert(isequal(sub_cat_nan,         find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)))));
                        assert(isequal(affordance_cat_nan,  find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)))));
                        assert(isequal(super_cat_name_nan,  find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)))));
                        assert(isequal(basic_cat_name_nan,  find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)))));
                        assert(isequal(sub_cat_name_nan,    find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)))));
                        assert(isequal(affordance_name_nan, find(isnan(curr_ses_table.condition_nr(strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side)))));
                        
                        
                        super_cat       = curr_ses_table.super_cat(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side);
                        basic_cat       = curr_ses_table.basic_cat(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side);
                        sub_cat         = curr_ses_table.sub_cat(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side);
                        affordance_cat  = curr_ses_table.affordance_cat(curr_ses_table.event_id==94 & strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side);
                        
                        delta_im      = curr_ses_table.stim2_delta(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side);
                        wm_test_im      = curr_ses_table.stim2_orient_dir(curr_ses_table.event_id==95 & strcmp(curr_ses_table.stim_class_name(:,side),'ns'),side);
                        assert(all(isnan(wm_test_im)));
                        
                        if side == 1
                            assert(all(ismember(super_cat,1:length(params.stim.ns.super_cat))))
                            assert(isequal(unique(super_cat)',1:length(params.stim.ns.super_cat)))
                            
                            assert(all(ismember(basic_cat,1:length(params.stim.ns.basic_cat{1}))))
                            assert(isequal(unique(basic_cat)',1:length(params.stim.ns.basic_cat{1})))
                            
                            assert(all(ismember(sub_cat,1:length(params.stim.ns.sub_cat{1}))))
                            assert(isequal(unique(sub_cat)',1:length(params.stim.ns.sub_cat{1})))
                            
                            assert(all(ismember(affordance_cat,1:4)))
                            assert(isequal(unique(affordance_cat)',1:4))
                            
                            assert(all(ismember(delta_im,params.stim.ns.change_im)))
                            assert(isequal(unique(delta_im)',params.stim.ns.change_im))
                            
                        elseif side == 2
                            assert(isequalwithequalnans(super_cat, NaN(length(super_cat),1)))
                            assert(isequalwithequalnans(basic_cat, NaN(length(basic_cat),1)))
                            assert(isequalwithequalnans(sub_cat, NaN(length(sub_cat),1)))
                            assert(isequalwithequalnans(affordance_cat, NaN(length(affordance_cat),1)))
                            assert(isequalwithequalnans(delta_im, NaN(length(delta_im),1)))
                        end
                end
            end
            
            
            for tc = 1:length(params.exp.taskclassnames)
                
                corr_rsp = curr_ses_table.correct_response(curr_ses_table.task_class==tc);
                cued_side = 1+mod(curr_ses_table.is_cued(curr_ses_table.task_class==tc & curr_ses_table.event_id==94)-1,2);
                
                switch tc
                    
                    case 1 % fix
                        block_offset_times = [];
                        assert(isequalwithequalnans(corr_rsp,NaN(length(corr_rsp),1)));
                        block_onset_idx     = (ismember(curr_ses_table.event_id,91) &  ~isnan(curr_ses_table.nr_fix_changes));
                        block_onset_times   = curr_ses_table.event_start(block_onset_idx);
                        block_nr            = unique(curr_ses_table.global_block_nr(block_onset_idx));
                        for bb = 1:length(block_nr)
                            block_offset_times(bb) = max(curr_ses_table.event_end(find(ismember(curr_ses_table.global_block_nr, block_nr(bb)))));
                        end
                        expected_presses = sum(floor((block_offset_times'-block_onset_times) ./params.stim.fix.dotmeanchange));
                        assert(isequal(expected_presses,sum(curr_ses_table.nr_fix_changes,'omitnan')));
                    case 2 % cd
                        nr_trials = params.exp.block.n_trials_single_epoch * length(unique(curr_ses_table.global_block_nr(curr_ses_table.task_class==tc)));
                        nr_expected_changes = round(nr_trials*params.exp.trial.cd.prob_change);
                        if strcmp(env_type,'BEHAVIOR')
                            assert(isequal(nr_expected_changes, sum(corr_rsp==1)));
                        elseif strcmp(env_type,'MRI')
                            assert(ismember(nr_expected_changes, sum(corr_rsp==1)+[-1,0,1])); % can be +/- 1 due to rounding error
                        end
                    case 3 % scc
                        tmp_name = curr_ses_table.stim_class_name(curr_ses_table.task_class==tc & curr_ses_table.event_id==94,:);
                        for pp = 1:length(cued_side)
                            tmp_name_cued(pp) = tmp_name(pp,cued_side(pp));
                        end
                        [~,sccnames] = ismember(tmp_name_cued,params.exp.stimclassnames([1,3,2,4])); % note that button presses for dot and rdk are flipped compared to stim class nr
                        assert(isequal(corr_rsp(~isnan(corr_rsp)), sccnames'));
                        n0 = histcounts(corr_rsp(~isnan(corr_rsp)),[1:5]);
                        assert( all(n0>=1)); % every stimulus class is pressed at least once
                        
                        % check difference between gabor/rdks and object/dot
                        % stimulus classes
                        count_diff = abs(diff(n0([1,3])))<=1 && abs(diff(n0([2,4])))<=1; % less than 2 counts difference in nr of button presses (1 for gabors/rdks and 1 for object/dot)
                        if strcmp(env_type,'MRI')
                            assert(count_diff==1)
                        elseif strcmp(env_type,'BEHAVIOR')
                            if (count_diff~=1)
                                warning('[%s]: SCC BEHAVIOR button press distribution is off by more than 1 when comparing gabors vs rdks AND objects vs dots',mfilename)
                            end
                        end
                        % check difference between each individual stimulus class
                        scc_ok = zeros(1,4);
                        for mm = 1:4
                            if all(abs(n0(mm)-n0(setdiff([1:4],mm)))<=1)
                                scc_ok(mm) = true;
                            else
                                scc_ok(mm) = false;
                            end
                        end
                        if strcmp(env_type,'MRI')
                            if is_wide
                                assert(sum(scc_ok)==0)
                            else
                                assert(sum(scc_ok)==length(scc_ok));
                            end
                        elseif strcmp(env_type,'BEHAVIOR')
                            if sum(scc_ok)~=length(scc_ok)
                                warning('[%s]: SCC BEHAVIOR button press distribution is off by more than 1 when comparing all four stimulus classes individually: gabors vs rdks vs objects vs dots!',mfilename)
                            end
                        end
                        
                    case 4 % pc
                        for stim = [1:4]
                            corr_rsp = curr_ses_table.correct_response(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim);
                            tmp_ori  = curr_ses_table.orient_dir(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim,:);
                            cued_side = 1+mod(curr_ses_table.is_cued(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim)-1,2);
                            cued_tmp_ori = [];
                            for cc = 1:length(cued_side)
                                cued_tmp_ori(cc) = tmp_ori(cc,cued_side(cc));
                            end
                            cued_tmp_ori=cued_tmp_ori';
                            if stim == 1
                                assert(isequal(sum(corr_rsp(~isnan(corr_rsp))==1),sum(ismember(cued_tmp_ori, params.stim.gabor.ori_deg(params.stim.gabor.ori_deg>45 & params.stim.gabor.ori_deg<135))))); % more horz
                                assert(isequal(sum(corr_rsp(~isnan(corr_rsp))==2),sum(ismember(cued_tmp_ori, params.stim.gabor.ori_deg(params.stim.gabor.ori_deg<45 | params.stim.gabor.ori_deg>135))))); % more vert
                            elseif stim == 2
                                more_horz = (params.stim.rdk.dots_direction>45 & params.stim.rdk.dots_direction<135) | (params.stim.rdk.dots_direction>225 & params.stim.rdk.dots_direction<315);
                                more_vert = (params.stim.rdk.dots_direction<45 | (params.stim.rdk.dots_direction>135 & params.stim.rdk.dots_direction<225) | params.stim.rdk.dots_direction>315);
                                assert(isequal(sum(corr_rsp(~isnan(corr_rsp))==1),sum(ismember(cued_tmp_ori, params.stim.rdk.dots_direction(more_horz))))); % more horz
                                assert(isequal(sum(corr_rsp(~isnan(corr_rsp))==2),sum(ismember(cued_tmp_ori, params.stim.rdk.dots_direction(more_vert))))); % more vert
                            elseif stim == 3
                                more_horz = (params.stim.dot.ang_deg>45 & params.stim.dot.ang_deg<135) | (params.stim.dot.ang_deg>225 & params.stim.dot.ang_deg<315);
                                more_vert = (params.stim.dot.ang_deg<45 | (params.stim.dot.ang_deg>135 & params.stim.dot.ang_deg<225) | params.stim.dot.ang_deg>315);
                                assert(isequal(sum(corr_rsp(~isnan(corr_rsp))==1),sum(ismember(cued_tmp_ori, params.stim.dot.ang_deg(more_horz))))); % more horz
                                assert(isequal(sum(corr_rsp(~isnan(corr_rsp))==2),sum(ismember(cued_tmp_ori, params.stim.dot.ang_deg(more_vert))))); % more vert
                            elseif stim == 4 % note that object catch trial use orientations that are different from core stimuli
                                assert(isequal(sum(corr_rsp(~isnan(corr_rsp))==1),sum(cued_tmp_ori > 45 & cued_tmp_ori < 135))); % forward
                                assert(isequal(sum(corr_rsp(~isnan(corr_rsp))==2),sum((cued_tmp_ori < 45 | cued_tmp_ori > 135)))); % sideways
                            end
                        end
                        
                        % NS indoor/outdoor im
                        corr_rsp = curr_ses_table.correct_response(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                        cued_tmp_ori = curr_ses_table.basic_cat(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5,1);
                        assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_ori, 1)))); % indoor
                        assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_ori, 2)))); % outdoor
                        
                        if strcmp(env_type, 'BEHAVIOR')
                            stim_nrs_pc = curr_ses_table.stim_nr_left(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                            basic_cat_inout  = curr_ses_table.basic_cat(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                            
                            nr_repeats = find(histcounts(stim_nrs_pc, [params.stim.ns.unique_im_nrs_core, params.stim.ns.unique_im_nrs_core(end)+1])>1);
                            repeated_im_nrs = params.stim.ns.unique_im_nrs_core(nr_repeats);
                            all_repeated_im_nrs_ns{tc} = repeated_im_nrs;
                            
                            for mm = 1:length(repeated_im_nrs)
                                repeated_im_nrs_idx = find(ismember(stim_nrs_pc,repeated_im_nrs(mm)));
                                repeated_im_nrs_resp(mm) = basic_cat_inout(repeated_im_nrs_idx(1));
                            end
                            
                            response_distr = histcounts(corr_rsp,[1:3]);
                            assert(all(abs(diff(response_distr-histcounts(repeated_im_nrs_resp,[1:3])))==0)); % equal nr of button presses (after account for uneven repeats
                        elseif strcmp(env_type, 'MRI')
                            if is_wide 
                                if st==1
                                    ns_pc = corr_rsp;
                                else
                                    corr_rsp = [ns_pc,corr_rsp];
                                    assert(diff(histcounts(corr_rsp,[1:3]))==0); % equal nr of button presses across two session types
                                end
                            else
                                assert(diff(histcounts(corr_rsp,[1:3]))==0); % equal nr of button presses
                            end
                        end
                        
                    case 5 % wm
                        for stim = [1:4]
                            tmp_ori1 = []; tmp_ori2 = [];
                            corr_rsp = curr_ses_table.correct_response(curr_ses_table.task_class==tc & curr_ses_table.event_id==95 & curr_ses_table.stim_class==stim);
                            cued_side = 1+mod(curr_ses_table.is_cued(curr_ses_table.task_class==tc & curr_ses_table.event_id==95 & curr_ses_table.stim_class==stim)-1,2);
                            tmp_ori1      = curr_ses_table.orient_dir(curr_ses_table.task_class==tc & curr_ses_table.event_id==95 & curr_ses_table.stim_class==stim,:);
                            tmp_ori2      = curr_ses_table.stim2_orient_dir(curr_ses_table.task_class==tc & curr_ses_table.event_id==95 & curr_ses_table.stim_class==stim,:);
                            tmp_delta     = curr_ses_table.stim2_delta(curr_ses_table.task_class==tc & curr_ses_table.event_id==95 & curr_ses_table.stim_class==stim,:);
                            cued_tmp_ori1 = []; cued_tmp_ori2 = []; cued_tmp_delta = [];
                            for cc = 1:length(cued_side)
                                cued_tmp_ori1(cc) = tmp_ori1(cc,cued_side(cc));
                                cued_tmp_ori2(cc) = tmp_ori2(cc,cued_side(cc));
                                cued_tmp_delta(cc) = tmp_delta(cc,cued_side(cc));
                            end
                            cued_tmp_ori1=cued_tmp_ori1';
                            cued_tmp_ori2=cued_tmp_ori2';
                            cued_tmp_delta=cued_tmp_delta';
                            
                            diff_ori = circulardiff(cued_tmp_ori2,cued_tmp_ori1,360);
                            assert(isequal(cued_tmp_delta,diff_ori));
                            
                            if stim == 1
                                assert(isequal(sum(corr_rsp==1),sum(diff_ori<0))); %  ccw
                                assert(isequal(sum(corr_rsp==2),sum(diff_ori>0))); %  cw
                            elseif stim == 2
                                assert(isequal(sum(corr_rsp==1),sum(diff_ori<0))); %  ccw
                                assert(isequal(sum(corr_rsp==2),sum(diff_ori>0))); %  cw
                            elseif stim == 3
                                assert(isequal(sum(corr_rsp==1),sum(diff_ori<0))); %  ccw
                                assert(isequal(sum(corr_rsp==2),sum(diff_ori>0))); %  cw
                            elseif stim == 4
                                assert(isequal(sum(corr_rsp==1),sum(diff_ori<0))); %  left
                                assert(isequal(sum(corr_rsp==2),sum(diff_ori>0))); %  right
                            end
                        end
                        
                        % NS change im
                        corr_rsp = curr_ses_table.correct_response(curr_ses_table.task_class==tc & curr_ses_table.event_id==95 & curr_ses_table.stim_class==5);
                        tmp_delta = curr_ses_table.stim2_delta(curr_ses_table.task_class==tc & curr_ses_table.event_id==95 & curr_ses_table.stim_class==5,1);
                        assert(isequal(sum(corr_rsp==1),sum(tmp_delta<0))); % remove
                        assert(isequal(sum(corr_rsp==2),sum(tmp_delta>0))); % add
                        
                        if strcmp(env_type, 'BEHAVIOR')
                            stim_nrs_wm = curr_ses_table.stim_nr_left(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                            change_cat_wm  = curr_ses_table.stim2_delta(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                            
                            nr_repeats = find(histcounts(stim_nrs_wm, [params.stim.ns.unique_im_nrs_core, params.stim.ns.unique_im_nrs_core(end)+1])>1);
                            repeated_im_nrs = params.stim.ns.unique_im_nrs_core(nr_repeats);
                            all_repeated_im_nrs_ns{tc} = repeated_im_nrs;
                            for mm = 1:length(repeated_im_nrs)
                                repeated_im_nrs_idx = find(ismember(stim_nrs_wm,repeated_im_nrs(mm)));
                                repeated_im_nrs_resp(mm) = change_cat_wm(repeated_im_nrs_idx(1));
                            end
                            response_distr = histcounts(corr_rsp,[1:3]);
                            if ~(diff(histcounts(response_distr,[1:3]))==0) % equal nr of button presses
                                assert(all(abs(diff(response_distr-histcounts(repeated_im_nrs_resp,[1:3])))==0)); % equal nr of button presses (after account for uneven repeats
                            end
                        elseif strcmp(env_type, 'MRI')
                            assert(diff(histcounts(corr_rsp,[1:3]))==0); % equal nr of button presses
                        end
                        
                    case 6 % ltm
                        % under construction
                    case 7 % img
                        % under construction
                    case 8 % what
                        
                        for stim = [4,5]
                            cued_tmp_cat = [];
                            corr_rsp = curr_ses_table.correct_response(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim);
                            if stim == 4
                                cued_side = 1+mod(curr_ses_table.is_cued(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim)-1,2);
                                tmp_cat = curr_ses_table.super_cat(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim,:);
                                for cc = 1:length(cued_side)
                                    cued_tmp_cat(cc) = tmp_cat(cc,cued_side(cc));
                                end
                                cued_tmp_cat=cued_tmp_cat';
                                
                            elseif stim == 5
                                cued_tmp_cat = curr_ses_table.super_cat(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim,1);
                                
                                if strcmp(env_type, 'BEHAVIOR')
                                    stim_nrs_what = curr_ses_table.stim_nr_left(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                                    nr_repeats = find(histcounts(stim_nrs_what, [params.stim.ns.unique_im_nrs_core, params.stim.ns.unique_im_nrs_core(end)+1])>1);
                                    repeated_im_nrs = params.stim.ns.unique_im_nrs_core(nr_repeats);
                                    all_repeated_im_nrs_ns{st,tc} = repeated_im_nrs;
                                end
                                
                            end
                            % distribution of button presses should match nr of supercategories
                            assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_cat,1)))); % human
                            assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_cat,2)))); % animal
                            assert(isequal(sum(corr_rsp==3),sum(ismember(cued_tmp_cat,[3,4])))); % object, food
                            assert(isequal(sum(corr_rsp==4),sum(ismember(cued_tmp_cat,5)))); % place
                        end
                        
                        % what task button presses merge objects and food
                        cued_tmp_cat2 = cued_tmp_cat;
                        cued_tmp_cat2(ismember(cued_tmp_cat,[3,4])) = 3;
                        cued_tmp_cat2(ismember(cued_tmp_cat,[5])) = 4;
                        
                        % distribution of button presses should match nr of supercategories
                        assert(isequal(histcounts(cued_tmp_cat2, [1:5]), histcounts(corr_rsp,[1:5]))); % distribution of button presses should match nr of supercategories
                        
                    case 9 % where
                        corr_rsp = curr_ses_table.correct_response(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                        cued_tmp_loc = curr_ses_table.sub_cat(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5,1);
                        assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_loc,1)))); % left
                        assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_loc,2)))); % center
                        assert(isequal(sum(corr_rsp==3),sum(ismember(cued_tmp_loc,3)))); % right
                        
                        if strcmp(env_type, 'BEHAVIOR')
                            stim_nrs_where = curr_ses_table.stim_nr_left(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                            sub_cat_where  = curr_ses_table.sub_cat(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                            
                            nr_repeats = find(histcounts(stim_nrs_where, [params.stim.ns.unique_im_nrs_core, params.stim.ns.unique_im_nrs_core(end)+1])>1);
                            repeated_im_nrs = params.stim.ns.unique_im_nrs_core(nr_repeats);
                            for mm = 1:length(repeated_im_nrs)
                                repeated_im_nrs_idx = find(ismember(stim_nrs_where,repeated_im_nrs(mm)));
                                repeated_im_nrs_resp(mm) = sub_cat_where(repeated_im_nrs_idx(1));
                            end
                            all_repeated_im_nrs_ns{st,tc} = repeated_im_nrs;
                            response_distr = histcounts(corr_rsp,[1:4]);
                            assert(all(abs(diff(response_distr-histcounts(repeated_im_nrs_resp,[1:4])))==0)); % equal nr of button presses (sfter account for uneven repeats
                        elseif sctrmp(env_type, 'MRI')
                            assert(diff(histcounts(corr_rsp,[1:4]))==0); % equal nr of button presses
                        end
                    case 10 % how
                        for stim = [4,5]
                            cued_tmp_afford = [];
                            corr_rsp = curr_ses_table.correct_response(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim);
                            if stim == 4
                                cued_side = 1+mod(curr_ses_table.is_cued(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim)-1,2);
                                tmp_cat = curr_ses_table.affordance_cat(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim,:);
                                for cc = 1:length(cued_side)
                                    cued_tmp_afford(cc) = tmp_cat(cc,cued_side(cc));
                                end
                                cued_tmp_afford=cued_tmp_afford';
                            elseif stim == 5
                                cued_tmp_afford = curr_ses_table.affordance_cat(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==stim,1);
                                if strcmp(env_type, 'BEHAVIOR')
                                    stim_nrs_how = curr_ses_table.stim_nr_left(curr_ses_table.task_class==tc & curr_ses_table.event_id==94 & curr_ses_table.stim_class==5);
                                    nr_repeats = find(histcounts(stim_nrs_how, [params.stim.ns.unique_im_nrs_core, params.stim.ns.unique_im_nrs_core(end)+1])>1);
                                    repeated_im_nrs = params.stim.ns.unique_im_nrs_core(nr_repeats);
                                    all_repeated_im_nrs_ns{st,tc} = repeated_im_nrs;
                                end
                            end
                            assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_afford,1)))); % greet
                            assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_afford,2)))); % grasp
                            assert(isequal(sum(corr_rsp==3),sum(ismember(cued_tmp_afford,3)))); % enter/walk
                            assert(isequal(sum(corr_rsp==4),sum(ismember(cued_tmp_afford,4)))); % observe/do nothing
                        end
                        assert(isequal(histcounts(cued_tmp_afford, [1:5]), histcounts(corr_rsp,[1:5]))); % distribution of button presses should match nr of supercategories
                end
            end
        end
        
        results.total_nr_of_conditions(st) = sum(length(results.condition_names{1})+length(results.condition_names{2}));
        
        correct_response_tally = zeros(st,length(params.exp.crossingnames),4);
        for cc = 1:length(params.exp.crossingnames)
            correct_response_tally(st,cc,:) = histcounts(curr_ses_table.correct_response(curr_ses_table.crossing_nr==cc), [1:5]);
        end
    end
end

results.trials_per_taskclass            = trials_per_taskclass;
results.trials_per_run                  = trials_per_run;
results.cues_per_block                  = cues_per_block;
results.nr_of_unique_conditions_per_run = nr_of_unique_conditions;
results.trials_per_block                = trials_per_block;
results.blocks_per_run                  = blocks_per_run;
results.total_nr_of_unique_conditions   = sum(nr_of_unique_conditions);
results.all_repeated_im_nrs_ns = all_repeated_im_nrs_ns;

% CD onset checks
total_cd_trials = length(find(time_table_master.task_class(time_table_master.event_id==94)==2));
cd_cond_names = time_table_master.condition_name(cd_start_not_nan,:);
cd_cond_nr = time_table_master.condition_nr(cd_start_not_nan,:);
cd_cond_names_no_nan = cell2mat(cellfun(@(x) isequalwithequalnans(x,NaN), cd_cond_names, 'UniformOutput',0));
side_cd_plus = []; cd_plus_cond_names = {};
for jj = 1:length(cd_start_not_nan)
    idx_cd_plus = find(~cellfun(@isempty, (regexp( cd_cond_names(jj,~cd_cond_names_no_nan(jj,:)), '+'))));
    assert(length(idx_cd_plus==1))
    side_cd_plus = cat(1,side_cd_plus,idx_cd_plus);
    cd_plus_cond_names(jj) = cd_cond_names(jj,idx_cd_plus);
    cd_plus_cond_nr(jj) = cd_cond_nr(jj,idx_cd_plus);
end
assert(isequal(all_cond_names(cd_plus_cond_nr),cd_plus_cond_names))

results.total_nr_of_CDplus_trials = length(cd_plus_cond_names);
results.total_nr_of_CD_trials = total_cd_trials;
% ltm stats
% assert(isequal(condition_nan, find(isnan(time_table_master.is_lure(:,side)))));

% cued stats
cue_ID        = time_table_master.is_cued(time_table_master.event_id==94);
is_NS         = ( time_table_master.stim_class(time_table_master.event_id==94)==5 );
is_FIX        = ( time_table_master.task_class(time_table_master.event_id==94)==1 );
is_NS_or_FIX  = ( time_table_master.task_class(time_table_master.event_id==94)==1 | time_table_master.stim_class(time_table_master.event_id==94)==5 );
is_SCC        = ( time_table_master.stim_class(time_table_master.event_id==94)==99 | time_table_master.task_class(time_table_master.event_id==94)==3 );
is_LTM        = ( time_table_master.stim_class(time_table_master.event_id==94)==99 | time_table_master.task_class(time_table_master.event_id==94)==6 );

c_names       = time_table_master.condition_name(time_table_master.event_id==94,:);
results.cued_count = histcounts(cue_ID,[1:4]);
assert(isequal(results.cued_count(1),results.cued_count(2))); % equal nr of left and right cues
assert(isequal(results.cued_count(3), sum(is_NS_or_FIX))); % equal nr of left and right cues

%%%%% Condition name checks:

loc1        = cellfun(@(x) strsplit(x,'-') ,c_names(:,1),'UniformOutput',0);
loc2        = cellfun(@(x) strsplit(x,'-') ,c_names(~isnan(time_table_master.condition_nr(time_table_master.event_id==94,2)),2),'UniformOutput',0);
is_FIX_lr   = cat(1,is_FIX,is_FIX(~isnan(time_table_master.condition_nr(time_table_master.event_id==94,2))));
is_NS_lr    = cat(1,is_NS,is_NS(~isnan(time_table_master.condition_nr(time_table_master.event_id==94,2))));
cue_ID_lr   = cat(1,cue_ID,cue_ID(~isnan(time_table_master.condition_nr(time_table_master.event_id==94,2))));
stimclass_name    = time_table_master.stim_class_name(time_table_master.event_id==94,:);
stimclass_name_lr = cat(1,stimclass_name(:,1), stimclass_name(~isnan(time_table_master.condition_nr(time_table_master.event_id==94,2)),2));
cond_names_lr     = cat(1,loc1,loc2);
cond_names_lr_idx = cat(1, ones(length(loc1),1), 2*ones(length(loc2),1));

extended_task_class_names = cat(2,upper(params.exp.taskclassnames),strcat(upper(params.exp.taskclassnames),'#'),{'PC+','CD+'});

for xx = 1:length(cond_names_lr)
    % get current condition name parts
    curr_stimclass   = cond_names_lr{xx}{1};
    curr_stimnr      = cond_names_lr{xx}{2};
    curr_loc         = cond_names_lr{xx}{3};
    curr_cued_status = cond_names_lr{xx}{4};
    curr_taskclass   = cond_names_lr{xx}{5};
    
    if strcmp(curr_taskclass,'PC+')
        curr_taskclass = 'PC';
    elseif  strcmp(curr_taskclass,'CD+')
        curr_taskclass = 'CD';
    end
    
    if cond_names_lr_idx(xx) == 1
        % check if left column location is as expected
        assert(ismember(curr_loc,{'L','C','X'}));
        
        % Middle letter should be L/C for first condition_name sub column, and cued
        % status should match is_cued column
        if cue_ID_lr(xx) == 1 && strcmp(curr_loc,'L')
            assert(strcmp(curr_cued_status,'CUED'))
        elseif cue_ID_lr(xx) == 2 && strcmp(curr_loc,'L')
            assert(strcmp(curr_cued_status,'UNCUED'))
        elseif cue_ID_lr(xx) == 3 && is_FIX_lr(xx)==1
            assert(strcmp(curr_cued_status,'NCUED'));
        elseif cue_ID_lr(xx) == 3 && is_NS_lr(xx)==1
            assert(strcmp(curr_cued_status,'NCUED'));
        elseif cue_ID_lr(xx) == 3 && strcmp(curr_loc,'C')
            assert(is_NS_lr(xx)==1)
            assert(strcmp(curr_cued_status,'NCUED'));
        elseif cue_ID_lr(xx) == 1 && strcmp(curr_loc,'X')
            assert(strcmp(curr_cued_status,'LCUE'));
        elseif cue_ID_lr(xx) == 2 && strcmp(curr_loc,'X')
            assert(strcmp(curr_cued_status,'RCUE'));
        elseif cue_ID_lr(xx) == 3 && strcmp(curr_loc,'X')
            assert(strcmp(curr_cued_status,'NCUE'));
        else
            error('[%s]: Condition name doesn''t match cued status of stimulus',mfilename);
        end
        
        
    elseif cond_names_lr_idx(xx) == 2
        % check if right column location is as expected
        assert(ismember(curr_loc,{'R','X'}));
        
        % Middle letter should be R for second condition_name sub column, and cued
        % status should match is_cued column
        if cue_ID_lr(xx) == 1 && strcmp(curr_loc,'R')
            assert(strcmp(curr_cued_status,'UNCUED'))
        elseif cue_ID_lr(xx) == 2 && strcmp(curr_loc,'R')
            assert(strcmp(curr_cued_status,'CUED'))
        elseif cue_ID_lr(xx) == 3 && is_FIX_lr(xx)==1
            assert(strcmp(curr_cued_status,'NCUED'));
        elseif cue_ID_lr(xx) == 1 && strcmp(curr_loc,'X')
            assert(strcmp(curr_cued_status,'LCUE'));
        elseif cue_ID_lr(xx) == 2 && strcmp(curr_loc,'X')
            assert(strcmp(curr_cued_status,'RCUE'));
        elseif cue_ID_lr(xx) == 3 && strcmp(curr_loc,'X')
            assert(strcmp(curr_cued_status,'NCUE'));
        else
            error('[%s]: Condition name doesn''t match cued status of stimulus.\n',mfilename);
        end
        
    else
        error('[%s]: Check indexing of condition names.\n',mfilename);
    end
    
    % Check stimulus class name and stimulus number correspondance
    if strcmp(curr_stimclass,'NS')
        assert(ismember(str2double(curr_stimnr),params.stim.ns.unique_im_nrs_core))
    elseif ismember(curr_stimclass,{'ALL','OBJ','DOT','RDK','GBR'})
        % we expect uneven stimulus numbers on the left, even stimulus
        % number on the right, and a stimulus nr of 0 for catch trials
        if strcmp(curr_loc,'L')
            if strcmp(stimclass_name_lr(xx),'dot') % dot stim have different order of stim nrs
                assert(ismember(str2double(curr_stimnr),params.stim.dot.unique_im_nrs_core(1:8)))
            else
                assert(mod(str2double(curr_stimnr),2)==1)
            end
        elseif strcmp(curr_loc,'R')
            if strcmp(stimclass_name_lr(xx),'dot') % dot stim have different order of stim nrs
                assert(ismember(str2double(curr_stimnr),params.stim.dot.unique_im_nrs_core(9:end)))
            else
                assert(mod(str2double(curr_stimnr),2)==0)
            end
        elseif strcmp(curr_loc,'X')
            assert(str2double(curr_stimnr)==0)
        else
            error('[%s]: Unrecognized stimulus location .\n',mfilename);
        end
        
        if strcmp(curr_stimclass,'OBJ')
            assert(ismember(str2double(curr_stimnr),params.stim.obj.unique_im_nrs_core))
        elseif strcmp(curr_stimclass,'DOT')
            assert(ismember(str2double(curr_stimnr),params.stim.dot.unique_im_nrs_core))
        elseif strcmp(curr_stimclass,'RDK')
            assert(ismember(str2double(curr_stimnr),params.stim.rdk.unique_im_nrs_core))
        elseif strcmp(curr_stimclass,'GBR')
            assert(ismember(str2double(curr_stimnr),params.stim.gabor.unique_im_nrs_core))
        end
    else
        error('[%s]: Condition name doesn''t match known stimulus classes.\n',mfilename);
    end
    % Check task class name
    assert(ismember(curr_taskclass,extended_task_class_names))
    
    if strcmp(curr_stimclass,'GBR')
        curr_stimclass = 'gabor';
    end
    % check if crossing exists in ses_blocks and in general
    if strcmp(curr_stimclass,'ALL')
        assert(ismember(curr_taskclass,{'SCC','LTM'}))
        assert(sum(reshape(ses_blocks(:,ismember(lower(curr_taskclass),params.exp.taskclassnames),:,:),1,[]))>0)
        if ismember(curr_taskclass,{'SCC'})
            assert(all(params.exp.crossings(1:4,strcmp(lower(curr_taskclass),params.exp.taskclassnames)))) %#ok<STCI>
        elseif ismember(curr_taskclass,{'LTM'})
            assert(all(params.exp.crossings(1:5,strcmp(lower(curr_taskclass),params.exp.taskclassnames)))) %#ok<STCI>
        end
    else
        assert(sum(reshape(ses_blocks(strcmp(lower(curr_stimclass),params.exp.stimclassnames),strcmp(lower(curr_taskclass),params.exp.taskclassnames),:,:),1,[]))>0) %#ok<STCI>
        assert(params.exp.crossings(strcmp(lower(curr_stimclass),params.exp.stimclassnames),strcmp(lower(curr_taskclass),params.exp.taskclassnames))==1) %#ok<STCI>
    end
end

% check if we equally cue left/right side within a classic stimulus block
equal_lr_cuing_count = []; neutral_cue_count = [];
for ses = 1:size(results.cues_per_block,1)
    for rr = 1:size(results.cues_per_block,2)
        for bb = 1:size(results.cues_per_block,3)
            if ~isempty(results.cues_per_block(ses,rr,bb,:))
                if results.cues_per_block(ses,rr,bb,1) ~= 0 && results.cues_per_block(ses,rr,bb,2) ~= 0
                    if isequal(results.cues_per_block(ses,rr,bb,1),results.cues_per_block(ses,rr,bb,2))
                        equal_lr_cuing_count = cat(1,equal_lr_cuing_count,1);
                    else
                        equal_lr_cuing_count = cat(1,equal_lr_cuing_count,0);
                    end
                elseif results.cues_per_block(ses,rr,bb,3)~=0
                    neutral_cue_count = cat(1,neutral_cue_count,1);
                end
            end
        end
    end
end
% nr of left and right cued stimuli should be equal per block
assert(isequal(sum(equal_lr_cuing_count),length(equal_lr_cuing_count)));
% sum of left/right and neutrally cued blocks should match the total nr of blocks
assert(isequal(sum([equal_lr_cuing_count;neutral_cue_count]),results.total_nr_of_blocks));
assert(isequal(sum(results.cues_per_block(:)),results.total_nr_of_trials))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we provided BEHAVIOR results, see if these results match with time_table_master
if ~isempty(behresults)
    fix_blocks = find(~cellfun(@isempty, regexp(params.exp.crossingnames,'fix')));
    
    sessions_behresults = unique(behresults(rr).trialinfo.session_nr);
    nr_runs_behresults = size(behresults,2);
    
    for ses = 1:length(sessions_behresults)
        curr_ses = sessions_behresults(ses);
        
        assert(isequal(results.total_nr_of_runs(ses,:),nr_runs_behresults));
        
        all_crossing_nr_beh = [];
        beh_runnumblocks  = zeros(results.total_nr_of_sessions, nr_runs_behresults,length(params.exp.crossingnames));
        tt_runnumblocks   = beh_runnumblocks;
        beh_runnumtrials  = beh_runnumblocks;
        tt_runnumtrials   = beh_runnumblocks;
        
        for rr = 1:nr_runs_behresults
            % find fixation blocks (they have pseudo trials)
            nr_fix_blocks = behresults(rr).summary.crossing_nr( ismember(behresults(rr).summary.crossing_nr,fix_blocks));
            fix_stim_task = params.exp.crossingnames(behresults(rr).summary.crossing_nr(ismember(behresults(rr).summary.crossing_nr,nr_fix_blocks)));
            
            % Get BEHAVIOR performance file condition nrs
            cond_nrs_beh = find(behresults(rr).conditioncount>0)';
            repeats = find(behresults(rr).conditioncount>1);
            if ~isempty(repeats)
                for cc = 1:length(repeats)
                    cond_nrs_beh = cat(1,cond_nrs_beh, repelem(repeats(cc), behresults(rr).conditioncount(repeats(cc))-1,1));
                end
            end
            cond_names_beh = all_cond_names(cond_nrs_beh)';
            % Get time table condition nrs
            cond_nrs_tt = time_table_master.condition_nr(time_table_master.session_nr == curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr,:);
            cond_nrs_tt = cond_nrs_tt(~isnan(cond_nrs_tt));
            % Get time table condition names
            cond_names_tt = time_table_master.condition_name(time_table_master.session_nr == curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr,:);
            cond_names_tt = cond_names_tt(:);
            cond_names_tt(cellfun(@(x) isequalwithequalnans(x,NaN), cond_names_tt)) = [];
            cond_nrs_tt_counts = histcounts(cond_nrs_tt,[1:length(all_cond_names)+1]);
            
            % Condition names
            assert(isequal(sort(cond_names_tt),sort(cond_names_beh)))
            
            % Condition nrs should match for every run
            assert(isequal(sort(cond_nrs_beh),sort(cond_nrs_tt)));
            
            % Condition nr tally should match for every run
            assert(isequal(cond_nrs_tt_counts,behresults(rr).conditioncount));
            
            % Nr of blocks for each crossing should match for every run
            beh_crossing_nrs           = cat(1,nr_fix_blocks, behresults(rr).trialinfo.crossing_nr(behresults(rr).trialinfo.trial_nr==1)); % we insert fixation blocks, because fix trial nrs are set to NaN
            beh_runnumblocks(ses,rr,:) = histcounts(beh_crossing_nrs,[1:(length(params.exp.crossingnames)+1)]);
            tt_runnumblocks(ses,rr,:)  = histcounts(time_table_master.crossing_nr(time_table_master.session_nr == curr_ses & time_table_master.trial_nr==1 & time_table_master.event_id==94 & time_table_master.run_nr==rr),[1:(length(params.exp.crossingnames)+1)]);
            assert(isequal(beh_runnumblocks,tt_runnumblocks));
            
            % Nr of trials for each crossing should match for every run
            beh_runnumtrials(ses,rr,:) = histcounts(behresults(rr).trialinfo.crossing_nr,[1:(length(params.exp.crossingnames)+1)]); % we insert fixation blocks, because fix trial nrs are set to NaN
            tt_runnumtrials(ses,rr,:) = histcounts(time_table_master.crossing_nr(time_table_master.session_nr == curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr),[1:(length(params.exp.crossingnames)+1)]);
            % insert nr of fixation trials
            for jj = 1:length(nr_fix_blocks)
                tmp_trials  = sum(isnan(behresults(rr).trialinfo.trial_nr(behresults(rr).trialinfo.crossing_nr==nr_fix_blocks(jj))));
                beh_runnumtrials(ses,rr,nr_fix_blocks(jj)) = tmp_trials;
                tt_runnumtrials(ses,rr,nr_fix_blocks(jj)) = sum(time_table_master.nr_fix_changes(time_table_master.session_nr == curr_ses & time_table_master.crossing_nr==nr_fix_blocks(jj) & time_table_master.run_nr==rr));
            end
            
            % Check if session nrs match
            % replace pseudo fixaton trials with NaN
            session_nr = behresults(rr).trialinfo.session_nr;
            run_nr = behresults(rr).trialinfo.run_nr;
            block_nr = behresults(rr).trialinfo.block_nr;
            crossing_nr = behresults(rr).trialinfo.crossing_nr;
            trial_nr = behresults(rr).trialinfo.trial_nr;
            is_catch = behresults(rr).trialinfo.is_catch;
            stim_cued = behresults(rr).trialinfo.stim_cued;
            if length(nr_fix_blocks)>0
                for ff = 1:length(nr_fix_blocks)
                    fix_stim_class = strsplit(fix_stim_task{ff},'-');
                    
                    fix_pseudo_trials = find(ismember(behresults(rr).trialinfo.crossing_nr,nr_fix_blocks(ff)));
                    session_nr(fix_pseudo_trials(9:end))= NaN;
                    run_nr(fix_pseudo_trials(9:end))= NaN;
                    block_nr(fix_pseudo_trials(9:end))= NaN;
                    crossing_nr(fix_pseudo_trials(9:end))= NaN;
                    trial_nr(fix_pseudo_trials(1:8)) = [1:8];
                    trial_nr(fix_pseudo_trials(9:end))= NaN;
                    is_catch(fix_pseudo_trials(9:end))= NaN;
                    stim_cued(fix_pseudo_trials(1:8)) = repmat(fix_stim_class(2),8,1);
                    stim_cued(fix_pseudo_trials(9:end))= {NaN};
                end
            end
            
            assert(isequal(session_nr(~isnan(session_nr)), time_table_master.session_nr(time_table_master.session_nr==curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
            
            % Check if run nrs match
            assert(isequal(run_nr(~isnan(run_nr)), time_table_master.run_nr(time_table_master.session_nr==curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
            
            % Check if block nrs match
            assert(isequal(block_nr(~isnan(block_nr)), time_table_master.block_nr(time_table_master.session_nr==curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
            
            % Check if crossing nrs match
            assert(isequal(crossing_nr(~isnan(crossing_nr)), time_table_master.crossing_nr(time_table_master.session_nr==curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
            
            % Check if trial nrs match
            assert(isequal(trial_nr(~isnan(trial_nr)), time_table_master.trial_nr(time_table_master.session_nr==curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
            
            % Check if catch trial nrs match
            assert(isequal(is_catch(~isnan(is_catch)), time_table_master.is_catch(time_table_master.session_nr==curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
            
            % Check if stim cued match
            stimclass = time_table_master.stim_class_name(time_table_master.session_nr==curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr,:);
            cued_side = 1+mod(time_table_master.is_cued(time_table_master.session_nr==curr_ses & time_table_master.event_id==94 & time_table_master.run_nr==rr)-1, 2);
            stimclass_cued = cell(size(stimclass,1),1);
            stimclass_cued(cued_side==1) = stimclass(cued_side==1,1);
            stimclass_cued(cued_side==2) = stimclass(cued_side==2,2);
            assert(isequal(stim_cued(cellfun(@isempty, cellfun(@(x) find(x==1), cellfun(@isnan, stim_cued,'UniformOutput',0),'UniformOutput',0))) ,stimclass_cued))
            
            all_crossing_nr_beh = cat(1,all_crossing_nr_beh, crossing_nr(~isnan(crossing_nr)));
        end
        assert(isequal(expected_crossing_nrs,unique(sort(all_crossing_nr_beh(:)))))
    end
    
    results.beh_crossing_count = beh_runnumblocks;
    results.tt_crossing_count = tt_runnumblocks;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRINT RESULTS %%%

fprintf('[%s]:Time_table_master statistics:\n',mfilename)
if strcmp(env_type,'BEHAVIOR')
    fprintf('*** Code assumes this is a behavioral session (!) ****\n');
end
fprintf('Total number of sessions: \t\t\t\t\t%d\n', results.total_nr_of_sessions);
fprintf('Total number of runs: \t\t\t\t\t\t%d\n', results.total_nr_of_runs);
fprintf('Total number of block: \t\t\t\t\t\t%d\n', results.total_nr_of_blocks);
fprintf('Total number of trials: \t\t\t\t\t%d\n', results.total_nr_of_trials);
fprintf('Total number of cues: \t\t\t\t\t\tleft: %d, right: %d, both: %d\n', results.cued_count(1),results.cued_count(2),results.cued_count(3));
fprintf('Total number of conditions (incl. repeats): \t\t\t%d \n', results.total_nr_of_conditions);
fprintf('Total number of unique condition across all sessions: \t\t%d \n', results.total_nr_of_unique_conditions)
fprintf('Total number of conditions repeated across all sessions: \t%d \n', results.total_nr_of_conditions-results.total_nr_of_unique_conditions)
fprintf('Total number of CD+ trials (out of %d CD trials): \t\t%d (%2.1f%%)\n', results.total_nr_of_CD_trials, results.total_nr_of_CDplus_trials, 100*(results.total_nr_of_CDplus_trials/results.total_nr_of_CD_trials))
fprintf('Total number of OBJ-PC+ trials (out of %d OBJ-PC trials): \t%d (%2.1f%%)\n', results.nr_of_objpc_trials, results.nr_of_objectcatch_trials, 100*(results.nr_of_objectcatch_trials/results.nr_of_objpc_trials))
fprintf('Total number of unique conditions for runs 1-%d: \t\t%s \n', results.total_nr_of_runs,num2str(results.nr_of_unique_conditions_per_run));
fprintf('Total number of blocks for runs 1-%d: \t\t\t\t%s \n', results.total_nr_of_runs,num2str(blocks_per_run))
fprintf('Total number of trials for runs 1-%d: \t\t\t\t%s \n', results.total_nr_of_runs,num2str(trials_per_run))
fprintf('Total number of trials per task class for runs 1-%d: \t\t%s \n', results.total_nr_of_runs, num2str(results.trials_per_taskclass))

if strcmp(env_type,'BEHAVIOR')
    
end

return






