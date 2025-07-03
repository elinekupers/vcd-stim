function [] = vcd_checkTimeTable(time_table_master,varargin)
% VCD bookkeeping function to check columns in the time table for any
% inconsistencies and deviations from expectations.
%   
%    [] = vcd_checkTimeTable(time_table_master,data_results)
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
% -Do we sample all at least 4 superordinate categories within an OBJ or NS block?
% -To what extent are button presses balanced for each crossing?
% -Does this hold across blocks/runs/the entire session?
% 
% Yoked columns:
% * Session_level (No NaN’s)
%     * session_nr
%     * session_type
%     * run_nr
%     * global_run_nr
%     * event_start
%     * event_dur
%     * event_end
%     * event_id
%     * event_name
% * Block level (NaN’s for IBIs/pre-/post-blank periods)
%     * block_nr
%     * global_block_nr
% * Semiblock level (NaN’s for ITIs/IBIs/pre-/post-blank periods)
%     * crossing_nr
%     * crossing_name
%     * nr_fix_changes
% * Trial level (NaN’s for ITIs)
%     * trial_nr
%     * global_trial_nr
%     * stim_class
%     * stim_class_name
%     * task_class
%     * task_class_name
%     * is_cued
%     * trial_type
% * Condition_level:
%     * condition_nr
%     * condition_name
%     * stim_nr_left
%     * stim_nr_right
%     * is_objectcatch
%     * correct_response
%     * cd_start
%     * orient_dir
%     * contrast
%     *  gbr_phase     
%     * rdk_coherence    
%     * super_cat          
%     * super_cat_name         
%     * basic_cat              
%     * basic_cat_name              
%     * sub_cat               
%     *  sub_cat_name               
%     * affordance_cat         
%     * affordance_name          
%     * stim2_im_nr    
%     * stim2_delta    
%     * stim2_orient_dir  
%     * repeat_nr
%     * is_catch
%     * is_special_core
%     * is_lure
%
% INPUT:
% * time_table_master    : (table) Table defining trials at the granularity 
%                         of sub-trial events. Table can be for behavioral 
%                         or MRI version of the VCD experiment. Rows 
%                         represent the individual trial events, columns  
%                         code different aspects of the experiment (e.g., 
%                         condition_nr, block_nr, etc.)
% OUTPUT:
% * None
%
% Example:
% subj_nr = 999; 
% dataDir = fullfile(vcd_rootPath,'data','BEHAVIOR',sprintf('vcd_subj%03d',subj_nr));
% d       = dir(fullfile(dataDir,subjID,sprintf('%s_time_table_master_PPROOM_EIZOFLEXSCAN*.mat',subjID)));
% a1      = load(fullfile(d(end).folder,d(end).name),'time_table_master')
% vcd_checkTimeTable(a1.time_table_master)
% 
%
% Written by E Kupers @ UMN 2025/06

% Check inputs
if nargin==1 && ~exist('behresults','var')
    behresults = [];
elseif nargin==2 && ~isempty(varargin{1})
    behresults = varargin{1};
end

% Get params
params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);
params.exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);

all_cond_names = vcd_getConditionNames;

% %%%% Bird's eye view stats %%%% 
results = struct();

% record
results.total_nr_of_sessions    = unique(time_table_master.session_nr);
results.total_nr_of_runs        = max(time_table_master.global_run_nr);
results.total_nr_of_blocks      = max(time_table_master.global_block_nr);
results.total_nr_of_trials      = max(time_table_master.global_trial_nr);

% check if nr of sessions is as we expect
assert(isequal(results.total_nr_of_sessions,[1:max(results.total_nr_of_sessions)]'))

% if we only have one session nr, tell the user 
if length(results.total_nr_of_sessions)==1 
    fprintf('[%s]: Only found 1 session number, assuming this is a behavioral experiment.\n',mfilename)
    assert(isequal(length(results.total_nr_of_sessions), params.exp.session.n_behavioral_sessions));
    env_type = 'BEHAVIORAL';
elseif length(results.total_nr_of_sessions)>1 
    fprintf('[%s]: Only found 1 session number, assuming this is the MRI experiment.\n',mfilename)
    assert(isequal(length(results.total_nr_of_sessions), params.exp.session.n_mri_sessions));
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
if strcmp(env_type,'BEHAVIORAL')
    ses_blocks = params.exp.session.behavior.ses_blocks; 
elseif strcmp(env_type,'MRI')
    ses_blocks = cat(2,params.exp.session.mri.wide.ses_blocks,params.exp.session.mri.deep.ses_blocks); % dim 1:stim, 2:task, 3:session, 4:session type (a/b)
end

% Compare crossing names
expected_crossing_names = {};
for kk = 1:size(ses_blocks,3)
    for ll = 1:size(ses_blocks,4)
        for mm = 1:size(ses_blocks,1)
            tc_names = params.exp.taskclassnames(ses_blocks(mm,:,kk,ll)>0);
            for nn = 1:length(tc_names)
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
% check if columns are yoked as expected

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

% prior and row relative to fixation block start and end should be either
% ibi or pre-/post-blank
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
% More run specific stats
total_nr_of_unique_conditions = length(unique(time_table_master.condition_nr(~isnan(time_table_master.condition_nr))));

% Stats per run:
[~,run_start] = unique(time_table_master.global_run_nr);
run_start = cat(1,run_start,size(time_table_master,1));

trials_per_run = zeros(1,results.total_nr_of_sessions);
for rr = 1:length(run_start)-1
    trials_per_run(rr) = length(time_table_master.trial_nr(run_start(rr):(run_start(rr+1)-1)));
end

trials_per_taskclass = zeros(1,results.total_nr_of_sessions);
for rr = 1:length(run_start)-1
    trials_per_taskclass(rr) = length(time_table_master.task_class(run_start(rr):(run_start(rr+1)-1)));
end

results.trials_per_run = trials_per_run;
results.trials_per_taskclass = trials_per_taskclass;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Condition specific stats

%%% Loop over left/right columns
for side = [1,2]
    
    % check conditon nrs and names
    condition_nan                   = find(isnan(time_table_master.condition_nr(:,side))); 
    results.condition_nrs{:,side}   = time_table_master.condition_nr(~isnan(time_table_master.condition_nr(:,side)),side);
    results.condition_names{:,side} = time_table_master.condition_name(~isnan(time_table_master.condition_nr(:,side)),side);
    
    % we expect conditon nrs and names columns to be yoked in terms of NaNs 
    assert(isequal(isnan(time_table_master.condition_nr(:,side)), cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.condition_name(:,side))))

    % we expect conditon nrs to be yoked in terms of NaNs for left/right
    % stimulus numbers
    if side==1
        assert(isequal(condition_nan,stim_nr_left_nan))
        filter_me = condition_nan;
        
    elseif side==2
        assert(isequal(condition_nan,stim_nr_right_nan))
        filter_me = ismember(time_table_master.stim_class,[1:5,99]);
    end
    
    % we should expect the same nr of NaNs as for contrast, repeat
    % nr and special core column
    assert(isequal(condition_nan, find(isnan(time_table_master.contrast(:,side)))));
    assert(isequal(condition_nan, find(isnan(time_table_master.repeat_nr(:,side)))));
    assert(isequal(condition_nan, find(isnan(time_table_master.is_special_core(:,side)))));
    
    % For WM cases, we check if this is true for those task-crossings only
    assert(isequal(isnan(time_table_master.condition_nr(time_table_master.task_class==5,side)), isnan(time_table_master.stim2_im_nr(time_table_master.task_class==5,side))))
    assert(isequal(isnan(time_table_master.condition_nr(time_table_master.task_class==5,side)), isnan(time_table_master.stim2_delta(time_table_master.task_class==5,side))))
    assert(isequal(isnan(time_table_master.condition_nr(time_table_master.task_class==5 & time_table_master.stim_class~=5,side)), ...
        isnan(time_table_master.stim2_orient_dir(time_table_master.task_class==5 & time_table_master.stim_class~=5,side)))); % exclude NS stim since they have no orientation/motion direction/angle change
    
    for sc = [1:length(params.exp.stimclassnames), 99] % 99 is 'all', for scc/ltm
       
        switch sc
            case 1 % Gabor
                gbr_phase_nan  = find(isnan(time_table_master.gbr_phase(strcmp(time_table_master.stim_class_name(:,side),'gabor'),side)));
                assert(isequal(gbr_phase_nan, find(isnan(time_table_master.condition_nr(strcmp(time_table_master.stim_class_name(:,side),'gabor'),side))))); 

            case 2 % RDK
                rdk_coherence_nan   = find(isnan(time_table_master.rdk_coherence(strcmp(time_table_master.stim_class_name(:,side),'rdk'),side)));
                assert(isequal(rdk_coherence_nan, find(isnan(time_table_master.condition_nr(strcmp(time_table_master.stim_class_name(:,side),'rdk'),side))))); 

            case 3 % dot
                % no special dot column
                
            case {4,5} % obj,ns
                obj_ns_idx = strcmp(time_table_master.stim_class_name(:,side),'obj') | strcmp(time_table_master.stim_class_name(:,side),'ns');
                super_cat_nan       = find(isnan(time_table_master.super_cat(obj_ns_idx,side)));
                basic_cat_nan       = find(isnan(time_table_master.basic_cat(obj_ns_idx,side)));
                sub_cat_nan         = find(isnan(time_table_master.sub_cat(obj_ns_idx,side)));
                affordance_cat_nan  = find(isnan(time_table_master.affordance_cat(obj_ns_idx,side)));

                super_cat_name_nan  = find(cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.super_cat_name(obj_ns_idx,side)));
                basic_cat_name_nan  = find(cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.basic_cat_name(obj_ns_idx,side)));
                sub_cat_name_nan    = find(cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.sub_cat_name(obj_ns_idx,side)));
                affordance_name_nan = find(cellfun(@(x) isequalwithequalnans(x,NaN), time_table_master.affordance_name(obj_ns_idx,side)));
                
                assert(isequal(super_cat_nan,       find(isnan(time_table_master.condition_nr(obj_ns_idx,side)))));
                assert(isequal(basic_cat_nan,       find(isnan(time_table_master.condition_nr(obj_ns_idx,side)))));
                assert(isequal(sub_cat_nan,         find(isnan(time_table_master.condition_nr(obj_ns_idx,side)))));
                assert(isequal(affordance_cat_nan,  find(isnan(time_table_master.condition_nr(obj_ns_idx,side)))));
                assert(isequal(super_cat_name_nan,  find(isnan(time_table_master.condition_nr(obj_ns_idx,side)))));
                assert(isequal(basic_cat_name_nan,  find(isnan(time_table_master.condition_nr(obj_ns_idx,side)))));
                assert(isequal(sub_cat_name_nan,    find(isnan(time_table_master.condition_nr(obj_ns_idx,side)))));
                assert(isequal(affordance_name_nan, find(isnan(time_table_master.condition_nr(obj_ns_idx,side)))));
        end
    end
    
    
    for tc = 1:length(params.exp.taskclassnames)
       
        corr_rsp = time_table_master.correct_response(time_table_master.task_class==tc);
        cued_side = 1+mod(time_table_master.is_cued(time_table_master.task_class==tc & time_table_master.event_id==94)-1,2);
        
        switch tc
            
            case 1 % fix
                block_offset_times = [];
                assert(isequalwithequalnans(corr_rsp,NaN(length(corr_rsp),1)));
                block_onset_idx     = (ismember(time_table_master.event_id,91) &  ~isnan(time_table_master.nr_fix_changes));
                block_onset_times   = time_table_master.event_start(block_onset_idx);
                block_nr            = unique(time_table_master.global_block_nr(block_onset_idx));
                for bb = 1:length(block_nr)
                    block_offset_times(bb) = max(time_table_master.event_end(find(ismember(time_table_master.global_block_nr, block_nr(bb)))));
                end                
                expected_presses = sum(floor((block_offset_times'-block_onset_times) ./params.stim.fix.dotmeanchange));
                assert(isequal(expected_presses,sum(time_table_master.nr_fix_changes,'omitnan')));
            case 2 % cd
                nr_trials = params.exp.block.n_trials_single_epoch * length(unique(time_table_master.global_block_nr(time_table_master.task_class==tc)));
                nr_expected_changes = round(nr_trials*params.exp.trial.cd.prob_change);
                assert(isequal(nr_expected_changes, sum(corr_rsp==1)));             
            case 3 % scc
                tmp_name = time_table_master.stim_class_name(time_table_master.task_class==tc & time_table_master.event_id==94,:);
                for pp = 1:length(cued_side)
                    tmp_name_cued(pp) = tmp_name(pp,cued_side(pp));
                end
                 [~,sccnames] = ismember(tmp_name_cued,params.exp.stimclassnames([1,3,2,4])); % note that button presses for dot and rdk are flipped compared to stim class nr
                assert(isequal(corr_rsp(~isnan(corr_rsp)), sccnames'));
                
                assert(all(abs(diff(histcounts(corr_rsp(~isnan(corr_rsp)),[1:5])))<=2)); % less than 2 counts difference in nr of button presses
                
            case 4 % pc
                for stim = [1:4]
                    tmp_ori = time_table_master.orient_dir(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim,:);
                    cued_side = 1+mod(time_table_master.is_cued(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim)-1,2);
                    for cc = 1:length(cued_side)
                        cued_tmp_ori(cc) = tmp_ori(cc,cued_side(cc));
                    end
                    cued_tmp_ori=cued_tmp_ori';
                    if stim == 1
                        assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_ori, params.stim.gabor.ori_deg([3,4,5,6]))))); % more horz
                        assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_ori, params.stim.gabor.ori_deg([1,2,7,8]))))); % more vert
                    elseif stim == 2
                        assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_ori, params.stim.rdk.dots_direction([3,4,5,6]))))); % more horz
                        assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_ori, params.stim.rdk.dots_direction([1,2,7,8]))))); % more vert
                    elseif stim == 3
                        assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_ori, params.stim.dot.ang_deg([3,4,5,6,11,12,13,14]))))); % more horz
                        assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_ori, params.stim.dot.ang_deg([1,2,7,8,9,10,15,16]))))); % more vert
                    elseif stim == 4
                        assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_ori, params.stim.obj.facing_dir_deg(params.stim.obj.facing_dir_deg < 45 | params.stim.obj.facing_dir_deg > 135))))); % more sideways
                        assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_ori, params.stim.obj.facing_dir_deg(params.stim.obj.facing_dir_deg > 45 & params.stim.obj.facing_dir_deg < 135))))); % more forward
                    end
                end
                
                % NS indoor/outdoor im
                cued_tmp_ori = time_table_master.basic_cat(time_table_master.task_class==tc & time_table_master.stim_class==5,1);
                assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_ori, 1)))); % indoor
                assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_ori, 2)))); % outdoor
                
                assert(all(abs(diff(histcounts(corr_rsp(~isnan(corr_rsp),[1:3])))==0))); % equal nr of button presses
                
            case 5 % wm
                for stim = [1:4]
                    cued_side = 1+mod(time_table_master.is_cued(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim)-1,2);
                    tmp_ori1      = time_table_master.orient_dir(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim,:);
                    tmp_ori2      = time_table_master.stim2_orient_dir(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim,:);
                    tmp_delta     = time_table_master.stim2_delta(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim,:);
                    
                    for cc = 1:length(cued_side)
                        cued_tmp_ori1(cc) = tmp_ori1(cc,cued_side(cc));
                        cued_tmp_ori2(cc) = tmp_ori2(cc,cued_side(cc));
                        cued_tmp_delta(cc) = tmp_delta(cc,cued_side(cc));
                    end
                    cued_tmp_ori1=cued_tmp_ori';
                    cued_tmp_ori2=cued_tmp_ori2';
                    cued_tmp_delta=cued_tmp_delta';
                    
                    diff_ori = cued_tmp_ori2-cued_tmp_ori1;
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
                tmp_delta = time_table_master.stim2_delta(time_table_master.task_class==tc & time_table_master.stim_class==5,1);
                assert(isequal(sum(corr_rsp==1),sum(tmp_delta<0))); % remove
                assert(isequal(sum(corr_rsp==2),sum(tmp_delta>0))); % add
                
                assert(all(abs(diff(histcounts(corr_rsp(~isnan(corr_rsp),[1:3])))==0))); % equal nr of button presses
                
            case 6 % ltm
                % under construction
            case 7 % img   
                % under construction
            case 8 % what
                all_cued_cat =[];
                for stim = [4,5]
                    if stim == 4
                        cued_side = 1+mod(time_table_master.is_cued(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim)-1,2);
                        tmp_cat = time_table_master.super_cat(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim,:);
                        for cc = 1:length(cued_side)
                            cued_tmp_cat(cc) = tmp_cat(cc,cued_side(cc));
                        end
                        cued_tmp_cat=cued_tmp_cat';
                                                
                    elseif stim == 5
                        cued_tmp_cat = time_table_master.super_cat(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim,1);
                    end
                    
                    assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_cat,1)))); % human
                    assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_cat,2)))); % animal
                    assert(isequal(sum(corr_rsp==3),sum(ismember(cued_tmp_cat,[3,4])))); % object, food
                    assert(isequal(sum(corr_rsp==4),sum(ismember(cued_tmp_cat,5)))); % place
                    
                    all_cued_cat = cat(1, all_cued_cat, cued_tmp_cat);
                end
                
                assert(isequal(histcounts(all_cued_cat, [1:6]), hiscounts(corr_rsp,[1:6]))); % distribution of button presses should match nr of supercategorie
                
            case 9 % where
                cued_tmp_loc = time_table_master.sub_cat(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==5,1);
                assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_loc,1)))); % left
                assert(isequal(sum(corr_rsp==3),sum(ismember(cued_tmp_loc,2)))); % center
                assert(isequal(sum(corr_rsp==4),sum(ismember(cued_tmp_loc,3)))); % right
                
                assert(diff(histcounts(corr_rsp,[1:4]))==0); % equal nr of button presses
                
            case 10 % how
                all_cued_aff =[];
                for stim = [4,5]
                    if stim == 4
                        cued_side = 1+mod(time_table_master.is_cued(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim)-1,2);
                        tmp_cat = time_table_master.super_cat(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim,:);
                        for cc = 1:length(cued_side)
                            cued_tmp_afford(cc) = tmp_cat(cc,cued_side(cc));
                        end
                        cued_tmp_afford=cued_tmp_afford';
                    elseif stim == 5
                        cued_tmp_afford = time_table_master.affordance_cat(time_table_master.task_class==tc & time_table_master.event_id==94 & time_table_master.stim_class==stim,1);
                    end
                    assert(isequal(sum(corr_rsp==1),sum(ismember(cued_tmp_afford,1)))); % greet
                    assert(isequal(sum(corr_rsp==2),sum(ismember(cued_tmp_afford,2)))); % grasp
                    assert(isequal(sum(corr_rsp==3),sum(ismember(cued_tmp_afford,3)))); % enter/walk
                    assert(isequal(sum(corr_rsp==4),sum(ismember(cued_tmp_afford,4)))); % observe/do nothing
                    
                    all_cued_aff = cat(1, all_cued_aff, cued_tmp_cat);
                end
                
                assert(isequal(histcounts(all_cued_aff, [1:5]), hiscounts(corr_rsp,[1:5]))); % distribution of button presses should match nr of supercategories  
        end
    end
end

% CD onset checks
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

% ltm stats
% assert(isequal(condition_nan, find(isnan(time_table_master.is_lure(:,side)))));

% cued stats 
results.cued_count = histcounts(time_table_master.is_cued(time_table_master.event_id==94),[1:4]);

% Middle letter should be L/C for first condition_name sub column
loc1 = cellfun(@(x) strsplit(x,'-') ,results.condition_names{1,1},'UniformOutput',0);
for xx = 1:length(loc1)
    assert(ismember(loc1{xx}{3},{'L','C'}));
end
% Middle letter should be R for second condition_name sub column
loc2 = cellfun(@(x) strsplit(x,'-') ,results.condition_names{1,2},'UniformOutput',0);
for xx = 1:length(loc2)
    assert(ismember(loc2{xx}{3},{'R'}));
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we provided behavioral results, see if these results match with time_table_master
if ~isempty(behresults)

    nr_runs_behresults = size(behresults,2);
    assert(isequal(total_nr_of_runs,nr_runs_behresults));
    
    all_crossing_nr_beh = [];
    
    for rr = 1:nr_runs_behresults
        cond_nrs_beh = find(behresults(rr).conditioncount>0)';
        repeats = find(behresults(rr).conditioncount>1);
        if ~isempty(repeats)
            for cc = 1:length(repeats)
                cond_nrs_beh = cat(1,cond_nrs_beh, repelem(repeats(cc), behresults(rr).conditioncount(repeats(cc))-1,1));
            end
        end
        cond_names_beh = all_cond_names(cond_nrs_beh)';
        
        
        cond_nrs_tt = time_table_master.condition_nr(time_table_master.event_id==94 & time_table_master.run_nr==rr,:);
        cond_nrs_tt = cond_nrs_tt(~isnan(cond_nrs_tt));

        cond_names_tt = time_table_master.condition_name(time_table_master.event_id==94 & time_table_master.run_nr==rr,:);
        cond_names_tt = cond_names_tt(:);
        cond_names_tt(cellfun(@(x) isequalwithequalnans(x,NaN), cond_names_tt)) = [];
        
        cond_nrs_tt_counts = histcounts(cond_nrs_tt,[1:length(all_cond_names)+1]);
        
        % Condition names
        assert(isequal(sort(cond_names_tt),sort(cond_names_beh)))
        
        % Condition nrs should match for every run
        assert(isequal(sort(cond_nrs_beh),sort(cond_nrs_tt)));
        
        % Condition nr tally should match for every run
        assert(isequal(cond_nrs_tt_counts,behresults(rr).conditioncount));
        
        % Check if session nrs match
        assert(isequal(behresults(rr).trialinfo.session_nr, time_table_master.session_nr(time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
        
        % Check if run nrs match
        assert(isequal(behresults(rr).trialinfo.run_nr, time_table_master.run_nr(time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));

        % Check if block nrs match
        assert(isequal(behresults(rr).trialinfo.block_nr, time_table_master.block_nr(time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
        
        % Check if crossing nrs match
        assert(isequal(behresults(rr).trialinfo.crossing_nr, time_table_master.crossing_nr(time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
        
        % Check if trial nrs match
        assert(isequal(behresults(rr).trialinfo.trial_nr, time_table_master.trial_nr(time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
        
        % Check if catch trial nrs match
        assert(isequal(behresults(rr).trialinfo.is_catch, time_table_master.is_catch(time_table_master.event_id==94 & time_table_master.run_nr==rr,:)));
        
        % Check if stim cued match
        stimclass = time_table_master.stim_class_name(time_table_master.event_id==94 & time_table_master.run_nr==rr,:);
        cued_side = 1+mod(time_table_master.is_cued(time_table_master.event_id==94 & time_table_master.run_nr==rr)-1, 2);
        stimclass_cued = cell(size(stimclass,1),1);
        stimclass_cued(cued_side==1) = stimclass(cued_side==1,1);
        stimclass_cued(cued_side==2) = stimclass(cued_side==2,2);
        assert(isequal(behresults(rr).trialinfo.stim_cued,stimclass_cued))
        
        all_crossing_nr_beh = cat(1,all_crossing_nr_beh, behresults(rr).summary.crossing_nr);
    end
    assert(isequal(expected_crossing_nrs,sort(all_crossing_nr_beh(:))))
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRINT RESULTS %%%

fprintf('[%s]: Stats of time_table_master:\n',mfilename) 
fprintf('Total number of sessions: \t\t%03d\n', total_nr_of_sessions);
fprintf('Total number of runs: \t\t\t%03d\n', total_nr_of_runs);
fprintf('Total number of block: \t\t\t%03d\n', total_nr_of_blocks);
fprintf('Total number of trials: \t\t%03d\n', total_nr_of_trials);
fprintf('Total number of unique conditions: \t%03d\n', total_nr_of_unique_conditions);
fprintf('Total number of trials for runs 1-%d: \t%s \n', total_nr_of_runs,num2str(trials_per_run))
fprintf('Total number of trials per task class for runs 1-%d: \t%s \n', total_nr_of_runs, num2str(trials_per_taskclass))
fprintf('Total number of conditions across all sessions: \t%d \n', length(results.condition_nrs(:)))
fprintf('Total number of condition repeats across all sessions: \t%d \n', length(unique(results.condition_nrs(:))))

return






