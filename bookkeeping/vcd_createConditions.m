function [params, condition_master, all_unique_im, all_cond] = vcd_createConditions(params,varargin)
% VCD function to define, shuffle, and organize trials into blocks:
%
%   [condition_master, all_unique_im, all_cond] = ...
%       vcd_createConditions(params,'load_params',<load_params>,'store_params',<store_params>)
%
% Purpose:
%   1. Define all possible stimulus-task combinations depending on session,
%   run, and subject.
% For every stimulus-task combination,
%   2. Define unique images (all images within a stimulus class, where
%       stimulus features of interest are fully-crossed.
%   3. Pseudo-randomly shuffle unique images according to prioritize
%       stimulus feature of interest.
%   4. Merge unique stimuli with left/right parafoveal stimulus locations
%       into a single trial, pseudo-randomly assign covert spatial
%       attention cue and fixation dot rim thickening direction.
%   5. Allocate trials to blocks
%
% INPUTS:
%  params         : (required) Struct with stimulus and session params
%  [load_params]  : (optional) Load prior stored parameters or not. Will
%                     search for a file name starting with trials_<dispname>*.mat
%                     in fullfile(vcd_rootPath,'workspaces','info').
%                     Default: true
%  [store_params] : (optional) Store generated parameters or not.
%                     Will store params in fullfile(vcd_rootPath,'workspaces','info');
%                     Default: true
%  [session_env] : (optional) Are we dealing with MRI or BEHAVIOR sessions
%                     Default: 'MRI'
%
% OUTPUT:
%  condition_master : N rows x M table that lists every unique image presented
%                       in VCD-core, for every task crossing separated by
%                       individual trial.
%                    Each row is an individual image.
%                    For Gabors/RDKs/Dots/ComplexObj, two subsequent rows
%                    define a unique trial. For NS, there is only 1 row per
%                    unique trial, as there is only one central image.
%
%   Columns contains the following information:
%   1: {'unique_im_nr'   } unique image nr"
%   2: {'stimloc'        } stimulus location relative to fixation dot: 1 = left, 2 = right, 3 = central
%   3: {'stimloc_name'   } same as column two but then in text 
%   4: {'orient_dir'     } orientation (gabor), motion direction (rdk), or angle (dot), or facing direction (obj) in deg
%   5: {'contrast'       } stimulus contrast (Michelson fraction)
%   6: {'gbr_phase'      } gabor stimulus phase  (NaN for non gabor stim)
%   7: {'rdk_coherence'  } rdk stimulus dot coherence (fraction of dots, NaN for non rdk stim)
%   8: {'super_cat'      } superordinate object category level (for obj and ns) 
%   9: {'basic_cat'      } basic object category level (for obj and ns)
%   10: {'sub_cat'        } subordinate object category level (for obj and ns)
%   11: {'super_cat_name' } same as 8, but then in text
%   12: {'basic_cat_name' } same as 9, but then in text
%   13: {'sub_cat_name'   } same as 10, but then in text
%   14: {'stim_class_name'} stimulus class name (gabor, rdk, dot, obj, ns)
%   15: {'stim_class'     } stim class nr (1:5)
%   16: {'task_class_name'} task class name (fix/cd/scc/pc/wm/ltm/img/what/where/how)
%   17: {'task_class'     } task class nr (1:10)
%   18: {'is_cued'         } cue status: 0 = uncued, 1 = cued
%   19: {'unique_trial_nr'} unique trial nr
%   20: {'thickening_dir' } thickening direction of fixation dot rim (for spatial attention cue) 1 = left, 2 = right, 3 = both/neutral
%   21: {'stim2_delta'    } for double epoch tasks, what predefined delta between stim 1 and stim 2 did we choose from p.stim.(<stim_class_name>).delta_ref
%   22: {'stim2'          } same as 21, but then update stim feature of stim2 (e.g., 80 deg orientation for a stim1: 95 - 15 deg)
%   23: {'ltm_stim_pair'  } for LTM task: each unique image nr is associated with
%                           another stimulus
%   24: {'islure'         } for LTM task: whether we use lure stimulus or not
%   25: {'repeat_nr'      } Keep track how many times has this unique image has
%                           been repeated thusfar in the experiment
%
%  There are two hierarchical structs, which are the building blocks of
%  the condition master table, listing the unique images and unique 
%  conditions for each stimulus class:
%  * all_unique_im  : matrix labeling the unique image features
%  * all_conds      : stimulus "condition master" matrix describing the
%                    pseudo-random shuffled order of trials for one
%                    stimulus class, across all sessions.

%
%
% Written by Eline Kupers November 2024 (kupers [at] umn [dot] edu)


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'        , @isstruct);
p0.addParameter('load_params'  , true, @islogical);
p0.addParameter('store_params' , true, @islogical);
p0.addParameter('session_env' , 'MRI', @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));

% Parse inputs
p0.parse(params,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

%% Load params if requested and we can find the file
if load_params
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('condition_master_%s_%s*.mat',params.disp.name,session_env)));
    fprintf('\n[%s]: Found %d condition_master .mat file(s)\n',mfilename,length(d));
    if ~isempty(d)
        if length(d) > 1
            warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
        end
        fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
        load(fullfile(d(end).folder,d(end).name));
    else
        error('[%s]: Can''t find trial params file!\n', mfilename)
    end
else % Recreate conditions and blocks and trials
    if ~isfield(params,'exp')
        warning('[%s]: params.exp doesn''t exist!)',mfilename)
        warning('[%s]: Will try to loading exp params from file.\n', mfilename);
        d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('exp_%s*.mat',params.disp.name)));
        if ~isempty(d)
            if params.verbose
                fprintf('\n[%s]: Found %d exp params .mat file(s)\n',mfilename,length(d));
                if length(d) > 1
                    warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
                end
                fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
            end
            tmp = load(fullfile(d(end).folder,d(end).name));
            params.exp = tmp.exp;
        else
            if params.verbose
                warning('[%s]: Can''t find exp session .mat files! Will run vcd_getSessionParams.m', mfilename);
            end
            params.exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
        end
    end
    if ~isfield(params,'stim')
        if params.verbose
            warning('[%s]: params.stim doesn''t exist!)',mfilename)
            warning('[%s]: Will try to loading stim params from file.\n', mfilename);
        end
        d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('stim_%s*.mat',params.disp.name)));
        if ~isempty(d)
            if params.verbose
                fprintf('\n[%s]: Found %d stim params .mat file(s)\n',mfilename,length(d));
                if length(d) > 1
                    warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
                end
                fprintf('[%s]: Loading stim params .mat file: %s\n', mfilename, d(end).name);
            end
            tmp = load(fullfile(d(end).folder,d(end).name));
            params.exp = tmp.exp;
        else
            if params.verbose
                warning('[%s]: Can''t find stim session .mat files! Will run vcd_getStimParams.m', mfilename);
            end
            params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);  
        end
    end
    
    %% Preallocate space and set up tables/structs
    if params.verbose
        tic
        fprintf('\n[%s]: Start creating conditions for.. \n',mfilename);
    end
    
    all_unique_im     = struct();  % details about unique core images present in VCD-core experiment. This information is also present in condition_master.
    all_cond          = struct();  % details about unique conditions (unique images x cueing condition) present in VCD-core experiment. This information is also present in condition_master.
    
    % master table with all conditions shown across all sessions. 
    % Note that block_nr counts continuously because we will shuffle the 
    % block order later for each subject. Same holds for run_nr, these
    % numbers do not represent the final nr of runs as they will change
    % depending on the allocation of blocks to each run within a session.
    condition_master  = table();   
    
    % Define block content for each stimulus class
    for stimClass_idx = 1:length(params.exp.stimclassnames)
        
        stim_table = table();
        
        % Get stimclass name
        stimClass_name = params.exp.stimclassnames{stimClass_idx};
        
        % Get corresponding task crossings
        bsc_idx = cell2mat(cellfun(@(x) strcmp(params.exp.stimclassnames,x), {stimClass_name}, 'UniformOutput', false));
        
        if strcmp(session_env,'MRI')
            task_crossings = find( params.exp.crossings(bsc_idx,:));
        elseif strcmp(session_env,'BEHAVIOR')
            task_crossings = find(params.exp.n_unique_trial_repeats_behavior(bsc_idx,:)>0);
        end
                
        % Loop over each task crossing for this stim class
        for curr_task = 1:length(task_crossings)

            % Define the unique images for Gabors
            [t_cond, n_unique_cases] = vcd_defineUniqueImages(params, stimClass_name);
            
            all_unique_im.(stimClass_name) = t_cond;
            
            % Add stimclass name and idx to temporary table
            t_cond.stim_class_name = repmat({stimClass_name}, size(t_cond,1),1);
            t_cond.stim_class      = repmat(stimClass_idx,    size(t_cond,1),1);
            
            % Add task name to temp table
            taskClass_name         = params.exp.taskclassnames{task_crossings(curr_task)};
            t_cond.task_class_name = repmat({taskClass_name},size(t_cond,1),1);
            t_cond.task_class      = repmat(task_crossings(curr_task),size(t_cond,1),1);
            
            % If IMG or LTM, we only want to use a subset of unique images
            if strcmp(taskClass_name,'ltm') || strcmp(taskClass_name,'img')
                t_cond = t_cond(t_cond.is_special_core,:);
            end
            
            % Get the number of trials and blocks, which depend on task
            if params.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                n_trials_per_block    = params.exp.block.n_trials_single_epoch;
            else
                n_trials_per_block    = params.exp.block.n_trials_double_epoch;
            end
            
            %% ---- IMPORTANT FUNCTION: Create condition master table ---- %%
            % Create condition master table, where unique images are
            % shuffled and distributed across blocks according to the
            % stimulus feature of interest that receive priority. For
            % example, we want to make sure we present all 3 gabor
            % contrasts within a block at least once.
            tbl = vcd_createConditionMaster(params, t_cond, n_trials_per_block,session_env);
            
            % --- Allocate trials to blocks ---
            tbl.stim_class_unique_block_nr = NaN(size(tbl,1),1);
            tbl.trial_nr = NaN(size(tbl,1),1);
            
            % assign block nr
            block_start_idx = 1:n_trials_per_block:size(tbl,1);
            for mm = block_start_idx
                selected_trials = mm:(mm+(n_trials_per_block-1));
                if mm==block_start_idx(end)
                    selected_trials = selected_trials(selected_trials<=size(tbl,1));
                end
                block_nr = find(mm==block_start_idx);
                tbl.stim_class_unique_block_nr(selected_trials) = block_nr;
                tbl.trial_nr(selected_trials) = 1:length(selected_trials);
            end
            
            % set Nans for cd onset (we will deal with it later)
            tbl.cd_start = NaN(size(tbl,1),1);
            
            % Concatete single stim-task crossing table to master table
            stim_table = cat(1,stim_table,tbl);
            
        end % task class idx
        
        % Accumulate condition master and info
        all_cond.(stimClass_name) = stim_table;
        
        condition_master = cat(1,condition_master,stim_table);
    end % stim idx
    
    
    % ---- IMPORTANT FUNCTION: Shuffle stimuli for SCC and LTM task ----
    condition_master = vcd_shuffleStimForTaskClass(params,  'scc', condition_master, params.exp.block.n_trials_single_epoch,session_env);
    if strcmp(session_env,'MRI')
        condition_master = vcd_shuffleStimForTaskClass(params,  'ltm', condition_master, params.exp.block.n_trials_double_epoch,session_env);
    end
    
    % ---- IMPORTANT FUNCTION: Add correct button press ----
    button_response = NaN(size(condition_master,1),1);
    for ii = 1:size(condition_master,1)
        if  condition_master.is_catch(ii)
            button_response(ii) = NaN;
        else
            button_response(ii) = vcd_getCorrectButtonResponse(params, condition_master(ii,:));
        end
    end
    condition_master.correct_response = button_response;

    % ---- bookkeeping: See if we can achieve roughly equally sample scc images within a block
    scc_trials0 = find(condition_master.stim_class==99);
    cued0 = condition_master.is_cued(scc_trials0);
    scc_trials = scc_trials0; % make a copy
    if ~isempty(scc_trials)
        while 1
            bb_bpress = reshape(condition_master.correct_response(scc_trials),params.exp.block.n_trials_single_epoch,[]);
            bb_cued   = reshape(cued0,params.exp.block.n_trials_single_epoch,[]);
            n0 = zeros(size(bb_bpress,2),4);
            m0 = zeros(size(bb_bpress,2),2);
            for bb = 1:size(bb_bpress,2)
                n0(bb,:) = histcounts(bb_bpress(:,bb),[1:5]);
                m0(bb,:) = histcounts(bb_cued(:,bb));
            end
            too_few = sum(n0==0,1);
            cued_ok = sum(m0==0,1);
            if all(too_few==false(1,4)) && all(cued_ok==false(1,2))
                break;
            else
                % shuffle trials while preserving 50:50 cue left/right
                cued00 = cued0; scc_trials00 = scc_trials;
                cueleft_idx0         = find(cued0==1);
                cueright_idx0        = find(cued0==2);
                cueleft_idx          = cueleft_idx0(randperm(length(cueleft_idx0),length(cueleft_idx0)));
                cueright_idx         = cueright_idx0(randperm(length(cueright_idx0),length(cueright_idx0)));
                cued0(cueleft_idx0)  = cued0(cueleft_idx); 
                cued0(cueright_idx0) = cued0(cueright_idx);
                assert(isequal(cued00,cued0)); % cued status should not change
                scc_trials(cueleft_idx0)  = scc_trials(cueleft_idx);
                scc_trials(cueright_idx0) = scc_trials(cueright_idx0);
            end
        end
        % reshuffle order
        condition_master(scc_trials0,:) = condition_master(scc_trials,:);
        % fix unique trial nr
        scc_trials = find(condition_master.stim_class==99);
        condition_master.unique_trial_nr(scc_trials) = [1:length(condition_master.unique_trial_nr(scc_trials))]';
        condition_master.trial_nr(scc_trials) = repmat([1:params.exp.block.n_trials_single_epoch]',length(condition_master.unique_trial_nr(scc_trials))/params.exp.block.n_trials_single_epoch,1);
        condition_master.stim_class_unique_block_nr(scc_trials) = repelem([1:length(condition_master.unique_trial_nr(scc_trials))/params.exp.block.n_trials_single_epoch],params.exp.block.n_trials_single_epoch)';
    end
    
    
    % ---- bookkeeping: Check if button presses are balanced to the extent possible
    curr_sc = unique(condition_master.stim_class);
    curr_tc = unique(condition_master.task_class);
    curr_tc = curr_tc(curr_tc~=1); % exclude fixation
    for jj = 1:length(curr_sc)
        for mm = 1:length(curr_tc)
            resp = condition_master.correct_response(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
            n = histcounts(resp);
            % For WHAT/HOW tasks we go by category info, because we combine
            % object and foods into one button press..
            if ismember(curr_sc(jj),[4,5]) && ismember(curr_tc(mm),8)
                supercat = condition_master.super_cat(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                if curr_sc(jj) == 4 % objects
                    cued0  = condition_master.is_cued(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
                    n0     = histcounts([supercat(cued0==1,1);supercat(cued0==2,2)]);
                elseif curr_sc(jj) == 5 % scenes, no left/right, only center
                    n0     = histcounts(supercat);
                end
                n1 = [n0([1,2]), n0(3)+n0(4), n0(5)];
                assert(all(n==n1))
            elseif ismember(curr_sc(jj),[4,5]) && ismember(curr_tc(mm),10)
                affordcat = condition_master.affordance_cat(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                if curr_sc(jj) == 4
                    cued0 = condition_master.is_cued(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
                    n1    = histcounts([affordcat(cued0==1,1);affordcat(cued0==2,2)]);
                elseif curr_sc(jj) == 5
                    n1     = histcounts(affordcat);    
                end
                assert(all(n==n1))
            elseif ismember(curr_sc(jj),99) && ismember(curr_tc(mm),3) % scc-all
                % Check if we sample all core images across trials
                stmclass0 = condition_master.stim_class_name(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                cued0     = condition_master.is_cued(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
                stmclass0_cued = [stmclass0(cued0==1,1);stmclass0(cued0==2,2)];
                [~,stmclass_cued_i] = ismember(stmclass0_cued,params.exp.stimclassnames([1,3,2,4]));
                n1 = histcounts(stmclass_cued_i);
                assert(all(n==n1))
            else
                assert(all(diff(n)==0))
            end
        end
    end
    
    %% ---- IMPORTANT FUNCTION: Allocate trials to blocks all unique trials and repeats of trials.
    % !!WARNING!! There is a randomization component involved in creating the
    % conditions (i.e., order of trials allocated to a block). If you don't
    % want this, set params.load_params = true to load an existing
    % condition_master.
    condition_master = vcd_allocateBlocksToRuns(params,condition_master,session_env);
    
    
    %% ---- bookkeeping: Add contrast decrement

    % shuffle [yes/no] responses to CD change such that cued trials have 
    % exactly 20% chance of a contrast decrement change to stimulus 
    % (defined by params.exp.trial.cd.prob_change).
    % For the uncued side, we do not change contrast.
    cd_blocks             = find(condition_master.task_class == 2); % get CD blocks
    cd_blocks_cue         = condition_master.is_cued(cd_blocks); % what side is cued for each CD trial?
    nr_cued_cd_trials     = length(cd_blocks_cue); % nr of cued trials across all CD blocks
    when                  = floor(linspace(1,nr_cued_cd_trials,nr_cued_cd_trials*params.exp.trial.cd.prob_change)); % when do we (roughly) expect a CD change?
    when                  = round(when+rand(1,8)); % add some jitter selected trials
    when(when>nr_cued_cd_trials) = nr_cued_cd_trials;
    assert(isequal(length(unique(when)),length(when)));
    
    expected_cue_cnt      = round(histcounts(cd_blocks_cue)*params.exp.trial.cd.prob_change); % round this probability as we can only deal with integer nr of trials.
    while 1 % Do some voodoo trying to balance the nr of cd changes across spatial cueing conditions..
        cue_cnt = histcounts(cd_blocks_cue(when));
        if any(cue_cnt < expected_cue_cnt)
            cue_type_to_add = find((cue_cnt - expected_cue_cnt)<0);
            [~,cue_type_to_remove] = max(cue_cnt);
            to_rm_idx  = find(cd_blocks_cue==cue_type_to_remove);
            to_add_idx = find(cd_blocks_cue==cue_type_to_add);
            
            rm_idx = find(ismember(when,to_rm_idx));
            already_added_idx = find(ismember(to_add_idx,when'));
            to_add_idx(already_added_idx)= [];
            idx_to_swap = rm_idx(randi(length(rm_idx),1));
            [~,idx_to_add] = min(abs(to_add_idx-when(idx_to_swap)));
            when(idx_to_swap) = to_add_idx(idx_to_add);
        else
            break
        end
        
    end
    
    % Get onset time frame of the cd modulation within a trial 
    % (the modulation function has 2 frames of high contrast, prior to
    % dip). The onset time refers to the time the actual contrast dip
    % happens (so after the second time frames finishes and the third
    % starts).
    cd_change_total = 2.*ones(1,nr_cued_cd_trials);
    cd_change_total(when)=1; % 1=yes change, 2=no change
    for rpt = 1:sum(cd_change_total==1)
        cd_start(rpt) = feval(params.stim.cd.cdsoafun);
    end
    % insert onsets into the condition_master table
    condition_master.cd_start(cd_blocks(cd_change_total==1),1) = cd_start';

    % check nr of changes across all cd trials
    assert(isequal(sum(cd_change_total==1),sum(~isnan(condition_master.cd_start))));
    assert(isequal(sum(cd_change_total==2),sum(isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued>0)))));
    % check nr of changes across left/right cueing conditions
    assert(isequal(sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==1))),...
                   sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==2)))));
               
    % check expected probability to nr of cd trials with a change
    assert(isequal(nr_cued_cd_trials*params.exp.trial.cd.prob_change,sum(~isnan(condition_master.cd_start))))
    
    

    %% Store condition_master if requested
    if params.store_params
        fprintf('[%s]:Storing condition_master..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir,sprintf('condition_master_%s_%s_%s.mat',params.disp.name,session_env,datestr(now,30))),'condition_master','all_unique_im','all_cond')
    end
    
    
    %% At last, do some checks:
    
    % Stim class nr should match with params
    assert(isequal(unique(condition_master.stim_class)',[1:length(params.exp.stimclassnames),99]))
    % Task class nr should match with params
    if strcmp(session_env,'MRI')
        assert(isequal(unique(condition_master.task_class)',1:length(params.exp.taskclassnames)));
    elseif strcmp(session_env,'BEHAVIOR')
        assert(isequal(unique(condition_master.task_class)', find(any(params.exp.n_unique_trial_repeats_behavior,1))))
    end
    % Cued vs uncued stimuli should be matched
    assert(isequal(sum(condition_master.is_cued==1),sum(condition_master.is_cued==2)))
    
    % Now dive into each stimulus class
    unique_im_from_table = [];
    N_tbl = NaN(5,30); 
    adjusted_crossings = params.exp.nr_unique_trials_per_crossing;
    
                
    all_sessions = vcd_getSessionEnvironmentParams(params, session_env);
    for ii = 1:length(params.exp.stimclassnames)
        
        % Check nr of unique stimuli per stimulus class
        if ii == 5
            tmp = condition_master(~condition_master.is_catch & strcmp(condition_master.stim_class_name(:,1), params.exp.stimclassnames{ii}),:);
        
            [N, N_edges] = histcounts(tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})));
        
            N_tbl(ii,1:length(N)) = N;
            assert(isequal(size(N,2), length(all_unique_im.(params.exp.stimclassnames{ii}).unique_im_nr)))

            unique_im_from_table(ii) = length(unique(tmp.stim_nr_left));

        
%             if strcmp(session_env,'MRI')
%                 % lower IMG/LTM block contribution (we have less unique images)
%                 unique_trial_repeats = params.exp.n_unique_trial_repeats_mri;
%             elseif strcmp(session_env,'BEHAVIOR') % NO LTM/IMG
%                 adjusted_crossings(ii,6) = 0;
%                 adjusted_crossings(ii,7) = 0;
%                 half_blocks = all_sessions>0 & all_sessions<1;
%                 if sum(half_blocks)>0
%                     unique_trial_repeats(half_blocks)=all_sessions(half_blocks);
%                 end
%             end
            expected_nr_of_trials = sum(all_sessions(ii,:).*params.exp.nr_trials_per_block(ii,:));
            empirical_nr_of_trials = size(tmp,1);
            assert(isequal(empirical_nr_of_trials, expected_nr_of_trials))
            
        else
            tmp = condition_master(~condition_master.is_catch & any(strcmp(condition_master.stim_class_name, params.exp.stimclassnames{ii}),2),:);
            
            [N, N_edges] = histcounts(cat(1,tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})),....
                                            tmp.stim_nr_right(strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))));
            
            N_tbl(ii,1:length(N)) = N;
            assert(isequal(size(N,2), length(all_unique_im.(params.exp.stimclassnames{ii}).unique_im_nr)))
            colsL = tmp.stim_nr_left(tmp.task_class ~=3 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii}));
            colsR = tmp.stim_nr_right(tmp.task_class ~=3 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}));
            colsLR = [colsL(:),colsR(:)];
            unique_im_from_table(ii) = length(unique(colsLR));
            colsLR_scc = [tmp.stim_nr_left(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); ...
                            tmp.stim_nr_right(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))];
            expected_nr_of_trials = sum(all_sessions(ii,:).*params.exp.nr_trials_per_block(ii,:));
            empirical_nr_of_trials = sum([size(colsLR,1),size(unique(colsLR_scc),1)/2]);
            assert((abs(floor(empirical_nr_of_trials)-expected_nr_of_trials))<=2); % we allow for 1-2 trials difference per stim class due to imbalanced nr of trials for SCC
        end
    end
    
    % Now dive into each task class
    M_tbl = NaN(length(params.exp.taskclassnames),length(params.stim.all_core_im_nrs)+1);
    for ii = unique(condition_master.task_class)'
        [M, edges_M] = histcounts([condition_master.stim_nr_left(~condition_master.is_catch & (condition_master.task_class==ii)); ...
                                    condition_master.stim_nr_left(~condition_master.is_catch & (condition_master.task_class==ii))],[0:1:length(params.stim.all_core_im_nrs)+1]);
        M_tbl(ii,1:length(M)) = M;
    end
    

    
    
    %% Plot figures to check condition master content
    if params.verbose

        vcd_visualizeMasterTable(condition_master, params.store_imgs,session_env);

        figure; set(gcf,'Position',[1,1,1200,300]);
        histogram([condition_master.stim_nr_left;condition_master.stim_nr_right],'numbins',length(params.stim.all_core_im_nrs))
        xlabel('unique condition nr'); ylabel('trial count')
         title('Total unique image nr')
        box off;
        if params.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_env));
            if ~exist(saveFigsFolder,'dir'); mkdir(saveFigsFolder); end
            filename = sprintf('vcd_totaluniqueim.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        figure; set(gcf,'Position',[1,1,1200,300]); imagesc(N_tbl)
        xlabel('unique im nr'); ylabel('Stim class');
        C = parula(20); colormap(C); set(gca,'CLim',[0 max(N_tbl(~isnan(N_tbl)))])
        cb = colorbar; cb.Label.String = 'unique im count';
        cb.Ticks = [0,unique(N_tbl(~isnan(N_tbl)))'];
        set(gca,'YTick',[1:5],'YTickLabel', params.exp.stimclassnames);  set(gca,'XTick',[1:30]); title('Sum of unique im per stimulus class')
        box off;
        if params.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_env));
            filename = sprintf('vcd_uniqueim_per_stimclass.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        
        figure; set(gcf,'Position',[1,1,1200,300]); imagesc(M_tbl)
        xlabel('unique condition nr'); ylabel('Task class');
        C = parula(20); colormap(C); set(gca,'CLim',[0 max(M_tbl(~isnan(M_tbl)))])
        cb = colorbar; cb.Label.String = 'unique im count';
        cb.Ticks = unique(M_tbl(~isnan(M_tbl)))';
        set(gca,'YTick',[1:10],'YTickLabel', params.exp.taskclassnames); 
        set(gca,'XTick',[0:25:110]); title('Sum of unique im per task class')
        box off;
        if params.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_env));
            filename = sprintf('vcd_uniqueim_per_taskclass.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        % Stim class distribution in SCC blocks
        foo = [condition_master.stim_nr_left((condition_master.task_class==3),:), ...
            condition_master.stim_nr_right((condition_master.task_class==3),:)]';
        figure; set(gcf,'Position',[1,1,1200,300])
        imagesc(foo)
        xlabel('Trial nr'); ylabel('left vs right stim');
        C = parula(max(foo(:))); colormap(C); set(gca,'CLim',[min(foo(:)),max(foo(:))])
        cb = colorbar; cb.Label.String = 'unique im nr';
        cb.Ticks = [min(foo(:)):10:max(foo(:))];
        set(gca,'YTick',[1:params.exp.block.n_trials_single_epoch]); title('SCC stimulus class distribution')
        box off;
        if params.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_env));
            filename = sprintf('vcd_scc_stimulusclass_trialdistr.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        if strcmp(session_env,'MRI')
            % Stim class distribution in LTM blocks
            foo = [condition_master.stim_nr_left((condition_master.task_class==6),:), ...
                condition_master.stim_nr_right((condition_master.task_class==6),:)]';
            figure; set(gcf,'Position',[1,1,1200,300])
            imagesc(foo)
            xlabel('Trial nr'); ylabel('left vs right stim');
            C = parula(max(foo(:))); colormap(C); set(gca,'CLim',[min(foo(:)),max(foo(:))])
            cb = colorbar; cb.Label.String = 'unique im nr';
            cb.Ticks = [min(foo(:)):10:max(foo(:))];
            set(gca,'YTick',[1:params.exp.block.n_trials_single_epoch]); title('LTM stimulus class distribution')
            box off;
            if params.store_imgs
                saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_env));
                filename = sprintf('vcd_ltm_stimulusclass_trialdistr.png');
                print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
            end
        end
    end
end % load params

return





