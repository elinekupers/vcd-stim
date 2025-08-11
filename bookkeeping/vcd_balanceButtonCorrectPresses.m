function condition_master = vcd_balanceButtonCorrectPresses(params, condition_master, env_type)
% VCD bookkeeping function to check if button presses are balanced across
% blocks and comply with our constraints. If not, then this function will
% reorder trials until they do comply.
%
%   condition_master = vcd_balanceButtonCorrectPresses(condition_master)
%
% This function will either use session nr/session type info when provided,
% or infer which trials belongs to which block and which session.
%
% Enforced constraints:
% 1. Each block within a session (and session type) but have 50/50
%    left/right cueing for classic stimuli.
% 2. We allow for max 2 repeats in a single-image block (which has 4 trials),
%    and 1 repeat within a double-image block (which has 8 trials)
% 3. Each block within a session must have equal distribution of correct
%    button presses, except for:
%       * WHAT/HOW tasks -- we go by category info, because we combine
%       object and foods into one button press.
%       * WHERE tasks -- we allow for one button press option to have +1
%       when trials cannot be equally divided across the three options
%       * OBJ-PC -- we allow for one button press option to have +1 when
%       trials contain object catch stimuli, which we ignore, resulting in
%       an uneven nr of trials within a block.
%
% Note that we ignore fixation and CD trials for now, as FIX has
% pseudotrials and we may have not selected CD+ trials yet (and will check
% button press distribution separately).
%
% INPUTS:
%   params                :  (struct) parameter struct needed to get subject
%                             nrs and run type params. REQUIRES params.trials
%                             to exist and contain condition_master v0.
%   condition_master      :  (table) general condition_master table with
%                             block and runs. If it does not have session
%                             nr or session type columns, then we will
%                             infer this from "ses_blocks",
%                             params.exp.session.mri.wide.ses_blocks, or
%                             params.exp.session.mri.deep.ses_blocks, or
%                             params.exp.session.behavior.ses_blocks.
%   env_type              :  (str) label to define what type of session we
%                             are defining: 'MRI' or 'BEHAVIOR'.
%
% OUTPUTS:
%  condition_master       : (table) condition_master table with reordered
%                            trials and updated columns for trial_nr,
%                            unique_trial_nr, and
%                            stim_class_unique_block_nr.
%
% Written by Eline Kupers @ UMN 2025/08

% check blocks per session type
[all_sessions,session_types] = vcd_getSessionEnvironmentParams(params, env_type);
% if we have a session number column, then use those to subselect trials
if any(ismember(condition_master.Properties.VariableNames,'session_nr'))
    session_nrs = unique(condition_master.session_nr);
    session_flag = true;
else
    % create temporary session/session type/block nr vector
    session_flag = false;
    session_nrs   = 1:size(all_sessions,3);
    t = [condition_master.stim_class, condition_master.task_class];
    allocated_trials = zeros(size(condition_master,1),1);
end

if strcmp(env_type, 'MRI')
    if params.is_wide
        total_catch_trials = params.exp.nr_catch_trials_wide;
    else
        total_catch_trials = params.exp.nr_catch_trials_deep;
    end
elseif strcmp(env_type, 'BEHAVIOR')
    total_catch_trials = zeros(5,10);
end

% Reset counters
all_stimtask_counter = zeros(length(params.exp.stimclassnames),length(params.exp.taskclassnames)); % 5 stim classes x 10 task classes

% Reset flag
reshuffle_me = false;

% Reset stim_class_unique_block_nr (used by vcd_allocateBlocksToRuns.m) and
% unique_trial_nr. We will update these numbers after reordering trials.
condition_master.stim_class_unique_block_nr = NaN(size(condition_master.stim_class_unique_block_nr));
condition_master.unique_trial_nr = NaN(size(condition_master.unique_trial_nr));

fprintf('[%s]: Check button press distribution:',mfilename);
for ses = session_nrs
    fprintf('\nSession nr %d.',ses);
    if session_flag % create trial index
        ses_idx       = condition_master.session_nr==ses;
        session_types = unique(condition_master.session_type(ses_idx));
    end
    
    for st = session_types
        if ~isnan(session_types(ses,st))
            fprintf('\nSession type %d.',st);
            if session_flag
                % index rows if we have this information in the condition master
                st_idx = (condition_master.session_nr==ses & condition_master.session_type==st);
            else
                stimtask_counter = all_sessions(:,:,ses,st);
                
                for sc = 1:size(stimtask_counter,1)
                    for tc = 2:size(stimtask_counter,2) % skip fixation
                        fprintf('.');
                        
                        % Allocated blocks and trials
                        nr_catch_trials = total_catch_trials(sc,tc);
                        % Get nr of blocks from the params struct
                        nr_blocks = stimtask_counter(sc,tc);
                        % Add catch trials
                        nr_trials = nr_catch_trials + (params.exp.nr_trials_per_block(sc,tc)*nr_blocks);
                        % Update nr of blocks after adding catch trials
                        nr_blocks = ceil(nr_trials/params.exp.nr_trials_per_block(sc,tc));
                        
                        if nr_trials > 0
                            if tc == 3
                                ti = find(t(:,1)==99 & t(:,2)==tc);
                            else
                                ti = find(t(:,1)==sc & t(:,2)==tc);
                            end
                            nr_unique_cond = params.exp.nr_unique_trials_per_crossing(sc,tc);
                            if nr_trials > nr_unique_cond
                                nr_unique_conds = ceil(nr_trials/nr_unique_cond);
                            else
                                nr_unique_conds = nr_unique_cond;
                            end
                            trials_to_be_allocated   = find(allocated_trials(ti,:)==0);
                            trial_to_condition_ratio = nr_trials/nr_unique_conds;
                            % If we have any left over trials from session
                            % type 1 (A), then we want to add a few more
                            % trials as possible optons to avoid getting
                            % stuck in a local minimum.
                            if nr_catch_trials > 0
                                if (st == 2) && (trial_to_condition_ratio >= 0.75 && trial_to_condition_ratio < 1)
                                    if length(ti) >= nr_unique_conds + params.exp.nr_trials_per_block(sc,tc)
                                        nr_unique_conds  = nr_unique_conds + params.exp.nr_trials_per_block(sc,tc);
                                    elseif length(ti) == nr_unique_conds + 1
                                        nr_unique_conds  = nr_unique_conds + 1;
                                    else
                                        error('not sure what to do here')
                                    end
                                end
                            end

                            while 1 % get trials for each stim-task crossing, separate for each session type
                                st_idx    = false(size(condition_master,1),1);
                                
                                % shuffle order within 24/16/30 unique
                                % conditions if button presses are not
                                % equally distributed (or within constraints)
                                while 1
                                    idx = trials_to_be_allocated(1:nr_unique_conds);
                                    if reshuffle_me
                                        if nr_unique_conds == 16 && trial_to_condition_ratio == 1
                                            if length(trials_to_be_allocated) > nr_unique_conds
                                                trials_to_add = choose(nr_unique_conds+4 <= length(trials_to_be_allocated),4,length(trials_to_be_allocated)-nr_unique_conds);
                                                idx = trials_to_be_allocated(1:nr_unique_conds+trials_to_add);
                                            end
                                            order = shuffle_concat(idx,1);
                                        else
                                            order = shuffle_concat(idx,1);
                                        end
                                    else
                                        order = idx;
                                    end
                                    ti(idx) = ti(order);
                                    % find trials
                                    potential_trials = find(allocated_trials(ti)==0);
                                    % get their catch status, object catch status and
                                    % cueing conditions
                                    catch_idx        = condition_master.is_catch(ti(potential_trials(1:nr_trials)))==1;
                                    objectcatch_idx  = condition_master.is_objectcatch(ti(potential_trials(1:nr_trials)))==1;
                                    cued_conds       = condition_master.is_cued(ti(potential_trials(1:nr_trials)));
                                    stim_nr_left     = condition_master.stim_nr_left(ti(potential_trials(1:nr_trials)));
                                    stim_nr_right    = condition_master.stim_nr_right(ti(potential_trials(1:nr_trials)));
                                    % check if (A) if cueing left/right is balanced within a block.
                                    % (B) if we have any stimulus number
                                    % repeats within a block (we allow for 2 repeats in a single-image block (which has 4 trials),
                                    % and 1 repeat within a double-image block (which has 8 trials)
                                    if all(cued_conds==3)
                                        stim_nr_repeats = zeros(nr_blocks,1);
                                        dt = 1:params.exp.nr_trials_per_block(sc,tc):nr_trials;
                                        for bb = 1:nr_blocks
                                            bbi = dt(bb):(dt(bb)+params.exp.nr_trials_per_block(sc,tc)-1);
                                            if length(stim_nr_left(bbi)) - length(unique(stim_nr_left(bbi))) <= 1
                                                stim_nr_repeats(bb) = 1;
                                            end
                                        end
                                        if sum(stim_nr_repeats)==nr_blocks
                                            break;
                                        else
                                            reshuffle_me = true;
                                        end
                                    else
                                        bb_ok = zeros(1,nr_blocks);
                                        stim_nr_repeats = zeros(nr_blocks,2);
                                        dt = 1:params.exp.nr_trials_per_block(sc,tc):nr_trials;
                                        for bb = 1:nr_blocks
                                            bbi = dt(bb):(dt(bb)+params.exp.nr_trials_per_block(sc,tc)-1);
                                            
                                            cued_counts = histcounts(cued_conds(bbi), [1:3]);
                                            if diff(cued_counts)==0
                                                bb_ok(bb) = 1;
                                            end
                                            % Check if we repeat the same stimulus number, if we have more than two repeats per block then we reshuffle
                                            if (params.exp.nr_trials_per_block(sc,tc) == 4) && (length(stim_nr_left(bbi)) - length(unique(stim_nr_left(bbi))) <= 1)
                                                stim_nr_repeats(bb,1) = 1;
                                            elseif (params.exp.nr_trials_per_block(sc,tc) == 8) && (length(stim_nr_left(bbi)) - length(unique(stim_nr_left(bbi))) <= 2)
                                                stim_nr_repeats(bb,1) = 1;
                                            end
                                            if (params.exp.nr_trials_per_block(sc,tc) == 4) && (length(stim_nr_right(bbi)) - length(unique(stim_nr_right(bbi))) <= 1)
                                                stim_nr_repeats(bb,2) = 1;
                                            elseif (params.exp.nr_trials_per_block(sc,tc) == 8) && (length(stim_nr_right(bbi)) - length(unique(stim_nr_right(bbi))) <= 2)
                                                stim_nr_repeats(bb,2) = 1;
                                            end
                                        end
                                        if (sum(stim_nr_repeats(:))/2) == nr_blocks && sum(bb_ok) == nr_blocks
                                            break;
                                            %                                         elseif trial_to_condition_ratio == 1 && sum(catch_idx)>0
                                            %                                             if nr_unique_conds < params.exp.nr_unique_trials_per_crossing(sc,tc)+sum(catch_idx)
                                            %                                                 nr_unique_conds = nr_unique_conds + 1;
                                            %                                             end
                                            %                                             reshuffle_me = true;
                                        else
                                            reshuffle_me = true;
                                        end
                                    end
                                    
                                end
                                assert(all(st_idx(ti(potential_trials(1:nr_trials)))==0))
                                st_idx = (ti(potential_trials(1:nr_trials)));
                                assert(isequal(length(st_idx),nr_trials));
                                if tc == 3, curr_sc = 99; else, curr_sc = sc; end
                                curr_tc = tc;
                                assert(all(ismember(condition_master.stim_class(st_idx), curr_sc)));
                                assert(all(ismember(condition_master.task_class(st_idx), curr_tc)));
                                
                                % check if we have such a crossing
                                if sum(ismember(condition_master.task_class(st_idx), curr_tc) & ismember(condition_master.stim_class(st_idx), curr_sc))>0
                                    
                                    % if so, then set the number of expected button presses
                                    % **** !!! careful this is hardcoded stuff !!! *****
                                    if ismember(curr_tc,[2,4,5,6,7]) % cd, pc, wm, ltm, img have 2 response options
                                        nr_responses = 2;
                                    elseif ismember(curr_tc,[3,8,10]) % scc, what, how have 4 response options
                                        nr_responses = 4;
                                    elseif ismember(curr_tc,[9]) % where has 3 response options
                                        nr_responses = 3;
                                    end
                                    
                                    resp = condition_master.correct_response(st_idx);
                                    resp = resp(~catch_idx);
                                    n = histcounts(resp,1:(nr_responses+1));
                                    
                                    % For WHAT/HOW tasks we go by category info, because we combine
                                    % object and foods into one button press..
                                    if ismember(curr_sc,[4,5]) && ismember(curr_tc,8) % **** !!! careful this is hardcoded stuff !!! *****
                                        supercat = condition_master.super_cat(st_idx,:);
                                        if curr_sc == 4 % objects
                                            cued0  = condition_master.is_cued(st_idx);
                                            n0     = histcounts([supercat(cued0==1,1);supercat(cued0==2,2)],1:6);
                                        elseif curr_sc == 5 % scenes, no left/right, only center
                                            n0     = histcounts(supercat,1:6);
                                        end
                                        n1 = [n0([1,2]), n0(3)+n0(4), n0(5)];  % we combine (3) objects and (4) food into a single button press (the ring finger)
                                        
                                        if (all(n==n1))
                                            % all responses are balanced
                                            reshuffle_me = false;
                                            break;
                                        elseif trial_to_condition_ratio == 1 && sum(catch_idx) > 0
                                            % if we sampled all unique conditions AND have catch trials,
                                            % then we will likely never reach a solution unless we allow for a small difference in response distribution.
                                            if (all(diff(n)<=sum(catch_idx)))
                                                reshuffle_me = false;
                                                break;
                                            else
                                                reshuffle_me = true;
                                            end
                                        else
                                            reshuffle_me = true;
                                        end
                                        
                                    elseif ismember(curr_sc,4) && ismember(curr_tc,4) % PC-OBJ
                                        resp2 = resp(~objectcatch_idx); % only check non-obj catch  trials
                                        n2 = histcounts(resp2,1:(nr_responses+1));
                                        
                                        if (diff(n2)==0)
                                            % ideally all responses are balanced
                                            reshuffle_me = false;
                                            break;
                                            
                                        elseif all(diff(n2)~=0) % but if not, we assume we have an uneven nr of trials due to objectcatch trials.
                                            if mod(length(resp2),2)==1 % if we have an uneven nr of trials after we removed objectcatch trials
                                                if trial_to_condition_ratio == 1 && sum(catch_idx) > 0
                                                    diff_allowed = sum(catch_idx) +1;
                                                else
                                                    diff_allowed = 1;
                                                end
                                                
                                                if ~(diff(n2)<=diff_allowed) % we allow for a difference of one, but not more
                                                    reshuffle_me = true;
                                                else
                                                    reshuffle_me = false;
                                                    break;
                                                end
                                            elseif trial_to_condition_ratio == 1 && sum(catch_idx) > 0
                                                % if we sampled all unique conditions AND have catch trials,
                                                % then we will likely never reach a solution unless we allow for a small difference in response distribution.
                                                if (all(diff(n)<=sum(catch_idx)))
                                                    reshuffle_me = false;
                                                    break;
                                                else
                                                    reshuffle_me = true;
                                                end
                                                
                                            else % assume we have an even nr of trials after we removed objectcatch trials
                                                if abs(diff(n2))<=1 % we don't allow for a difference more than one
                                                    reshuffle_me = false;
                                                    break;
                                                else
                                                    reshuffle_me = true;
                                                end
                                            end
                                        end
                                        
                                    elseif ismember(curr_sc,5) && ismember(curr_tc,9) % NS-WHERE
                                        
                                        if all(diff(n)==0) % we want equal distribution of L/C/R
                                            reshuffle_me = false;
                                            break % hurray
                                        elseif rem(length(st_idx),3)>0
                                            % but if that is not possible,
                                            % we can be 1 off.
                                            if sum(abs(diff(n)))==1
                                                reshuffle_me = false;
                                                break % hurray
                                            elseif trial_to_condition_ratio == 1 && sum(catch_idx) > 0
                                                % if we sampled all unique conditions AND have catch trials,
                                                % then we will likely never reach a solution unless we allow for a small difference in response distribution.
                                                if (all(diff(n)<=sum(catch_idx)))
                                                    reshuffle_me = false;
                                                    break;
                                                else
                                                    reshuffle_me = true;
                                                end
                                            else
                                                reshuffle_me = true;
                                            end
                                        end
                                    elseif ismember(curr_sc,[4,5]) && ismember(curr_tc,10) % OBJ-HOW & NS-HOW
                                        affordcat = condition_master.affordance_cat(st_idx,:);
                                        affordcat = affordcat(~catch_idx,:);
                                        if curr_sc == 4 % OBJ-HOW
                                            cued0 = condition_master.is_cued(st_idx);
                                            cued0 = cued0(~catch_idx);
                                            n1    = histcounts([affordcat(cued0==1,1);affordcat(cued0==2,2)],1:(nr_responses+1));
                                            
                                        elseif curr_sc == 5 % NS-HOW
                                            n1     = histcounts(affordcat,1:(nr_responses+1));
                                        end
                                        
                                        if (all(n==n1))
                                            reshuffle_me = false;
                                            break % hurray
                                        elseif trial_to_condition_ratio == 1 && sum(catch_idx) > 0
                                            % if we sampled all unique conditions AND have catch trials,
                                            % then we will likely never reach a solution unless we allow for a small difference in response distribution.
                                            if (all(diff(n)<=sum(catch_idx)))
                                                reshuffle_me = false;
                                                break;
                                            else
                                                reshuffle_me = true;
                                            end
                                        else
                                            reshuffle_me = true;
                                        end
                                        
                                    elseif ismember(curr_sc,99) || ismember(curr_tc,3) % scc-all
                                        % Check if we sample all core images across trials
                                        stmclass0 = condition_master.stim_class_name(st_idx,:);
                                        stmclass0 = stmclass0(~catch_idx,:);
                                        cued0     = condition_master.is_cued(st_idx);
                                        cued0     = cued0(~catch_idx);
                                        stmclass0_cued      = [stmclass0(cued0==1,1);stmclass0(cued0==2,2)];
                                        [~,stmclass_cued_i] = ismember(stmclass0_cued,params.exp.stimclassnames([1,3,2,4]));
                                        
                                        n1 = histcounts(stmclass_cued_i,1:(nr_responses+1));
                                        
                                        if (all(n==n1))
                                            reshuffle_me = false;
                                            break % hurray
                                        elseif trial_to_condition_ratio == 1 && sum(catch_idx) > 0
                                            % if we sampled all unique conditions AND have catch trials,
                                            % then we will likely never reach a solution unless we allow for a small difference in response distribution.
                                            if (all(diff(n)<=sum(catch_idx)))
                                                reshuffle_me = false;
                                                break;
                                            else
                                                reshuffle_me = true;
                                            end
                                        else
                                            reshuffle_me = true;
                                        end
                                        
                                    elseif all(n==0) && curr_tc~=2 && ~ismember(curr_tc,[8,9,10])
                                        error('[%s]: No correct responses found?!',mfilename);
                                    else
                                        if (all(diff(n)==0)) % check if we have balanced button presses
                                            reshuffle_me = false;
                                            break % hurray
                                        elseif trial_to_condition_ratio == 1 && sum(catch_idx) > 0
                                            % if we sampled all unique conditions AND have catch trials,
                                            % then we will likely never reach a solution unless we allow for a small difference in response distribution.
                                            if (all(diff(n)<=sum(catch_idx)))
                                                reshuffle_me = false;
                                                break;
                                            else
                                                reshuffle_me = true;
                                            end
                                        else
                                            reshuffle_me = true;
                                        end
                                        
                                    end
                                end
                            end
                            % Updated allocated trials to this session
                            allocated_trials(ti(potential_trials(1:nr_trials))) = 1;
                            
                            % Check if any stimulus numbers are the same
                            if curr_tc==2
                                assert(isequalwithequalnans(condition_master.stim_nr_left(st_idx,:),condition_master.stim_nr_left(ti(potential_trials(1:nr_trials)),:)))
                                assert(isequalwithequalnans(condition_master.stim_nr_right(st_idx,:),condition_master.stim_nr_right(ti(potential_trials(1:nr_trials)),:)))
                            else
                                assert(isequal(condition_master.stim_nr_left(st_idx,:),condition_master.stim_nr_left(ti(potential_trials(1:nr_trials)),:)))
                                if curr_sc~=5
                                    assert(isequal(condition_master.stim_nr_right(st_idx,:),condition_master.stim_nr_right(ti(potential_trials(1:nr_trials)),:)))
                                end
                            end
                            % Check if nr of trials are the same
                            assert(isequal(st_idx,ti(potential_trials(1:nr_trials))));
                            
                            % Check if left/right cueing is balanced
                            cue_tmp = condition_master.is_cued(ti(potential_trials(1:nr_trials)));
                            for xx = 1:nr_blocks
                                xi = ((xx-1)*params.exp.nr_trials_per_block(sc,tc)) + [1:params.exp.nr_trials_per_block(sc,tc)];
                                if ~all(cue_tmp==3)
                                    assert(diff(histcounts(cue_tmp(xi),[1:3]))==0)
                                end
                            end
                            
                            % Apply new trial order
                            condition_master(st_idx,:) = condition_master(ti(potential_trials(1:nr_trials)),:);
                            
                            % Update unique block nr (as trials may have been shuffled)
                            if tc ==3
                                unique_block_nr = sum(all_stimtask_counter(:,tc)) + (1:nr_blocks);
                            else
                                unique_block_nr = all_stimtask_counter(sc,tc) + (1:nr_blocks);
                            end
                            condition_master.stim_class_unique_block_nr(st_idx,:) = repelem(unique_block_nr, params.exp.nr_trials_per_block(sc,tc))';
                            condition_master.unique_trial_nr(st_idx,:)         = potential_trials(1:nr_trials);
                            condition_master.trial_nr(st_idx,:)                = repmat(1:params.exp.nr_trials_per_block(sc,tc), 1, nr_blocks)';
                            
                            % Update crossing counter
                            all_stimtask_counter(sc,tc) = all_stimtask_counter(sc,tc) + nr_blocks;
                        end
                    end
                    
                end
            end
        end
    end
end

% update unused trials
unique_stim_classes = unique(condition_master.stim_class);
unique_task_classes = unique(condition_master.task_class);

leftover_trials = (isnan(condition_master.unique_trial_nr));

for curr_sc = 1:(length(unique_stim_classes)-1)
    for curr_tc = 1:length(unique_task_classes)
        if unique_task_classes(curr_tc) == 3 && unique_stim_classes(curr_sc) < 5 % SCC
            all_trials = (condition_master.stim_class==99 & condition_master.task_class==unique_task_classes(curr_tc));
        else
            all_trials = (condition_master.stim_class==unique_stim_classes(curr_sc) & condition_master.task_class==unique_task_classes(curr_tc));
        end
        unused_idx    = find(all_trials & leftover_trials);
        allocated_idx = find(all_trials & ~leftover_trials);
        
        if curr_tc == 1 % FIX
            assert(isequal(unused_idx,find(all_trials)))
            
            % update stim_class_unique_block_nr
            nr_blocks = length(unused_idx)/params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc));
            
            if ~isnan(nr_blocks)
                if mod(nr_blocks,1)~=0
                    tmp = repelem(1:nr_blocks,params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc)));
                    tmp = cat(2,tmp, repmat(tmp(end)+1, 1, mod(nr_blocks,1)*params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc))));
                else
                    tmp = repelem(1:nr_blocks,params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc)));
                end
                
                if length(tmp) > length(unused_idx) % trim if needed
                    tmp = tmp(1:length(unused_idx))';
                elseif length(unused_idx) > length(tmp)
                    error('wtf')
                end
                condition_master.stim_class_unique_block_nr(unused_idx) = tmp;
                
                % update unique_trial_nr
                nr_unique_condition_repeats = length(unused_idx)/params.exp.nr_unique_trials_per_crossing(unique_stim_classes(curr_sc),unique_task_classes(curr_tc));
                tmp = repmat(1:params.exp.nr_unique_trials_per_crossing(unique_stim_classes(curr_sc),unique_task_classes(curr_tc)),1,ceil(nr_unique_condition_repeats));
                if length(tmp) > length(unused_idx) % trim if needed
                    tmp = tmp(1:length(unused_idx))';
                elseif length(unused_idx) > length(tmp)
                    error('wtf')
                end
                condition_master.unique_trial_nr(unused_idx) = tmp;
                
                % update trial_nr
                tmp = repmat(1:params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc)),1,ceil(nr_blocks));
                if length(tmp) > length(unused_idx) % trim if needed
                    tmp = tmp(1:length(unused_idx))';
                elseif length(unused_idx) > length(tmp)
                    error('wtf')
                end
                condition_master.trial_nr(unused_idx) = tmp;
            end
        else
            
            start_block_nr = max(condition_master.stim_class_unique_block_nr(allocated_idx)) + 1;
            
            % update stim_class_unique_block_nr
            nr_blocks = length(unused_idx)/params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc));
            
            if ~isnan(nr_blocks)
                if mod(nr_blocks,1)~=0
                    tmp = repelem(start_block_nr:(start_block_nr+nr_blocks-1),params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc)));
                    tmp = cat(2,tmp, repmat(tmp(end)+1, 1, mod(nr_blocks,1)*params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc))));
                else
                    tmp = repelem(start_block_nr:(start_block_nr+nr_blocks-1),params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc)));
                end
                
                if length(tmp) > length(unused_idx) % trim if needed
                    tmp = tmp(1:length(unused_idx))';
                elseif length(unused_idx) > length(tmp)
                    error('wtf')
                end
                condition_master.stim_class_unique_block_nr(unused_idx) = tmp;
                
                % update unique_trial_nr
                nr_unique_condition_repeats = length(unused_idx)/params.exp.nr_unique_trials_per_crossing(unique_stim_classes(curr_sc),unique_task_classes(curr_tc));
                tmp = repmat(1:params.exp.nr_unique_trials_per_crossing(unique_stim_classes(curr_sc),unique_task_classes(curr_tc)),1,ceil(nr_unique_condition_repeats));
                if length(tmp) > length(unused_idx) % trim if needed
                    tmp = tmp(1:length(unused_idx))';
                elseif length(unused_idx) > length(tmp)
                    error('wtf')
                end
                condition_master.unique_trial_nr(unused_idx) = tmp;
                
                % update trial_nr
                tmp = repmat(1:params.exp.nr_trials_per_block(unique_stim_classes(curr_sc),unique_task_classes(curr_tc)),1,ceil(nr_blocks));
                if length(tmp) > length(unused_idx) % trim if needed
                    tmp = tmp(1:length(unused_idx))';
                elseif length(unused_idx) > length(tmp)
                    error('wtf')
                end
                condition_master.trial_nr(unused_idx) = tmp;
            end
        end
    end
end

fprintf('Done!\n')

return