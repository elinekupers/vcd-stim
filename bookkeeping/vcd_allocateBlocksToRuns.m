function condition_master = vcd_allocateBlocksToRuns(params,condition_master_in,env_type, varargin)
% VCD function to allocate the unique trials and blocks to given runs.
%
%   condition_master = vcd_allocateBlocksToRuns(params)
%
% INPUTS:
%   params                :  (struct) parameter struct needed to get subject
%                             nrs and run type params. REQUIRES params.trials
%                             to exist and contain condition_master v0.
%   condition_master_in   :  (struct) input condition_master table with
%                             block and runs for each subject and session.
%   env_type              :  (str) label to define what type of session we
%                             are defining: 'MRI' or 'BEHAVIOR'.
%
% OUTPUTS:
%   condition_master      :  (struct) updated condition master table with
%                             block and runs for each subject and session.

%% Parse inputs
p0 = inputParser;
p0.addRequired('params'                 , @isstruct);
p0.addRequired('condition_master'       , @istable);
p0.addRequired('env_type'               , @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));
p0.addParameter('load_params'           , false, @islogical);
p0.addParameter('store_params'          , true, @islogical);
p0.addParameter('store_imgs'            , false, @islogical);
p0.addParameter('verbose'               , false, @islogical);

% Parse inputs
p0.parse(params,condition_master_in,env_type,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

%% copy condition master
condition_master = condition_master_in;

% Add extra columns to master table
if sum(strcmp(condition_master.Properties.VariableNames,'session_nr'))==0
    t_session = table(...
        NaN(size(condition_master,1),1), ... session_nr
        NaN(size(condition_master,1),1), ... session_type
        NaN(size(condition_master,1),1), ... run_nr
        NaN(size(condition_master,1),1), ... block_nr
        NaN(size(condition_master,1),1), ... global_run_nr
        NaN(size(condition_master,1),1), ... global_block_nr (resets per subject)
        NaN(size(condition_master,1),1), ... global_trial_nr (resets per subject)
        NaN(size(condition_master,1),2), ... condition_nr
        num2cell(NaN(size(condition_master,1),2)), ... condition_name
        NaN(size(condition_master,1),1), ... crossing_nr
        num2cell(NaN(size(condition_master,1),1)), ... crossing_name
        NaN(size(condition_master,1),1)); % trial_type

    t_session.Properties.VariableNames = {...
        'session_nr',...
        'session_type',...
        'run_nr',...
        'block_nr',...
        'global_run_nr',...
        'global_block_nr', ...
        'global_trial_nr', ...
        'condition_nr', ...
        'condition_name', ...
        'crossing_nr',...
        'crossing_name',...
        'trial_type',...
        };
    
    condition_master = cat(2, t_session, condition_master);
end


%% Preallocate counters/matrices

% Keep track of allocated blocks for each stim-task crossing and continue
% counting to avoid overwriting of trials
stimtask_tracker_global = ones(length(params.exp.stimclassnames), length(params.exp.taskclassnames));

% change stim/task labels to upper case and abbreviated version
stim_class_abbr = upper(params.exp.stimclassnames);
stim_class_abbr{1} = 'GBR';
task_class_abbr = upper(params.exp.taskclassnames);

% check session type
[all_sessions,session_types] = vcd_getSessionEnvironmentParams(params, env_type);

% create separate counters for scc and ltm task crossings, because trials
% are spread across stimulus classes
scc_bb = 0; ltm_bb = 0;

for ses = 1:size(all_sessions,3)
    for st = 1:size(all_sessions,4) % session types
        if ~isnan(session_types(ses,st))
            %% First we check how many repeats of stim-task crossings we want to run within a session
            if verbose
               fprintf('%sSESSION %03d, session_type %02d..\n',choose(params.is_demo, 'DEMO ',''), ses,session_types(ses,st));
            end
            
            % also keep track of local stim-task block allocation
            stimtask_tracker_local = ones(length(params.exp.stimclassnames), length(params.exp.taskclassnames));
            
            % Nr of blocks for each stim-task crossing. We use the nr of blocks,
            % same trials, and same unique image for each subject
            block_distr = all_sessions(:,:,ses,st);

            % Reset block counter within a session
            % Note: we continue counting scc and ltm blocks instead of relying on stimtask_tracker_global
            bb = 1;  half_block_counter = 0;
            run_nr = 1;
            
            for sc = 1:length(params.exp.stimclassnames)
                
                % check how many half blocks we have for this stim class?
                half_block_count = sum(mod(block_distr(sc,:),1)<1 & mod(block_distr(sc,:),1)>0);
                full_block_count = sum(mod(block_distr(sc,:),1)==0);
                
                for tc = 1:length(params.exp.taskclassnames)
                    
                    if strcmp(env_type,'MRI')
                        task_start = params.exp.session.mri.task_start(tc);
                        
                    elseif strcmp(env_type,'BEHAVIOR')
                        task_start = params.exp.session.behavior.task_start(tc);
                    end
                    if  ses >= task_start
                        % get nr of blocks for this specific crossing
                        n_blocks = block_distr(sc,tc);
                        
                            
                        % if it isn't zero, then add it to the table
                        if (n_blocks > 0)
                            
                            % get nr of blocks we want to allocate
                            if n_blocks < 1
                                curr_blocks = [stimtask_tracker_global(sc,tc)+n_blocks-1];
                            else 
                                curr_blocks = [stimtask_tracker_global(sc,tc) : (stimtask_tracker_global(sc,tc)+n_blocks-1)];
                            end
                            
                            for ii = 1:length(curr_blocks)
                                
                                if ~isint(curr_blocks(ii))
                                    half_block_flag = true;
                                else
                                    half_block_flag = false;
                                end

                                if ismember(params.exp.taskclassnames{tc},'scc')
                                    % We can't sort on stim class nr for SCC
                                    % otherwise stim will be allocated to different
                                    % blocks. Hence we use an if statement and call
                                    % the crossing_name 'scc-all'
                                    idx = ((condition_master.stim_class==99) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr== scc_bb + find(curr_blocks==curr_blocks(ii))));
                                    condition_master.crossing_name(idx) = {sprintf('%s-all',params.exp.taskclassnames{tc})};
                                    stim_class_tmp_name = 'ALL';
                                elseif ismember(params.exp.taskclassnames{tc},'ltm')
                                    % We can't sort on stim class nr for LTM either (they are all defined as 99)
                                    idx = ((condition_master.stim_class==99) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr == ltm_bb + find(curr_blocks==curr_blocks(ii))));
                                    condition_master.crossing_name(idx) = {sprintf('%s-all',params.exp.taskclassnames{tc})};
                                    stim_class_tmp_name = 'ALL';
                                else
                                    idx = ((condition_master.stim_class==sc) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr == curr_blocks(ii)));
                                    condition_master.crossing_name(idx) = {sprintf('%s-%s',params.exp.taskclassnames{tc}, params.exp.stimclassnames{sc})};
                                    stim_class_tmp_name = stim_class_abbr{sc};
                                end
                                if verbose
                                    if ii == 1
                                        fprintf('\n %s   \t: %3.1f block(s)\n',cell2mat(unique(condition_master.crossing_name(find(idx)))),n_blocks)
                                    end
                                    fprintf(' %s   \t: Found %d trials, total stim-task trials in condition master = %d, %d have already been allocated\n', ...
                                        cell2mat(unique(condition_master.crossing_name(find(idx)))), sum(idx), length(condition_master.stim_class_unique_block_nr( (condition_master.stim_class==choose(strcmp(stim_class_tmp_name,'ALL'),99,sc)) & (condition_master.task_class==tc))), ...
                                        length(condition_master.stim_class_unique_block_nr(~isnan(condition_master.session_nr) & (condition_master.stim_class==choose(strcmp(stim_class_tmp_name,'ALL'),99,sc)) & (condition_master.task_class==tc))));
                                end
                                sub_idx   = find(idx);
                                nr_trials = length(sub_idx);
                                
                                if half_block_flag
                                    % ensure we sample trials equally wrt
                                    % left/right cue and stimclasses
                                    cue_types  = unique(condition_master.is_cued(sub_idx));
                                    stim_clsss = condition_master.stim_class_name(sub_idx,:);
                                    nr_trials_half_block = (nr_trials*curr_blocks(ii));
                                    
                                        if all(ismember(cue_types,[1,2]))
                                            cue_l_idx = find(condition_master.is_cued(sub_idx)==1);
                                            cue_r_idx = find(condition_master.is_cued(sub_idx)==2);
                                            
                                            if tc == 3 % scc trials
                                                
                                                [~,stm_idx] = ismember(stim_clsss, params.exp.stimclassnames);
                                                nn_stim_clsss = histcounts(stm_idx);
                                                sub_idx0 = [];
                                                for mm = 1:4
                                                    stmclss_idx(:,:,mm) = (stm_idx==mm);
                                                end
                                            
                                                [~,mm_order] = sort(nn_stim_clsss,'ascend');
                                                
                                                cue_counter = 1;
                                                while cue_counter <= nr_trials_half_block
                                                    % alternate stim location (left or right)
                                                    tmp_loc = mod(cue_counter-1,2)+1;  % either 1:left or 2:right
                                                    
                                                    % find trial that is most diverse in stimulusclass
                                                    tmp_trial_unique = {}; chosen_trial = [];
                                                    for mmm = mm_order
                                                        chosen_trial0 = find(stmclss_idx(:,tmp_loc,mmm));
                                                        if ~isempty(chosen_trial0) % if we get no trials, go to the next stim class..
                                                            
                                                            if length(chosen_trial0)>1 % if we happen to get multipe trials, then find the trial with 2 stim classes
                                                                tmp_trial_unique = abs(diff(stm_idx(stmclss_idx(:,tmp_loc,mmm),:),[],2))>0;
                                                                if ~isempty(tmp_trial_unique)
                                                                    [~,tmp_trial] = max(tmp_trial_unique);
                                                                    if length(tmp_trial)>1
                                                                        chosen_trial_idx = max(tmp_trial);
                                                                    else
                                                                        chosen_trial_idx = tmp_trial;
                                                                    end
                                                                else
                                                                    chosen_trial_idx = NaN;
                                                                end
                                                            else
                                                                chosen_trial_idx = 1;
                                                            end

                                                            if ~isnan(chosen_trial_idx) && chosen_trial_idx~=0
                                                                chosen_trial = chosen_trial0(chosen_trial_idx);
                                                                break;
                                                            end
                                                        end
                                                    end
                                                    
                                                    if ~isempty(chosen_trial) % if we found one, add it to the list
                                                        sub_idx0 = cat(1,sub_idx0,chosen_trial);
                                                    elseif isempty(chosen_trial) % if we can't find any, then just take another trial from a stimclass we haven't sampled yet
                                                        leftover_trials = setdiff([1:4]',sub_idx0);
                                                        sub_idx0 = cat(1,sub_idx0,leftover_trials(1));
                                                    end
                                                    % remove that trial from list, to avoid double sampling
                                                    stmclss_idx(sub_idx0(end),:,:) = false;
                                                    
                                                    % get new priority list of stim classes to sample from
                                                    mm_order = setdiff([1:4], reshape(stm_idx(sub_idx0,:),1,[]));
                                                    if isempty(mm_order)
                                                        [~,mm_order] = sort(histcounts(reshape(stm_idx(sub_idx0,:),[],1)),'ascend');
                                                    end
                                                    cue_counter = cue_counter+1;
                                                end
                                            end
                                            % Otherwise we just go by cue condition (equally sampling
                                            % left and right cueing conditons)
                                            sub_idx0 = [cue_l_idx(1:(nr_trials_half_block/2)); cue_r_idx(1:(nr_trials_half_block/2))];
                                        else
                                            cue_c_idx = find(condition_master.is_cued(sub_idx)==3);
                                            sub_idx0 = cue_c_idx(1:(nr_trials_half_block/2));
                                        end
                                    
                                    % find trial indices for selected trials
                                    sub_idx1 = sub_idx(sort(sub_idx0));
                                    sub_idx  = sub_idx1; clear sub_idx1;
                                    idx(setdiff(find(idx),sub_idx))=false;
                                    nr_trials = length(sub_idx);
                                    assert(isequal(nr_trials,nr_trials_half_block));
                                end
                                
                                % Check if these rows in the condition master are already occupied, we want them to be empty (otherwise we would overwrite the trial)!
                                assert(all(isnan(condition_master.session_nr(idx))))
                                assert(all(isnan(condition_master.session_type(idx))))
                                assert(all(isnan(condition_master.block_nr(idx))))
                                assert(all(isnan(condition_master.global_block_nr(idx))))
                                assert(isequal(sum(condition_master.is_cued(idx)==1),sum(condition_master.is_cued(idx)==2)))
                                
                                % get block numbers for this stim-task crossing
                                condition_master.session_nr(idx)        = ses;
                                condition_master.session_type(idx)      = session_types(ses,st);
                                blck_name                               = condition_master.crossing_name(idx);
                                condition_master.crossing_name(idx)     = blck_name(1);
                                condition_master.crossing_nr(idx)       = find(strcmp(blck_name(1),params.exp.crossingnames));
                                if half_block_flag
                                    condition_master.block_nr(idx)      = curr_blocks(ii);
                                else
                                    condition_master.block_nr(idx)      = bb;
                                end
                                condition_master.run_nr(idx)            = run_nr; % put in temporary run nr for now.
                                
                                % trial type (single/double-stim
                                % presentation)
                                if params.exp.trial.double_epoch_tasks(tc)
                                    condition_master.trial_type(idx) = 2; % double epoch;
                                else
                                    condition_master.trial_type(idx) = 1; % double epoch;
                                end
                                
                                % get condition label for each trial and stim loc
                                
                                for tt = 1:nr_trials
                                    
                                    if condition_master.is_catch(sub_idx(tt))==1
                                            % Left-side general catch trials
                                        if condition_master.is_cued(sub_idx(tt))==1
                                            cue_label = {sprintf('%s 0000 X LCUED %s#',stim_class_tmp_name,task_class_abbr{tc}), ...
                                                sprintf('%s 0000 X LCUED %s#',stim_class_tmp_name,task_class_abbr{tc})};
                                            % Right-side general catch trials
                                        elseif condition_master.is_cued(sub_idx(tt))==2
                                            cue_label = {sprintf('%s 0000 X RCUED %s#',stim_class_tmp_name,task_class_abbr{tc}), ...
                                                sprintf('%s 0000 X RCUED %s#',stim_class_tmp_name,task_class_abbr{tc})};
                                            % Center or FIX general catch trials
                                        elseif condition_master.is_cued(sub_idx(tt))==3
                                            if strcmp(stim_class_tmp_name,'NS')
                                                cue_label = {sprintf('%s 0000 X NCUED %s#',stim_class_tmp_name,task_class_abbr{tc}), NaN};
                                            else
                                                cue_label = {sprintf('%s 0000 X NCUED %s#',stim_class_tmp_name,task_class_abbr{tc}), ...
                                                    sprintf('%s 0000 X NCUED %s#',stim_class_tmp_name,task_class_abbr{tc})};
                                            end
                                        else
                                            error('[%s]: Can''t tell what type of trial this is..',mfilename)
                                        end
                                    elseif condition_master.is_cued(sub_idx(tt))==1
                                        % Left-side object catch trials
                                        if ~isnan(condition_master.is_objectcatch(sub_idx(tt))) && condition_master.is_objectcatch(sub_idx(tt))>0
                                            cue_label = {sprintf('%s %04d L CUED %s+',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}), ...
                                                sprintf('%s %04d R UNCUED %s',stim_class_tmp_name,condition_master.stim_nr_right(sub_idx(tt)),task_class_abbr{tc})};
                                            condition_master.stim_nr_left(sub_idx(tt)) = condition_master.is_objectcatch(sub_idx(tt)); % replace image nr with unique im nr associated with objectcatch
                                            condition_master.is_objectcatch(sub_idx(tt)) = 1; % Now we reset is_objectcatch to 1
                                            % Left-side NON-object catch trials
                                        elseif isnan(condition_master.is_objectcatch(sub_idx(tt))) || condition_master.is_objectcatch(sub_idx(tt))==0
                                            cue_label = {sprintf('%s %04d L CUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}), ...
                                                sprintf('%s %04d R UNCUED %s',stim_class_tmp_name,condition_master.stim_nr_right(sub_idx(tt)),task_class_abbr{tc})};
                                        else
                                            error('[%s]: Can''t tell what type of trial this is..',mfilename)
                                        end
                                       
                                        
                                    elseif condition_master.is_cued(sub_idx(tt))==2
                                           % Right-side object catch trials
                                        if ~isnan(condition_master.is_objectcatch(sub_idx(tt))) && condition_master.is_objectcatch(sub_idx(tt))>0
                                            cue_label = {sprintf('%s %04d L UNCUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}),...
                                            sprintf('%s %04d R CUED %s+',stim_class_tmp_name,condition_master.stim_nr_right(sub_idx(tt)),task_class_abbr{tc})};
                                            condition_master.stim_nr_right(sub_idx(tt)) = condition_master.is_objectcatch(sub_idx(tt)); % replace image nr with unique im nr associated with objectcatch
                                            condition_master.is_objectcatch(sub_idx(tt)) = 1; % Now we reset is_objectcatch to 1
                                            % Right-side NON-object catch trials
                                        elseif isnan(condition_master.is_objectcatch(sub_idx(tt)))  || condition_master.is_objectcatch(sub_idx(tt))==0
                                            cue_label = {sprintf('%s %04d L UNCUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}),...
                                            sprintf('%s %04d R CUED %s',stim_class_tmp_name,condition_master.stim_nr_right(sub_idx(tt)),task_class_abbr{tc})};
                                        else
                                            error('[%s]: Can''t tell what type of trial this is..',mfilename)
                                        end
                                        
                                    elseif condition_master.is_cued(sub_idx(tt))==3
                                        
                                        % NS trial - left stim / right is nan
                                        if any(ismember(condition_master.stim_nr_left(sub_idx(tt)),params.stim.ns.unique_im_nrs_core)) && isnan(condition_master.stim_nr_right(sub_idx(tt)))
                                            cue_label = {sprintf('%s %04d C NCUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}), NaN};
                                            
                                            % NS catch trial
                                        elseif condition_master.stim_nr_left(sub_idx(tt))==0 && strcmp(stim_class_tmp_name,'NS')
                                            cue_label = {sprintf('%s %04 C NCUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}), NaN};
                                            
                                        else
                                            % fix trial clssic left/right stim
                                            cue_label = {sprintf('%s %04d L NCUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}), ...
                                                sprintf('%s %04d R NCUED %s',stim_class_tmp_name,condition_master.stim_nr_right(sub_idx(tt)),task_class_abbr{tc})};
                                        end
                                    else
                                        error('[%s]: Can''t tell cue type this is..',mfilename)
                                    end
                                    
                                    % add dash here, otherwise sprintf interprets dash as left align
                                    for cl = 1:length(cue_label)
                                        if isnan(cue_label{cl})
                                            if strcmp(stim_class_tmp_name,'NS') && cl==2
                                                condition_master.condition_name(sub_idx(tt),cl) = {NaN};
                                                condition_master.condition_nr(sub_idx(tt),cl) = NaN;
                                            else
                                                error('[%s]: Check condition name!')
                                            end
                                        else
                                            condition_master.condition_name(sub_idx(tt),cl) = strrep(cue_label(cl),' ', '-');
                                            condition_master.condition_nr(sub_idx(tt),cl)   = vcd_conditionName2Number(condition_master.condition_name(sub_idx(tt),cl));
                                        end
                                    end
                                end
                                
                                % update counters 
                                if half_block_flag
                                    stimtask_tracker_local(sc,tc)  = stimtask_tracker_local(sc,tc) + mod(curr_blocks(ii),1);
                                    stimtask_tracker_global(sc,tc) = stimtask_tracker_global(sc,tc) + mod(curr_blocks(ii),1);
                                    half_block_counter = half_block_counter + curr_blocks(ii);
                                    if half_block_counter == 1
                                        bb = bb + 1;
                                        half_block_counter = 0;
                                    elseif half_block_count == 1 && full_block_count == 1
                                        bb = bb + 1;
                                    end
                                else
                                    stimtask_tracker_local(sc,tc)  = stimtask_tracker_local(sc,tc)+1;
                                    stimtask_tracker_global(sc,tc) = stimtask_tracker_global(sc,tc)+1;
                                    bb = bb + 1;
                                end
                                
                                if mod(bb,8)==0 % every 7 blocks we add a run 
                                    run_nr = run_nr + 1; 
                                end
                            end
                            
                            if strcmp(params.exp.taskclassnames{tc},'scc')
                                scc_bb = scc_bb + length(curr_blocks);
                            elseif strcmp(params.exp.taskclassnames{tc},'ltm')
                                ltm_bb = ltm_bb + length(curr_blocks);
                            end
                        end
                    end
                end
            end
            
            % Merge half blocks
            if any(condition_master.block_nr < 1)
                single_stim_half_blocks = find(~isnan(condition_master.block_nr) & (condition_master.block_nr < 1) & (condition_master.trial_type==1));
                double_stim_half_blocks = find(~isnan(condition_master.block_nr) & (condition_master.block_nr < 1) & (condition_master.trial_type==2));

                nr_comb_single_blocks = length(single_stim_half_blocks)/params.exp.block.n_trials_single_epoch;
                nr_comb_double_blocks = length(double_stim_half_blocks)/params.exp.block.n_trials_double_epoch;
                nr_of_uncomb_single_blocks = length(unique(condition_master.run_nr(single_stim_half_blocks)));
                nr_of_uncomb_double_blocks = length(unique(condition_master.run_nr(double_stim_half_blocks)));

                if nr_comb_single_blocks > 0
                    for nn = 1:params.exp.block.n_trials_single_epoch:length(single_stim_half_blocks)
                        trial_idx = single_stim_half_blocks(nn:(nn+(params.exp.block.n_trials_single_epoch-1)));
                        condition_master.run_nr(trial_idx) = condition_master.run_nr(trial_idx(1));
                        all_block_nrs = unique(condition_master.block_nr(...
                                                condition_master.session_nr==condition_master.session_nr(trial_idx(1)) & ...
                                                condition_master.session_type==condition_master.session_type(trial_idx(1)) & ...
                                                condition_master.run_nr==condition_master.run_nr(trial_idx(1))));
                        all_block_nrs = setdiff(all_block_nrs, condition_master.block_nr(trial_idx));
                        if all(all_block_nrs~=1) && all(all_block_nrs<8)
                            allocated_block_nr = setdiff([1:max(all_block_nrs)],all_block_nrs);
                        else
                            if isequal([min(all_block_nrs):max(all_block_nrs)]',all_block_nrs)
                                unused_block_nrs = setdiff([1:max(condition_master.block_nr(~isnan(condition_master.block_nr)))]', unique(sort(condition_master.block_nr(~isnan(condition_master.block_nr)))));
                                mi0 = min(abs(all_block_nrs-unused_block_nrs'));
                                [~,mi1] = min(mi0);
                                allocated_block_nr = unused_block_nrs(mi1);
                            else
                                allocated_block_nr = setdiff([min(all_block_nrs):max(all_block_nrs)],all_block_nrs);
                            end
                        end
                        if length(allocated_block_nr)>1
                            allocated_block_nr = allocated_block_nr(1);
                        end
                        condition_master.block_nr(trial_idx) = repmat(allocated_block_nr,size(trial_idx,1),1);
                        condition_master.trial_nr(trial_idx) = [1:params.exp.block.n_trials_single_epoch]';
                    end
                    shave_me_off1 = nr_of_uncomb_single_blocks - nr_comb_single_blocks;
                    bb = bb + shave_me_off1;
                end
                
                if nr_comb_double_blocks > 0
                    for nn = 1:params.exp.block.n_trials_double_epoch:length(double_stim_half_blocks)
                        trial_idx = double_stim_half_blocks(nn+(params.exp.block.n_trials_double_epoch-1));
                        condition_master.run_nr(trial_idx)   = condition_master.run_nr(trial_idx(1));
                        condition_master.block_nr(trial_idx) = ceil(condition_master.block_nr(trial_idx(1)));
                        condition_master.trial_nr(trial_idx) = 1:params.exp.block.n_trials_single_epoch;
                    end
                    shave_me_off2 = nr_of_uncomb_double_blocks - nr_comb_double_blocks;
                    bb = bb + shave_me_off2;  % correct nr of blocks
                end
               
                % check if the block nrs are ascending
                unused_block_nrs = setdiff([1:max(condition_master.block_nr(~isnan(condition_master.block_nr)))], ...
                    sort(unique(condition_master.block_nr(~isnan(condition_master.block_nr))))');
                assert(isempty(unused_block_nrs));
            end
            
            if verbose
                fprintf('\nTOTAL MINIBLOCKS: %d for %sSESSION %03d, session_type %02d\n',bb, choose(params.is_demo,'DEMO ',''), ses,session_types(ses,st))
            end
            % the number of blocks in the table should be equal to the number of
            % tracked stim-task-crossings allocated, as well as the sum vs vector length
            assert(isequal(stimtask_tracker_local-1,block_distr))
            tmp = max(condition_master.block_nr(~isnan(condition_master.block_nr) & condition_master.session_nr==ses & condition_master.session_type==st));
            
            assert(isequal(tmp,sum(sum(stimtask_tracker_local-1))))
        end
    end
end

% Remove obsolete columns
condition_master(isnan(condition_master.session_nr),:) = [];
condition_master.stim_class_unique_block_nr = [];
condition_master.unique_trial_nr = [];

% Reorder runs
[~,run_idx] = sort(condition_master.run_nr);

% We keep and reorder block nr for now (we will shuffle blocks later for each subject's session)
[~,block_idx] = sort(condition_master.block_nr(run_idx));

condition_master = condition_master(run_idx(block_idx),:);
    
% Define column order by their names (yes, this is tricky//not ideal)
newOrderNames = {'session_nr','session_type','run_nr','block_nr','trial_nr',...
   'global_run_nr', 'global_block_nr','global_trial_nr',...
    'condition_nr','condition_name', ...
    'stim_class','stim_class_name', ...
    'task_class','task_class_name', ...
    'crossing_nr','crossing_name',...
    'stim_nr_left', 'stim_nr_right', ...
    'is_cued', 'is_catch', 'is_objectcatch',...
    'correct_response', ...
    'orient_dir', ...
    'contrast', ...
    'gbr_phase', ...
    'rdk_coherence', ...
    'super_cat', 'super_cat_name', ...
    'basic_cat', 'basic_cat_name', ...
    'sub_cat', 'sub_cat_name', ...
    'affordance_cat','affordance_name', ...
    'cd_start', ...
    'stim2_delta',...
    'stim2_im_nr',...
    'stim2_orient_dir',...
    'is_special_core', 'is_lure', ...
    'repeat_nr', 'trial_type'};

% check if we missed any names:
assert(isequal(sort(newOrderNames),sort(condition_master.Properties.VariableNames)));

% Find the right column-name match
newOrder_idx = [];
for fn = 1:length(newOrderNames)
    [~,col_idx] = ismember(condition_master.Properties.VariableNames, newOrderNames{fn});
    newOrder_idx(fn) = find(col_idx);
end

% apply reodering!
condition_master = condition_master(:,newOrder_idx);


% IF WE DEAL WITH A DEMO SESSION, DELETE 50% OF TRIALS
if params.is_demo
    curr_sessions = unique(condition_master.session_nr);
    for ss = 1:length(curr_sessions)
        curr_runs = unique(condition_master.run_nr(condition_master.session_nr == curr_sessions(ss)));
        
        for rr = 1:length(curr_runs)
            delete_me = [];
            session_trials = (condition_master.session_nr== curr_sessions(ss)  &  condition_master.run_nr==curr_runs(rr)); % find trials in one run, for a given session
            session_sub    = find(session_trials);
            blocks         = unique(condition_master.block_nr(session_trials)); % how many unique blocks do we have a given session?
            for bb = 1:length(blocks)
                block_idx = condition_master.block_nr(session_trials) == blocks(bb);  % how many trials do we have in a given blovk?
                block_sub = find(block_idx);
                delete_me = cat(1, delete_me, block_sub((round(length(block_sub)/2)+1):length(block_sub))); % keep the first half of trials.
            end
            condition_master(session_sub(delete_me),:) = [];
        end
    end
end
  
% Update trial repeat nr
condition_master = vcd_getTrialRepeatNr(condition_master);

% add global_trial_nr
condition_master = vcd_updateGlobalCounters(params, condition_master, env_type);
% condition_master.global_trial_nr = [1:size(condition_master,1)]';
% condition_master.block_nr = condition_master.global_block_nr;

%Visualize results
if verbose
    figure; set(gcf,'Position',[1 1 1200 600]);
    clims = [0, max([max(condition_master.block_nr),max(condition_master.run_nr)])];
    cmap = cmapturbo(max(clims)+1);                
    close all; makeprettyfigures
    for ses = 1:size(all_sessions,3)
        for st = 1:size(all_sessions,4) % session types
            
            if ~isnan(session_types(ses,st))
                
                % visualize blocks and trials, now after shuffle
                ses_idx = (condition_master.session_nr==ses & condition_master.session_type==st);
                ses_data = condition_master(ses_idx,:);
                [~,new_run_line] = unique(ses_data.run_nr,'stable');
                dataToPlot = [ses_data.session_nr, ses_data.run_nr, ses_data.block_nr]';
                imagesc(dataToPlot);
                colormap(cmap);
                set(gca,'CLim',clims);
                cb = colorbar;

                hold on;
                for ff = 1:length(new_run_line)
                    plot([new_run_line(ff),new_run_line(ff)],[0,6],'k','linewidth',4)
                end
                
                set(gca,'YTick',[1:params.exp.total_subjects+2],'YTickLabel',{'session nr','run nr', 'block nr'});
                title(sprintf('SESSION %02d %s OVERVIEW',ses, choose(st==1, 'A','B')),'FontSize',20);
                xlabel('BLOCK NR')
                
                
                if store_imgs
                    saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master1_%s',env_type));
                    if ~exist(saveFigsFolder,'dir'); mkdir(saveFigsFolder);end
                    filename = sprintf('vcd_session%02d_%s_subjblocks_post_shuffle.png',ses,choose(st==1, 'A','B'));
                    print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                end
            end
        end
    end
end

return

%     % visualize blocks and trials
%     cmap = cmapturbo(500);
%     figure; set(gcf,'Position',[1 1 1200 300]);
%     colormap(cmap); sgtitle(sprintf('Session %02d',ses))
%     new_block_line = find(abs(diff(condition_master.block_nr(~isnan(condition_master.block_nr))))>1);
%     subplot(311); cla;
%     imagesc(condition_master.block_nr(~isnan(condition_master.block_nr))');
% %         hold on;
% %     for ff = 1:length(new_block_line)
% %         plot([new_block_line(ff),new_block_line(ff)]+0.05,[0,4],'k','linewidth', 0.01)
% %     end
%     set(gca,'YTick',1,'YTickLabel',{'trial nr'})
%     subplot(312);cla;
%     imagesc(condition_master.block_nr(~isnan(condition_master.block_nr))');
%
% %         hold on;
% %     for ff = 1:length(new_block_line)
% %         plot([new_block_line(ff),new_block_line(ff)]+0.05,[0,4],'k','linewidth', 0.01)
% %     end
%     set(gca,'YTick',1,'YTickLabel',{'local block nr'})
%     subplot(313)
%     imagesc(condition_master.block_nr(~isnan(condition_master.block_nr))');
% %     hold on;
% %     for ff = 1:length(new_block_line)
% %         plot([new_block_line(ff),new_block_line(ff)]+0.05,[0,4],'k','linewidth', 0.01)
% %     end
%     set(gca,'YTick',1,'YTickLabel',{' block nr'})
%
%     if store_imgs
%         saveFigsFolder = fullfile(vcd_rootPath,'figs');
%         filename = sprintf('vcd_session%02d_blocks_pre_shuffle.png', ses);
%         print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
%     end

