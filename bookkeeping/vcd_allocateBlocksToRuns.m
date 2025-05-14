function condition_master = vcd_allocateBlocksToRuns(params,condition_master_in,session_type)
% VCD function to allocate the unique trials and blocks to given runs.
%
%   condition_master = vcd_allocateBlocksToRuns(params)
%
% INPUTS:
%   params              : 	(struct) parameter struct needed to get subject
%                           nrs and run type params. REQUIRES params.trials
%                           to exist and contain condition_master v0.
%   condition_master    :  (struct) condition master table with
%                           block and runs for each subject and session.
%   session_type        :   (str) label to define what type of session we
%                           are defining: 'MRI' or 'BEHAVIOR'.
%
% OUTPUTS:
%   condition_master    :  (struct) updated condition master table with
%                           block and runs for each subject and session.

%%

% copy condition master
condition_master =condition_master_in;

% Add extra columns to master table
if sum(strcmp(condition_master.Properties.VariableNames,'session_nr'))==0
    t_session = table( ...
        NaN(size(condition_master,1),1), ... session_nr
        NaN(size(condition_master,1),1), ... session_type
        NaN(size(condition_master,1),1), ... run_nr
        NaN(size(condition_master,1),1), ... block_nr
        NaN(size(condition_master,1),1), ...
        num2cell(NaN(size(condition_master,1),2)), ... condition_name
        NaN(size(condition_master,1),2), ... condition_nr
        NaN(size(condition_master,1),1), ... crossing_nr
        num2cell(NaN(size(condition_master,1),1)), ... crossing_name
        NaN(size(condition_master,1),1), ... trial_type
        NaN(size(condition_master,1),1), ... %  global_session_nr (resets per subject)
        NaN(size(condition_master,1),1), ... %  global_block_nr (resets per subject)
        NaN(size(condition_master,1),1));    %  global_trial_nr (resets per subject)
    t_session.Properties.VariableNames = {...
        'session_nr',...
        'session_type',...
        'run_nr',...
        'block_nr',...
        'subj_block_nr', ...
        'condition_name', ...
        'condition_nr', ...
        'crossing_nr',...
        'crossing_name',...
        'trial_type',...
        'global_session_nr', ...
        'global_block_nr', ...
        'global_trial_nr', ...
        };
    
    condition_master = cat(2, t_session(:,[1:7]), condition_master,t_session(:,[8,9]));
end


%% Preallocate counters/matrices

% Keep track of allocated blocks for each stim-task crossing and continue
% counting to avoid overwriting of trials
stimtask_tracker_global = ones(length(params.exp.stimclassnames), length(params.exp.taskclassnames));

% keep track of all the sessions
global_session_counter = 0;
global_block_counter   = 0;
global_trial_counter   = 0;

% change stim/task labels to upper case and abbreviated version
stim_class_abbr = upper(params.exp.stimclassnames);
stim_class_abbr{1} = 'GBR';
task_class_abbr = upper(params.exp.taskclassnames);

% check session type
if strcmp(session_type,'MRI')
    all_sessions     = cat(3, params.exp.session.wide.ses_blocks, params.exp.session.deep.ses_blocks);
    runs_per_session = cat(1, size(params.exp.session.wide.ses_blocks,3),size(params.exp.session.deep.ses_blocks,3));
    session_types    = cat(1, params.exp.session.mri.wide.session_types,params.exp.session.mri.deep.session_types);
    
elseif strcmp(session_type,'BEHAVIOR')
    all_sessions         = params.exp.session.behavior.ses_blocks;
    runs_per_session     = params.exp.session.behavior.n_runs_per_session;
    session_types        = params.exp.session.behavior.session_types;
end



for ses = 1:size(all_sessions,3)
    
    global_session_counter = global_session_counter+1;
    
    for st = 1:size(all_sessions,4) % session types
        
        if ~isnan(session_types(ses,st))
            %% First we check how many repeats of stim-task crossings we want to run within a session
            fprintf('SESSION %03d, session_type %02d..\n',ses,session_types(ses,st))
            
            % also keep track of local stim-task block allocation
            stimtask_tracker_local = ones(length(params.exp.stimclassnames), length(params.exp.taskclassnames));
            
            % Nr of blocks for each stim-task crossing. We use the nr of blocks,
            % same trials, and same unique image for each subject
            block_distr = all_sessions(:,:,ses,st);
            
            bb = 1; % block tracker
            scc_bb = 0; ltm_bb = 0; % separate counters for scc and ltm task crossings, because it is spread across stimulus classes
            
            for sc = 1:length(params.exp.stimclassnames)
                for tc = 1:length(params.exp.taskclassnames)
                    
                    if strcmp(session_type,'MRI')
                        task_start = params.exp.session.mri.task_start(tc);
                        
                    elseif strcmp(session_type,'BEHAVIOR')
                        task_start = params.exp.session.behavior.task_start(tc);
                    end
                    if  ses >= task_start
                        % get nr of blocks for this specific crossing
                        n_blocks = block_distr(sc,tc);
                        
                        % if it isn't zero, then add it to the table
                        if (n_blocks > 0)
                            
                            % get nr of blocks we want to allocate
                            curr_blocks = [stimtask_tracker_global(sc,tc) : (stimtask_tracker_global(sc,tc)+n_blocks-1)];
                            
                            for ii = 1:length(curr_blocks)
                                
                                
                                if ismember(params.exp.taskclassnames{tc},'scc')
                                    % We can't sort on stim class nr for SCC
                                    % otherwise stim will be allocated to different
                                    % blocks. Hence we use an if statement and call
                                    % the crossing_name 'scc-all'
                                    idx = ((condition_master.stim_class==99) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr== scc_bb + curr_blocks(ii)));
                                    condition_master.crossing_name(idx) = {sprintf('%s-all',params.exp.taskclassnames{tc})};
                                    stim_class_tmp_name = 'ALL';
                                elseif ismember(params.exp.taskclassnames{tc},'ltm')
                                    % We can't sort on stim class nr for LTM either (they are all defined as 99)
                                    idx = ((condition_master.stim_class==99) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr == ltm_bb + curr_blocks(ii)));
                                    condition_master.crossing_name(idx) = {sprintf('%s-all',params.exp.taskclassnames{tc})};
                                    stim_class_tmp_name = 'ALL';
                                else
                                    idx = ((condition_master.stim_class==sc) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr==curr_blocks(ii)));
                                    condition_master.crossing_name(idx)          = {sprintf('%s-%s',params.exp.taskclassnames{tc}, params.exp.stimclassnames{sc})};
                                    stim_class_tmp_name = stim_class_abbr{sc};
                                end
                                
                                % get block numbers for this stim-task crossing
                                condition_master.session_nr(idx)        = ses;
                                condition_master.session_type(idx)      = session_types(ses,st);
                                blck_name                               = condition_master.crossing_name(idx);
                                condition_master.crossing_name(idx)     = blck_name(1);
                                condition_master.crossing_nr(idx)       = find(strcmp(blck_name(1),params.exp.crossingnames));
                                condition_master.block_nr(idx)          = bb;
                                condition_master.global_session_nr(idx) = global_session_counter;
                                

                                
                                % get condition label for each trial and stim loc
                                sub_idx = find(idx);
                                for tt = 1:length(sub_idx)
                                    
                                    if condition_master.is_catch(sub_idx(tt))==1
                                        if condition_master.is_cued(sub_idx(tt))==1
                                            cue_label = {sprintf('%s 0000 X LCUED %s',stim_class_tmp_name,task_class_abbr{tc}), ...
                                                sprintf('%s 0000 X LCUED %s',stim_class_tmp_name,task_class_abbr{tc})};
                                        elseif condition_master.is_cued(sub_idx(tt))==2
                                            cue_label = {sprintf('%s 000 X RCUED %s',stim_class_tmp_name,task_class_abbr{tc}), ...
                                                sprintf('%s 0000 X RCUED %s',stim_class_tmp_name,task_class_abbr{tc})};
                                        elseif condition_master.is_cued(sub_idx(tt))==3
                                            if strcmp(stim_class_tmp_name,'NS')
                                                cue_label = {sprintf('%s 0000 X NCUED %s',stim_class_tmp_name,task_class_abbr{tc}), NaN};
                                            else
                                                cue_label = {sprintf('%s 0000 X NCUED %s',stim_class_tmp_name,task_class_abbr{tc}), ...
                                                    sprintf('%s 0000 X NCUED %s',stim_class_tmp_name,task_class_abbr{tc})};
                                            end
                                        end
                                    elseif condition_master.is_cued(sub_idx(tt))==1
                                        cue_label = {sprintf('%s %04d L CUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}), ...
                                            sprintf('%s %04d R UNCUED %s',stim_class_tmp_name,condition_master.stim_nr_right(sub_idx(tt)),task_class_abbr{tc})};
                                    elseif condition_master.is_cued(sub_idx(tt))==2
                                        cue_label = {sprintf('%s %04d L UNCUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}),...
                                            sprintf('%s %04d R CUED %s',stim_class_tmp_name,condition_master.stim_nr_right(sub_idx(tt)),task_class_abbr{tc})};
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
                                        
                                    end
                                    
                                    % add dash here, otherwise sprintf interprets dash as left align
                                    for cl = 1:length(cue_label)
                                        if isnan(cue_label{cl})
                                            condition_master.condition_name(sub_idx(tt),cl) = {NaN};
                                            condition_master.condition_nr(sub_idx(tt),cl) = NaN;
                                        else
                                            condition_master.condition_name(sub_idx(tt),cl) = strrep(cue_label(cl),' ', '-');
                                            condition_master.condition_nr(sub_idx(tt),cl) = vcd_conditionName2Number(condition_master.condition_name(sub_idx(tt),cl));
                                        end
                                    end
                                end
                                
                                if params.exp.trial.double_epoch_tasks(tc)
                                    condition_master.trial_type(idx) = 2; % double epoch;
                                else
                                    condition_master.trial_type(idx) = 1; % double epoch;
                                end
                                
                                bb = bb+1;
                                
                                stimtask_tracker_local(sc,tc) = stimtask_tracker_local(sc,tc)+1;
                                stimtask_tracker_global(sc,tc) = stimtask_tracker_global(sc,tc)+1;
                                
                            end
                            
                            if strcmp(params.exp.taskclassnames{tc},'scc')
                                scc_bb = scc_bb+max(curr_blocks);
                            elseif strcmp(params.exp.taskclassnames{tc},'ltm')
                                ltm_bb = ltm_bb+max(curr_blocks);
                            end
                            fprintf('\n %s   \t:\t %d block(s)',blck_name{1},n_blocks)
                        end
                    end
                end
            end
            fprintf('\nTOTAL MINIBLOCKS: %d',bb-1)
            
            % the number of blocks in the table should be equal to the number of
            % tracked stim-task-crossings allocated, as well as the sum vs vector length
            assert(isequal(stimtask_tracker_local-1,block_distr))
            tmp = max(condition_master.block_nr(~isnan(condition_master.block_nr) & condition_master.session_nr==ses));
            assert(isequal(tmp,sum(sum(stimtask_tracker_local-1))))
            
            %     % visualize blocks and trials
            %     cmap = cmapturbo(500);
            %     figure; set(gcf,'Position',[1 1 1200 300]);
            %     colormap(cmap); sgtitle(sprintf('Session %02d',ses))
            %     new_block_line = find(abs(diff(condition_master.block_local_trial_nr(~isnan(condition_master.global_block_nr))))>1);
            %     subplot(311); cla;
            %     imagesc(condition_master.block_local_trial_nr(~isnan(condition_master.global_block_nr))');
            % %         hold on;
            % %     for ff = 1:length(new_block_line)
            % %         plot([new_block_line(ff),new_block_line(ff)]+0.05,[0,4],'k','linewidth', 0.01)
            % %     end
            %     set(gca,'YTick',1,'YTickLabel',{'trial nr'})
            %     subplot(312);cla;
            %     imagesc(condition_master.block_nr(~isnan(condition_master.global_block_nr))');
            %
            % %         hold on;
            % %     for ff = 1:length(new_block_line)
            % %         plot([new_block_line(ff),new_block_line(ff)]+0.05,[0,4],'k','linewidth', 0.01)
            % %     end
            %     set(gca,'YTick',1,'YTickLabel',{'local block nr'})
            %     subplot(313)
            %     imagesc(condition_master.global_block_nr(~isnan(condition_master.global_block_nr))');
            % %     hold on;
            % %     for ff = 1:length(new_block_line)
            % %         plot([new_block_line(ff),new_block_line(ff)]+0.05,[0,4],'k','linewidth', 0.01)
            % %     end
            %     set(gca,'YTick',1,'YTickLabel',{' block nr'})
            %
            %     if params.store_imgs
            %         saveFigsFolder = fullfile(vcd_rootPath,'figs');
            %         filename = sprintf('vcd_session%02d_blocks_pre_shuffle.png', ses);
            %         print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
            %     end
            
            %% SHUFFLE RUN AND BLOCK ORDER WITHIN RUNS
            
            block_vec        = condition_master.block_nr;
            trialtype_vec    = condition_master.trial_type;
            crossing_vec     = condition_master.crossing_nr;
            ses_idx          = (condition_master.session_nr==ses);
            ses_sub          = find(ses_idx);
            ses_blocks       = block_vec(ses_idx);
            ses_trialtype    = trialtype_vec(ses_idx);
            ses_crossing_vec = crossing_vec(ses_idx);
            
            % get unique block nrs and associated trial types
            [unique_blocks, ia] = unique(ses_blocks);
            crossings_unique_blocks  = ses_crossing_vec(ia);
            unique_trialtypes = ses_trialtype(ia);
            n_trialtype1 = sum(unique_trialtypes==1); % single
            n_trialtype2 = sum(unique_trialtypes==2); % double
            
            if strcmp(session_type,'MRI')
                runtype1 = cat(1,params.exp.session.wide.nr_of_type1_runs,params.exp.session.deep.nr_of_type1_runs);
                runtype2 = cat(1,params.exp.session.wide.nr_of_type2_runs,params.exp.session.deep.nr_of_type2_runs);
                runtype3 = cat(1,params.exp.session.wide.nr_of_type3_runs,params.exp.session.deep.nr_of_type3_runs);
                runtype4 = cat(1,params.exp.session.wide.nr_of_type4_runs,params.exp.session.deep.nr_of_type4_runs);
            elseif strcmp(session_type,'BEHAVIOR')
                runtype1 = params.exp.session.behavior.nr_of_type1_runs;
                runtype2 = params.exp.session.behavior.nr_of_type2_runs;
                runtype3 = params.exp.session.behavior.nr_of_type3_runs;
                runtype4 = params.exp.session.behavior.nr_of_type4_runs;
            end
            
            block_type_allocation  = [params.exp.run.run_type1,params.exp.run.run_type2,params.exp.run.run_type3,params.exp.run.run_type4];
            ses_run_types          = [runtype1(ses,st),runtype2(ses,st),runtype3(ses,st),runtype4(ses,st)];
            nr_blocks_per_run_type = ses_run_types * block_type_allocation';

            assert(isequal(unique_blocks,[1:length(unique_blocks)]'))
            if strcmp(session_type,'MRI')
                if ses==27
                    assert(length(unique_blocks) < sum(nr_blocks_per_run_type))
                    assert(n_trialtype1 < nr_blocks_per_run_type(1))
                    assert(n_trialtype2 < nr_blocks_per_run_type(2))
                else
                    assert(isequal(length(unique_blocks), sum(nr_blocks_per_run_type)))
                    assert(isequal(n_trialtype1 , nr_blocks_per_run_type(1)))
                    assert(isequal(n_trialtype2 , nr_blocks_per_run_type(2)))
                end
            elseif strcmp(session_type,'BEHAVIOR')
                % we have one block too many for the allocated runs, so we shove if that block into a run of runtype3
                assert(length(unique_blocks) == sum(nr_blocks_per_run_type)+1)
                assert(n_trialtype1 == nr_blocks_per_run_type(1)+1)
                assert(n_trialtype2 == nr_blocks_per_run_type(2))
            end
            
            
            runs_to_fill = [ones(1, runtype1(ses)), 2*ones(1,runtype2(ses)),3*ones(1,runtype3(ses)),4*ones(1,runtype4(ses))];
            
            curr_runs = zeros(length(runs_to_fill),7); % max 7 possible slots
            
            % separate single and double stim blocks:
            s_blocks = unique_blocks(find(unique_trialtypes==1));
            d_blocks = unique_blocks(find(unique_trialtypes==2));

            assert(isequal(length(s_blocks)+length(d_blocks),length(unique_blocks)))
            
            % Shuffle block order within a session
            while 1
                s_blocks_shuffled = s_blocks(randperm(length(s_blocks),length(s_blocks)));
                d_blocks_shuffled = d_blocks(randperm(length(d_blocks),length(d_blocks)));
                
                shuffled_s_block_crossings = crossings_unique_blocks(s_blocks_shuffled);
                shuffled_d_block_crossings = crossings_unique_blocks(d_blocks_shuffled);
                
                for nn = 1:4:length(shuffled_s_block_crossings)
                    if (nn+3) <= length(shuffled_s_block_crossings)
                        tmp_blocks = shuffled_s_block_crossings(nn:(nn+3));
                        
                        if length(unique(tmp_blocks)) >= 4
                            shuffle_ok = true;
                        else
                            shuffle_ok = false;
                            break;
                        end
                    end
                end
                
                if shuffle_ok
                    for nn = 1:4:length(shuffled_d_block_crossings)
                        if (nn+3) <= length(shuffled_d_block_crossings)
                            tmp_blocks = shuffled_s_block_crossings(nn:(nn+3));
                            
                            if length(unique(tmp_blocks)) >= 2
                                shuffle_ok = true;
                            else
                                shuffle_ok = false;
                                break;
                            end
                        end
                    end
                end
                
                if shuffle_ok,
                    break;
                end
            end
                    
            
            % reset counters
            bb_s = 0; bb_d = 0;
            
            % Double check if we have at least two double-epoch trial types per run
            for rr = 1:length(runs_to_fill)
                
                blocks_to_fill = block_type_allocation(:,runs_to_fill(rr))';
                
                s_block_idx = (bb_s+1):(bb_s+blocks_to_fill(1));
                d_block_idx = (bb_d+1):(bb_d+blocks_to_fill(2));
                
                % if we can still index blocks
                if all(s_block_idx <= length(s_blocks_shuffled)) && all(d_block_idx <= length(d_blocks_shuffled))
                    
                    if ~isempty(s_block_idx) && ~isempty(d_block_idx)
                        
                        tmp = cat(1,s_blocks_shuffled(s_block_idx), d_blocks_shuffled(d_block_idx))';
                        curr_runs(rr,1:length(tmp)) = shuffle_concat(tmp,1);
                        
                    elseif isempty(s_block_idx) && ~isempty(d_block_idx)
                        
                        tmp = cat(1,d_blocks_shuffled(d_block_idx))';
                        curr_runs(rr,1:length(tmp)) = shuffle_concat(tmp,1);
                        
                    elseif ~isempty(s_block_idx) && isempty(d_block_idx)
                        
                        tmp = cat(1,s_blocks_shuffled(s_block_idx))';
                        curr_runs(rr,1:length(tmp)) = shuffle_concat(tmp,1);
                        
                    end
                    
                    
                    % if are at the end of available stim presentation blocks
                elseif any(s_block_idx <= length(s_blocks_shuffled)) || any(d_block_idx <= length(d_blocks_shuffled))
                    d_block_idx = d_block_idx(d_block_idx<= length(d_blocks_shuffled));
                    if ~isempty(d_block_idx)
                        tmp = cat(1,d_blocks_shuffled(d_block_idx))';
                        curr_runs(rr,1:length(tmp)) = shuffle_concat(tmp,1);
                    end
                    s_block_idx = s_block_idx(s_block_idx<= length(s_blocks_shuffled));
                    if ~isempty(s_block_idx)
                        s_block_idx = s_block_idx(s_block_idx <= length(s_blocks_shuffled));
                        tmp = cat(1,s_blocks_shuffled(s_block_idx), s_blocks_shuffled(s_block_idx))';
                        curr_runs(rr,1:length(tmp)) = shuffle_concat(tmp,1);
                    end
                end
                
                
                if ~isempty(s_block_idx)
                    bb_s = s_block_idx(end);
                end
                if ~isempty(d_block_idx)
                    bb_d = d_block_idx(end);
                end
            end
            
            if strcmp(session_type,'BEHAVIOR')
                if length(s_blocks_shuffled) > s_block_idx(end)
                    empty_slot = find(curr_runs(end,:)==0);
                    empty_slot_idx = find(curr_runs(:,empty_slot(1))==0);
                    curr_runs(empty_slot_idx(1), empty_slot(1)) = s_blocks_shuffled(end);
                end
                
                if length(d_blocks_shuffled) > d_block_idx(end)
                    empty_slot = find(curr_runs(end,:)==0);
                    empty_slot_idx = find(curr_runs(:,empty_slot(1))==0);
                    curr_runs(empty_slot_idx(1), empty_slot(1)) = d_blocks_shuffled(end);
                end
            end
            
            %% Shuffle run order (runs x blocks)
            clear curr_runs_shuffled;
            for rr = 1:size(curr_runs,1)
                tmp_run_blocks =  curr_runs(rr,curr_runs(rr,:)>0);
                curr_runs_shuffled(rr,1:length(tmp_run_blocks)) = shuffle_concat(tmp_run_blocks,1);
            end
            
            run_shuffle = shuffle_concat(1:size(curr_runs_shuffled,1),1);
            
            for run_idx = 1:length(run_shuffle)
                %
                block_in_this_run = curr_runs_shuffled(run_shuffle(run_idx),:);
                block_in_this_run = block_in_this_run(block_in_this_run>0);
                
                fprintf('[%s]: Run %02d, block order for all subjects:',mfilename, run_idx)
                disp(block_in_this_run)
                
                for iii = 1:length(block_in_this_run)
                    
                    global_block_counter = global_block_counter + 1;
                    
                    % find the OG trials that correspond to the block from the shuffled list
                    idx0 = find(ses_blocks==block_in_this_run(iii));
                    idx1 = ses_sub(idx0);
                    
                    nr_trials = length(idx1);
                    
                    % get the new row nr for these trials
                    condition_master.subj_block_nr(idx1) = iii;
                    condition_master.run_nr(idx1)   = run_idx;
                    
                    % add global trial and block nr
                    condition_master.global_trial_nr(idx1)   = global_trial_counter + (1:nr_trials);
                    condition_master.global_block_nr(idx1)   = global_block_counter;
                    
                    global_trial_counter = global_trial_counter + nr_trials ;
                    
                    assert(unique(condition_master.session_nr(idx1))==ses)
                end
                
            end
     
            % sort by run_nr, then block_nr
            all_trial_idx = find(condition_master.session_nr==ses & condition_master.session_type==st);
            [~,sort_idx0]  = sort(condition_master.run_nr(all_trial_idx));
            run_idx_sorted = all_trial_idx(sort_idx0);
            condition_master(all_trial_idx,:) = condition_master(run_idx_sorted,:);
            
            all_trial_idx = find(condition_master.session_nr==ses & condition_master.session_type==st);

            for zz = 1:length(unique(condition_master.run_nr(all_trial_idx)))
                run_block_idx = condition_master.run_nr(all_trial_idx) == zz;
                single_run_blocks = condition_master.subj_block_nr(all_trial_idx(run_block_idx));
                [~,sort_idx1]     = sort(single_run_blocks);
                tmp_idx = all_trial_idx(run_block_idx);
                run_block_idx_sorted = tmp_idx(sort_idx1);
                condition_master(tmp_idx,:) = condition_master(run_block_idx_sorted,:);
            end
            
        end % session type
    end
end % session nr

% Remove unused trials and obsolete columns
condition_master(isnan(condition_master.session_nr),:) = [];
condition_master.block_nr = [];
condition_master.block_nr = condition_master.subj_block_nr;
condition_master.subj_block_nr = [];
condition_master.stim_class_unique_block_nr = [];
condition_master.unique_trial_nr = [];

% resort columns
newOrderNames = {'session_nr','session_type','run_nr','block_nr','trial_nr',...
    'global_session_nr', 'global_block_nr','global_trial_nr',...
    'condition_nr','condition_name', ...
    'stim_class','stim_class_name', ...
    'task_class','task_class_name', ...
    'crossing_nr','crossing_name',...
    'stim_nr_left', 'stim_nr_right', ...
    'is_cued', 'is_catch', ...
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
newOrder_idx = [];
for fn = 1:length(newOrderNames)
    [~,col_idx] = ismember(condition_master.Properties.VariableNames, newOrderNames{fn});
    newOrder_idx(fn) = find(col_idx);
end
condition_master = condition_master(:,newOrder_idx);


% Update repeat nr
repeat_nr = zeros(size(condition_master,1),2);
for cols = [1,2]
    [uniquerow, rowidx] = unique(condition_master.condition_nr(:,cols),'stable');
    nr_of_occurrences = accumarray(rowidx, 1);
    
    for nrow = 1:size(condition_master,1)
        
        curr_cond = condition_master.condition_nr(nrow,cols);
        if ~isnan(curr_cond)
            [~,comb_nr] = intersect(uniquerow,curr_cond);
            repeat_nr(comb_nr, cols) = repeat_nr(comb_nr,cols)+1;
        end
    end
end
% insert nans for zeros
repeat_nr(repeat_nr(:,1)==0,1) = NaN;
repeat_nr(repeat_nr(:,2)==0,2) = NaN;

% update repeat_nr col
condition_master.repeat_nr = repeat_nr;

if params.store_params
    fprintf('\n[%s]:Storing updated condition_master data..\n',mfilename)
    saveDir = fullfile(vcd_rootPath,'workspaces','info');
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir,sprintf('condition_master_with_blocks_%s_%s_%s.mat',params.disp.name, session_type, datestr(now,30))),'condition_master')
end

% Visualize results
if params.verbose
    close all; makeprettyfigures
    for ses = 1:size(all_sessions,3)
        % visualize blocks and trials, now after shuffle
        figure; set(gcf,'Position',[1 1 1200 600]);
        ses_data = condition_master(~isnan(condition_master.session_nr),:);
        x1 = (ses_data.session_nr==ses);
        [yy,x2] = sort(ses_data.run_nr(x1));
        x3 = find(x1);
        x3 = x3(x2);
        [~,new_run_line] = unique(yy,'stable');
        dataToPlot = [ses_data.session_nr(x3), ...
            ses_data.run_nr(x3), ...
            ses_data.subj_block_nr(x3,:)]';
        imagesc(dataToPlot);
        cb = colorbar;
        cmap = cmapturbo(max(dataToPlot(:))+1);
        colormap(cmap);
        set(gca,'CLim',[min(dataToPlot(:)) max(dataToPlot(:))])
        hold on;
        for ff = 1:length(new_run_line)
            plot([new_run_line(ff),new_run_line(ff)],[0,6],'k','linewidth',4)
        end
        for ii = 1:params.exp.total_subjects, label{ii} = sprintf('subj %d',ii); end
        set(gca,'YTick',[1:params.exp.total_subjects+2],'YTickLabel',{'session nr','run nr', label{:}});
        title(sprintf('SESSION %02d OVERVIEW',ses),'FontSize',20);
        xlabel('BLOCK NR')
        
        
        if params.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master1_%s',session_type));
            if ~exist(saveFigsFolder,'dir'); mkdir(saveFigsFolder);end
            filename = sprintf('vcd_session%02d_subjblocks_post_shuffle.png',ses);
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
    end
end
return