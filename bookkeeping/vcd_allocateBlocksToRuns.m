function condition_master = vcd_allocateBlocksToRuns(p)

condition_master = p.trials;

if sum(strcmp(condition_master.Properties.VariableNames,'session_nr'))==0
    t_session = table( NaN(size(condition_master,1),1), NaN(size(condition_master,1),1),...
        NaN(size(condition_master,1),1), NaN(size(condition_master,1),p.exp.total_subjects), ...
        NaN(size(condition_master,1),1), num2cell(NaN(size(condition_master,1),1)), ...
        NaN(size(condition_master,1),1), NaN(size(condition_master,1),1));
    t_session.Properties.VariableNames = {'session_nr','run_nr',...
        'within_ses_block_nr','subj_within_run_block_nr',...
        'block_ID','block_name','trial_type','across_session_block_nr'};
    
    condition_master = cat(2,condition_master, t_session);
end


%%
stimtask_tracker = ones(length(p.exp.stimClassLabels), length(p.exp.taskClassLabels));

global_block_counter = 1;

for ses = 1:p.exp.n_sessions
    
    %% First we check how many repeats of stim-task crossings we want to run within a session
    fprintf('\nSESSION %d:',ses)
    
    % Same blocks / same trials / same retinal image per subject
    block_distr = p.exp.ses_blocks(:,:,ses);
    
    bb = 1; % block tracker
    scc_bb = 1; % separate counter for scc task, because it is spread across stimulus classes
    
    for sc = 1:length(p.exp.stimClassLabels)
        for tc = 1:length(p.exp.taskClassLabels)
            
            if ses >= p.exp.session.task_start(tc)
                
                n_blocks = block_distr(sc,tc);
                
                if (n_blocks > 0)
                    
                    curr_blocks = [stimtask_tracker(sc,tc) : (stimtask_tracker(sc,tc)+n_blocks-1)];
                    
                    for ii = 1:length(curr_blocks)

                        if strcmp(p.exp.taskClassLabels{tc},'scc')
                            % we can't sort on stim class name for SCC,
                            % otherwise stim will be allocated to different
                            % blocks
                            idx = ((condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr==scc_bb));
                            condition_master.session_nr(idx)        = ses;
                            condition_master.block_name(idx)        = {sprintf('%s-all',p.exp.taskClassLabels{tc})};
                            cond_name                               = condition_master.block_name(idx); 
                            cond_name = cond_name(1);
                            condition_master.block_ID(idx)      = find(strcmp(cond_name,p.exp.stimTaskLabels));
                            condition_master.within_ses_block_nr(idx) = bb;
                            condition_master.across_session_block_nr(idx) = global_block_counter;
                            scc_bb = scc_bb+1;
                        else
                            idx = ((condition_master.stim_class==sc) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr==curr_blocks(ii)));
                            condition_master.session_nr(idx)          = ses;
                            condition_master.block_name(idx)          = {sprintf('%s-%s',p.exp.taskClassLabels{tc},p.exp.stimClassLabels{sc})};
                            cond_name                                 = condition_master.block_name(idx); 
                            cond_name                                 = cond_name(1);
                            condition_master.block_ID(idx)            = find(strcmp(cond_name,p.exp.stimTaskLabels));
                            condition_master.within_ses_block_nr(idx) = bb;
                            condition_master.across_session_block_nr(idx) = global_block_counter;
                        end
                        
%                         condition_master.within_session_repeat(idx) = curr_blocks(ii);

                        if p.exp.trial.double_epoch_tasks(tc)
                           condition_master.trial_type(idx) = 2; % double epoch;
                        else
                           condition_master.trial_type(idx) = 1; % double epoch;
                        end
                        
                        bb = bb+1;
                        global_block_counter = global_block_counter+1;
                        stimtask_tracker(sc,tc) = stimtask_tracker(sc,tc)+1;
                    end
                    
                    fprintf('\n%s   \t:\t %d block(s)',cond_name{:},n_blocks)
                end
            end
        end
    end
    fprintf('\nTOTAL MINIBLOCKS: %d',bb-1)
    
    % the number of blocks in the table should be equal to the number of
    % tracked stim-task-crossings allocated
    tmp = max(condition_master.across_session_block_nr(~isnan(condition_master.across_session_block_nr)));
    assert(isequal(tmp,sum(sum(stimtask_tracker-1))))
    
    %% Check runs for block order
    % We want to make sure that runs have at least 2 blocks that use 
    % double-epoch trials, otherwise some runs will be much longer than others..
    ses_blocks = condition_master.within_ses_block_nr((condition_master.session_nr==ses),:);
    trial_type = condition_master.trial_type((condition_master.session_nr==ses),:);
    assert(isequal(unique(ses_blocks),[1:length(unique(ses_blocks))]'))
    
%     blocks_nrs = condition_master.stim_class_unique_block_nr((condition_master.session_nr==ses),:);
%     block_IDs = condition_master.block_ID((condition_master.session_nr==ses),:);
    
    while 1
        run_ok = 0;
        run_not_ok = false;

        % Shuffle blocks
        shuffle_idx = randperm(length(unique(ses_blocks)),length(unique(ses_blocks)));
        ses_blocks_new_order = [];

        % Find the corresponding blocks for each shuffle idx and
        % reorder blocks defined for this sesesion
        for ii = 1:length(shuffle_idx)
            ses_block_idx = find(ses_blocks==shuffle_idx(ii));
            idx = (length(ses_blocks_new_order)+1):((length(ses_blocks_new_order)+1)+length(ses_block_idx)-1);
            ses_blocks_new_order(idx,1) = ses_blocks(ses_block_idx);
            ses_blocks_new_order(idx,2) = trial_type(ses_block_idx);
        end

        % Divide single vector of session blocks into runs
        [bo, bi] = unique(ses_blocks_new_order(:,1),'stable');
        rz_ses_blocks1 = reshape(bo, p.exp.run.blocks_per_run, []);
        rz_ses_blocks2 = reshape(ses_blocks_new_order(bi,2), p.exp.run.blocks_per_run, []);

        for jj = 1:size(rz_ses_blocks2,2)
                if sum(rz_ses_blocks2(:,jj)==2)>2
                    run_not_ok = true;
                else
                    run_ok = run_ok + 1;
                end

                if run_not_ok
                    break;
                end
        end
        if run_ok == p.exp.n_runs_per_session
            break;
        end
    end

    assert(sum(sum(rz_ses_blocks2==2,1)>0) >= (size(rz_ses_blocks2,2)-2))
    %% EK CHECK BLOCK NR ALLOCATION ACROSS SUBJECTS
    % Once we are happy, we want to shuffle the order of blocks within 
    % a run for each subject
    for run_idx = 1:size(rz_ses_blocks1,2)
        
        for sj = 1:p.exp.total_subjects
            
            fprintf('\n[%s]: Subject %03d, run %02d, block order: ',mfilename, sj, run_idx)
            subj_ses_block_shuffle = shuffle_concat(rz_ses_blocks1(:,run_idx),1);
        
            for block_idx = 1:length(subj_ses_block_shuffle)
            
                image_idx = find(condition_master.session_nr==ses & condition_master.within_ses_block_nr == subj_ses_block_shuffle(block_idx));
                
                condition_master.subj_within_run_block_nr(image_idx,sj) = block_idx;
                condition_master.run_nr(image_idx) = run_idx;
                fprintf('%d ',subj_ses_block_shuffle(block_idx))
            end
        end % subject 
    end
    fprintf('\n')
end % session





return