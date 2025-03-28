function condition_master = vcd_allocateBlocksToRuns(p)

condition_master = p.trials;

% Add extra columns to master table
if sum(strcmp(condition_master.Properties.VariableNames,'session_nr'))==0
    t_session = table( NaN(size(condition_master,1),1), NaN(size(condition_master,1),1),...
        NaN(size(condition_master,1),1), NaN(size(condition_master,1),p.exp.total_subjects), ...
        NaN(size(condition_master,1),1), num2cell(NaN(size(condition_master,1),1)), ...
        NaN(size(condition_master,1),1), NaN(size(condition_master,1),1));
    t_session.Properties.VariableNames = {'session_nr','run_nr',...
        'local_block_nr_within_ses',...
        'subj_within_run_block_nr',...
        'block_ID',...
        'block_name',...
        'trial_type',...
        'global_block_nr_across_ses'};
    
    condition_master = cat(2,condition_master, t_session);
end


%%
% Keep track of allocated blocks for each stim-task crossing
stimtask_tracker = ones(length(p.exp.stimClassLabels), length(p.exp.taskClassLabels));

global_block_counter = 1; % keep track of all the blocks allocated across sessions

% Find the corresponding table row for each block in this session
ses_blocks_new_order = NaN(size(condition_master,1),p.exp.total_subjects);

for ses = 1:p.exp.n_sessions
    
    %% First we check how many repeats of stim-task crossings we want to run within a session
    fprintf('\nSESSION %d:',ses)
    
    % Nr of blocks for each stim-task crossing. We use the nr of blocks,
    % same trials, and same unique image for each subject
    block_distr = p.exp.ses_blocks(:,:,ses);
    
    bb = 1; % block tracker
    scc_bb = 1; % separate counter for scc task, because it is spread across stimulus classes
    
    for sc = 1:length(p.exp.stimClassLabels)
        for tc = 1:length(p.exp.taskClassLabels)
            
            if ses >= p.exp.session.task_start(tc)
                % get nr of blocks for this specific crossing
                n_blocks = block_distr(sc,tc);
                
                % if it isn't zero, then add it to the table
                if (n_blocks > 0)
                    
                    % get nr of blocks we want to allocate 
                    curr_blocks = [stimtask_tracker(sc,tc) : (stimtask_tracker(sc,tc)+n_blocks-1)];
                    
                    for ii = 1:length(curr_blocks)

                        if strcmp(p.exp.taskClassLabels{tc},'scc')
                            % We can't sort on stim class name for SCC,
                            % otherwise stim will be allocated to different
                            % blocks. Hence we use an if statement and call
                            % the block_name 'scc-all'
                            idx = ((condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr==scc_bb));
                            condition_master.session_nr(idx)        = ses;
                            condition_master.block_name(idx)        = {sprintf('%s-all',p.exp.taskClassLabels{tc})};
                            cond_name                               = condition_master.block_name(idx); 
                            cond_name = cond_name(1);
                            condition_master.block_ID(idx)      = find(strcmp(cond_name,p.exp.stimTaskLabels));
                            condition_master.local_block_nr_within_ses(idx) = bb;
                            condition_master.global_block_nr_across_ses(idx) = global_block_counter;
                            scc_bb = scc_bb+1;
                        else
                            % find block numbers for this stim-task
                            % crossing
                            idx = ((condition_master.stim_class==sc) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr==curr_blocks(ii)));
                            condition_master.session_nr(idx)          = ses;
                            condition_master.block_name(idx)          = {sprintf('%s-%s',p.exp.taskClassLabels{tc},p.exp.stimClassLabels{sc})};
                            cond_name                                 = condition_master.block_name(idx); 
                            cond_name                                 = cond_name(1);
                            condition_master.block_ID(idx)            = find(strcmp(cond_name,p.exp.stimTaskLabels));
                            condition_master.local_block_nr_within_ses(idx)  = bb;
                            condition_master.global_block_nr_across_ses(idx) = global_block_counter;
                        end
                        
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
    tmp = max(condition_master.global_block_nr_across_ses(~isnan(condition_master.global_block_nr_across_ses)));
    assert(isequal(tmp,sum(sum(stimtask_tracker-1))))
    
    % visualize blocks and trials
    cmap = cmapturbo(500);
    figure; set(gcf,'Position',[1 1 1200 300]); 
    colormap(cmap); sgtitle(sprintf('Session %02d',ses))
    new_block_line = find(abs(diff(condition_master.block_local_trial_nr(~isnan(condition_master.global_block_nr_across_ses))))>1);
    subplot(311)
    imagesc(condition_master.block_local_trial_nr(~isnan(condition_master.global_block_nr_across_ses))');
        hold on;
    for ff = 1:length(new_block_line)
        plot([new_block_line(ff),new_block_line(ff)]+0.05,[0,4],'k')
    end
    set(gca,'YTick',1,'YTickLabel',{'trial nr'})
    subplot(312)
    imagesc(condition_master.local_block_nr_within_ses(~isnan(condition_master.global_block_nr_across_ses))');

        hold on;
    for ff = 1:length(new_block_line)
        plot([new_block_line(ff),new_block_line(ff)]+0.05,[0,4],'k')
    end
    set(gca,'YTick',1,'YTickLabel',{'local block nr'})
    subplot(313)
    imagesc(condition_master.global_block_nr_across_ses(~isnan(condition_master.global_block_nr_across_ses))');
    hold on;
    for ff = 1:length(new_block_line)
        plot([new_block_line(ff),new_block_line(ff)]+0.05,[0,4],'k')
    end
    set(gca,'YTick',1,'YTickLabel',{'global block nr'})
    
    if p.store_imgs
        saveFigsFolder = fullfile(vcd_rootPath,'figs');
        filename = sprintf('vcd_session%02d_blocks_pre_shuffle.png', ses);
        print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
    end
    
    %% Check runs for block order
    % We want to make sure that runs have at least 2 blocks that use 
    % double-epoch trials, otherwise some runs will be much longer than others..
    block_vec      = condition_master.local_block_nr_within_ses;
    trialtype_vec  = condition_master.trial_type;
    ses_idx        = (condition_master.session_nr==ses);
    ses_sub        = find(ses_idx);
    ses_blocks     = block_vec(ses_idx);
    ses_trialtype  = trialtype_vec(ses_idx);

    [unique_blocks, ia] = unique(ses_blocks);
    unique_trialtypes = ses_trialtype(ia);
    
    assert(isequal(unique_blocks,[1:length(unique_blocks)]'))
    assert(isequal(length(unique_blocks),p.exp.n_runs_per_session*p.exp.run.blocks_per_run))
    
    while 1
        run_ok = 0;
        run_not_ok = false;

        % Shuffle block order within a session
        shuffle_idx = randperm(length(unique_blocks),length(unique_blocks));
        
        unique_blocks_rnd = unique_blocks(shuffle_idx);
        unique_trialtypes_rnd = unique_trialtypes(shuffle_idx);
        
        % reshape blocks into runs
        rz_ses_blocks1 = reshape(unique_blocks_rnd, p.exp.run.blocks_per_run, []);
        rz_ses_blocks2 = reshape(unique_trialtypes_rnd, p.exp.run.blocks_per_run, []);
        
        % Double check if we have at least two double-epoch trial types per run
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


    
    for run_idx = 1:size(rz_ses_blocks1,2)
        
        reordered_blocks = rz_ses_blocks1(:,run_idx);
        
        %% BLOCK NR ALLOCATION PER SUBJECTS
        % Once we are happy with the blocks allocated to each run, we want to
        % shuffle the order of blocks within a run for each subject
        for sj = 1:p.exp.total_subjects
            subj_ses_block_shuffle = shuffle_concat(reordered_blocks,1);

            fprintf('\n[%s]: Run %02d, block order for subject %03d:',mfilename, run_idx,sj)
            disp(subj_ses_block_shuffle')
            
            for ii = 1:length(subj_ses_block_shuffle)
        
                % find the OG trials that correspond to the block from the shuffled list
                idx0 = find(ses_blocks==subj_ses_block_shuffle(ii));
                idx1 = ses_sub(idx0);
        
                nr_trials = length(idx1);
        
                % get the new row nr for these trials
                ses_blocks_new_order(idx1,sj) = repmat(ii,nr_trials,1);
                condition_master.run_nr(idx1) = run_idx;
                assert(unique(condition_master.session_nr(idx1))==ses)
            end
        end
    end
end % session

condition_master.subj_within_run_block_nr = ses_blocks_new_order;

for ses = 1:p.exp.n_sessions
    % visualize blocks and trials, now after shuffle
    cmap = cmapturbo(8);
    figure; set(gcf,'Position',[1 1 1200 600]);
    colormap(cmap); 
    x1 = (condition_master.session_nr==ses);
    [yy,x2] = sort(condition_master.run_nr(x1));
    x3 = find(x1);
    x3 = x3(x2);
    [~,new_run_line] = unique(yy,'stable');
    imagesc([condition_master.session_nr(x3), ...
            condition_master.run_nr(x3), ...
            condition_master.subj_within_run_block_nr(x3,:)]');
    hold on;
    for ff = 1:length(new_run_line)
        plot([new_run_line(ff),new_run_line(ff)],[0,6],'k')
    end
    for ii = 1:p.exp.total_subjects, label{ii} = sprintf('subj %d',ii); end
    set(gca,'YTick',[1:p.exp.total_subjects+2],'YTickLabel',{'session nr','run nr', label{:}}); title('subj block nr')

    if p.store_imgs
        saveFigsFolder = fullfile(vcd_rootPath,'figs');
        filename = sprintf('vcd_session%02d_subjblocks_post_shuffle.png',ses);
        print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
    end
end
return