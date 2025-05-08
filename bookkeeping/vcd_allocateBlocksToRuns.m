function condition_master = vcd_allocateBlocksToRuns(params,session_type)
% VCD function to allocate the unique trials and blocks to given runs.
% 
%   condition_master = vcd_allocateBlocksToRuns(params)
%
% INPUTS:
%   params              : 	(struct) parameter struct needed to get subject
%                           nrs and run type params. REQUIRES params.trials 
%                           to exist and contain condition_master v0. 
%   session_type        :   (str) label to define what type of session we
%                           are defining: 'MRI' or 'BEHAVIOR'.
%   
% OUTPUTS:
%   condition_master    :  (struct) updated condition master table with
%                           block and runs for each subject and session.

%% 

% copy condition master
condition_master = params.trials;

% Add extra columns to master table
if sum(strcmp(condition_master.Properties.VariableNames,'session_nr'))==0
    t_session = table( ...
        NaN(size(condition_master,1),1), ... session_nr
        num2cell(NaN(size(condition_master,1),1)), ... session_name
        NaN(size(condition_master,1),1), ... run_nr
        NaN(size(condition_master,1),1), ... block_nr
        NaN(size(condition_master,1), params.exp.total_subjects), ...
        num2cell(NaN(size(condition_master,1),2)), ... cond_name
        NaN(size(condition_master,1),1), ... block_ID
        num2cell(NaN(size(condition_master,1),1)), ... block_name
        NaN(size(condition_master,1),1), ... trial_type
        NaN(size(condition_master,1),1)); %  global_block_nr
    t_session.Properties.VariableNames = {...
        'session_nr',...
        'session_name', ...
        'run_nr',...
        'block_nr',...
        'subj_block_nr', ...
        'cond_name', ...
        'block_ID',...
        'block_name',...
        'trial_type',...
        'global_block_nr', ...
        };
    
    condition_master = cat(2, t_session(:,[1:7]), condition_master,t_session(:,[8,9]));
end


%% Preallocate counters/matrices

% Keep track of allocated blocks for each stim-task crossing and continue
% counting to avoid overwriting of trials
stimtask_tracker_global = ones(length(params.exp.stimclassnames), length(params.exp.taskclassnames));
    
global_block_counter = 1; % keep track of all the blocks allocated across sessions

% change stim/task labels to upper case and abbreviated version
stim_class_abbr = upper(params.exp.stimclassnames);
stim_class_abbr{1} = 'GBR';
task_class_abbr = upper(params.exp.taskclassnames);

% check session type
if strcmp(session_type,'MRI')
    session_names = {'WIDE1A', 'WIDE1B'};
    all_sessions = cat(3, params.exp.session.wide.ses_blocks, params.exp.session.deep.ses_blocks);
    for ii = params.exp.session.deep.session_nrs, session_names = cat(2,session_names{:},{sprintf('DEEP%02d',ii)}); end
    session_names{length(session_names)} = 'DEEP26A';
     session_names{length(session_names)+1} = 'DEEP26B';
elseif strcmp(session_type,'BEHAVIOR')
    session_names = {'BEHAVIOR01'};
    all_sessions = params.exp.session.behavior.ses_blocks;
end


for ses = 1:size(all_sessions,3)
    
    %% First we check how many repeats of stim-task crossings we want to run within a session
    fprintf('\nSESSION %s:',session_names{ses})
    
    % also keep track of local stim-task block allocation
    stimtask_tracker_local = ones(length(params.exp.stimclassnames), length(params.exp.taskclassnames));
    
    % Nr of blocks for each stim-task crossing. We use the nr of blocks,
    % same trials, and same unique image for each subject
    block_distr = all_sessions(:,:,ses);
    
    bb = 1; % block tracker
    scc_bb = 0; ltm_bb = 0; % separate counters for scc and ltm task crossings, because it is spread across stimulus classes
    
    for sc = 1:length(params.exp.stimclassnames)
        for tc = 1:length(params.exp.taskclassnames)
            
            if strcmp(session_type,'MRI')
                if ismember(ses, params.exp.session.wide.session_nrs)
                    task_start = params.exp.session.wide.task_start(tc);
                else
                    task_start = params.exp.session.deep.task_start(tc);
                end
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
                            % the block_name 'scc-all'
                            idx = ((condition_master.stim_class==99) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr== scc_bb + curr_blocks(ii)));
                            condition_master.block_name(idx) = {sprintf('%s-all',params.exp.taskclassnames{tc})};
                            stim_class_tmp_name = 'ALL';
                        elseif ismember(params.exp.taskclassnames{tc},'ltm')
                            % We can't sort on stim class nr for LTM either (they are all defined as 99)
                            idx = ((condition_master.stim_class==99) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr == ltm_bb + curr_blocks(ii)));
                            condition_master.block_name(idx) = {sprintf('%s-all',params.exp.taskclassnames{tc})};
                            stim_class_tmp_name = 'ALL';
                        else
                            idx = ((condition_master.stim_class==sc) & (condition_master.task_class==tc) & (condition_master.stim_class_unique_block_nr==curr_blocks(ii)));
                            condition_master.block_name(idx)          = {sprintf('%s-%s',params.exp.taskclassnames{tc}, params.exp.stimclassnames{sc})};
                            stim_class_tmp_name = stim_class_abbr{sc};
                        end
                        % get block numbers for this stim-task crossing
                        condition_master.session_nr(idx)        = ses;
                        condition_master.session_name(idx)      = session_names(ses);
                        blck_name                               = condition_master.block_name(idx);
                        blck_name                               = blck_name(1);
                        condition_master.block_ID(idx)          = find(strcmp(blck_name,params.exp.crossingnames));
                        condition_master.block_nr(idx)          = bb;
                        condition_master.global_block_nr(idx)   = global_block_counter;
                        
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
                                    cue_label = {sprintf('%s 000 X NCUED %s',stim_class_tmp_name,task_class_abbr{tc}), ...
                                        sprintf('%s 0000 X NCUED %s',stim_class_tmp_name,task_class_abbr{tc})};
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
                                    cue_label = {sprintf('%s %04d C NCUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}), ''};
                                
                                % NS catch trial
                                elseif condition_master.stim_nr_left(sub_idx(tt))==0 && strcmp(stim_class_tmp_name,'NS')
                                    cue_label = {sprintf('%s %04 C NCUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}), ''};
                                
                                else
                                    % fix trial clssic left/right stim
                                    cue_label = {sprintf('%s %04d L NCUED %s',stim_class_tmp_name,condition_master.stim_nr_left(sub_idx(tt)),task_class_abbr{tc}), ...
                                        sprintf('%s %04d R NCUED %s',stim_class_tmp_name,condition_master.stim_nr_right(sub_idx(tt)),task_class_abbr{tc})};
                                end
                            end
                            
                            % add dash here, otherwise sprintf interprets dash as left align
                            condition_master.cond_name(sub_idx(tt),:) = strrep(cue_label,' ', '-'); 
                        end
                        
                        if params.exp.trial.double_epoch_tasks(tc)
                            condition_master.trial_type(idx) = 2; % double epoch;
                        else
                            condition_master.trial_type(idx) = 1; % double epoch;
                        end
                        
                        bb = bb+1;
                        global_block_counter = global_block_counter+1;
                        stimtask_tracker_local(sc,tc) = stimtask_tracker_local(sc,tc)+1;
                        stimtask_tracker_global(sc,tc) = stimtask_tracker_global(sc,tc)+1;
                    end
                    
                    if strcmp(params.exp.taskclassnames{tc},'scc')
                        scc_bb = scc_bb+max(curr_blocks);
                    elseif strcmp(params.exp.taskclassnames{tc},'ltm')
                        ltm_bb = ltm_bb+max(curr_blocks);
                    end
                    fprintf('\n%s   \t:\t %d block(s)',blck_name{:},n_blocks)
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
    %     set(gca,'YTick',1,'YTickLabel',{'global block nr'})
    %
    %     if params.store_imgs
    %         saveFigsFolder = fullfile(vcd_rootPath,'figs');
    %         filename = sprintf('vcd_session%02d_blocks_pre_shuffle.png', ses);
    %         print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
    %     end
    
    %% SHUFFLE RUN AND BLOCK ORDER WITHIN RUNS

    block_vec      = condition_master.block_nr;
    trialtype_vec  = condition_master.trial_type;
    ses_idx        = (condition_master.session_nr==ses);
    ses_sub        = find(ses_idx);
    ses_blocks     = block_vec(ses_idx);
    ses_trialtype  = trialtype_vec(ses_idx);
    
    % get unique block nrs and associated trial types
    [unique_blocks, ia] = unique(ses_blocks);
    unique_trialtypes = ses_trialtype(ia);
    n_trialtype1 = sum(unique_trialtypes==1); % single
    n_trialtype2 = sum(unique_trialtypes==2); % double
    
    if strcmp(session_type,'MRI')
        runtype1 = cat(2,params.exp.session.wide.nr_of_type1_runs,params.exp.session.deep.nr_of_type1_runs);
        runtype2 = cat(2,params.exp.session.wide.nr_of_type2_runs,params.exp.session.deep.nr_of_type2_runs);
        runtype3 = cat(2,params.exp.session.wide.nr_of_type3_runs,params.exp.session.deep.nr_of_type3_runs);
        runtype4 = cat(2,params.exp.session.wide.nr_of_type4_runs,params.exp.session.deep.nr_of_type4_runs);
    elseif strcmp(session_type,'BEHAVIOR')
        runtype1 = params.exp.session.behavior.nr_of_type1_runs;
        runtype2 = params.exp.session.behavior.nr_of_type2_runs;
        runtype3 = params.exp.session.behavior.nr_of_type3_runs;
        runtype4 = params.exp.session.behavior.nr_of_type4_runs;
    end
    
    block_type_allocation  = [params.exp.run.run_type1,params.exp.run.run_type2,params.exp.run.run_type3,params.exp.run.run_type4];
    ses_run_types          = [runtype1(ses),runtype2(ses),runtype3(ses),runtype4(ses)];
    nr_blocks_per_run_type = ses_run_types * block_type_allocation';
    
    assert(isequal(unique_blocks,[1:length(unique_blocks)]'))
    if strcmp(session_type,'MRI')
        assert(isequal(length(unique_blocks), sum(nr_blocks_per_run_type)))
        assert(isequal(n_trialtype1 , nr_blocks_per_run_type(1)))
   	    assert(isequal(n_trialtype2 , nr_blocks_per_run_type(2)))
    elseif strcmp(session_type,'BEHAVIOR') % we don't fill all runs in behavioral session, so as long as there are less blocks allocated than possible across runs we are good
        assert(length(unique_blocks) < sum(nr_blocks_per_run_type))
        assert(n_trialtype1 < nr_blocks_per_run_type(1))
   	    assert(n_trialtype2 < nr_blocks_per_run_type(2))
    end
    
    
    runs_to_fill = [ones(1, runtype1(ses)), 2*ones(1,runtype2(ses)),3*ones(1,runtype3(ses)),3*ones(1,runtype4(ses))];
  
    curr_runs = zeros(length(runs_to_fill),7); % max 7 possible slots

    % separate single and double stim blocks:
    s_blocks = unique_blocks(find(unique_trialtypes==1));
    d_blocks = unique_blocks(find(unique_trialtypes==2));
    
    assert(isequal(length(s_blocks)+length(d_blocks),length(unique_blocks)))
    
    % Shuffle block order within a session
    s_blocks_shuffled = s_blocks(randperm(length(s_blocks),length(s_blocks)));
    d_blocks_shuffled = d_blocks(randperm(length(d_blocks),length(d_blocks)));
    
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
            
        % if we ran out of double stim presentation blocks, but not single 
        elseif all(s_block_idx <= length(s_blocks_shuffled)) && ~all(d_block_idx <= length(d_blocks_shuffled))
            if ~isempty(s_block_idx)
                tmp = cat(1,s_blocks_shuffled(s_block_idx), d_blocks_shuffled(d_block_idx))';
                curr_runs(rr,1:length(tmp)) = shuffle_concat(tmp,1);
            end
            if ~isempty(d_block_idx)
                d_block_idx = d_block_idx(d_block_idx <= length(d_blocks_shuffled));
                tmp = cat(1,d_blocks_shuffled(d_block_idx))';
                curr_runs(rr,1:length(tmp)) = shuffle_concat(tmp,1);
            end
             
        % if we ran out of single stim presentation blocks, but not double
        elseif ~all(s_block_idx <= length(s_blocks_shuffled)) && all(d_block_idx <= length(d_blocks_shuffled))
             if ~isempty(d_block_idx)
                tmp = cat(1,d_blocks_shuffled(d_block_idx))';
                curr_runs(rr,1:length(tmp)) = shuffle_concat(tmp,1);
             end
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
    
    %% Shuffle run order
    curr_runs_shuffled = curr_runs(randperm(size(curr_runs,1),size(curr_runs,1)),:);
    
    
    for run_idx = 1:size(curr_runs_shuffled,1)
        
        %% BLOCK NR ALLOCATION PER SUBJECTS
        % Once blocks are allocated to each run, we want to
        % shuffle the order of blocks within a run for each subject
        
        for sj = 1:params.exp.total_subjects
            
            this_run = curr_runs_shuffled(run_idx,:);
            this_run = this_run(this_run>0);
            
            subj_ses_block_shuffle = this_run(randperm(length(this_run),length(this_run)));
            
            fprintf('\n[%s]: Run %02d, block order for subject %03d:',mfilename, run_idx,sj)
            disp(subj_ses_block_shuffle)
            
            for iii = 1:length(subj_ses_block_shuffle)
                
                % find the OG trials that correspond to the block from the shuffled list
                idx0 = find(ses_blocks==subj_ses_block_shuffle(iii));
                idx1 = ses_sub(idx0);
                
                nr_trials = length(idx1);
                
                % get the new row nr for these trials
                condition_master.subj_block_nr(idx1,sj) = iii;
                condition_master.run_nr(idx1) = run_idx;
                assert(unique(condition_master.session_nr(idx1))==ses)
            end
        end
    end
end % session

if params.store_params
    fprintf('\n[%s]:Storing updated condition_master data..\n',mfilename)
    saveDir = fullfile(vcd_rootPath,'workspaces','info');
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir,sprintf('condition_master_w_subject_blocks_%s_%s_%s.mat',params.disp.name, session_type, datestr(now,30))),'condition_master')
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