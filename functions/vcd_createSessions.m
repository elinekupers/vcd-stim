function subject_sessions = vcd_createSessions(p)

if p.load_params
    
    % check if trial struct is already defined
    d = dir(fullfile(vcd_rootPath,'workspaces','info','trials.mat'));
    
    if length(d) > 1
        error('[%s]: Multiple trial.mat files! Please check', mfilename);
    end
    
    if ~isempty(d.name)
        load(fullfile(d.folder,d.name));
    else
        error('[%s]: Can''t find trial.mat files! Please check', mfilename);
    end
    
else
    all_trials = vcd_makeTrials(p);
end

if ~isfield(p,'exp')
    p.exp = vcd_getSessionParams;
end

%% per task/stim crossing: get nr of trials per miniblock

% Different trial order per block (already accomplished in vcd_makeTrials.m)
% Across a single session, each subject will experience the same miniblocks and unique images.



session = [];
miniblock = [];

for ses = 1:p.exp.n_sessions
    fprintf('\nSESSION %d:',ses)
    
    t_total = [];
    trial_order_timing = num2cell(zeros(1,5));
    
    % Same blocks / same trials / same retinal image per subject
    miniblock_distr = p.exp.ses_blocks(:,:,ses);
    
    stimtask_tracker = ones(size(miniblock_distr));
    
    bb = 1;
    
    for sc = 1:length(p.exp.stimClassLabels)
        for tc = 1:length(p.exp.taskClassLabels)
            
            if ses >= p.exp.session.task_start(tc)
                
                n_blocks = miniblock_distr(sc,tc);
                
                if (n_blocks > 0)
                    
                    fprintf('\n%s\t x %s\t:\t %d block(s)',p.exp.stimClassLabels{sc},p.exp.taskClassLabels{tc},n_blocks)
                    
                    curr_miniblocks = [stimtask_tracker(sc,tc) : (stimtask_tracker(sc,tc)+n_blocks-1)];
                    
                    for ii = 1:length(curr_miniblocks)
                        miniblock(bb).name        = sprintf('%s-%s',p.exp.taskClassLabels{tc},p.exp.stimClassLabels{sc});
                        miniblock(bb).ID          = find(strcmp(miniblock(bb).name,p.exp.stimTaskLabels));
                        miniblock(bb).st_block    = curr_miniblocks(ii);
                        miniblock(bb).trial       = all_trials.stim(sc).tasks(tc).miniblock(curr_miniblocks(ii)).trial;
                        
                        if p.exp.trial.double_epoch_tasks(tc)
                            miniblock(bb).trial_type = 2; % double epoch;
                        else
                            miniblock(bb).trial_type = 1; % double epoch;
                        end
                        bb = bb+1;
                        stimtask_tracker(sc,tc) = stimtask_tracker(sc,tc)+1;
                    end
                end
                
            end
        end
    end
    
    % Shuffle blocks
    shuffle_idx = randperm(length(miniblock),length(miniblock));
    miniblock_shuffled = miniblock(shuffle_idx);
    itis = shuffle_concat(p.exp.trial.ITI,1);
    
    fudge_factor = 20; % seconds
    
    % Divide blocks into runs
    rr = 1;
    total_block_i = 1;
    
    while rr < p.exp.n_runs_per_session && ...
            (total_block_i <= length(miniblock_shuffled))
        
        % reset counter
        total_run_time = 0;
        bb = 1;
        
        if total_run_time == 0
            if rr > 1 
                run_t0 = num2cell(zeros(1,5));
                 trial_order_timing = [trial_order_timing; run_t0];
            end
            % Add pre-blank period
            total_run_time = p.exp.run.pre_blank_dur;
            t_preblank = num2cell([0, NaN, NaN, p.exp.run.pre_blank_dur, total_run_time]);
            trial_order_timing = [trial_order_timing; t_preblank];
        end
        
        if isempty(itis)
            itis = shuffle_concat(p.exp.trial.ITI,3);
        end
        
        
        while (total_run_time < (p.exp.run.total_run_dur - p.exp.run.post_blank_dur - fudge_factor)) && ...
            (total_block_i <= length(miniblock_shuffled))
            
            if size(miniblock_shuffled(total_block_i).trial,1) > size(miniblock_shuffled(total_block_i).trial,2)
                n_trials_per_block = size(miniblock_shuffled(total_block_i).trial,1);
            else
                n_trials_per_block = size(miniblock_shuffled(total_block_i).trial,2);
            end
            
            
            if isempty(itis) || (length(itis) < n_trials_per_block)
                itis = cat(2, itis, shuffle_concat(p.exp.trial.ITI,2));
            end
            
            if miniblock_shuffled(total_block_i).trial_type == 1 % block with single stim epoch trials
                
                time_passed = (n_trials_per_block*p.exp.trial.single_epoch_dur) + ...
                    sum(itis(1:n_trials_per_block));
            elseif miniblock_shuffled(total_block_i).trial_type == 2 % block with double stim epoch trials
                time_passed = (n_trials_per_block*p.exp.trial.double_epoch_dur) + ...
                    sum(itis(1:n_trials_per_block));
            end
            
            % Task cue
            t_taskcue = num2cell([length(p.exp.stimTaskLabels)+1, NaN, NaN, p.exp.miniblock.task_cue_dur,  total_run_time+p.exp.miniblock.task_cue_dur]);
            
            % Trial IDs
            fn = fieldnames(miniblock_shuffled(total_block_i).trial);
            tmp = struct2cell(miniblock_shuffled(total_block_i).trial);
            unique_im_nr_idx = strcmp(fn,{'unique_im_nr'});
            t_id = cell(n_trials_per_block*2,1);
            tmp2 = squeeze(tmp(unique_im_nr_idx,:,:));
            t_id(1:2:end,:) = mat2cell(cell2mat(tmp2),ones(1,n_trials_per_block));
            t_id(2:2:end,:) = {[0]};
            
            % Thickening dir
            thicking_nr_idx = strcmp(fn,{'covert_att_loc_cue'});
            t_spat_cue =  cell(n_trials_per_block*2,1);
            t_spat_cue(1:2:end,:) = squeeze(tmp(thicking_nr_idx,:,1));
            t_spat_cue(2:2:end,:) = {[0]};
            
            % Get block ID
            b_id = num2cell(repmat(miniblock_shuffled(total_block_i).ID,2*n_trials_per_block,1));
            
            % get timing
            if miniblock_shuffled(total_block_i).trial_type == 1
                t_t = [repmat(p.exp.trial.single_epoch_dur,1,n_trials_per_block);itis(1:n_trials_per_block)];
            else
                t_t = [repmat(p.exp.trial.double_epoch_dur,1,n_trials_per_block);itis(1:n_trials_per_block)];
            end
            t_t = t_t(:);
            t_c = total_run_time + cumsum(t_t);
            
            % concat columns
            t_total = [b_id, t_id, t_spat_cue, num2cell(t_t), num2cell(t_c)];
            
            % Add inter-block-interval (1 TR + whatever is left)
            IBI = p.exp.TR + mod(t_total{end,end},p.exp.TR);
            t_ibi = {0,NaN,NaN,IBI,t_total{end,end} + IBI};
            
            % concat more columns
            t_total = [t_taskcue; t_total; t_ibi];
            
            % keep track of order and timing of trials
            trial_order_timing = cat(1,trial_order_timing, t_total);
            
            % add block to run
            run(rr,bb) = miniblock_shuffled(total_block_i);
            
            % Update total run time
            total_run_time = t_total{end,end};
            
            % Update counters
            itis(1:n_trials_per_block)=[];
            bb = bb + 1;
            total_block_i = total_block_i + 1;
            
            % Add blank/task cue between blocks
            
        end % total run time
        
       
        % Add pre-blank period
        t_postblank = num2cell([0, NaN, NaN, p.exp.run.post_blank_dur, total_run_time+p.exp.run.post_blank_dur]);
        trial_order_timing = [trial_order_timing; t_postblank];
        
        % start a new run
        rr = rr + 1;
    end % run check
    
    session(ses).run = run;
    session(ses).trial_order_timing = trial_order_timing;
    
end % session

% s = struct2cell(session(1).run);
% s = permute(s,[2,3,1]); % change to [runs, blocks, fieldnames]
%
% run_cond_order = {};
%
% for jj = 1:size(s,1)
%     run_cond_order{jj} = repelem(cell2mat(s(jj,:,2)), 8);
% end
%
for ses = 1:length(session)
   

    newRun_idx = [find(cell2mat(session(ses).trial_order_timing(:,5))==0); length(session(ses).trial_order_timing(:,5))];
    for ii = 1:(length(newRun_idx)-1)
        figure(ii); clf; set(gcf,'Position',[2663,297,665,517])
        timepoints = cell2mat(session(ses).trial_order_timing(newRun_idx(ii): newRun_idx(ii+1),5));
        conditions = cell2mat(session(ses).trial_order_timing(newRun_idx(ii): newRun_idx(ii+1),1));
        
%         subplot(size(session(ses).run,1),1,ii)
        stem(timepoints,conditions); hold on;
        xlabel('time (s)')
        ylabel('stim-task crossing')
        set(gca,'YTick',[0:length(p.exp.stimTaskLabels)+1],'YTickLabel',[{'blank'};p.exp.stimTaskLabels;{'taskcue'}])
        set(gca,'TickDir', 'out')
        title(sprintf('run %d',ii))
    end
end


%% For each subject, we randomize miniblock order within a run,
subject_sessions = [];

for ses = 1:size(session,2)
    for sj = 1:p.exp.total_subjects
        
        for rr = 1:size(session(ses).run,1)
            num_blocks_per_run = length(session(ses).run(rr,:));
            new_block_order    = randperm(num_blocks_per_run,num_blocks_per_run);
            
            subject_sessions(ses,sj).run(rr,1:num_blocks_per_run) = session(ses).run(rr,new_block_order);
        end
    end
end

save(fullfile(vcd_rootPath,'workspaces','info',sprintf('subject_sessions_%s.mat',datestr(now,30))),'subject_sessions')

return
