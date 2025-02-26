function [subject_sessions,session_master] = vcd_createSessions(p,load_params,store_params)

if ~exist('store_params','var') || isempty(store_params)
    store_params = true;
end

if ~exist('load_params','var') || isempty(load_params)
    load_params = true;
end

if load_params
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('subject_sessions*.mat')));
    if ~isempty(d)
        load(fullfile(d(end).folder,d(end).name),'subject_sessions');
    else
        error('[%s]: Can''t find subject sessions file!', mfilename)
    end
else
    % Create subject sessions

    % check if trial struct is already defined
    if isempty(p.trials) || ~isfield(p,'trials')
        if p.load_params
            
            % load trial info
            d = dir(fullfile(vcd_rootPath,'workspaces','info','trials*.mat'));
            
            if ~isempty(d(end).name) && length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
                load(fullfile(d(end).folder,d(end).name),'all_trials');
                p.trials = all_trials;
            else
                error('[%s]: Can''t find trial.mat files! Please check', mfilename);
            end
            
        else % or create it
            p.trials = vcd_makeTrials(p,load_params,store_params);
        end
    end
    
    if isempty(p.exp) ||  ~isfield(p,'exp')
        if p.load_params
            
            % load trial info
            d = dir(fullfile(vcd_rootPath,'workspaces','info','exp_session*.mat'));
            
            if ~isempty(d(end).name) && length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
                load(fullfile(d(end).folder,d(end).name),'exp_session');
                p.exp = exp_session;
            else
                error('[%s]: Can''t find exp session .mat files! Please check', mfilename);
            end
            
        else % or create it
            p.exp = vcd_getSessionParams(load_params,store_params);
        end
    end
    
    %% per task/stim crossing: get nr of trials per miniblock
    
    % Different trial order per block (already accomplished in vcd_makeTrials.m)
    % Across a single session, each subject will experience the same miniblocks and unique images.
    tolerance = 1.^-10;
    nearZero = @(x,tol) abs(x) < tol;

    
    % reset RNG
    rng('shuffle','twister');
    
    % Time we shave off between end of run, to check if we can fit another run
    fudge_factor = 10; % seconds
    
    %% %%%%%%%%%%%%%%%%%% EK TEMPORARY HACK
    p.exp.n_runs_per_session = 8;
    %%%%%%%%%%%%%%%%%%
    
    session_master = [];
    miniblock = [];
    
    for ses = 1:p.exp.n_sessions
        fprintf('\nSESSION %d:',ses)
        
        t_total = [];
        % trial_order_timing:
        % 1: run nr,
        % 2: block nr,
        % 3: blockID,
        % 4: unique image nr,
        % 5: spatial cue dir,
        % 6: event duration,
        % 7: total run duration
        trial_order_timing = {}; %cell(1,7);
        
        % Same blocks / same trials / same retinal image per subject
        miniblock_distr = p.exp.ses_blocks(:,:,ses);
        
        stimtask_tracker = ones(size(miniblock_distr));
        
        bb = 1;
        
        for sc = 1:length(p.exp.stimClassLabels)
            for tc = 1:length(p.exp.taskClassLabels)
                
                if ses >= p.exp.session.task_start(tc)
                    
                    n_blocks = miniblock_distr(sc,tc);
                    
                    if (n_blocks > 0)
                        
                        fprintf('\n%s\t x %s\t\t:\t %d block(s)',p.exp.stimClassLabels{sc},p.exp.taskClassLabels{tc},n_blocks)
                        
                        curr_miniblocks = [stimtask_tracker(sc,tc) : (stimtask_tracker(sc,tc)+n_blocks-1)];
                        
                        for ii = 1:length(curr_miniblocks)
                            miniblock(bb).name                      = sprintf('%s-%s',p.exp.taskClassLabels{tc},p.exp.stimClassLabels{sc});
                            miniblock(bb).ID                        = find(strcmp(miniblock(bb).name,p.exp.stimTaskLabels));
                            miniblock(bb).within_session_repeat     = sprintf('%d-out-of-%d',curr_miniblocks(ii),n_blocks);
                            miniblock(bb).trial                     = p.trials.stim(sc).tasks(tc).miniblock(curr_miniblocks(ii)).trial;
                            
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
        
        
        % Make sure that runs have at least 2 miniblocks that use double-epoch
        % trials, otherwise some runs will be much longer than others..
        while 1
            run_ok = 0;
            run_not_ok = false;
            % Shuffle blocks
            shuffle_idx = randperm(length(miniblock),length(miniblock));
            miniblock_shuffled = miniblock(shuffle_idx);
            
            tt = struct2cell(miniblock_shuffled);
            tm = cell2mat(tt(end,:));
            for jj = 1:6:length(tm)-1
                if sum(tm(jj:(jj+5))==2)>2
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
        
        
        itis = shuffle_concat(p.exp.trial.ITI,1);
        
        % Divide blocks into runs
        rr = 1;
        total_block_i = 1;
        
        while rr <= p.exp.n_runs_per_session && ...
                (total_block_i <= length(miniblock_shuffled))
            
            % reset counter
            total_run_time = 0;
            bb = 1;
            
            if total_run_time == 0 % Add pre-blank period

                total_run_time = p.exp.run.pre_blank_dur;
                t_preblank = [num2cell([rr, zeros(1,7)]); ...
                    num2cell([rr, 0, 0, 0, NaN, 0, p.exp.run.pre_blank_dur, total_run_time])];
                
                trial_order_timing = [trial_order_timing; t_preblank];
                
                run(rr).block(bb).name       = 'pre-blank';
                run(rr).block(bb).ID         = 0;
                run(rr).block(bb).within_session_repeat   = NaN;
                run(rr).block(bb).trial      = [];
                run(rr).block(bb).trial_type = []; % single or double epoch
                run(rr).block(bb).timing     = cell2table(t_preblank, 'VariableNames',{'run','block','stimtaskID','unique_im','spatial_cue','onset_time','event_dur','run_time'});
                
                bb=bb+1;
            end
            
            if isempty(itis)
                itis = shuffle_concat(p.exp.trial.ITI,3);
            end
            
            
            while (total_run_time < (p.exp.run.total_run_dur - p.exp.run.post_blank_dur - fudge_factor)) && ...
                    (total_block_i <= length(miniblock_shuffled)) && bb <= p.exp.run.miniblocks_per_run+1
                
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
                t_taskcue = num2cell([rr, bb-1, miniblock_shuffled(total_block_i).ID, p.exp.miniblock.task_cue_ID, NaN, total_run_time, p.exp.miniblock.task_cue_dur,  total_run_time+p.exp.miniblock.task_cue_dur]);
                
                
                % Image IDs
                fn = fieldnames(miniblock_shuffled(total_block_i).trial);
                tmp = struct2cell(miniblock_shuffled(total_block_i).trial);
                unique_im_nr_idx = strcmp(fn,{'unique_im_nr'});
                t_im_id = cell(n_trials_per_block*2,1);
                tmp2 = squeeze(tmp(unique_im_nr_idx,:,:));
                t_im_id(1:2:end,:) = mat2cell(cell2mat(tmp2),ones(1,n_trials_per_block));
                t_im_id(2:2:end,:) = {p.exp.miniblock.ITI_ID};
                
                % Thickening dir
                thicking_nr_idx = strcmp(fn,{'covert_att_loc_cue'});
                t_spat_cue =  cell(n_trials_per_block*2,1);
                t_spat_cue(1:2:end,:) = squeeze(tmp(thicking_nr_idx,:,1));
                t_spat_cue(2:2:end,:) = {NaN};
                
                % Assign stimtask ID
                t_st_id = num2cell(repmat(miniblock_shuffled(total_block_i).ID,2*n_trials_per_block,1));
                t_st_id(2:2:end,1) = {0};
                
                % get timing
                if miniblock_shuffled(total_block_i).trial_type == 1
                    t_trial_time = [repmat(p.exp.trial.single_epoch_dur,1,n_trials_per_block); itis(1:n_trials_per_block)];
                else
                    t_trial_time = [repmat(p.exp.trial.double_epoch_dur,1,n_trials_per_block); itis(1:n_trials_per_block)];
                end
                
                t_trial_time = t_trial_time(:);
                t_trial_time(end,end) = 0; % we don't need the final ITI for the last trial of the block, because we have an IBI
                t_cumul_time = t_taskcue{end} + cumsum(t_trial_time);
                t_onset_time = t_cumul_time - t_trial_time;
                
                assert( all(nearZero(mod((t_onset_time./p.stim.fps),1),tolerance)))
                assert( all(nearZero(mod((t_cumul_time./p.stim.fps),1),tolerance)))

                % concat columns
                t_total = [num2cell(repmat(rr,2*n_trials_per_block,1)), num2cell(repmat(bb-1,2*n_trials_per_block,1)), ...
                    t_st_id, t_im_id, t_spat_cue, num2cell(t_onset_time), num2cell(t_trial_time), num2cell(t_cumul_time)];
                
                
                % Add inter-block-interval (1 TR + whatever is left)
                sampled_IBI = p.exp.miniblock.IBI(randi(length(p.exp.miniblock.IBI),1));
                IBI =  sampled_IBI + (1 - mod(t_total{end,end}+sampled_IBI,1));
                
                assert( all(nearZero(mod((IBI./p.stim.fps),1),tolerance)))
                
                t_ibi = {rr, bb-1, 0, p.exp.miniblock.IBI_ID, NaN, t_total{end,end}, IBI, t_total{end,end} + IBI};
                
                
                % concat more columns
                t_total = [t_taskcue; t_total; t_ibi];
                
                % keep track of order and timing of trials
                trial_order_timing = cat(1,trial_order_timing, t_total);
                
                % add timing to miniblock struct
                miniblock_shuffled(total_block_i).timing = cell2table(t_total, ...
                    'VariableNames',{'run','block','stimtaskID','unique_im','spatial_cue','onset_time','event_dur','run_time'});
                
                % add block to run
                run(rr).block(bb) = miniblock_shuffled(total_block_i);
                
                % Update total run time
                total_run_time = trial_order_timing{end,end};
                
                % Update counters
                itis(1:(n_trials_per_block-1))=[];
                bb = bb + 1;
                total_block_i = total_block_i + 1;
                
            end % total run time
            
            % get last IBI and change to round out block to a full second +
            % postblank period
            bad_IBI = t_total{end,end-1};
            total_run_time = total_run_time - bad_IBI; % revert total run time
            
            % round out post blank period to finish on full TR
            total_run_time2 = total_run_time+p.exp.run.post_blank_dur;
            round_me_out = (ceil(total_run_time2/p.exp.TR) - (total_run_time2/p.exp.TR)).*p.exp.TR;
            
            postblank_dur_fullTR = p.exp.run.post_blank_dur + round_me_out;
            
            % Add post-blank period and update total duration
            run(rr).block(bb-1).timing.event_dur(end) = 0; % reset postblank_dur to 0
            assert(nearZero(run(rr).block(bb-1).timing.run_time(end-1)-total_run_time, tolerance));
            
            run(rr).block(bb-1).timing.run_time(end) = total_run_time; % restate previous trial total_run_time;
            
            % do the same for trial_order_timing matrix
            trial_order_timing{end,end-1} = 0;
            trial_order_timing{end,end} = total_run_time;
            t_postblank = num2cell([rr, bb-1, 0, 0, NaN, total_run_time, postblank_dur_fullTR, total_run_time+postblank_dur_fullTR]);
            
            trial_order_timing = [trial_order_timing; t_postblank];
            
            run(rr).block(bb).name       = 'post-blank';
            run(rr).block(bb).ID         = 0;
            run(rr).block(bb).within_session_repeat   = NaN;
            run(rr).block(bb).trial      = [];
            run(rr).block(bb).trial_type = []; % single or double epoch
            run(rr).block(bb).timing     = cell2table(t_postblank, 'VariableNames',{'run','block','stimtaskID','unique_im','spatial_cue','onset_time','event_dur','run_time'});
            
            assert(nearZero(mod(run(rr).block(bb).timing.run_time(end),p.exp.TR),tolerance))
            
            % go to the next run
            rr = rr + 1;
            
        end % run nr check
        
        session_master(ses).run = run;
        session_master(ses).trial_order_timing = cell2table(trial_order_timing, 'VariableNames',{'run','block','stimtaskID','unique_im','spatial_cue','onset_time','event_dur','run_time'});
        
    end % session
    
    %
    % foo = {};
    % for xx= 1:8
    %     for yy = 1:8
    %         foo = cat(1, foo, table2cell(session_master(ses).run(xx).block(yy).timing));
    %     end
    % end
    %
    % assert(isequal(foo,table2cell(session_master(ses).trial_order_timing)))
    
    for ses = 1:length(session_master)
        
        figure(ses); clf; set(gcf,'Position',[1,1,1920,1080])
        for rr = 1:length(session_master(ses).run)
            timepoints =[]; conditions=[]; imID = [];
            for ii = 1:length(session_master(ses).run(rr).block)
                timepoints = [timepoints;session_master(ses).run(rr).block(ii).timing.run_time];
                conditions = [conditions;session_master(ses).run(rr).block(ii).timing.stimtaskID];
                
                if iscell(session_master(ses).run(rr).block(ii).timing.unique_im(1))
                    yy = session_master(ses).run(rr).block(ii).timing.unique_im{1};
                else
                    yy = session_master(ses).run(rr).block(ii).timing.unique_im(1);
                end
                xx = session_master(ses).run(rr).block(ii).timing.run_time(1);
                imID = cat(1, imID, [xx,yy]);
            end
            imID(2:end-1,2) = 40;
            subplot(ceil(size(session_master(ses).run,2)/2),2,rr)
            stem(timepoints,conditions); hold all;
            stem(imID(:,1),imID(:,2),'r');
            xlabel('time (s)')
            ylabel('stim-task crossing')
            set(gca,'YTick',[0:length(p.exp.stimTaskLabels)+2],'YTickLabel',[{'blank'}; p.exp.stimTaskLabels;{'taskcue'};{'IBI'}])
            set(gca,'TickDir', 'out')
            title(sprintf('run %d',ii))
        end
        if p.store_imgs
            print(gcf,'-dpng','-r300',fullfile(vcd_rootPath,'figs',sprintf('ses%02d_run_order_master',ses)));
        end
    end
    
    
    
    %% For each subject, we randomize miniblock order within a run,
    subject_sessions = [];
    
    for ses = 1:size(session_master,2)
        for sj = 1:p.exp.total_subjects
            
            for rr = 1:size(session_master(ses).run,2)
                num_blocks_per_run = length(session_master(ses).run(rr).block)-2;
                new_block_order    = [1, shuffle_concat(2:(1+num_blocks_per_run),1), length(session_master(ses).run(rr).block)];
                
                %             blockstart_idx = find(session_master(ses).trial_order_timing.stimtaskID == p.exp.miniblock.task_cue_ID);
                %             IBIstart_idx   = find(session_master(ses).trial_order_timing.stimtaskID == p.exp.miniblock.IBI_ID);
                %             runstart_idx = [2, find(diff(session_master(ses).trial_order_timing.run)==1)', size(session_master(ses).trial_order_timing.run,1)-1];
                
                % isolate final IBI
                newidx_lastblack = find(new_block_order==(1+num_blocks_per_run));
                final_IBI_dur = session_master(ses).run(rr).block(end-1).timing.event_dur(end);
                
                if newidx_lastblack ~= (1+num_blocks_per_run)
                    replace_IBI_dur = session_master(ses).run(rr).block(newidx_lastblack).timing.event_dur(end);
                else
                    replace_IBI_dur = final_IBI_dur;
                end
                
                % Block "1" stays the same (pre-blank)
                subject_sessions(ses,sj).run(rr).block(1) = session_master(ses).run(rr).block(new_block_order(1));
                
                % Block "2:7" are corrected (actual experimental trials)
                for nn = 2:(1+num_blocks_per_run)
                    subject_sessions(ses,sj).run(rr).block(nn) = session_master(ses).run(rr).block(new_block_order(nn));
                    
                    % update block nr
                    subject_sessions(ses,sj).run(rr).block(nn).timing.block = (nn-1).*ones(size(subject_sessions(ses,sj).run(rr).block(nn).timing.block));
                    
                    % swap IBI
                    if nn == newidx_lastblack
                        subject_sessions(ses,sj).run(rr).block(nn).timing.event_dur(end) = replace_IBI_dur;
                    elseif nn == (1+num_blocks_per_run)
                        subject_sessions(ses,sj).run(rr).block(nn).timing.event_dur(end) = final_IBI_dur;
                    end
                    
                    % update onset time
                    subject_sessions(ses,sj).run(rr).block(nn).timing.onset_time = [subject_sessions(ses,sj).run(rr).block(nn-1).timing.run_time(end);
                        subject_sessions(ses,sj).run(rr).block(nn-1).timing.run_time(end) + ...
                        cumsum(subject_sessions(ses,sj).run(rr).block(nn).timing.event_dur(1:end-1))];
                    % update total time
                    subject_sessions(ses,sj).run(rr).block(nn).timing.run_time = subject_sessions(ses,sj).run(rr).block(nn).timing.onset_time + ...
                        subject_sessions(ses,sj).run(rr).block(nn).timing.event_dur;
                    
                end
                
                % Last Block stays the same (post-blank period)
                subject_sessions(ses,sj).run(rr).block(new_block_order(end)) = session_master(ses).run(rr).block(new_block_order(end));
                
                subject_sessions(ses,sj).run(rr).block(new_block_order(end)).timing.onset_time = ...
                        subject_sessions(ses,sj).run(rr).block(nn).timing.run_time(end); 
                    
                    % update total time
                    tmptime = [subject_sessions(ses,sj).run(rr).block(end).timing.onset_time + p.exp.run.post_blank_dur];
                    round_me_out = (ceil(tmptime/p.exp.TR) - (tmptime/p.exp.TR)).*p.exp.TR;
            
                    tmptime2 = tmptime + round_me_out;
                    subject_sessions(ses,sj).run(rr).block(end).timing.run_time = tmptime2;

            end
            
              
            
            figure; clf; set(gcf,'Position',[1,1,1920,1080])
            for rr = 1:length(subject_sessions(ses,sj).run)
                timepoints =[]; conditions=[]; imID = [];
                for ii = 1:length(subject_sessions(ses,sj).run(rr).block)
                    timepoints = [timepoints;subject_sessions(ses,sj).run(rr).block(ii).timing.run_time];
                    conditions = [conditions;subject_sessions(ses,sj).run(rr).block(ii).timing.stimtaskID];
                    
                    if iscell(subject_sessions(ses,sj).run(rr).block(ii).timing.unique_im(1))
                        yy = subject_sessions(ses,sj).run(rr).block(ii).timing.unique_im{1};
                    else
                        yy = subject_sessions(ses,sj).run(rr).block(ii).timing.unique_im(1);
                    end
                    xx = subject_sessions(ses,sj).run(rr).block(ii).timing.run_time(1);
                    imID = cat(1, imID, [xx,yy]);
                end
                imID(2:end-1,2) = 40;
                subplot(ceil(size(subject_sessions(ses,sj).run,2)/2),2,rr)
                stem(timepoints,conditions); hold all;
                stem(imID(:,1),imID(:,2),'r');
                xlabel('time (s)')
                ylabel('stim-task crossing')
                set(gca,'YTick',[0:length(p.exp.stimTaskLabels)+2],'YTickLabel',[{'blank'}; p.exp.stimTaskLabels;{'taskcue'};{'IBI'}])
                set(gca,'TickDir', 'out')
                title(sprintf('run %d',rr))
                sgtitle(sprintf('Subject %02d, session %02d',sj, ses));
            end
            if p.store_imgs
                print(gcf,'-dpng','-r300',fullfile(vcd_rootPath,'figs',sprintf('subj%02d_ses%02d_run_order',sj,ses)));
            end
        end
    end

    
    if store_params
        fprintf('\n[%s]:Storing subject session data..\n',mfilename)
        saveDir = fileparts(fullfile(p.stim.fix.infofile));
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir,sprintf('subject_sessions_%s_%s.mat',p.disp.name, datestr(now,30))),'subject_sessions')
    end
end
return
