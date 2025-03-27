function time_table_master = vcd_createRunTimeTables(p)

% Define event ID within trial for single and double epoch trial types
trial_ID_single_epoch = [p.exp.miniblock.trial_start_ID, p.exp.miniblock.spatial_cue_ID, p.exp.miniblock.stim_epoch1_ID, p.exp.miniblock.response_ID];
trial_ID_double_epoch = [p.exp.miniblock.trial_start_ID, p.exp.miniblock.spatial_cue_ID, p.exp.miniblock.stim_epoch1_ID, p.exp.miniblock.delay_ID, p.exp.miniblock.stim_epoch2_ID, p.exp.miniblock.response_ID];

% Preallocate space
time_table_master = [];
tbl_nrows = 1000; % pick an arbitrary large number of rows (otherwise table uses 0 for missing values and we don't want that)
tr_in_frames = (p.exp.TR*p.stim.presentationrate_hz);
for sj = 1:p.exp.total_subjects
    fprintf('\nSUBJECT %03d..',sj)
    
    session_time_table = [];
    
    for ses = 1:p.exp.n_sessions
        fprintf('SESSION %02d..\n',ses)
        
        subj_time_table = [];
        
        % Loop over runs
        for rr = 1:p.exp.n_runs_per_session
            
            % copy condition order table headers and scrub content
            sz = [tbl_nrows size(p.trials(1,:),2)];
            varTypes = varfun(@class,p.trials(1,:),'OutputFormat','cell');
            time_table = table('Size',sz,'VariableTypes',varTypes, 'VariableNames',p.trials.Properties.VariableNames);
            for vt = 1:length(varTypes)
                if strcmp(varTypes(vt),'double')
                    time_table.(p.trials.Properties.VariableNames{vt}) = NaN(tbl_nrows,1);
                end
            end
            
            % add column for event timing and id
            time_table.event_start = NaN(tbl_nrows,1);
            time_table.event_dur   = NaN(tbl_nrows,1);
            time_table.event_end   = NaN(tbl_nrows,1);
            time_table.event_id    = NaN(tbl_nrows,1);
            time_table.event_name  = repmat({''},tbl_nrows,1);
            time_table.subj_nr     = sj.*ones(tbl_nrows,1);
            time_table.subj_within_run_block_nr = NaN(tbl_nrows,p.exp.total_subjects);
            
            
            % Find corresponding trials (rows in table) and block nrs for
            % this particular subject, run, and session. The t_trial table
            % will be inserted into a bigger time table that includes other
            % event types, like iti, cues, etc.
            idx0      = find((p.trials.session_nr==ses) & (p.trials.run_nr==rr));
            t_trial   = p.trials(idx0,:);
            block_nrs = p.trials.subj_within_run_block_nr(idx0,sj);
            
            % shuffle ITIs prior to allocating them after a trial
            itis      = shuffle_concat(p.exp.trial.ITI,1);
            
            % reset counter
            total_run_frames = 0;
            run_finished     = 0;
            
            % Add pre-blank period
            if total_run_frames == 0
                table_idx = 1;
                time_table.event_start(table_idx)    = total_run_frames;
                time_table.event_dur(table_idx)      = p.exp.run.pre_blank_dur;
                time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                time_table.event_id(table_idx)       = 0;
                time_table.event_name(table_idx)     = {'pre-blank'};
                
                time_table.run_nr(table_idx)         = rr;
                time_table.session_nr(table_idx)     = ses;
                total_run_frames = total_run_frames + time_table.event_end(table_idx);
                
                % update time table idx
                table_idx = table_idx+1;
            end
            
            % reset block vector counter
            block_vec_idx = 1;
            
            while 1
                
                if isempty(t_trial) || run_finished
                    run_finished = 1;
                    break;
                end
                
                % restock ITIs if needed
                if isempty(itis) || length(itis)<2
                    itis = cat(2,itis,shuffle_concat(p.exp.trial.ITI,1));
                end
                
                % check if we use single or 2 stimulus patches
                if strcmp(t_trial.stim_class_name{1},'ns')
                    single_stim_flag = 1;
                else
                    single_stim_flag = 0;
                end
                
                % get current miniblock nr
                curr_block = block_nrs(block_vec_idx);
                curr_across_session_miniblock_nr = t_trial.across_session_miniblock_nr(1);
                curr_within_ses_block_nr = t_trial.within_ses_block_nr(1);
                curr_subj_within_run_block_nr = t_trial.subj_within_run_block_nr(1,sj);
                
                % first trial of the block has a task cue
                if t_trial.miniblock_local_trial_nr(1) == 1 || ...
                        (block_vec_idx == 1) || (curr_block ~= block_nrs(block_vec_idx-1)) % new block
                    time_table.event_start(table_idx)    = total_run_frames;
                    time_table.event_dur(table_idx)      = p.exp.miniblock.task_cue_dur;
                    time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)       = p.exp.miniblock.task_cue_ID;
                    time_table.event_name(table_idx)     = {'task-cue'};
                    
                    % add ses, run nr
                    time_table.session_nr(table_idx) = ses;
                    time_table.run_nr(table_idx) = rr;
                    total_run_frames = time_table.event_end(table_idx);
                    time_table.across_session_miniblock_nr(table_idx) = curr_across_session_miniblock_nr;
                    time_table.within_ses_block_nr(table_idx)         = curr_within_ses_block_nr;
                    time_table.subj_within_run_block_nr(table_idx,sj) = curr_subj_within_run_block_nr;
                    
                    table_idx = table_idx+1;
                end
                
                % Get trial IDs of current trial
                if t_trial.trial_type(1) ==1
                    trial_IDs = trial_ID_single_epoch;
                else
                    trial_IDs = trial_ID_double_epoch;
                end
                
                % Loop over events within trial
                for id = 1:length(trial_IDs)
                    
                    % event start is the same for all IDs
                    time_table.event_start(table_idx) = total_run_frames;
                    
                    % add ses, run nr
                    time_table.session_nr(table_idx) = ses; 
                    time_table.run_nr(table_idx) = rr;
                    time_table.across_session_miniblock_nr(table_idx) = curr_across_session_miniblock_nr;
                    time_table.within_ses_block_nr(table_idx)         = curr_within_ses_block_nr;
                    time_table.subj_within_run_block_nr(table_idx,sj) = curr_subj_within_run_block_nr;
                    
                    % Add individual trial events
                    switch trial_IDs(id)
                        
                        case p.exp.miniblock.trial_start_ID % trial start (dot thickening)
                            time_table.event_dur(table_idx)      = p.exp.trial.start_cue_dur;
                            time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                            time_table.event_id(table_idx)       = trial_IDs(id);
                            time_table.event_name(table_idx)     = {'trial-start'};
                            total_run_frames = time_table.event_end(table_idx);
                            table_idx   = table_idx+1;
                            
                        case p.exp.miniblock.spatial_cue_ID  % spatial cue
                            time_table.event_dur(table_idx)      = p.exp.trial.spatial_cue_dur;
                            time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                            time_table.event_id(table_idx)       = trial_IDs(id);
                            time_table.event_name(table_idx)     = {'spatial-cue'};
                            total_run_frames = time_table.event_end(table_idx);
                            table_idx   = table_idx+1;
                            
                        case p.exp.miniblock.stim_epoch1_ID % stim epoch 1 - left&right
                            if ~single_stim_flag
                                tmp1 = t_trial([1,2],:);
                                tmp1.event_start    = repmat(total_run_frames,2,1);
                                tmp1.event_dur      = repmat(p.exp.trial.stim_array_dur,2,1);
                                tmp1.event_end      = repmat(tmp1.event_start(1) + tmp1.event_dur(1),2,1);
                                tmp1.event_id       = repmat(trial_IDs(id),2,1);
                                tmp1.event_name     = repmat({'stim1'},2,1);
                                tmp1.subj_nr        = repmat(sj,2,1);
                                
                                % replace time table with trial table
                                time_table([table_idx,table_idx+1],:) = tmp1;
                                
                                total_run_frames = time_table.event_end(table_idx+1);
                                table_idx   = table_idx+2;
                                
                                if t_trial.trial_type(1) == 1 % if there is a second stim epoch, we wait with updating trial table
                                    t_trial([1,2],:) = []; % otherwise, we remove trial from trial table
                                end
                                
                            else % stim epoch 1 - center
                                tmp1 = t_trial(1,:);
                                tmp1.event_start    = total_run_frames(end);
                                tmp1.event_dur      = p.exp.trial.stim_array_dur;
                                tmp1.event_end      = tmp1.event_start(1) + tmp1.event_dur(1);
                                tmp1.event_id       = trial_IDs(id);
                                tmp1.event_name     = {'stim1'};
                                tmp1.subj_nr        = sj;
                                
                                time_table(table_idx,:) = tmp1;
                                total_run_frames = time_table.event_end(table_idx);
                                
                                table_idx   = table_idx+1;
                                
                                if t_trial.trial_type(1) == 1 % if there is a second stim epoch, we wait with updating trial table
                                    t_trial(1,:) = [];
                                end
                            end
                            
                        case p.exp.miniblock.stim_epoch2_ID % stim epoch 2 (after delay) - left&right
                            
                            if ~single_stim_flag
                                tmp2 = tmp1;
                                
                                tmp2.event_start    = repmat(total_run_frames(end),2,1);
                                tmp2.event_dur      = repmat(p.exp.trial.stim_array_dur,2,1);
                                tmp2.event_end      = repmat(tmp2.event_start(1) + tmp2.event_dur(1),2,1);
                                tmp2.event_id       = repmat(trial_IDs(id),2,1);
                                tmp2.event_name     = repmat({'stim2'},2,1);
                                tmp2.unique_im_nr   = NaN(2,1);
                                tmp2.orient_dir     = cell2mat(tmp2.stim2_delta);
                                tmp2.stim2_delta    = cell(2,1);
                                
                                time_table([table_idx, table_idx+1],:) = tmp2;
                                
                                total_run_frames = time_table.event_end(table_idx+1);
                                
                                table_idx = table_idx+2;
                                
                                % Now we update trial
                                t_trial([1,2],:) = [];
                                
                            else % stim epoch 2 - center
                                tmp2 = tmp1; %time_table(t_idx-2,:);
                                tmp2.event_start    = total_run_frames(end);
                                tmp2.event_dur      = p.exp.trial.stim_array_dur;
                                tmp2.event_end      = tmp2.event_start(1) + tmp2.event_dur(1);
                                tmp2.event_id       = trial_IDs(id);
                                tmp2.event_name     = {'stim2'};
                                tmp2.unique_im_nr   = NaN;
                                tmp2.stim2_delta     = tmp2.stim2_delta{1};

                                time_table(table_idx,:) = tmp2;
                                total_run_frames = time_table.event_end(table_idx);
                                
                                table_idx = table_idx+1;
                                
                                % Now we update trial
                                t_trial(1,:) = [];
                            end
                            
                        case p.exp.miniblock.response_ID % response
                            time_table.event_dur(table_idx)      = p.exp.trial.response_win_dur;
                            time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                            time_table.event_id(table_idx)       = p.exp.miniblock.response_ID;
                            time_table.event_name(table_idx)     = {'response-cue'};
                            total_run_frames = time_table.event_end(table_idx);
                            
                            table_idx = table_idx+1;
                            
                        case p.exp.miniblock.delay_ID % delay
                            time_table.event_dur(table_idx)      = p.exp.trial.response_win_dur;
                            time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                            time_table.event_id(table_idx)       = trial_IDs(id);
                            time_table.event_name(table_idx)     = {'delay'};
                            total_run_frames = time_table.event_end(table_idx);
                            
                            table_idx = table_idx+1;
                    end
                    
                end % trial IDs
                
                %% Check for inter-trial, inter-block
                if single_stim_flag && (block_vec_idx == length(block_nrs)) % one more trial to go
                    % add ITI_ID
                    ITI =  itis(1); % grab the first one from the shuffled list
                    itis(1) = []; % and remove it
                    time_table.event_start(table_idx)  = total_run_frames;
                    time_table.event_dur(table_idx)    = ITI;
                    time_table.event_name(table_idx)   = {'ITI'};
                    time_table.event_id(table_idx)     = p.exp.miniblock.ITI_ID;
                    time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    
                    % add ses, run nr
                    time_table.session_nr(table_idx) = ses;
                    time_table.run_nr(table_idx) = rr;
                    
                    total_run_frames = time_table.event_end(table_idx);
                    
                    table_idx = table_idx+1;
                    
                    % last trial of the last block    
                elseif ~single_stim_flag && ((block_vec_idx+1) == length(block_nrs))
                    
                    run_finished = 1;
                    break;
                    
                    % last trial of the last block    
                elseif single_stim_flag && (block_vec_idx == length(block_nrs))
                    run_finished = 1;
                    break;
                
                    % more trials to go
                elseif block_vec_idx < length(block_nrs) 
                    
                    if ~single_stim_flag
                        next_block_id = block_nrs(block_vec_idx+2);
                    else
                        next_block_id = block_nrs(block_vec_idx+1);
                    end
                    
                    if curr_block ~= next_block_id  % next trial is from different block (tt is already updated), so this is the last trial of the current block
                        
                        % add IBI
                        time_table.event_start(table_idx)  = total_run_frames;
                        
                        sampled_IBI = p.exp.miniblock.IBI(randi(length(p.exp.miniblock.IBI),1));
                        IBI =  sampled_IBI + (1 - mod(time_table.event_start(table_idx)+sampled_IBI,1)); % round out to a full second
                        
                        time_table.event_start(table_idx)  = total_run_frames;
                        time_table.event_dur(table_idx)    = IBI;
                        time_table.event_name(table_idx)   = {'IBI'};
                        time_table.event_id(table_idx)     = p.exp.miniblock.IBI_ID;
                        time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                        % add ses, run nr
                        time_table.session_nr(table_idx) = ses;
                        time_table.run_nr(table_idx) = rr;

                        total_run_frames = time_table.event_end(table_idx);
                        
                        table_idx = table_idx+1;
                        
                    elseif curr_block == next_block_id % % next trial is from the same block (tt is already updated), so not last trial of the block
                        % add ITI
                        ITI =  itis(1); % grab the first one from the shuffled list
                        itis(1) = []; % and remove it
                        time_table.event_start(table_idx)  = total_run_frames;
                        time_table.event_dur(table_idx)    = ITI;
                        time_table.event_name(table_idx)   = {'ITI'};
                        time_table.event_id(table_idx)     = p.exp.miniblock.ITI_ID;
                        time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                        % add ses, run nr
                        time_table.session_nr(table_idx) = ses;
                        time_table.run_nr(table_idx) = rr;
                        time_table.subj_nr(table_idx) = sj;
                    
                        total_run_frames = time_table.event_end(table_idx);
                        
                        table_idx = table_idx+1;
                    end
                    
                    if ~single_stim_flag
                        block_vec_idx = block_vec_idx+2;
                    else
                        block_vec_idx = block_vec_idx+1;
                    end
                end % ITI/IBI if statement
                
            end % while loop
            
            % add post-run blank
            if run_finished || block_vec_idx > length(block_nrs) || (block_vec_idx+1)==length(block_nrs) % last trial of block and last block
                
                % remove last ITI
                if strcmp(time_table.event_name(table_idx-1),'ITI')
                    time_table(table_idx-1,:) = [];
                    table_idx = table_idx-1;
                end
                
                % round out post blank period to finish on full TR
                total_run_time2 = total_run_frames + p.exp.run.post_blank_dur;
                round_me_out = (ceil(total_run_time2/ tr_in_frames) - (total_run_time2/tr_in_frames)).*tr_in_frames;
                postblank_dur_fullTR = p.exp.run.post_blank_dur + round_me_out;
                
                % add post_blank
                time_table.event_start(table_idx)  = total_run_frames;
                time_table.event_dur(table_idx)    = postblank_dur_fullTR;
                time_table.event_name(table_idx)   = {'post-blank'};
                time_table.event_id(table_idx)     = 0;
                time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                % add ses, run nr
                time_table.session_nr(table_idx) = ses;
                time_table.run_nr(table_idx) = rr;
                time_table.subj_nr(table_idx) = sj;
                
                total_run_frames = time_table.event_end(table_idx);

                assert(nearZero(mod(total_run_frames,p.exp.TR)))
                run_finished = 0; % reset flag

                % trim time table
                time_table2 = time_table;
                if strcmp(time_table2.event_name(table_idx),{'post-blank'})
                    time_table2((table_idx+1):end,:) = [];
                end
            
                % add run to master
                subj_time_table = cat(1, subj_time_table, time_table2);
                clear time_table2
            end
        end
        
        session_time_table = cat(1, session_time_table, subj_time_table);
    end
    time_table_master = cat(1,time_table_master,session_time_table);
end

if p.store_params
    fprintf('\n[%s]:Storing subject session data..\n',mfilename)
    saveDir = fileparts(fullfile(p.stim.fix.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir,sprintf('subject_sessions_%s_%s.mat',p.disp.name, datestr(now,30))),'time_table_master')
end

return