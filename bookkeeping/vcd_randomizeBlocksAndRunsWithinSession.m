function [condition_master, randomization_filename] = vcd_randomizeBlocksAndRunsWithinSession(params, condition_master, session_env)

%% For every session/session type: sort by run_nr, then block_nr
unique_sessions      = unique(condition_master.session_nr);

[~,session_types] = vcd_getSessionEnvironmentParams(params, session_env);


for ses = 1:length(unique_sessions)
    for st = 1:length(session_types)
        if ~isnan(session_types(ses,st))            
           
            
            %% SHUFFLE RUN AND BLOCK ORDER WITHIN RUNS
            ses_idx           = (condition_master.session_nr==ses & condition_master.session_type==st);
            ses_blocks        = condition_master.block_nr(ses_idx);
            ses_trials        = condition_master.trial_nr(ses_idx);
            ses_trialtype     = condition_master.trial_type(ses_idx);
            ses_crossing_vec  = condition_master.crossing_nr(ses_idx);

            % get unique block nrs and associated trial types
            [unique_blocks, block_start_idx] = unique(ses_blocks);
            crossings_unique_blocks          = ses_crossing_vec(block_start_idx);
            unique_trialtypes                = ses_trialtype(block_start_idx);
            n_trialtype1 = sum(unique_trialtypes==1); % single
            n_trialtype2 = sum(unique_trialtypes==2); % double

            % Check if we have the expect nr of single and double stim
            % presentation blocks
            assert(isequal(n_trialtype1+n_trialtype2,length(unique_blocks)))
            
            % set up block structure
            runs_to_fill = runs_per_session(ses,st); 
            
            % separate single and double stim blocks:
            s_idx = (unique_trialtypes==1);
            d_idx = (unique_trialtypes==2);
            s_blocks = unique_blocks(s_idx);
            d_blocks = unique_blocks(d_idx);

            assert(isequal(length(s_blocks)+length(d_blocks),length(unique_blocks)))
            
            % Shuffle block order within a session
            while 1
                
                while 1
                    bb_rnd = randperm(length(unique_blocks),length(unique_blocks));
                    crossings_shuffled = crossings_unique_blocks(bb_rnd);
                    s_blocks_shuffled   = unique_blocks(bb_rnd(s_idx));
                    d_blocks_shuffled   = unique_blocks(bb_rnd(d_idx));
                    block_start_idx_shuffled = block_start_idx(bb_rnd);
                    trialtypes_shuffled = unique_trialtypes(bb_rnd);
                    shuffle_ok = [];
                    % for every 4 executive single stimulus presentation blocks,
                    %  we check their crossing, if they have the same crossing
                    %  we shuffle the order of all single stimulus presentation
                    %  blocks
                    block_cnt = 1:4:length(crossings_shuffled);
                    for nn = 1:length(block_cnt)
                        if (block_cnt(nn)+3) <= length(crossings_shuffled)
                            tmp_blocks = crossings_shuffled(block_cnt(nn):(block_cnt(nn)+3));
                            else
                            tmp_blocks = crossings_shuffled(block_cnt(nn):end);
                        end
                        nr_s_blocks = sum(ismember(tmp_blocks,s_blocks));
                        nr_d_blocks = sum(ismember(tmp_blocks,d_blocks));
                        if length(unique(tmp_blocks)) >= length(tmp_blocks)-1
                            shuffle_ok(nn) = true;
                        else
                            break;
                        end
                        if nr_d_blocks > 0
                            shuffle_ok(nn) = true;
                        else
                            break;
                        end

                    end
                    if sum(shuffle_ok) == length(block_cnt)-1
                        break;
                    end
                end
                    
            
                % Add blocks to runs within this session
                run_matrix = zeros(runs_to_fill,7); % max 7 possible slots
                
                % reset counters
                bb_glbl_cnt = 1;
                bb_cnt = 1;
                rr_cnt = 1;
                run_dur = run_dur_min + sum(IBIs);
                block_dur = [params.exp.block.total_single_epoch_dur, params.exp.block.total_double_epoch_dur];
                total_run_dur = [];
                slack = 10;
                while 1
                    
                    if bb_cnt > size(run_matrix,2)
                        total_run_dur(rr_cnt) = run_dur;
                        rr_cnt  = rr_cnt +1;   % count next run
                        bb_cnt  = 1;           % reset within run block counter when moving to the next run
                        run_dur = run_dur_min + sum(IBIs); % reset run duration when moving to the next run
                    end
                    
                    run_dur_tmp = run_dur + block_dur(trialtypes_shuffled(bb_glbl_cnt));
                    
                    if run_dur_tmp-slack > run_dur_max
                        total_run_dur(rr_cnt) = run_dur;
                        % go to next run
                        rr_cnt  = rr_cnt +1;   % count next run
                        bb_cnt  = 1;           % reset within run block counter when moving to the next run
                        run_dur = run_dur_min + sum(IBIs); % reset run duration when moving to the next run
                    else
                        run_dur = run_dur_tmp;
                    end
                    
                    % add block to run
                    run_matrix(rr_cnt, bb_cnt) = crossings_shuffled(bb_glbl_cnt);
                    
                    % update counters
                    bb_glbl_cnt = bb_glbl_cnt+1;
                    bb_cnt      = bb_cnt+1;
                    if bb_glbl_cnt > length(crossings_shuffled)
                        break;
                    end
                    
                end
            
                % if the run is 30 s shorter than expected, we try again.
                run_too_short = (run_dur_max-total_run_dur) > (30*params.stim.presentationrate_hz);
                
                if sum(run_too_short)==0
                    break;
                end
            
                % ensure we added all the blocks
                assert(length(run_matrix(run_matrix>0))==length(crossings_shuffled))
                assert(length(unique(run_matrix(run_matrix>0)))==length(unique(crossings_shuffled)))
            
            end
            
            % Shuffle trials within each block
            block_cnt = 0; 
            trial_order_local_shuffled  = cell(size(run_matrix));
            trial_order_global_shuffled = cell(size(run_matrix));
            fprintf('[%s]: \n',mfilename)
            for run_idx = 1:size(run_matrix,1)
                % Show user the block order within a run
                curr_blocks = run_matrix(run_idx,:);
                curr_blocks = curr_blocks(curr_blocks>0);
                fprintf('Run %02d, block order:',run_idx)
                disp(curr_blocks)

                for bb_idx = 1:length(curr_blocks)
                    block_cnt = block_cnt +1;
                    
                    if trialtypes_shuffled(block_cnt)==1
                        nr_trials = params.exp.block.n_trials_single_epoch;
                    elseif trialtypes_shuffled(block_cnt)==2
                        nr_trials = params.exp.block.n_trials_double_epoch;
                    end
                    
                    trial_order = ses_trials(block_start_idx_shuffled):(ses_trials(block_start_idx_shuffled)+nr_trials-1);
                    trial_order_local_shuffled{run_idx,bb_idx} = shuffle_concat(1:nr_trials,1);
                    trial_order_global_shuffled{run_idx,bb_idx} = shuffle_concat(trial_order,1);
                    % Show user the trial order within a block
                    fprintf('\nBlock %02d, trial order:',curr_blocks(bb_idx))
                    disp(trial_order_local_shuffled{run_idx,bb_idx})
                end
            end

        end % ~isnan session type
    end % session type
end % session nr

% get a new condition master with shuffled blocks? 
% or just save the file and use those in the vcd_createRunTimeTableMaster.m

% Sort block order by ascending nr
condition_master0 = condition_master;
for ses = unique_sessions
    for st = unique_session_types
        if ~isnan(session_types(ses,st))
            all_runs_bool = (condition_master0.session_nr==ses & condition_master0.session_type==st);
            all_runs_idx = find(all_runs_bool);
            [~,sort_idx0]  = sort(condition_master0.run_nr(all_runs_idx));
            condition_master0(all_runs_bool,:) = condition_master0(all_runs_idx(sort_idx0),:);
            
            unique_runs = unique(condition_master0.run_nr(all_runs_bool,:));
            
            for run_idx = 1:length(unique_runs)                
                all_trial_idx2    = find(condition_master0.session_nr==ses & condition_master0.session_type==st & condition_master0.run_nr==run_idx);
                single_run_blocks = condition_master0.block_nr(all_trial_idx2);
                [~,sort_idx1]     = sort(single_run_blocks);
                condition_master0(all_trial_idx2,:) = condition_master0(all_trial_idx2(sort_idx1),:);
            end
        end
    end
end
