function [condition_master_shuffled, condition_master_shuffle_idx, session_crossing_matrix, session_block_matrix] ...
                = vcd_randomizeBlocksAndRunsWithinSession(params, condition_master, session_env, varargin)
% VCD function to shuffle stimulus blocks and trials within a block across
% all sessions/session types. And then resave individual subject's runs
% of shuffled condition_masters that has sorted the new order of run nrs,
% block_nrs and trial nrs.
%
%  [condition_master_shuffled, run_matrix,trial_order_local_shuffled, trial_order_global_shuffled] ...
%    = vcd_randomizeBlocksAndRunsWithinSession(params, condition_master, session_env)
%
% INPUTS:
%
%
% OUTPUTS:
%
%
% Written by E Kupers @ UMN 2025/06

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'           , @isstruct);
p0.addRequired('condition_master' , @istable);
p0.addRequired('session_env'      , @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));
p0.addParameter('prefix'          , 'subjXXX_', @ischar);

% Parse inputs
p0.parse(params,condition_master,session_env,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0



%% Define params

% User can allow for slack in the run duration threshold before moving on
% to the next run. We do this to allow for adding one more block and avoid
% ending ending upwith long blank at the end.
%
% To use run_dur as a constraint, we overestimate the total duration for 
% the allocated IBIs by summing the longest 7 IBIs (3780/60 = 63 s).
% For reference: difference between the sum of the shortest 7 and the sum
% of longest 7 IBIs is 1050/60 = 17.5 s. Hence we want to also introduce
% some slack, to add avoid adding too few blocks to a run.
slack                       = 10*params.stim.presentationrate_hz; % time frames we extract from run_dur when allocating blocks. 
max_run_deviation           = 40*params.stim.presentationrate_hz; % time frames we use as cutoff to start over (if run is 40 s shorter than expect run duration, then we reshuffle blocks). 
max_allowable_block_repeats = 2;  % nr of blocks we allow for repeat in their crossing (so finding A-A and later on B-B)
block_dur     = [params.exp.block.total_single_epoch_dur, params.exp.block.total_double_epoch_dur];

% Check inputs
unique_sessions  = unique(condition_master.session_nr);

% Get session parameters depending on whether this is the Behavioral
% experiment or MRI experiment.
[~,session_types,runs_per_session,run_dur_min,run_dur_max,~,IBIs, ~, ~, ~, ...
    ~,nr_blocks_per_run] = vcd_getSessionEnvironmentParams(params, session_env);

% Preallocate space
condition_master_shuffled    = [];
condition_master_shuffle_idx = cell(length(unique_sessions),length(session_types));
session_crossing_matrix      = cell(length(unique_sessions),length(session_types));
session_block_matrix         = cell(length(unique_sessions),length(session_types));

tic;

%% Loop over sessions..
for ses = 1:length(unique_sessions)
    for st = 1:length(session_types)
        if ~isnan(session_types(ses,st))
            blocks_per_run = nr_blocks_per_run(ses,st);
            
            
            % Get session information about blocks and trials.
            ses_idx           = (condition_master.session_nr==ses & condition_master.session_type==st);
            ses_blocks        = condition_master.global_block_nr(ses_idx);
            ses_trials        = condition_master.trial_nr(ses_idx);
            ses_trialtype     = condition_master.trial_type(ses_idx);
            ses_crossing_vec  = condition_master.crossing_nr(ses_idx);
            assert(isequal([1:length(unique(ses_blocks))]',sort(unique(ses_blocks))))
            
            % Get unique block nrs and associated trial types
            [unique_blocks, block_start_idx] = unique(ses_blocks);
            crossings_unique_blocks          = ses_crossing_vec(block_start_idx);
            unique_trialtypes                = ses_trialtype(block_start_idx);
            n_trialtype1 = sum(unique_trialtypes==1); % single
            n_trialtype2 = sum(unique_trialtypes==2); % double
            
            % Check if we have the expect nr of single and double stim
            % presentation blocks
            assert(isequal(n_trialtype1+n_trialtype2,length(unique_blocks)))
            
            % Get the nr of runs for this session (either 11 for behavior
            % or for MRI 10 runs for wide and all but last deep sessions.
            % Last deep session has 5 runs.
            runs_to_fill = runs_per_session(ses,st);
            
            % Separate single and double stim blocks:
            s_idx = (unique_trialtypes==1);
            d_idx = (unique_trialtypes==2);
            s_blocks = unique_blocks(s_idx);
            d_blocks = unique_blocks(d_idx);
            
            % Check if separation still results in total nr of blocks
            % across the session.
            assert(isequal(length(s_blocks)+length(d_blocks),length(unique_blocks)))
            
            % Shuffle block order within a session:
            % we only break this while loop when we are happy with the
            % order of blocks: that is, we want every 4 blocks to be from a
            % different crossing nr
            while 1
                reshuffle_me = false;
                attempts = 0;
                
                fprintf('[%s]: (Re-)Start shuffle attempts..', mfilename)
                while 1
                    attempts = attempts + 1;
                    fprintf('.');

                    % Create a indexing vector that randomizes the block order
                    % across all runs within the session.
                    bb_rnd = randperm(length(unique_blocks),length(unique_blocks));

                    % Apply the randomization vector to reorder blocks,
                    % crossings, trial types accordingly..
                    block_start_idx_shuffled = block_start_idx(bb_rnd);
                    crossings_shuffled       = crossings_unique_blocks(bb_rnd);
                    trialtypes_shuffled      = unique_trialtypes(bb_rnd);

                    % check if we have more than two crossings repeat back to back
                    % which may happen by chance. If so, we repeat shuffle of
                    % blocks.
                    repeat_blocks = sum(diff(crossings_shuffled)==0);
                    if repeat_blocks > max_allowable_block_repeats
                        reshuffle_me = true;
                    else
                        reshuffle_me = false;
                    end   
                    if reshuffle_me
                        continue;
                    end
                    % For every group of 7 single stimulus presentation blocks,
                    %  we check their crossing, if they have the same crossing
                    %  we shuffle the order of all single stimulus presentation
                    %  blocks. We use 7 because that is the max nr of blocks
                    %  we expect within a run.
                    block_cnt = 1:blocks_per_run:length(crossings_shuffled);
                    for nn = 1:length(block_cnt)
                        if (block_cnt(nn)+(blocks_per_run-1)) <= length(crossings_shuffled)
                            tmp_blocks = crossings_shuffled(block_cnt(nn):(block_cnt(nn)+(nr_blocks_per_run-1)));
                        else
                            tmp_blocks = crossings_shuffled(block_cnt(nn):end);
                        end
                        d_block_idx = ismember(tmp_blocks,d_blocks);
                        s_block_idx = ismember(tmp_blocks,s_blocks);
                        nr_s_blocks = sum(d_block_idx);
                        nr_d_blocks = sum(s_block_idx);
                        
                        if nr_d_blocks>=4 % if we have too many double-stim blocks, we restart
                            shuffle_ok(nn) = false;
                        end
%                         if nr_s_blocks==7 % if we have only single-stim blocks in a row, we restart
%                             shuffle_ok(nn) = false;
%                         end
                        if length(unique(tmp_blocks)) >= length(tmp_blocks)-1 % if we have 6 or 7 unique crossings in a row, we are happy
                            shuffle_ok(nn) = true;
                        else % if we don't have at least 6 unique crossings in a row, we break out of loop and restart shuffling
                            shuffle_ok(nn) = false;
                            continue; 
                        end
                        if nr_d_blocks > 0  % && nr_s_blocks > 0 we want at least one double stim presentation block per run
                            shuffle_ok(nn) = true;
                        else
                            shuffle_ok(nn) = false;
                            continue; % break out of brick laying loop and restart?
                        end

                    end

                    if (sum(shuffle_ok) == length(block_cnt)) && ~reshuffle_me
                        fprintf('Done!\n')
                        break;
                    end
                end % first smaller while loop to check randomization/block repeats
            
                fprintf('\n[%s]: Total of %d attempts\n', mfilename, attempts);
            
                % Add blocks to runs within this session
                run_matrix    = zeros(runs_to_fill,blocks_per_run); % max 7 possible slots
                run_crossings = run_matrix;

                % reset counters
                bb_glbl_cnt   = 0;
                bb_cnt        = 0;
                rr_cnt        = 1;
                run_dur       = run_dur_min + sum(max(IBIs)*blocks_per_run);
                total_run_dur = [];
                run_too_short = [];
                runs_ok       = true;
                
                while 1
                    % update counters
                    bb_glbl_cnt = bb_glbl_cnt+1;
                    bb_cnt      = bb_cnt+1;
                    if bb_glbl_cnt > length(block_start_idx_shuffled)
                        if sum(run_matrix(end,:)>0)<5
                            runs_ok = false;
                        end
                        break; % if we used all the blocks, we continue
                    end
                    
                    if bb_cnt > size(run_matrix,2)
                        total_run_dur(rr_cnt) = run_dur;
                        rr_cnt  = rr_cnt +1;   % count next run
                        bb_cnt  = 1;           % reset within run block counter when moving to the next run
                        run_dur = run_dur_min + sum(max(IBIs)*blocks_per_run); % reset run duration when moving to the next run
                    end
                    
                    run_dur_tmp = run_dur + block_dur(trialtypes_shuffled(bb_glbl_cnt));
                    
                    if run_dur_tmp-slack > run_dur_max
                        total_run_dur(rr_cnt) = run_dur;
                        
                        % if the run is 30 s shorter than expected, we try again.
                        run_too_short = (run_dur_max-total_run_dur) > max_run_deviation;
                        
                        % go to next run
                        rr_cnt  = rr_cnt +1;   % count next run
                        bb_cnt  = 1;           % reset within run block counter when moving to the next run
                        run_dur = run_dur_min + sum(IBIs); % reset run duration when moving to the next run
                    else
                        run_dur = run_dur_tmp;
                    end
                    
                    if ~isempty(run_too_short) && sum(run_too_short)>0
                        runs_ok = false;
                        break;
                    end
                    
                    % add block to run
                    run_crossings(rr_cnt, bb_cnt) = crossings_shuffled(bb_glbl_cnt);
                    run_matrix(rr_cnt, bb_cnt) = block_start_idx_shuffled(bb_glbl_cnt);
                end % second smaller while loop to check run duration

                if runs_ok
                    break;
                end
                
            end % big while loop
            % --- ONCE SHUFFLING AND ALLOCATING IS DONE --- 
            
            
            % ensure we added all the blocks
            assert(length(run_matrix(run_matrix>0))==length(block_start_idx_shuffled))
            assert(length(unique(run_crossings(run_crossings>0)))==length(unique(crossings_shuffled)))
            
            
            
            % Shuffle trials within each block
            block_cnt = 0;
            trial_order_local_shuffled  = cell(size(run_matrix));
            trial_order_global_shuffled = cell(size(run_matrix));
            fprintf('[%s]: \n',mfilename)
            for run_idx = 1:size(run_matrix,1)
                % Show user the block order within a run
                curr_block_starts = run_matrix(run_idx,:);
                curr_block_starts = curr_block_starts(curr_block_starts>0);
                fprintf('Run %02d, block start:',run_idx)
                disp(curr_block_starts)
                
                for bb_idx = 1:length(curr_block_starts)
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
                    fprintf('\nBlock start %02d, trial order:',curr_block_starts(bb_idx))
                    disp(trial_order_local_shuffled{run_idx,bb_idx})
                end
            end
            
            
            % For convenience, we make a copy of the condition master with shuffled blocks
            condition_master0 = condition_master;
            
            % reset the temporary run and block nrs
            condition_master0.run_nr   = NaN(size(condition_master0.run_nr));
            condition_master0.block_nr = NaN(size(condition_master0.block_nr));
            
            % Sort block order by ascending nr
            for rr = 1:size(run_matrix,1)
                
                % get blocks for a given run
                curr_block_start = run_matrix(rr,:);
                curr_block_end   = curr_block_start+cellfun(@length, trial_order_local_shuffled(rr,:));
                curr_block_start = curr_block_start(curr_block_start>0);
                curr_block_end   = curr_block_end(1:length(curr_block_start))-1;
                is_run = [];
                for dd = 1:length(curr_block_start)
                    is_run = cat(1,is_run, [curr_block_start(dd):curr_block_end(dd)]');
                end
                assert(all(ismember(run_crossings(rr,run_crossings(rr,:)>0),condition_master0.crossing_nr(is_run))'))
                condition_master0.run_nr(is_run)=rr;
                
                for cc = 1:length(curr_block_start)
                    is_block = curr_block_start(cc):curr_block_end(cc);
                    condition_master0.block_nr(is_block) = cc;
                    condition_master0.trial_nr(is_block) = trial_order_local_shuffled{rr,cc}';
                    
                    % sort trial order
                    [~,ti] = sort(condition_master0.trial_nr(is_block));
                    condition_master0.trial_nr(is_block) = condition_master0.trial_nr(is_block(ti));
                end
                
            end
            
            % sort blocks
            [~,bi] = sort(condition_master0.block_nr);
            
            % sort runs
            [~,ri] = sort(condition_master0.run_nr(bi));
            
            % sort both blocks and runs at once
            condition_master0 = condition_master0(bi(ri),:);
            
            % Update trial repeat nr
            condition_master0 = vcd_getTrialRepeatNr(condition_master0);
            
            % Get order of shuffled trials within a session
            [~,global_trial_row_idx] = intersect(condition_master0.global_trial_nr,condition_master.global_trial_nr);
            assert(isequal(condition_master0.global_trial_nr(global_trial_row_idx),condition_master.global_trial_nr));
            
            

            % concatenate condition_master, run_matrix, and global trial
            % indices
            condition_master_shuffled = cat(1,condition_master_shuffled,condition_master0);
            condition_master_shuffle_idx{ses,st} = global_trial_row_idx;
            session_block_matrix{ses,st}         = run_matrix;
            session_crossing_matrix{ses,st}      = run_crossings;
        end
    end
end

toc;

% Update global_block_nr 
condition_master_shuffled = vcd_updateGlobalCounters(params, condition_master_shuffled, session_env);

% Store randomization file and shuffled condition master
if params.store_params
    fprintf('\n[%s]:Storing shuffled condition master and randomization file..\n',mfilename)
    saveDir = fileparts(fullfile(params.stim.fix.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir,sprintf('%scondition_master_shuffled_%s_%s_%s.mat',prefix,params.disp.name, session_env, datestr(now,30))),'condition_master_shuffled')
    save(fullfile(saveDir,sprintf('%scondition_master_randomization_%s_%s_%s.mat',prefix,params.disp.name, session_env, datestr(now,30))),'condition_master_shuffle_idx','session_block_matrix','session_crossing_matrix','params')
end

