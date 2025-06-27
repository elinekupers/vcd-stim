function [condition_master_shuffled, fname, condition_master_shuffle_idx, ...
          session_crossing_matrix, session_block_matrix]  = ...
            vcd_randomizeBlocksAndRunsWithinSession(params, condition_master, env_type, varargin)
% VCD function to shuffle stimulus blocks and trials within a block across
% all sessions/session types. Once we find the order of stimulus blocks that 
% adheres to our contraints, we resave the condition table with the 
% individual subject nr in the name as "condition_master_shuffled", a table 
% that has sorted the new order of run nrs, block_nrs, and trial nrs for 
% all the sessions in the MRI or BEHAVIORAL experiment.
%
%  [condition_master_shuffled, run_matrix,trial_order_local_shuffled, trial_order_global_shuffled] ...
%    = vcd_randomizeBlocksAndRunsWithinSession(params, condition_master, env_type)
%
% This code will do the following:
% 1. Grab the crossing numbers for each stimulus block and separate them
%    into two vectors, according to trial type (1: single-stimulus 
%    presentation, 2: double-stimulus presentaton).
% 2. Shuffle the order of single- and double-stim blocks separately.
% 3. Check if there are any repeats for each shuffle. If there are more
%    than X blocks (defined by "max_block_repeats") with the same
%    crossing nr happen back to back (what we call a "block repeat"), which
%    may happen by chance, we repeat shuffle of blocks for both trial types.
%    Note that block repeats are not necessarily an issue, as the code will
%    try to mix single- and double-stim blocks when allocating blocks to 
%    runs. In other words, if we have a repeat of a double-stim block 
%    (11-11), they may actually get separated by a single-stim block during 
%    a run (11-31-11).
% 4. Add blocks to runs within this session:
%    * We set how many possible block slots each run can have 
%    (“blocks_per_run”) and how many runs a session can have 
%    (“runs_to_fill”).
%    * Randomly pick the number [1] or [2], which determines which list we 
%    will use first to add block to a run (1 = try single-stim block first,
%    2 = try double-stim block first).
%    * Check how many blocks do we have in the run right now? If we already
%    filled up the run with max blocks, we move the potential blocks to the
%    next run. If we have block slots left within a run AND adding this
%    block to the run will not exceed the run's max duration, then we will
%    the first potential block option (e.g., [2]). If we have block slots 
%    left, BUT the potential block option makes the run exceed its max 
%    duration, then we will try adding the other potential block option (in
%    this example, [1]). If both potential block options make the run too
%    long, we will move to the next run and try adding block option 1.
%    * Once a block is added, we remove it from its block list.
%    * Once a list of single-/double-stim blocks is empty, just use the  
%    other block list to fill the runs, until we run out.
%    Note: Users can constrain the combinations of single-/double-stim   
%    blockswithin a run, using the parameter "allowed_block_combinations". 
%    For example, if you only want runs with either 7 single-stim and no
%    double-stim blocks, or 5 double-stim and no single-stim blocks you can
%    set "allowed_block_combinations" to [7 0; 0 5].
% 5. Once the blocks are allocated to runs, we will do the following "HARD"
%    constraints checks: If the allocation of blocks does not adhere to 
%    constraints, we start over and reshuffle all blocks within a session.
%  * Run duration. The total sum of all stimulus blocks + overhead in a run 
%    (eye tracking block, pre/post rest periods, minimum duration of IBIs)
%    cannot exceed the number of time frames defined by:
%       params.exp.run.total_run_dur_[MRI/BEHAVIOR].
%    Runs cannot be too short in their duration either, here defined as the 
%    max duration - "max_run_deviation". If runs are shorter than expected, 
%    we check if adding max allowable IBIs would help reach the threshold 
%    set by "max_run_deviation".
%  * Nr of blocks per run. The number of blocks within a run cannot exceed 
%    the number defined by "blocks_per_run", which is defined by:
%       params.exp.session.mri.[wide/deep/behavior].nr_blocks_per_run.
%    and the number cannot be less or equal to ceil(blocks_per_run/2).
%  % Nr of runs per session. The number of runs within a sessions cannot  
%    exceed the number defined by: 
%       params.exp.session.mri.[wide/deep/behavior].n_runs_per_session.
%  % Block allocation. We must assign all the blocks to the session.
% If we shuffled all blocks more than 10,000 times, code will throw an error.
%
% Additional (tweakable) constraints:
% User can set some additional constraints when to stop inserting blocks 
% into a run and/or reshuffle the block order:
% * "slack" allows user to extend the max run duration threshold before
% moving on to the next run. One may want to do this to allow for adding
% one more block to a run (because the variable IBI could technically allow
% for that) and avoid ending up with long blank at the end.
% * "max_run_deviation" is another way to avoid ending up with really long
% blanks at the end of the run as it defines how much shorter a run can be
% before we will reshuffle the block order within a run again.
% * "max_block_repeats" specifies how many blocks do we allow to
% repeat back to back in the long sequence of shuffled session blocks.
% * "allowed_block_combinations" is for users that know the specific nr of
% single- and double-stim presentation blocks they want in a given run.
% This constraint can reduce time searching for the optimal solution.
%
% For reference, the optimal combinations of single- and double-stimulus 
% presentation blocks per 376 s run (6'16" or 235 TRs) are :
%
% single    double   min min   max min   mean min  unaccount time (s) from max run time
% 7.0000         0    6.2667    6.6667    6.4667  -12.0000
% 4.0000    2.0000    5.9750    6.3083    6.1417    7.5000
% 3.0000    3.0000    6.2417    6.5750    6.4083   -8.5000
%      0    5.0000    5.9500    6.2167    6.0833    9.0000
%
% IBIs: For VCD, the possible IBI ranges are between 5-9 seconds, which
% sums to a total IBI duration of [min-max] 25-54 seconds or [5*300=2400 to
% 6*540=3240] time frames.
%
% INPUTS:
% * params  
% * condition_master            : (table) master table with all the stimuli 
%                                 for each trial, in each session/run/block
% * env_type                 : (char) session environment for this session: 
%                                 choose from 'BEHAVIORAL' or 'MRI'.
% * [subj_id]                    : (char) fixed name added to beginning of
%                                 filename that stores the 
%                                 condition_master_shuffled file. Default
%                                 is 'vcd_subj000'
% * [slack]                     : (int) time frames per run to extend the 
%                                 run duration threshold before moving on 
%                                 to the next run to add new blocks in a 
%                                 particular session. Default: 0 seconds 
%                                 (assuming 60 Hz presentation timing).
% * [max_run_deviation]         : (int) time frames per run to check if a 
%                                 run is X seconds too short or too long. 
%                                 Default: 9 seconds (assuming 60 Hz 
%                                 presentation timing).
% * [max_block_repeats]         : (int) nr of blocks we allow for repeat in
%                                 their crossing (so finding A-A and later 
%                                 on B-B). Default: 10.
% * [allowed_block_combinations]: (int) vector of possible [single;double]
%                                 -stimulus blocks we allow within a run.
%                                 If left empty ([]), code will not look
%                                 for specific combination of 
%                                 single-/double-stim blocks in a run and
%                                 assigned blocks randomly with the given
%                                 constraints. Default: [].
%                                 Example: [6 4 2 1; 1 2 4 5]. 
% * [saveDir]                   : (char) folder where subject's shuffled
%                                  condition_master and randomization files will be stored.
%
% OUTPUTS:
% * condition_master_shuffled    : (table) new version of condition_master
%                                  where the stimulus blocks are shuffled 
%                                  across runs within a session and trials 
%                                  are shuffled within a stimulus blocks.
% * fname                        : (char) filename of stored shuffled
%                                   condition master.
% * condition_master_shuffle_idx : (cell) indexing order of all trials
%                                   within a single session (dim 1) and
%                                   session_type (dim 2) of the original to 
%                                   condition_master to get the trial order 
%                                   in condition_master_shuffled.
% * session_crossing_matrix      : (cell) indexing order of the block's 
%                                   crossing nrs within a single shuffled  
%                                   session (dim 1) and session_type (dim
%                                   2). Each cell is matrix with dimensions
%                                   [runs,blocks].
% * session_block_matrix         : (cell) indexing order of the block's 
%                                   crossing nrs within a single shuffled  
%                                   session (dim 1) and session_type (dim
%                                   2). Each cell is matrix with dimensions
%                                   [runs,blocks].
%
% Written by E Kupers @ UMN 2025/06

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'                         , @isstruct);
p0.addRequired('condition_master'               , @istable);
p0.addRequired('env_type'                       , @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));
p0.addParameter('subj_id'                       , 'vcd_subj000', @ischar);
p0.addParameter('slack'                         , 0*60 , @isnumeric); % time frames
p0.addParameter('max_run_deviation'             , 15*60, @isnumeric); % time frames
p0.addParameter('max_block_repeats'             , 10   , @isnumeric);
p0.addParameter('allowed_block_combinations'    , []   , @isnumeric); % first row is max single-stim block nrs / run, second row is max double-stim block nrs / run (e.g., [6 4 2 3; 1 2 4 3]) 
p0.addParameter('saveDir'                       , ''   , @ischar); % where to store the shuffled condition_master
p0.addParameter('load_params'                   , false, @islogical);
p0.addParameter('store_params'                  , true, @islogical);
p0.addParameter('store_imgs'                    , false, @islogical);
p0.addParameter('verbose'                       , false, @islogical);

% Parse inputs
p0.parse(params,condition_master,env_type,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

%% Infer other inputs
unique_sessions  = unique(condition_master.session_nr);

% Get session parameters depending on whether this is the Behavioral
% experiment or MRI experiment.
[~,session_types,runs_per_session,run_dur_min, ...
 run_dur_max,~, IBIs, ~, ~, ~, ~, nr_blocks_per_run,~,all_block_dur] = ...
    vcd_getSessionEnvironmentParams(params, env_type);

% Get additional IBI duration (in case we deal with a short run)
additional_IBIs = max(IBIs)-min(IBIs);

% Preallocate space
condition_master_shuffled    = [];
condition_master_shuffle_idx = cell(length(unique_sessions),find(sum(~isnan(session_types))));
session_crossing_matrix      = cell(length(unique_sessions),find(sum(~isnan(session_types))));
session_block_matrix         = cell(length(unique_sessions),find(sum(~isnan(session_types))));

tic;

%% Loop over sessions..
for ses = 1:length(unique_sessions)
    for st = 1:size(session_types,2)
        if ~isnan(session_types(ses,st))
            blocks_per_run = nr_blocks_per_run(ses,st);
            
            % Get session information about blocks and trials.
            ses_idx           = (condition_master.session_nr==ses & condition_master.session_type==st);
            ses_blocks        = condition_master.global_block_nr(ses_idx);
            ses_trials        = condition_master.trial_nr(ses_idx);
            ses_trialtype     = condition_master.trial_type(ses_idx);
            ses_crossing_vec  = condition_master.crossing_nr(ses_idx);
            assert(isequal(length(ses_trials),length(ses_blocks)))
            if params.is_demo
                assert(isequal(1+sum(abs(diff(ses_crossing_vec))>0),numel(unique(ses_blocks))))
            end
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
            all_blocks{1} = unique_blocks(s_idx);
            all_blocks{2} = unique_blocks(d_idx);
            
            % Check if separation still results in total nr of blocks
            % across the session.
            assert(isequal(length(all_blocks{1})+length(all_blocks{2}),length(unique_blocks)))
            
            % Shuffle block order within a session:
            % we only break this while loop when we are happy with the
            % order of blocks: that is, we want every 4 blocks to be from a
            % different crossing nr
            fprintf('[%s]: Start shuffle attempts..', mfilename)
            attempts = 0;
            while 1
                reshuffle_me = false(1,2);

                while 1
                    attempts = attempts + 1;
                    
                    if mod(attempts,100)==0
                        fprintf('.');
                    end
                    
                    block_start_idx_shuffled = cell(1,2);
                    crossings_shuffled       = cell(1,2);
                    trialtypes_shuffled      = cell(1,2);
                    bb_rnd                   = cell(1,2);
                    for ii = [1,2] % single, double
                        % Create a indexing vector that randomizes the block order
                        % of single-stim and double-stim blocks separately.
                        bb_rnd{ii} = randperm(length(all_blocks{ii}),length(all_blocks{ii}));
                        
                        % Apply the randomization vector to reorder blocks,
                        % crossings, trial types accordingly..
                        block_start_idx_shuffled{ii} = block_start_idx(all_blocks{ii}(bb_rnd{ii}));
                        crossings_shuffled{ii}       = crossings_unique_blocks(all_blocks{ii}(bb_rnd{ii}));
                        trialtypes_shuffled{ii}      = unique_trialtypes(all_blocks{ii}(bb_rnd{ii}));
                        
                        % check if we have more than X crossings repeat back to back
                        % which may happen by chance. If so, we repeat shuffle of
                        % blocks.
                        repeat_blocks = sum(diff(block_start_idx_shuffled{ii})==0);
                        if repeat_blocks > max_block_repeats
                            reshuffle_me(ii) = true;
                        else
                            reshuffle_me(ii) = false;
                        end
                    end

                    if ~reshuffle_me
                        break;
                    end
                end % first smaller while loop to check randomization/block repeats
            
                % Add blocks to runs within this session
                run_matrix                = zeros(runs_to_fill,blocks_per_run); % max 7 possible slots
                run_crossings             = run_matrix;
                run_trial_types           = run_matrix;
                run_dur                   = [];
                trialtypes_shuffled0      = trialtypes_shuffled;
                block_start_idx_shuffled0 = block_start_idx_shuffled;
                crossings_shuffled0       = crossings_shuffled;
                
                run_types         = cell(size(allowed_block_combinations));
                run_types_counter = zeros(size(allowed_block_combinations));
                check_blocks      = true; % we don't check blocks if we divide the final ones across the left over runs
                curr_run_dur      = 0; % duration tracker
                rr_cnt            = 1; % run counter
                
                while 1
                    % randomly pick a starting point that determines if we 
                    % try to add a single-stim block first or a double-stim block first?
                    bbi = randi([1,2],[1,1]);
                    
                    % how many blocks do we have in the run right now?
                    curr_nr_blocks = sum(run_matrix(rr_cnt,:)>0);
                    
                    % if double-stim list is empty, use single-stim list
                    if ~isempty(trialtypes_shuffled0{1}) && isempty(trialtypes_shuffled0{2})
                        
                        total_nr_blocks_left = length(trialtypes_shuffled0{1});
                        nr_incomplete_runs   = runs_to_fill-rr_cnt;
                        
                        potential_block_dur(1)      = all_block_dur(trialtypes_shuffled0{1}(1));
                        potential_block_crossing(1) = crossings_shuffled0{1}(1);
                        potential_block_idx(1)      = block_start_idx_shuffled0{1}(1);
                        potential_block_dur(2)      = [NaN];
                        potential_block_crossing(2) = [NaN];
                        potential_block_idx(2)      = [NaN];
                        
                        % check if accidentally picked the empty side,
                        % switch block index to the non empty side
                        if isnan(potential_block_idx(bbi))
                            bbi = setdiff([1,2],bbi);
                            assert(~isnan(potential_block_idx(bbi)))
                        end
                        
                        if (ibi + curr_run_dur + run_dur_min + max_run_deviation) < run_dur_max
                            % if this run needs another block, then let's continue as usual
                            check_blocks = true;
                            
                        elseif ceil(total_nr_blocks_left/nr_incomplete_runs) <= blocks_per_run
                            % record last run duration
                            run_dur(rr_cnt) = curr_run_dur + run_dur_min;
                            curr_run_dur    = 0;
                            
                            % split single-stim blocks evenly across the last runs
                            nr_blocks_incomplete_runs = repmat(floor(total_nr_blocks_left/nr_incomplete_runs),1,nr_incomplete_runs);
                            
                            remainder_blocks = rem(total_nr_blocks_left,nr_incomplete_runs);
                            if ~isempty(remainder_blocks) && remainder_blocks~=0
                                tmp_idx = find(nr_blocks_incomplete_runs<blocks_per_run);
                                for tt = 1:length(tmp_idx)
                                    if remainder_blocks~=0
                                        nr_blocks_incomplete_runs(tmp_idx(tt)) = nr_blocks_incomplete_runs(tmp_idx(tt)) + 1;
                                        remainder_blocks = remainder_blocks - 1;
                                    end
                                end
                            end
                            % make sure we use the right ibi duration
                            ibi = min(IBIs);

                            for rmb = 1:nr_incomplete_runs
                                rr_cnt = rr_cnt+1;
                                final_trialtypes      = trialtypes_shuffled0{1}(1:nr_blocks_incomplete_runs(rmb))';
                                final_block_start_idx = block_start_idx_shuffled0{1}(1:nr_blocks_incomplete_runs(rmb))';
                                final_crossings       = crossings_shuffled0{1}(1:nr_blocks_incomplete_runs(rmb))';
                                
                                run_matrix(rr_cnt, 1:nr_blocks_incomplete_runs(rmb)) = final_block_start_idx;
                                run_crossings(rr_cnt, 1:nr_blocks_incomplete_runs(rmb)) = final_crossings;
                                run_trial_types(rr_cnt, 1:nr_blocks_incomplete_runs(rmb)) = final_trialtypes;
                                curr_run_dur = (ibi*(nr_blocks_incomplete_runs(rmb)-1)) + sum(all_block_dur(final_trialtypes));
                            
                                % log run duration and reset curr_run_dur
                                run_dur(rr_cnt) = curr_run_dur + run_dur_min;
                                curr_run_dur = 0;
                                
                                % remove block from lists
                                crossings_shuffled0{1}(1:nr_blocks_incomplete_runs(rmb)) = [];
                                block_start_idx_shuffled0{1}(1:nr_blocks_incomplete_runs(rmb)) = [];
                                trialtypes_shuffled0{1}(1:nr_blocks_incomplete_runs(rmb)) = [];
                                
                            end
                            
                            check_blocks = false;
                        end
                       % if no more single-stim blocks, then only use double-stim blocks 
                    elseif isempty(trialtypes_shuffled0{1}) && ~isempty(trialtypes_shuffled0{2})
                           potential_block_dur(1)      = [NaN];
                           potential_block_crossing(1) = [NaN];
                           potential_block_idx(1)      = [NaN];
                           potential_block_dur(2)      = all_block_dur(trialtypes_shuffled0{2}(1));
                           potential_block_crossing(2) = crossings_shuffled0{2}(1);
                           potential_block_idx(2)      = block_start_idx_shuffled0{2}(1);
                           
                           total_nr_blocks_left = length(trialtypes_shuffled0{2});
                           nr_incomplete_runs   = runs_to_fill-rr_cnt;
                           
                           if (ibi + curr_run_dur + run_dur_min + max_run_deviation) < run_dur_max
                               % if this run needs that block, then let's continue as usual
                               check_blocks = true;
                               
                           elseif ceil(total_nr_blocks_left/nr_incomplete_runs) <= blocks_per_run
                               
                               % record last run duration
                               run_dur(rr_cnt) = curr_run_dur + run_dur_min;
                               curr_run_dur    = 0;
                               
                               % split double-stim blocks evenly across the last runs
                               nr_blocks_incomplete_runs = repmat(floor(total_nr_blocks_left/nr_incomplete_runs),1,nr_incomplete_runs);
                               remainder_blocks = rem(total_nr_blocks_left,nr_incomplete_runs);
                               if ~isempty(remainder_blocks) && remainder_blocks~=0
                                   tmp_idx = find(nr_blocks_incomplete_runs<blocks_per_run);
                                   for tt = 1:length(tmp_idx)
                                       if remainder_blocks~=0
                                           nr_blocks_incomplete_runs(tmp_idx(tt)) = nr_blocks_incomplete_runs(tmp_idx(tt)) + 1;
                                           remainder_blocks = remainder_blocks - 1;
                                       end
                                   end
                               end
                               
                               % make sure we use the right ibi duration
                               ibi = min(IBIs);
                               
                               for rmb = 1:nr_incomplete_runs
                                    rr_cnt = rr_cnt+1;
                                    final_trialtypes      = trialtypes_shuffled0{2}(1:nr_blocks_incomplete_runs(rmb))';
                                    final_block_start_idx = block_start_idx_shuffled0{2}(1:nr_blocks_incomplete_runs(rmb))';
                                    final_crossings       = crossings_shuffled0{2}(1:nr_blocks_incomplete_runs(rmb))';
                                   
                                    run_matrix(rr_cnt, 1:nr_blocks_incomplete_runs(rmb)) = final_block_start_idx;
                                    run_crossings(rr_cnt, 1:nr_blocks_incomplete_runs(rmb)) = final_crossings;
                                    run_trial_types(rr_cnt, 1:nr_blocks_incomplete_runs(rmb)) = final_trialtypes;
                                    curr_run_dur = (ibi*(nr_blocks_incomplete_runs(rmb)-1)) + sum(all_block_dur(final_trialtypes));
                                   
                                    % log run duration and reset curr_run_dur
                                    run_dur(rr_cnt) = curr_run_dur + run_dur_min;
                                    curr_run_dur = 0;
                                    
                                    % remove block from lists
                                    crossings_shuffled0{2}(1:nr_blocks_incomplete_runs(rmb)) = [];
                                    block_start_idx_shuffled0{2}(1:nr_blocks_incomplete_runs(rmb)) = [];
                                    trialtypes_shuffled0{2}(1:nr_blocks_incomplete_runs(rmb)) = [];
                               end
                               
                               check_blocks = false;
                           end
                    else
                        % if both single and double blocks exists, we
                        % consider them both as potential blocks (1:
                        % single; 2: double)
                        potential_block_dur(1)      = all_block_dur(trialtypes_shuffled0{1}(1));
                        potential_block_dur(2)      = all_block_dur(trialtypes_shuffled0{2}(1));
                        potential_block_crossing(1) = crossings_shuffled0{1}(1);
                        potential_block_crossing(2) = crossings_shuffled0{2}(1);
                        potential_block_idx(1)      = block_start_idx_shuffled0{1}(1);
                        potential_block_idx(2)      = block_start_idx_shuffled0{2}(1);
                        
                        check_blocks = true;
                    end
                    % if we are ready with potential block options
                    if check_blocks
                    
                        % see if we can add another block..
                        if (curr_nr_blocks < blocks_per_run)

                            if curr_nr_blocks == 0
                                ibi = 0;
                            else
                                ibi = min(IBIs);
                            end

                            % if we have specified allowed combinations, use those
                            if ~isempty(allowed_block_combinations)
                                if potential_block_dur(bbi) == all_block_dur(1)

                                    for rt = 1:size(allowed_block_combinations,2)
                                        if run_types_counter(1,rt) < allowed_block_combinations(1,rt)
                                            run_types{1,rt} = potential_block_idx(bbi);
                                            run_types_counter(1,rt)= run_types_counter(1,rt)+1;
                                            break;
                                        end
                                    end

                                elseif potential_block_dur(bbi) == all_block_dur(2)

                                    for rt = 1:size(allowed_block_combinations,2)
                                        if run_types_counter(2,rt) < allowed_block_combinations(2,rt)
                                            run_types{2,rt} = potential_block_idx(bbi);
                                            run_types_counter(2,rt)= run_types_counter(2,rt)+1;
                                            break;
                                        end
                                    end
                                end
                                
                                curr_run_dur = curr_run_dur + ibi + potential_block_dur(bbi);

                            else % we just add to the available run
                                % if we have enough time left in the run, add the first potential block
                                if (ibi + curr_run_dur + potential_block_dur(bbi)) < (run_dur_max-run_dur_min)

                                    run_matrix(rr_cnt, curr_nr_blocks+1)     = potential_block_idx(bbi);
                                    run_crossings(rr_cnt, curr_nr_blocks+1)  = potential_block_crossing(bbi);
                                    run_trial_types(rr_cnt, curr_nr_blocks+1) = bbi;
                                    curr_run_dur = curr_run_dur + ibi + potential_block_dur(bbi);

                                    % remove block from lists
                                    crossings_shuffled0{bbi}(1) = [];
                                    block_start_idx_shuffled0{bbi}(1) = [];
                                    trialtypes_shuffled0{bbi}(1) = [];

                                % if not, try adding the other stim block..
                                elseif (curr_run_dur + ibi + potential_block_dur(setdiff([1,2],bbi))) < (run_dur_max-run_dur_min)

                                    run_matrix(rr_cnt,curr_nr_blocks+1)       = potential_block_idx(setdiff([1,2],bbi));
                                    run_crossings(rr_cnt, curr_nr_blocks+1)   = potential_block_crossing(setdiff([1,2],bbi));
                                    run_trial_types(rr_cnt, curr_nr_blocks+1) = setdiff([1,2],bbi);
                                    curr_run_dur = curr_run_dur + ibi + potential_block_dur(setdiff([1,2],bbi));

                                    % remove block from list
                                    crossings_shuffled0{setdiff([1,2],bbi)}(1) = [];
                                    block_start_idx_shuffled0{setdiff([1,2],bbi)}(1) = [];
                                    trialtypes_shuffled0{setdiff([1,2],bbi)}(1) = [];
                                    
                                else % if both blocks make the run too long, we move potential blocks to the next run
                                    run_dur(rr_cnt) = curr_run_dur + run_dur_min;
                                    rr_cnt = rr_cnt+1;
                                    curr_run_dur = 0;
                                end
                            end

                        else
                            % if we already filled up the run with max blocks, we move the potential blocks to the next run
                            run_dur(rr_cnt) = curr_run_dur + run_dur_min;
                            rr_cnt = rr_cnt+1;
                            curr_run_dur = 0;
                        end
                    end
                    
                    if rr_cnt > runs_to_fill
                        % check if we added last run duration
                        if length(run_dur) < sum(run_matrix(:,1)>0)
                            run_dur(rr_cnt) = curr_run_dur + run_dur_min;
                        end
                        break;
                    end
                    
                    if isempty(trialtypes_shuffled0{1}) && isempty(trialtypes_shuffled0{2})
                        % check if we added last run duration
                        if length(run_dur) < sum(run_matrix(:,1)>0)
                            run_dur(rr_cnt) = curr_run_dur + run_dur_min;
                        end
                        break;
                    end
                end % second smaller while loop to add blocks to runs
                

                % --- Do some checks --- 
                runs_ok = true;
                clear tmp_blocks_per_run;
                for xx = 1:size(run_matrix,1); tmp_blocks_per_run(xx) = sum(run_matrix(xx,:)>0); end
                        
                % if runs are longer than expected, then start over
                if sum(run_dur > run_dur_max) > 0
                    runs_ok = false;
                end
                
                % if we have more runs than expected, then start over
                if runs_ok
                    if size(run_matrix,1) > runs_to_fill 
                        runs_ok = false;
                    end
                end
                
                % if runs are shorter than expected, then check if additional IBIs would help
                if runs_ok
                    short_runs = find( (run_dur + max_run_deviation) < run_dur_max);
                    
                    % for debug purposes, print run duration...
                    %fprintf('[%s]: Run(s) %s are on the short side..\n',mfilename, num2str(short_runs));
                    %fprintf('[%s]: Run duration(s): %s time frames (%s seconds)\n', mfilename, num2str(run_dur(short_runs)), num2str(run_dur(short_runs)./params.stim.presentationrate_hz));
                    
                    if ~isempty(short_runs)
                        short_runs_more_IBIs = run_dur(short_runs) + [(tmp_blocks_per_run(short_runs)-1)*additional_IBIs]; % we add the difference between min and max IBI to simulate the run dur when we would use all max IBIs
                        too_short = find( (short_runs_more_IBIs + max_run_deviation) < run_dur_max );
                        if length(too_short)>0
                            % for debug purposes print outcome of adding extra IBI 
                            % %fprintf('[%s]: Even with extra-long IBIs, run(s) %s are still considered too short: %s time frames\n', mfilename, num2str(short_runs(too_short)), num2str(short_runs_more_IBIs(too_short)));
                            runs_ok = false; % if not, we start over
                        end
                    end
                end
                
                % if have more than 7 blocks per run, then start over
                if runs_ok
                    if size(run_matrix,2) > blocks_per_run 
                        runs_ok = false;
                    end
                end
                
                % if runs have 4 blocks or less, then start over
                if runs_ok
                    if any(tmp_blocks_per_run <= ceil(blocks_per_run/2))
                        runs_ok = false;
                    end
                end
                
                % if we haven't used all the blocks, then start over
                if runs_ok
                    if ~isempty(trialtypes_shuffled0{1}) || ~isempty(trialtypes_shuffled0{2})
                        runs_ok = false;
                    end
                end
                
                % if this is a demo, and there are not 7 blocks in one run, then start over
                if runs_ok && params.is_demo
                    if tmp_blocks_per_run ~= params.exp.session.demo.nr_blocks_per_run(ses,st)
                        runs_ok = false;
                    end
                end

                % Did we pass all the checks???
                if runs_ok
                    break; % if we passed all checks, we break the while loop and move on
                end
                
                % If we shuffled all blocks more than 10,000 times, we will throw an error.
                if attempts>10000
                    error('[%s]: Can''t seem to find a solution after 10000 attempts! Try running the same code again!\n',mfilename)
                end
                
            end % big while loop
            
            fprintf('\n[%s]: Total of %d attempts. \n', mfilename, attempts);
            % --- ONCE SHUFFLING AND ALLOCATING ARE FINISHED --- 
            
            % ensure we added all the blocks
            assert(isequal(length(run_matrix(run_matrix>0)),(length(block_start_idx_shuffled{1})+length(block_start_idx_shuffled{2}))))
            assert(isequal(length(unique(run_crossings(run_crossings>0))),length(unique(cat(1,crossings_shuffled{1},crossings_shuffled{2})))))
            sorted_crossing_start_tmp = sort(cat(1,crossings_shuffled{1},crossings_shuffled{2}));
            sorted_block_start_tmp    = sort(cat(1,block_start_idx_shuffled{1},block_start_idx_shuffled{2}));
            sorted_trial_types_tmp    = sort(cat(1,trialtypes_shuffled{1},trialtypes_shuffled{2}));
            if size(sorted_crossing_start_tmp,1) < size(sorted_crossing_start_tmp,2)
                sorted_crossing_start_tmp = sorted_crossing_start_tmp;
            end
            if size(sorted_block_start_tmp,1) < size(sorted_block_start_tmp,2)
                sorted_block_start_tmp = sorted_block_start_tmp';
            end
            if size(sorted_trial_types_tmp,1) < size(sorted_trial_types_tmp,2)
                sorted_trial_types_tmp = sorted_trial_types_tmp';
            end
            assert(isequal(sort(reshape(run_crossings(run_crossings>0),[],1)),sorted_crossing_start_tmp))
            assert(isequal(sort(reshape(run_matrix(run_matrix>0),[],1)),sorted_block_start_tmp));
            assert(isequal(sort(reshape(run_trial_types(run_trial_types>0),[],1)),sorted_trial_types_tmp))
            
            % Shuffle trials within each block
            glbl_block_cnt = 0;
            trial_order_local_shuffled  = cell(size(run_matrix));
            trial_order_global_shuffled = cell(size(run_matrix));
            
            % First, we make a copy of the condition master with shuffled blocks
            condition_master0 = condition_master(condition_master.session_nr==ses & condition_master.session_type==st,:);
            
            % reset the temporary run and block nrs
            condition_master0.run_nr   = NaN(size(condition_master0.run_nr));
            condition_master0.block_nr = NaN(size(condition_master0.block_nr));
            
            % shuffle the order of runs, to avoid that the last runs
            % of a session have one stimulus type.
            run_order = shuffle_concat(1:size(run_matrix,1),1);
            
            % Now go through each run and shuffle the trial order for each
            % allocated block, as well as adding the new block/run order to
            % a copy of the condition_master
            for run_idx = 1:length(run_order)
                
                % get blocks for a given run
                curr_block_start = run_matrix(run_idx,:);
                curr_block_start = curr_block_start(curr_block_start>0);
 
                for bb_idx = 1:length(curr_block_start)
                    glbl_block_cnt = glbl_block_cnt +1;
                    % What type of block are we dealing with
                    if run_trial_types(run_idx,bb_idx)==1
                        if params.is_demo
                            nr_trials = params.exp.block.demo.n_trials_single_epoch;
                        else
                            nr_trials = params.exp.block.n_trials_single_epoch;
                        end
                    elseif run_trial_types(run_idx,bb_idx)==2
                        if params.is_demo
                            nr_trials = params.exp.block.demo.n_trials_double_epoch;
                        else
                            nr_trials = params.exp.block.n_trials_double_epoch;
                        end
                    end
                    
                    trial_order_not_ok = true;
                    while trial_order_not_ok
                        % Shuffle local trial order (e.g.: [1,2,3,4] --> [2,4,1,3])
                        local_trial_order = shuffle_concat(1:nr_trials,1);
                        % get trial order within the run
                        trial_order       = curr_block_start(bb_idx):(curr_block_start(bb_idx)+nr_trials-1);
                        trial_order       = trial_order(local_trial_order);

                        % check order of image nrs, if we have any repeats, we shuffle
                        % the order of stimuli in a block
                        if all(diff(condition_master0.stim_nr_left(trial_order))~=0) || ...
                                all(diff(condition_master0.stim_nr_right(trial_order))~=0)
                            trial_order_not_ok = false;
                        end
                    end
                    
                    trial_order_global_shuffled{run_idx,bb_idx} = trial_order;
                    trial_order_local_shuffled{run_idx,bb_idx}  = local_trial_order;

                    % For debug purposes, you can print trial order within a block
                    % %fprintf('\nRun %02d, block row start %02d in condition_master, trial order: %s',run_idx,curr_block_start(bb_idx),num2str(local_trial_order))
                    
                    % Add the new run/block/trial order to condition_master
                    condition_master0.run_nr(trial_order)   = run_order(run_idx);
                    condition_master0.block_nr(trial_order) = bb_idx;
                    condition_master0.trial_nr(trial_order) = local_trial_order;
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
            [~,global_trial_row_idx] = intersect(condition_master0.global_trial_nr,condition_master.global_trial_nr(condition_master.session_nr==ses & condition_master.session_type==st));
            assert(isequal(condition_master0.global_trial_nr(global_trial_row_idx),condition_master.global_trial_nr(condition_master.session_nr==ses & condition_master.session_type==st)));

            % concatenate condition_master, run_matrix, and global trial
            % indices
            condition_master_shuffled            = cat(1, condition_master_shuffled,condition_master0);
            condition_master_shuffle_idx{ses,st} = global_trial_row_idx;
            session_block_matrix{ses,st}         = run_matrix;
            session_crossing_matrix{ses,st}      = run_crossings;
        end
    end
end

fprintf('\n')
toc;

% Update global_block_nr 
condition_master_shuffled = vcd_updateGlobalCounters(params, condition_master_shuffled, env_type);

% Store randomization file and shuffled condition master
if store_params
    randomization_params = struct();
    randomization_params.slack = slack;
    randomization_params.max_run_deviation = max_run_deviation;
    randomization_params.max_block_repeats = max_block_repeats;
    randomization_params.allowed_block_combinations = allowed_block_combinations;
    
    if isempty(saveDir)
        saveDir = fullfile(vcd_rootPath,'data',env_type,subj_id);
    end
    if ~exist(saveDir,'dir'), mkdir(saveDir); end    
    
    fname = sprintf('%scondition_master_%s%s_%s.mat',[subj_id '_'],choose(params.is_demo,'demo_',''),params.disp.name, datestr(now,30));
    fprintf('\n[%s]: Storing shuffled condition master and randomization file here:\n',mfilename)
    fprintf('\t%s', fullfile(saveDir,fname))
    save(fullfile(saveDir,fname),'condition_master_shuffled','condition_master_shuffle_idx','session_block_matrix','session_crossing_matrix','randomization_params');
end


