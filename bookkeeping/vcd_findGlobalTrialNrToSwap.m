function [A1, A2] = vcd_findGlobalTrialNrToSwap(condition_master,condition_master0,trial_idx)
% VCD function to navigate stimulus blocks where trials have the repeated stimulus numbers.
% This typically happens for WM blocks with 4 trials, for OBJ/Gabor/RDKs.
% To fix this, you can run this function to find a trial in a nearby stimulus block with
% the same crossing number, cued side and response number. After that, you
% can decide if you want to swap the trials as shown in the example below.
% 
% Note 1: This function is meant to be used in
%           vcd_randomizeBlocksAndRunsWithinSession.m, and not written to
%           be stand-alone function.
% Note 2: This function does not take into account session type and will
% search for similar trials across both session types (A/B). 
%
% INPUTS: 
%   condition_master    : (table) condition_master with all session and all trials
%   condition_master0   : (table) subset of condition_master with all trials indexed for one session nr and session type.
%   trial_idx           : (int) integers indexing the trials in one block of the condition_master0 that have multiple stimulus repeats.
%
% OUTPUTS:
%   A1                  : (int) global trial number for the new matching trial we want to swap bad trial with
%   A2                  : (int) global trial number for the bad trial (with repeat) 
%
% % Example how to use this function:
% [A1, A2] = vcd_findGlobalTrialNrToSwap(condition_master,condition_master0,trial_idx);
% A1_idx = condition_master.global_trial_nr==A1; % get condition master table index
% A2_idx = condition_master.global_trial_nr==A2; % get condition master table index
%
% %% Define trials to swap
% NOTE: We start from column 9 and leave column 1:8 with session_nr/session_type/block_nr/trial_nr/etc. 
% We only want to swap the stimulus-related trial info.
% tmp1 = condition_master(A1_idx,9:end); 
% tmp2 = condition_master(A2_idx,9:end);
%
% %% Display trials to swap, and the trials in the blocks
% tmp1
% tmp2
% condition_master(condition_master.global_block_nr==condition_master.global_block_nr(A1_idx),:)
% condition_master(condition_master.global_block_nr==condition_master.global_block_nr(A2_idx),:)
%
% %% Resave large condition_master
% condition_master(A2_idx,9:end) = tmp1; condition_master(A1_idx,9:end) = tmp2;
% fname = sprintf('condition_master_%s%s%s_%s.mat',choose(params.is_wide,'wide_','deep_'), choose(params.is_demo,'demo_',''),params.disp.name,datestr(now,30));
% saveDir = fullfile(vcd_rootPath, 'workspaces','info');
% save(fullfile(saveDir,fname),'condition_master','all_unique_im','all_cond');

%%

% Get stimulus numbers for left and right position in the flagged trials.
stim_nrs_right = condition_master0.stim_nr_right(trial_idx);
stim_nrs_right = stim_nrs_right(stim_nrs_right>0);
stim_nrs_left = condition_master0.stim_nr_left(trial_idx);
stim_nrs_left = stim_nrs_left(stim_nrs_left>0);

% Find the repeated stimulus numbers and their trial position in the stimulus block on either side
repeated_stim_nrs_right = find(accumarray(stim_nrs_right,1)>1);
repeated_stim_nrs_left  = find(accumarray(stim_nrs_left,1)>1);

if ~isempty(repeated_stim_nrs_right)
    rep_idx = find(condition_master0.stim_nr_right(trial_idx)==repeated_stim_nrs_right);
elseif ~isempty(repeated_stim_nrs_left)
    rep_idx = find(condition_master0.stim_nr_left(trial_idx)==repeated_stim_nrs_left);    
end

% For the trial with the repeated stimulus, get the corresponding global block nr, 
% global trial number, crossing number (we want to match this!), 
% cued side (we want to match this!), and correct response (we want to match this!)
global_block_nr = unique(condition_master0.global_block_nr(trial_idx(rep_idx)));
global_trial_nr = condition_master0.global_trial_nr(trial_idx(rep_idx));
rep_cued_side   = condition_master0.is_cued(trial_idx(rep_idx));          % cued side 
rep_corr_resp   = condition_master0.correct_response(trial_idx(rep_idx)); % correct response
cr_nr           = unique(condition_master0.crossing_nr(trial_idx));

% Find stimulus blocks with the same crossing
block_nrs_with_same_crossing = unique(condition_master0.global_block_nr(condition_master0.crossing_nr==cr_nr,:));
block_nrs_with_same_crossing = block_nrs_with_same_crossing(~isnan(block_nrs_with_same_crossing));
block_nrs_with_same_crossing2 = setdiff(block_nrs_with_same_crossing,global_block_nr);

% Find trials within the same session with the crossing, cued side, and
% correct response. We first search nearby, to minimally mess up the distribution of stimuli across sessions. 
if length(block_nrs_with_same_crossing2)>=1
    for ii = 1:length(block_nrs_with_same_crossing2)
        % find blocks with similar crossings
        trial_idx2  = find(condition_master0.global_block_nr == block_nrs_with_same_crossing2(ii));
        for jj = 1:length(rep_cued_side)
            % see if trial matches cued side and correct response
            match_trial = find( ismember(condition_master0.is_cued(trial_idx2),rep_cued_side(jj)) & ismember(condition_master0.correct_response(trial_idx2),rep_corr_resp(jj)) );
            % if we found one, we get out of the for loop
            if ~isempty(match_trial)
                break;
            end
        end
        % Pick the first trial if we have multiple matches
        if length(match_trial)>1
            match_trial = match_trial(1);
        end
        % get the global trial number for the bad trial (with repeat) and new matching trial we want to swap it with 
        A1 = condition_master0.global_trial_nr(trial_idx2(match_trial));
        A2 = global_trial_nr(jj);
        if ~isempty(A1)
            break;
        end
            
    end
end

% If we didn't find any trials within the session, find trials in nearby
% sessions. 
if isempty(block_nrs_with_same_crossing2) || isempty(A1)
    
    % check the other runs
    block_nrs_with_same_crossing = unique(condition_master.global_block_nr(condition_master.crossing_nr==cr_nr,:));
    block_nrs_with_same_crossing = block_nrs_with_same_crossing(~isnan(block_nrs_with_same_crossing));
    block_nrs_with_same_crossing2 = setdiff(block_nrs_with_same_crossing,global_block_nr);
    
    % with the blocks that are the closest to the block with the trial we
    % want to swap
    closest_blocks = abs(block_nrs_with_same_crossing2-global_block_nr);
    [~,ri] = sort(closest_blocks,'ascend');
    block_nrs_with_same_crossing3 = block_nrs_with_same_crossing2(ri);
    
    for ii = 1:length(block_nrs_with_same_crossing3)
        % find blocks with similar crossings
        trial_idx2  = find(condition_master.global_block_nr == block_nrs_with_same_crossing3(ii));
        if ~isempty(trial_idx2)
            for jj = 1:length(rep_cued_side)
                % see if trial matches cued side and correct response
                match_trial = find( ismember(condition_master.is_cued(trial_idx2),rep_cued_side(jj)) & ismember(condition_master.correct_response(trial_idx2),rep_corr_resp(jj)) );
                % if we found one, we get out of the for loop
                if ~isempty(match_trial)
                    break;
                end
            end
            % Pick the first trial if we have multiple matches
            if length(match_trial)>1
                match_trial = match_trial(1);
            end
            % get the global trial number for the bad trial (with repeat) and new matching trial we want to swap it with 
            A1 = condition_master.global_trial_nr(trial_idx2(match_trial));
            A2 = global_trial_nr(jj);
            if ~isempty(A1)
                break;
            end
        end
    end
end

return

