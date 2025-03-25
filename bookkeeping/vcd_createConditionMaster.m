function cond_master = vcd_createConditionMaster(p, unique_im, n_unique_cases, n_trials_per_block, stimClass, taskClass)
% VCD function to create condition master table.
%
%  cond_master = vcd_createConditionMaster(p, unique_im, stimClass, use_fix_flag)
%
% Purpose:
% This function creates a table with N rows (trials) by M columns (stimulus
% class conditions) for the requested stimulus class. The trials are
% shuffled pseudorandomly, such that the subject sees all K unique
% number of cases occur every K trials, and where prioritized stimulus
% features of interest (e.g., 3 gabor contrast levels) are presented within
% a miniblock.
%
%
% INPUTS:
%  p                    : (struct) stimulus and experimental session parameters
%  unique_im            : (matrix) unique stimulus images x M conditions
%  n_unique_cases       : (int) number of unique cases. NOTE: this
%                           number is NOT the same as size(unique_im,1), as
%                           unique cases only counts the fully-crossed
%                           stim features of interest.
%  n_trials_per_block   : (int) number of trials per miniblock (4 or 8)
%  stimClass            : (str) stimulus super class
%  taskClass            : (str) task super class
%
% OUTPUTS:
%  cond_master          : (matrix) master table with N rows for each trial
%                          and M columns containing stimulus/task condition
%                          information. Trials are shuffled pseudorandomly
%                          to ensure that (1) we show prioritized stimulus
%                          features  within a miniblock and (2) all unique
%                          stimulus images and cueing cases within the
%                          fewest number of repeats.
%
% -- MORE INFO ABOUT CONDITION MASTER TABLE --
% Rows contain number of trials, where trials are defined by the unique
% cases x nr of repeats.
%   > Unique cases are defined by the number of unique images (e.g., 8
%   gabors orientations x 3 contrast levels) x number of unique cuing
%   scenarios (i.e., being cued or not cued).
%
% For Gabor/RDK/Dot/Complex Objects crossed with NON-FIX tasks (i.e.,
% cd/scc/pc/wm/ltm/img/what/how), each unique  image will have a set
% stimulus location (either left or right from fixation) and will be cued
% or uncued. This results in conds_master with dims N rows x M cols.
%   If the task is FIX (i.e., Did fixation dot luminance become brighter or
% darker?), we show each unique image both on the left and right side of
% fixation, resulting in conds_master dims: 2xN rows x M cols.
%
% For Natural scenes, there is only one central stimulus, hence the number
% of unique images is identical for each task.
%
% Columns contain stimulus conditions, where the first 5 and the last 3 are the
% same for each stimulus class:
%   1: trial nr       (for parafoveal stimulus patches: 1 trial occupies 2 rows)
%   2: thickening_dir (1 = left, 2 = right, 3 = both/neutral)
%   3: unique im nr   (1-K, where K depends on stimClass)
%   4: stim loc       (1 = left, 2 = right, 3 = central)
%   5: cue status     (0 = uncued, 1 = cued)
%  ...
%  10: paired_stim    (for LTM task: each unique image nr is associated with
%                       another stimulus)
%  11: lure           (for LTM task: XX of the trials we show a lure stimulus)
%  12: repeat nr      (Keep track how many times has this unique image has
%                       been repeated thusfar in the experiment)
%
% The other columns contain stimulus info related to each class.
% For Gabors:
%   col 6: angle
%   col 7: contrast
%   col 8: phase
%   col 9: delta ref (WM query stim)
%
% For RDK:
%   col 6: motion direction
%   col 7: coherence
%   col 8: delta ref (WM query stim)
%
% For Dot:
%   col 6: angle
%   col 9: delta ref (WM query stim)
%
% For Complex Objects (cobj):
%   col 6: super_cat
%   col 7: basic_cat
%   col 8: sub_cat
%   col 9: facing direction
%   col 10: delta ref (WM query stim)

% Natural scenes (NS):
%   col 6: super_cat
%   col 7: basic_cat
%   col 8: sub_cat
%   col 9: change_blindness (WM query stim)
%   col 10: lure nr
%
% Written by Eline Kupers Feb 2025 @ UMN


%% Create variable
cond_master = [];

switch stimClass
    
    case 'gabor'
        
        % Get gabor stim manipulations
        n_contrasts    = length(p.stim.gabor.contrast);
        assert(isequal(length(unique(unique_im(:,4))),n_contrasts));
        
        n_ori_bins     = length(p.stim.gabor.ori_deg);
        assert(isequal(length(unique(unique_im(:,3))),n_ori_bins));
        
        if strcmp(taskClass,'fix')
            stimloc_cues   = NaN; % {1:cued,0:uncued, NaN:no cue} or 3??
            n_stimloc_cues = length(stimloc_cues);
        else
            stimloc_cues   = [1,0]; %{1:cued,0:uncued, NaN:no cue, 3:both sides cued?}
            n_stimloc_cues = length(stimloc_cues);
        end
        
        for rep = 1:p.exp.n_unique_trial_repeats
            
            conds_single_rep = [];
            
            cue_vec = []; im_cue_vec_loc1 = [];
            for loc = 1:n_stimloc_cues % cued vs uncued
                
                % shuffle 8 orientations within every 8 unique images
                left_ori = unique(unique_im(unique_im(:,2)==1,3));
                right_ori = unique(unique_im(unique_im(:,2)==2,3));
                
                shuffle_ori = shuffle_concat([1:n_ori_bins],n_contrasts);
                shuffle_ori = shuffle_ori' + repelem([0:n_ori_bins:(n_unique_cases-1)],n_ori_bins)';
                assert(isequal(unique(shuffle_ori),[1:(n_ori_bins*n_contrasts)]'))
                
                % create a contrast case vector where im 1-8 is low contrast, 9-17 is mid and 18-24 is high
                contrast_vec0 = reshape(1:n_unique_cases,[],n_contrasts);
                contrast_vec0 = contrast_vec0';
                contrast_vec0 = contrast_vec0(:); % we want low,mid,high followed by low,mid,high, etc.
                
                % shuffle contrasts every 3 trials
                shuffle_c  = shuffle_concat(1:n_contrasts,n_unique_cases/n_contrasts);
                shuffle_c  = contrast_vec0(shuffle_c + repelem([0:n_contrasts:(n_unique_cases-1)],n_contrasts));
                
                assert(isequal(unique(shuffle_c),[1:n_unique_cases]'))
                
                % SHUFFLE ORDER OF UNIQUE IMAGES
                conds_shuffle0 = unique_im(shuffle_ori,:); % <-- first shuffle based on unique nr of orientations
                conds_shuffle0 = conds_shuffle0(shuffle_c,:); % <-- then shuffle based on unique nr of contrasts, such that new trial order prioritizes contrast order
                
                % get covert spatial attention cuing direction vector
                [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, unique_im, stimloc_cues, taskClass);
                
                % Insert cue vec as 3rd column
                conds_shuffle1_c = [conds_shuffle0, cue_vec];
                newOrder = [1,2,6,3,4,5];  % unique im nr, stim loc, cue stat, orient, contrast, phase
                conds_shuffle1 = conds_shuffle1_c(:,newOrder); clear conds_shuffle1_c;
                
                % Do some checks:
                assert(isequal(sort(conds_shuffle1(:,1),'ascend'),[1:n_unique_cases]')); % Check unique image nr
                assert(isequal(sum(conds_shuffle1(:,2)==1),n_unique_cases/2));  % Check left stim loc
                assert(isequal(sum(conds_shuffle1(:,2)==2),n_unique_cases/2));  % Check right stim loc
                assert(isequal(sort(conds_shuffle1(:,4),'ascend'),repelem(p.stim.gabor.ori_deg,n_contrasts)')); % Check orientations
                assert(isequal(sort(conds_shuffle1(:,5),'ascend'),repelem(p.stim.gabor.contrast,n_ori_bins)')); % Check contrast
                assert(isequal(sort(conds_shuffle1(:,6),'ascend'),repelem(p.stim.gabor.ph_deg,n_unique_cases/length(p.stim.gabor.ph_deg))')); % Check phase
                
                assert(isequal(sort(reshape(conds_shuffle1(:,5),n_contrasts,[])),repmat(p.stim.gabor.contrast',1,n_ori_bins))); % Check contrast order
                
                % Vertcat single repeat to create this master table
                conds_single_rep = [conds_single_rep;conds_shuffle1];
            end
            
            % Merge trials and add fix cue thickening direction
            conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_single_rep, taskClass);
            
            % Add WM change
            if strcmp(taskClass,'wm')
                n_deltas = length(p.stim.gabor.delta_from_ref);
                shuffle_delta = shuffle_concat(1:length(p.stim.gabor.delta_from_ref), (n_unique_cases/(n_deltas/2)));
                delta_vec = p.stim.gabor.delta_from_ref(shuffle_delta)';
                orient2_vec = conds_single_rep_merged(:,6) + delta_vec;
            else
                delta_vec = NaN(size(conds_single_rep_merged,1),1);
                orient2_vec = NaN(size(conds_single_rep_merged,1),1);
            end
            
            %% TODO Preallocate space for LTM pairing (we'll do this later)
            pair_vec = NaN(size(conds_single_rep_merged,1),1);
            lure_vec = pair_vec; % should be boolean
            
            % Keep track of the times a unique condition has been repeated
            rep_vec = rep.*ones(size(conds_single_rep_merged,1),1);
            
            conds_single_rep_merged2 = [conds_single_rep_merged, delta_vec, orient2_vec, pair_vec, lure_vec, rep_vec];
            
            %% Accummulate
            cond_master = [cond_master; conds_single_rep_merged2];
            clear conds_single_rep_merged2 conds_master_single_rep
            
            % thickening direction doesn't have to match
            % between left and right, from my
            % understanding...
            if ~strcmp(taskClass,'fix')
                assert(isequal(sum(cond_master(:,2)==1),sum(cond_master(:,2)==2)))
            end
        end
        
    case 'rdk'
        
        % Get stim manipulations
        n_coh       = length(p.stim.rdk.dots_coherence); % levels of dot coherence
        n_motdir    = length(p.stim.rdk.dots_direction); % number of motion direction bins per hemifield
        
        if strcmp(taskClass,'fix')
            stimloc_cues   = NaN; % {1:cued,0:uncued, NaN:neutral/no cue} or 3??
            n_stimloc_cues = length(stimloc_cues);
        else
            stimloc_cues = [1,0]; % {1:cued,0:uncued, NaN:neutral/no cue}
            n_stimloc_cues = length(stimloc_cues);
        end
        
        for rep = 1:p.exp.n_unique_trial_repeats
            
            conds_master_single_rep = [];
            
            cue_vec = []; im_cue_vec_loc1 = [];
            for loc = 1:n_stimloc_cues % cued vs uncued
                
                % shuffle orientation every 8 trials
                shuffle_motdir = shuffle_concat([1:n_motdir],n_coh);
                shuffle_motdir = shuffle_motdir' + repelem([0:n_motdir:(n_unique_cases-1)],n_motdir)';
                
                % shuffle coherence every 3 trials
                shuffle_c = shuffle_concat(1:n_coh,n_unique_cases/n_coh);
                case_vec = reshape(1:n_unique_cases,[],3);
                
                % Reshape number of unique image cases, such that we can
                % allocate unique image nrs across combinations of coh and
                % motdir
                case_vec = case_vec';
                case_vec = case_vec(:);
                shuffled_c = case_vec(shuffle_c + repelem([0:n_coh:(n_unique_cases-1)],n_coh),:);
                
                % SHUFFLE ORDER OF UNIQUE IMAGES
                conds_shuffle0 = unique_im(shuffle_motdir,:); % <-- first shuffle based on unique nr of orientations
                conds_shuffle0 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of coh, such that new trial order prioritized coh order
                
                % get covert spatial attention cuing direction vector
                [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, unique_im, stimloc_cues, taskClass);
                
                % Insert cue vec as 3rd column
                conds_shuffle1_c = [conds_shuffle0, cue_vec];
                newOrder         = [1,2,5,3,4];  % unique im nr, stim loc, cue stat, motdir, coh
                conds_shuffle1 = conds_shuffle1_c(:,newOrder); clear conds_shuffle1_c;
                
                % Do some checks:
                assert(isequal(sort(conds_shuffle1(:,1),'ascend'),[1:n_unique_cases]')); % Check unique image nr
                assert(isequal(sum(conds_shuffle1(:,2)==1),n_unique_cases/2));  % Check left stim loc
                assert(isequal(sum(conds_shuffle1(:,2)==2),n_unique_cases/2));  % Check right stim loc
                assert(isequal(sort(conds_shuffle1(:,4),'ascend'),repelem(p.stim.rdk.dots_direction,n_coh)')); % Check motdir
                assert(isequal(sort(conds_shuffle1(:,5),'ascend'),repelem(p.stim.rdk.dots_coherence,n_motdir)')); % Check coherence
                
                assert(isequal(sort(reshape(conds_shuffle1(:,5),n_coh,[])),repmat(p.stim.rdk.dots_coherence',1,n_motdir))); % Check coh order
                
                conds_master_single_rep = [conds_master_single_rep;conds_shuffle1];
            end
            
            clear conds_shuffle0 conds_shuffle1
            
            % Merge unique im into trials and add fix cue thickening direction.
            conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep, taskClass);
            
            % Add WM change.
            if strcmp(taskClass,'wm')
                n_deltas = length(p.stim.rdk.delta_from_ref);
                shuffle_delta = shuffle_concat(1:length(p.stim.rdk.delta_from_ref), (n_unique_cases/(n_deltas/2)));
                delta_vec = p.stim.rdk.delta_from_ref(shuffle_delta)';
                motdir2 = conds_single_rep_merged(:,6) + delta_vec;
            else
                delta_vec = NaN(size(conds_single_rep_merged,1),1);
                motdir2 = NaN(size(conds_single_rep_merged,1),1);
            end
           
            
            % Add ltm pair.
            %                         if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
            % %                             pair_match = [];
            % %                             delta_vec = shuffle_concat(p.stim.rdk.ltm_pairs, (n_unique_cases/(n_deltas/2)));
            % %                             pair_vec = conds_master_single_rep3(:,7) + delta_vec';
            %                         else
            pair_vec = NaN(size(conds_single_rep_merged,1),1);
            lure_vec = pair_vec; % should be boolean
            %                         end
            
            % Keep track of repeat
            rep_vec = rep.*ones(size(conds_single_rep_merged,1),1);
            conds_single_rep_merged2 = [conds_single_rep_merged, delta_vec, motdir2, pair_vec, lure_vec, rep_vec];
            
            % Accummulate
            cond_master = [cond_master; conds_single_rep_merged2];
            clear conds_single_rep_merged2 conds_master_single_rep
        end
        
        % thickening direction doesn't have to match
        % between left and right, or do they??...
        if ~strcmp(taskClass,'fix')
            assert(isequal(sum(cond_master(:,2)==1),sum(cond_master(:,2)==2)))
        end
        
        
    case 'dot'
        
        % Get stim manipulations
        n_dot_loc  = size(p.stim.dot.loc_deg,2);
        
        if strcmp(taskClass,'fix')
            stimloc_cues = NaN;  % {1:cued,0:uncued, NaN:neutral/no cue} or 3??
            n_stimloc_cues = length(stimloc_cues);
        else
            stimloc_cues = [1,0];  % {1:cued,0:uncued, NaN:neutral/no cue}
            n_stimloc_cues = length(stimloc_cues);
        end
        
        n_unique_cases = n_dot_loc;
        
        for rep = 1:n_repeats
            
            conds_master_single_rep = [];
            
            cue_vec = []; im_cue_vec_loc1 = [];
            for loc = 1:n_stimloc_cues % cued vs uncued
                
                % shuffle dot angle every 8 trials
                shuffle_loc = shuffle_concat([1:n_dot_loc],1);
                
                % SHUFFLE TRIALS
                conds_shuffle0 = unique_im(shuffle_loc,:); % <-- shuffle based on unique nr of orientations
                
                % get covert spatial attention cuing direction vector
                [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, unique_im, stimloc_cues, taskClass);
                
                % Insert cue vec as 3rd column
                conds_shuffle1_c = [conds_shuffle0, cue_vec];
                newOrder         = [1,2,4,3];  % unique im nr, stim loc, cue stat, dot_angle
                conds_shuffle1 = conds_shuffle1_c(:,newOrder); clear conds_shuffle1_c;
                
                % Do some checks:
                assert(isequal(sort(conds_shuffle1(:,1),'ascend'),[1:n_unique_cases]')); % Check unique image nr
                assert(isequal(sum(conds_shuffle1(:,2)==1),n_unique_cases/2));  % Check left stim loc
                assert(isequal(sum(conds_shuffle1(:,2)==2),n_unique_cases/2));  % Check right stim loc
                assert(isequal(sort(conds_shuffle1(:,4),'ascend'),repelem(p.stim.dot.ang_deg,1)')); % Check dot angle
                
                conds_master_single_rep = [conds_master_single_rep;conds_shuffle1];
                
            end
            clear conds_shuffle0 conds_shuffle1
            
            % Merge trials and add fix cue thickening direction
            conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep,taskClass);
            
            % Add WM change.
            if strcmp(taskClass,'wm')
                n_deltas = length(p.stim.dot.delta_from_ref);
                shuffle_deta = shuffle_concat(1:length(p.stim.dot.delta_from_ref), (n_unique_cases/(n_deltas/2)));
                delta_vec = p.stim.dot.delta_from_ref(shuffle_deta)';
                loc_deg2 = conds_single_rep_merged(:,6) + delta_vec;
            else
                delta_vec = NaN(size(conds_single_rep_merged,1),1);
                loc_deg2 = NaN(size(conds_single_rep_merged,1),1);
            end
            
            % Add ltm pair.
            %                         if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
            % %                             pair_match = [];
            % %                             delta_vec = shuffle_concat(p.stim.cobj.delta_from_ref, (n_unique_cases/(n_deltas/2)));
            % %                             pair_vec = conds_master_single_rep3(:,7) + delta_vec';
            %                         else
            pair_vec = NaN(size(conds_single_rep_merged,1),1);
            lure_vec = pair_vec; % should be boolean
            %                         end
            
            % Keep track of repeat
            rep_vec = rep.*ones(size(conds_single_rep_merged,1),1);
            conds_single_rep_merged2 = [conds_single_rep_merged, delta_vec, loc_deg2, pair_vec, lure_vec, rep_vec];
            
            % Accummulate
            cond_master = [cond_master; conds_single_rep_merged];
            clear conds_single_rep_merged2 conds_master_single_rep
        end
        % thickening direction doesn't have to match
        % between left and right, or do they??...
        if ~strcmp(taskClass,'fix')
            assert(isequal(sum(cond_master(:,2)==1),sum(cond_master(:,2)==2)))
        end
        
        
    case 'cobj'
        
        % Get stim manipulations
        n_super_cat = length(p.stim.cobj.super_cat);
        stimloc_vec = [1,2]; %{'left','right'}
        
        % UNIQUE COBJ array dims: 8 basic categories
        basic_cat_vec = [];
        for ni = 1:n_super_cat
            tmp = unique(p.stim.cobj.basic_cat{ni}, 'stable');
            n_basic_cat(ni) = length(tmp);
            
            basic_tmp = [];
            for nj = 1:length(tmp)
                basic_tmp = cat(1,basic_tmp,nj.*(arrayfun(@(x) strcmp(x, tmp(nj)),p.stim.cobj.basic_cat{ni})));
            end
            basic_cat_vec = cat(2, basic_cat_vec, sum(basic_tmp,1));
        end
        
        % Get super and sub category info
        super_cat_vec = []; sub_cat_vec = [];
        for ii = 1:length(n_basic_cat)
            n_sub_cat(ii) = length(p.stim.cobj.sub_cat{ii});
            super_cat_vec = cat(2, super_cat_vec, repelem(ii,n_sub_cat(ii)));
            sub_cat_vec = cat(2, sub_cat_vec, 1:n_sub_cat(ii));
        end
        
        if strcmp(taskClass,'fix')
            stimloc_cues = NaN;  % {1:cued,0:uncued, NaN:neutral/no cue} or 3??
            n_stimloc_cues = length(stimloc_cues);
        else
            stimloc_cues = [1,0]; % {1:cued,0:uncued, NaN:neutral/no cue}
            n_stimloc_cues = length(stimloc_cues);
        end
        
        % Shuffle trials for each repeat
        for rep = 1:n_repeats
            
            basic_cat_shuffle_idx = [];
            basic_cat_start = find(sub_cat_vec==1);
            for bi = 1:length(basic_cat_start)
                if bi < length(basic_cat_start)
                    
                    basic_cat_shuffle_idx = cat(2, basic_cat_shuffle_idx, ...
                        shuffle_concat(basic_cat_start(bi):(basic_cat_start(bi+1)-1),1));
                else
                    basic_cat_shuffle_idx = cat(2, basic_cat_shuffle_idx, ...
                        shuffle_concat(basic_cat_start(bi):length(sub_cat_vec),1));
                end
            end
            assert(isequal(1:length(stimloc_vec),unique(basic_cat_shuffle_idx)))
            
            
            conds_master_single_rep = [];
            
            cue_vec = []; im_cue_vec_loc1 = [];
            for loc = 1:n_stimloc_cues % cued vs uncued
                
                % Make sure every 8 trials, subject sees images from 4 super
                % classes
                super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                shuffled_super_cat = super_cat_vec(super_cat_shuffle_idx);
                for jj = 1:n_trials_per_block:(length(shuffled_super_cat))
                    
                    while ~isempty(setdiff([1:n_super_cat],unique(shuffled_super_cat(jj:(jj+n_trials_per_block-1)),'legacy')))
                        super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                        shuffled_super_cat = super_cat_vec(super_cat_shuffle_idx);
                        
                    end
                end
                % Check if n_sub_cat still holds
                assert(isequal(histcounts(shuffled_super_cat),n_sub_cat))
                
                % NOW SHUFFLE!
                conds_shuffle0 = unique_im(basic_cat_shuffle_idx,:); % <-- shuffle based on basic category
                conds_shuffle1 = conds_shuffle0(super_cat_shuffle_idx,:); % <-- shuffle based on super category, last one to give it priority
                
                % get covert spatial attention cuing direction vector
                [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, unique_im, stimloc_cues, taskClass);
                
                % Insert cue vec as 3rd column
                conds_shuffle1_c = [conds_shuffle0, cue_vec];
                newOrder         = [1,2,7,3,4,5,6];  % unique im nr, stim loc, cue stat, super category, basic category,sub category,facing direction
                conds_shuffle1 = conds_shuffle1_c(:,newOrder); clear conds_shuffle1_c;
                
                % Do some checks:
                assert(isequal(sort(conds_shuffle1(:,1),'ascend'),[1:n_unique_cases]')); % Check unique image nr
                assert(isequal(sum(conds_shuffle1(:,2)==1),n_unique_cases/2));  % Check left stim loc
                assert(isequal(sum(conds_shuffle1(:,2)==2),n_unique_cases/2));  % Check right stim loc
                assert(isequal(sort(conds_shuffle1(:,4),'ascend'),repelem(p.stim.cobj.super_cat,n_sub_cat)')); % Check supercat
                
                conds_master_single_rep2 = [conds_master_single_rep;conds_shuffle1];
            end
            
            clear conds_shuffle0 conds_shuffle1
            
            % Merge trials and add fix cue thickening direction
            conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep,taskClass);
            
            % Add WM change.
            if strcmp(taskClass,'wm')
                n_deltas = length(p.stim.cobj.delta_from_ref);
                shuffle_delta = shuffle_concat(1:length(p.stim.cobj.delta_from_ref), (n_unique_cases/(n_deltas/2)));
                delta_vec = p.stim.cobj.delta_from_ref(shuffle_delta)';
                facing_dir2 = conds_single_rep_merged(:,7) + delta_vec;
            else
                delta_vec = NaN(size(conds_single_rep_merged,1),1);
                facing_dir2 = NaN(size(conds_single_rep_merged,1),1);
            end
            
            % Add ltm pair.
            %                         if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
            % %                             pair_match = [];
            % %                             delta_vec = shuffle_concat(p.stim.cobj.delta_from_ref, (n_unique_cases/(n_deltas/2)));
            % %                             pair_vec = conds_master_single_rep3(:,7) + delta_vec';
            %                         else
            pair_vec = NaN(size(conds_single_rep_merged,1),1);
            lure_vec = pair_vec;
            %                         end
            
            % Keep track of repeat
            rep_vec = rep.*ones(size(conds_single_rep_merged,1),1);
            conds_single_rep_merged2 = [conds_single_rep_merged, delta_vec, facing_dir2, pair_vec, lure_vec, rep_vec];
            
            % Accummulate
            cond_master = [cond_master; conds_single_rep_merged];
        end
        
        % thickening direction doesn't have to match
        % between left and right, or do they??...
        if ~strcmp(taskClass,'fix')
            assert(isequal(sum(cond_master(:,2)==1),sum(cond_master(:,2)==2)))
        end
        
    case 'ns'
        
        
        stimloc_cues   = 3; %{1:'cued',0:'uncued',3:'bothcued'}; % should we make bothcued=2?
        n_stimloc_cues = length(stimloc_cues);
        
        for rep = 1:n_repeats
            
            conds_master_single_rep = [];
            
            % Make sure every 6 trials, subject sees images from 4 super
            % classes
            super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
            shuffled_super_cat = super_cat_vec(super_cat_shuffle_idx);
            for jj = 1:6:(length(shuffled_super_cat))
                
                while ~isempty(setdiff([1:n_super_cat],unique(shuffled_super_cat(jj:(jj+6-1)),'legacy')))
                    super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                    shuffled_super_cat = super_cat_vec(super_cat_shuffle_idx);
                    
                end
            end
            % Check if n_sub_cat still holds
            assert(isequal(histcounts(shuffled_super_cat),sum(n_sub_cat,2)'))
            
            % NOW SHUFFLE!
            conds_shuffle0 = unique_im(basic_cat_shuffle_idx,:); % <-- shuffle based on basic category
            conds_shuffle1 = conds_shuffle0(super_cat_shuffle_idx,:); % <-- shuffle based on super category, last one to give it priority
            
            % No need for cued vs uncued status
            cue_vec = NaN(n_unique_cases,1);
            
            % Single stim loc
            stim_loc = ones(size(conds_master_single_rep2,1));
            
            % No need to merge unique im into trials because there is only one
            % stim.
            trial_vec = [1:size(conds_master_single_rep2,1)]';
            
            % Both sides will be thickened given single central image
            thickening_dir = repmat(stimloc_cues,1,size(conds_master_single_rep2,1)/n_stimloc_cues)';
            
            % unique im nr, stim loc, cue stat, unique im nr, stim loc, thickening, super category, basic category,sub category
            conds_master_single_rep = [trial_vec, cue_vec, conds_shuffle1(:,1), stim_loc, thickening_dir, conds_shuffle1(:,2:end)];
           
            % add WM change
            if strcmp(taskClass,'wm')
                n_changes = length(p.stim.ns.change_im);
                change_blindness_vec = shuffle_concat(1:length(p.stim.ns.change_im), ceil(size(trial_vec,1)/n_changes))';
                change_blindness_vec = change_blindness_vec(1:size(trial_vec,1));
            else
                change_blindness_vec = NaN(size(conds_master_single_rep,1),1);
            end
            
            % add ltm change
            %                         if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
            %                             n_stim_pairs = length(p.stim.ns.stim_pairs_ltm);
            %                             pair_vec = shuffle_concat(1:length(p.stim.ns.stim_pairs_ltm), ceil(size(trial_vec,1)/n_changes))';
            %                             pair_vec = pair_vec(1:size(trial_vec,1));
            %
            %                             n_lures = length(p.stim.ns.lure_im);
            %                             lure_vec = (pair_vec==;
            %                         else
            pair_vec = NaN(size(conds_master_single_rep,1),1);
            lure_vec = pair_vec;
            %                         end
            
            % Keep track of repeat
            rep_vec = rep.*ones(size(conds_master_single_rep,1),1);
            conds_master_single_rep2 = [conds_master_single_rep, change_blindness_vec, pair_vec, lure_vec, rep_vec]; % change blindness: 1 = easy, 2 = hard
            
            % Accummulate
            cond_master = [cond_master; conds_master_single_rep2];
            
        end
        
end

return