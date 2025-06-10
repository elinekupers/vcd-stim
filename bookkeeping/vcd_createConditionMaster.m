function cond_master = vcd_createConditionMaster(params, cond_table, n_trials_per_block, session_env)
% VCD function to create condition master table.
%
%  cond_master = vcd_createConditionMaster(params, cond_table, n_trials_per_block)
%
% Purpose:
% This function creates a table with N rows (trials) by M columns (stimulus
% class conditions) for the requested stimulus class. The trials are
% shuffled pseudorandomly, such that the subject sees all K unique
% number of cases occur every K trials, and where prioritized stimulus
% features of interest (e.g., 3 gabor contrast levels) are presented within
% a block.
%
%
% INPUTS:
%  params               : (struct) stimulus and experimental session parameters
%  cond_table           : (table) table with core stimuli information for a
%                          particular stimulus class
%  n_trials_per_block   : (int) number of trials per block (4 or 8)
%  session_env          : (char) character string defining if this is for 
%                          the "MRI" or "BEHAVIORAL" version of the 
%                          vcd-core experiment? Choose: "MRI" or "BEHAVIORAL"
%
% OUTPUTS:
%  cond_master          : (matrix) master table with N rows for each trial
%                          and M columns containing stimulus/task condition
%                          information. Trials are shuffled pseudorandomly
%                          to ensure that (1) we show prioritized stimulus
%                          features  within a block and (2) all unique
%                          stimulus images and cueing cases within the
%                          fewest number of repeats.
%
% -- MORE INFO ABOUT CONDITION MASTER TABLE --
% Rows contains a single trial, with either 1 or 2 unique stimuli,
% then repeated for the n_repeat times.
%
% Unique trials are defined by the number of unique images (e.g., 8
%   gabors orientations x 3 contrast levels) x number of unique cuing
%   scenarios (i.e., being cued or not cued).
%
% For Gabor/RDK/Dot/Complex Objects crossed with each task (i.e., fix/
% cd/scc/pc/wm/ltm/img/what/how), each unique  image will have a fixed
% stimulus location (either left or right from fixation) and will be cued
% or uncued. This results in conds_master with dims N rows x M cols.
% For Natural scenes, there is only one central stimulus, hence the number
% of unique images is identical for each task.
%
% Columns contains the following information:
%   1: {'unique_stim_nr'   } unique images (parafoveal stimulus patches: 1 trial occupies 2 rows)
%   2: {'stimloc'        } stimulus location relative to fixation dot: 1 = left, 2 = right, 3 = central
%   3: {'stimloc_name'   } same as column two but then in text 
%   4: {'orient_dir'     } orientation (gabor), motion direction (rdk), or angle (dot), or facing direction (obj) in deg
%   5: {'contrast'       } stimulus contrast (Michelson fraction)
%   6: {'gbr_phase'      } gabor stimulus phase  (NaN for non gabor stim)
%   7: {'rdk_coherence'  } rdk stimulus dot coherence (fraction of dots, NaN for non rdk stim)
%   8: {'super_cat'      } superordinate object category level (for obj and ns) 
%   9: {'basic_cat'      } basic object category level (for obj and ns)
%   10: {'sub_cat'        } subordinate object category level (for obj and ns)
%   11: {'super_cat_name' } same as 8, but then in text
%   12: {'basic_cat_name' } same as 9, but then in text
%   13: {'sub_cat_name'   } same as 10, but then in text
%   14: {'stim_class_name'} stimulus class name (gabor, rdk, dot, obj, ns)
%   15: {'stim_class'     } stim class nr (1:5)
%   16: {'task_class_name'} task class name (fix/cd/scc/pc/wm/ltm/img/what/where/how)
%   17: {'task_class'     } task class nr (1:10)
%   18: {'is_cued'        } what spatial location did we cue with rim of
%                           fixation circle: 1 = left, 2 = right, 3 = both sides/neutral cue
%   19: {'unique_trial_nr'} unique trial nr
%   20: {'stim2_delta'    } for  WM/IMG/LTM task: what delta gabor tilt/
%                           rdk motion dir / dot angle / object rotation /
%                           scene was used for test stim in second stimulus array (after the delay) (see
%                           p.stim.(xx).delta_from_ref).
%   21: {'stim2_im_nr'    } for  WM/IMG/LTM task: unique test image nr in stimulus array 2
%                           (after the delay). Defined for double epoch
%                           tasks only.
%   22: {'stim2_orient_dir'} for  WM/IMG/LTM task: absolute gabor tilt/rdk motion direction/dot
%                           angle/obj rotation of test stimulus (so after
%                           delta is applied to core stimulus feature)
%   23: {'is_lure'        } for LTM task: is stim2_im_nr a lure stimulus or not.
%   24: {'ltm_stim_pair'  } for LTM task: each unique image nr is associated with
%                           another stimulus
%   25: {'repeat_nr'      } Keep track how many times has this unique image has
%                           been repeated thusfar in the experiment
%
% Written by Eline Kupers Feb 2025 @ UMN


%% Get variables from input table


if (nargin < 4 && ~exist('session_env','var')) || isempty(session_env)
    session_env = 'MRI';
end
    

stimClass = unique(cond_table.stim_class_name);
taskClass = unique(cond_table.task_class_name);
n_unique_cases = length(unique(cond_table.unique_im_nr));

assert(length(stimClass)==1)
assert(length(taskClass)==1)

% Prepare cond_master
cond_master = [];

% get nr of repetitions (depends on session type: mri or behavior)
[~,~,~,~, ~,~, ~, ~, ~, ~, unique_trial_repeats, ~,catch_trial_flag] = vcd_getSessionEnvironmentParams(params, session_env);
nr_reps = ceil(unique_trial_repeats(strcmp(stimClass, params.exp.stimclassnames),strcmp(taskClass, params.exp.taskclassnames)));

% only insert one catch trial every other block, given that special core is 50% or less stimuli per class
if catch_trial_flag && (strcmp(taskClass{:},'ltm') || strcmp(taskClass{:},'img'))
    catch_trial_occurance = shuffle_concat([0,1],ceil(nr_reps/2));
end

if nr_reps > 0

    switch stimClass{:}

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'gabor'

            % Get gabor stim manipulations
            n_contrasts  = length(unique(cond_table.contrast));
            n_ori_bins   = length(unique(cond_table.orient_dir));
            if strcmp(taskClass{:},'ltm') || strcmp(taskClass{:},'img')
                assert(isequal(n_contrasts,1));
                assert(isequal(n_ori_bins,length(params.stim.gabor.ori_deg)));
            else
                assert(isequal(n_contrasts,length(params.stim.gabor.contrast)));
                assert(isequal(n_ori_bins,length(params.stim.gabor.ori_deg)));
            end

            if strcmp(taskClass{:},'fix')
                stimloc_cues   = NaN; % {1:cued,0:uncued, NaN:no cue} or 3??
                n_cued_stimlocs = length(stimloc_cues);
            else
                stimloc_cues   = [1,0]; %{1:cued,0:uncued, NaN:no cue, 3:both sides cued?}
                n_cued_stimlocs = length(stimloc_cues);
            end


            % loop over repeats of unique trials
            for rep = 1:nr_reps

                % allocate space for single repetition
                conds_single_rep              = [];
                cue_vec                       = []; 
                im_cue_vec_loc1               = [];
                prior_cue_dir_for_catch_trial = [];

                % Loop over stimulus locations being cued 
                for loc = 1:n_cued_stimlocs % cued vs uncued

                    % shuffle 8 orientations within every 8 unique images
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
                    conds_shuffle0 = cond_table(shuffle_ori,:); % <-- first shuffle based on unique nr of orientations
                    conds_shuffle0 = conds_shuffle0(shuffle_c,:); % <-- then shuffle based on unique nr of contrasts, such that new trial order prioritizes contrast order

                    % Add catch trial vector
                    conds_shuffle0.is_catch = false(size(conds_shuffle0,1),1);

                    % get covert spatial attention cuing direction vector
                    [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, stimloc_cues, taskClass);

                    % Add cue vec column
                    conds_shuffle0.is_cued = cue_vec;

                    % Do some checks on unique trial condition master table:
                    assert(isequal(length(conds_shuffle0.unique_im_nr),n_unique_cases)); % Check unique image nr
                    assert(isequal(sum(conds_shuffle0.stimloc==1),size(conds_shuffle0,1)/2));  % Check left stim loc
                    assert(isequal(sum(conds_shuffle0.stimloc==2),size(conds_shuffle0,1)/2));  % Check right stim loc
                    assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),repelem(params.stim.gabor.ori_deg,n_contrasts)')); % Check orientations
                    assert(isequal(sort(conds_shuffle0.gbr_phase,'ascend'),repelem(params.stim.gabor.ph_deg,n_unique_cases/length(params.stim.gabor.ph_deg))')); % Check phase
                    if ~ismember(taskClass{:}, {'ltm','img'})
                        assert(isequal(sort(conds_shuffle0.contrast,'ascend'),repelem(params.stim.gabor.contrast,n_ori_bins)')); % Check contrast levels
                        assert(isequal(sort(reshape(conds_shuffle0.contrast,n_contrasts,[])),repmat(params.stim.gabor.contrast',1,n_ori_bins))); % Check contrast order
                    else
                        assert(isequal(conds_shuffle0.contrast,repelem(params.stim.gabor.contrast(3),n_ori_bins)')); % Check contrast
                    end


                    % Add catch trials
                    if sum(isnan(cue_vec)) == size(conds_shuffle0.stimloc,1)
                        cue_dir_for_catch_trial = [];
                    elseif sum(conds_shuffle0.stimloc==3)==size(conds_shuffle0.stimloc,1)
                        cue_dir_for_catch_trial = [];
                    elseif sum(conds_shuffle0.stimloc==1 & conds_shuffle0.is_cued==1) > sum(conds_shuffle0.stimloc==2 & conds_shuffle0.is_cued==1)  % if we have more lefts than rights
                        cue_dir_for_catch_trial = 2; % catch trial will cue right
                    elseif sum(conds_shuffle0.stimloc==2 & conds_shuffle0.is_cued==1) > sum(conds_shuffle0.stimloc==1 & conds_shuffle0.is_cued==1)   % if we have more rights than lefts
                        cue_dir_for_catch_trial = 1; % catch trial will cue right
                    elseif sum(conds_shuffle0.stimloc==2 & conds_shuffle0.is_cued==1) == sum(conds_shuffle0.stimloc==1 & conds_shuffle0.is_cued==1) 
                        if ~isempty(prior_cue_dir_for_catch_trial)
                            cue_dir_for_catch_trial = setdiff([1,2],prior_cue_dir_for_catch_trial); % do the opposite of previous trial
                        else
                            cue_dir_for_catch_trial = shuffle_concat([1,2],1); % Define catch trial stim loc cue dir
                            cue_dir_for_catch_trial = cue_dir_for_catch_trial(1);
                        end
                    end

                    % Check if we want to skip a catch trial
                    if catch_trial_flag
                        if ismember(taskClass{:}, {'ltm','img'})
                            n_catch_trials = catch_trial_occurance(rep); 
                        else
                            n_catch_trials = 1;
                        end
                    else
                        n_catch_trials = 0;
                    end
                    
                    if n_catch_trials > 0 
                        conds_shuffle1 = vcd_addCatchTrials(conds_shuffle0, n_catch_trials, cue_dir_for_catch_trial);
                        prior_cue_dir_for_catch_trial = cue_dir_for_catch_trial;
                    else
                        conds_shuffle1 = conds_shuffle0;
                    end
                    % Vertcat single repeat to create this master table
                    conds_single_rep = [conds_single_rep;conds_shuffle1];

                    clear conds_shuffle0 conds_shuffle1
                end

                % Merge trials and add fix cue thickening direction
                conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_single_rep);

                clear conds_master_single_rep

                % Add WM change
                if strcmp(taskClass{:},'wm')
                    % create offset vector
                    n_deltas            = length(params.stim.gabor.delta_from_ref);
                    delta_vec           = NaN(size(conds_single_rep_merged.orient_dir));
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    delta_vec(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1,1)  = shuffle_concat(repelem(params.stim.gabor.delta_from_ref, n_unique_cases/2/n_deltas),1); % left cued
                    delta_vec(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1),1) = shuffle_concat(repelem(params.stim.gabor.delta_from_ref, n_unique_cases/2/n_deltas),1); % left uncued
                    delta_vec(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2,2)  = shuffle_concat(repelem(params.stim.gabor.delta_from_ref, n_unique_cases/2/n_deltas),1); % right cued
                    delta_vec(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2),2) = shuffle_concat(repelem(params.stim.gabor.delta_from_ref, n_unique_cases/2/n_deltas),1); % right uncued
                    
                    % calculate absolute orientation of stim 2
                    orient_dir2 = conds_single_rep_merged.orient_dir + delta_vec;
                    
                    % find unique WM test im nr
                    idx_wm      = reshape(params.stim.gabor.unique_im_nrs_wm_test,4,[]);
                    [~,idx_l]   = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.gabor.unique_im_nrs_core);
                    [~,idx_r]   = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.gabor.unique_im_nrs_core);
                    assert(isequal(length(idx_l),length(idx_r)))
                    
                    stim2_im_nr = NaN(size(conds_single_rep_merged,1),2);
                    [~,shuffle_delta] = ismember(delta_vec,params.stim.gabor.delta_from_ref);
                    clear tmp;
                    
                    non_nan_trials = find(noncatch_trials_idx);
                    for kk = 1:length(idx_l)
                        [~,idx_ori_l] = find(conds_single_rep_merged.orient_dir(non_nan_trials(kk),1)==params.stim.gabor.ori_deg);
                        [~,idx_ori_r] = find(conds_single_rep_merged.orient_dir(non_nan_trials(kk),2)==params.stim.gabor.ori_deg);
                        tmp(kk,:) = [idx_wm(shuffle_delta(kk,1),idx_ori_l); idx_wm(shuffle_delta(kk,2),idx_ori_r)]';
                    end
                    stim2_im_nr(noncatch_trials_idx,:)  = tmp;
                     % add to table
                     conds_single_rep_merged.stim2_delta       = delta_vec;
                     conds_single_rep_merged.stim2_im_nr       = stim2_im_nr;
                     conds_single_rep_merged.stim2_orient_dir  = orient_dir2;
                     conds_single_rep_merged.is_lure          = false(size(conds_single_rep_merged,1),2); % should be boolean                
                
                elseif strcmp(taskClass{:},'img')
                    % Add IMG text prompt and quiz dots
                    n_quiz_images       = length(unique(params.stim.gabor.imagery_quiz_images));
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    shuffle_delta       = [shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images)); ...
                                            shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images))]';
                    
                    % create yes/no overlap vector 
                    delta_vec                        = NaN(size(conds_single_rep_merged,1),2);
                    delta_vec(noncatch_trials_idx,:) = shuffle_delta;
                    
                    % find unique IMG test im nr                    
                    [~,idx_l] = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.gabor.unique_im_nrs_specialcore);
                    [~,idx_r] = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.gabor.unique_im_nrs_specialcore);
                    idx_img   = reshape(params.stim.gabor.unique_im_nrs_img_test,length(params.stim.gabor.imagery_quiz_images),[]); % 10 yes + 10 no overlap test images per unique image
                    idx_img_yes = idx_img(params.stim.gabor.imagery_quiz_images==1,:);
                    idx_img_no  = idx_img(params.stim.gabor.imagery_quiz_images==2,:);
                    tmp = NaN(size(delta_vec(noncatch_trials_idx,:)));
                    for ll = unique(idx_l)'
                        tmp( (delta_vec(noncatch_trials_idx,1)==1 & idx_l==ll),1) = idx_img_no(randi(size(idx_img_no,1),[sum(delta_vec(noncatch_trials_idx,1)==1 & idx_l==ll),1]),ll);
                        tmp( (delta_vec(noncatch_trials_idx,1)==2 & idx_l==ll),1) = idx_img_yes(randi(size(idx_img_yes,1),[sum(delta_vec(noncatch_trials_idx,1)==2 & idx_l==ll),1]),ll);
                    end
                    for ll = unique(idx_r)'
                        tmp( (delta_vec(noncatch_trials_idx,2)==1 & idx_r==ll),2) = idx_img_no(randi(size(idx_img_no,1),[sum(delta_vec(noncatch_trials_idx,2)==1 & idx_r==ll),1]),ll);
                        tmp( (delta_vec(noncatch_trials_idx,2)==2 & idx_r==ll),2) = idx_img_yes(randi(size(idx_img_yes,1),[sum(delta_vec(noncatch_trials_idx,2)==2 & idx_r==ll),1]),ll);
                    end
                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx,:)==2),idx_img_yes)));
                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx,:)==1),idx_img_no)));
                    stim2_im_nr = NaN(size(conds_single_rep_merged.orient_dir,1),2);
                    stim2_im_nr(noncatch_trials_idx,:) = tmp;
                    
                    conds_single_rep_merged.stim2_delta    = delta_vec;
                    conds_single_rep_merged.stim2_im_nr    = stim2_im_nr;
                    conds_single_rep_merged.stim2_orient_dir = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure          = false(size(conds_single_rep_merged,1),2); % should be boolean
                    
                % Add LTM pairs
                elseif strcmp(taskClass{:}, 'ltm')
                    conds_single_rep_merged.stim2_im_nr    = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure        = false(size(conds_single_rep_merged,1),2); % should be boolean
                else
                    conds_single_rep_merged.stim2_delta      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_im_nr      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure          = false(size(conds_single_rep_merged,1),2); % should be boolean
                end
                
                % Keep track of the times a unique condition has been repeated
                rep_vec = rep.*ones(size(conds_single_rep_merged,1),1);
                conds_single_rep_merged.repeat_nr       = rep_vec;
                unique_trial_nr                         = conds_single_rep_merged.unique_trial_nr + ((rep-1)*max(conds_single_rep_merged.unique_trial_nr));
                conds_single_rep_merged.unique_trial_nr = unique_trial_nr;

                %% Accummulate
                cond_master = [cond_master; conds_single_rep_merged];

                % thickening direction doesn't have to match
                % between left and right, from my
                % understanding...
                if ~strcmp(taskClass{:},'fix')
                    assert(isequal(sum(cond_master.is_cued==1),sum(cond_master.is_cued==2)))
                end
            end

            
            
            
            
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'rdk'

            % Get stim manipulations
            n_coh       = length(unique(cond_table.rdk_coherence)); % levels of dot coherence
            n_motdir    = length(unique(cond_table.orient_dir)); % number of motion direction bins per hemifield

            if ismember(taskClass{:},{'ltm','img'})
                assert(isequal(n_coh,1));
                assert(isequal(n_motdir,length(params.stim.rdk.dots_direction)));
            else
                assert(isequal(n_coh,length(params.stim.rdk.dots_coherence)));
                assert(isequal(n_motdir,length(params.stim.rdk.dots_direction)));
            end

            if strcmp(taskClass{:},'fix')
                stimloc_cues   = NaN; % {1:cued,0:uncued, NaN:neutral/no cue} or 3??
                n_cued_stimlocs = length(stimloc_cues);
            else
                stimloc_cues = [1,0]; % {1:cued,0:uncued, NaN:neutral/no cue}
                n_cued_stimlocs = length(stimloc_cues);
            end

            % Keep track of stim loc cue loc for catch trials.
            prior_cue_dir_for_catch_trial = [];

            for rep = 1:nr_reps

                conds_master_single_rep = [];

                cue_vec = []; im_cue_vec_loc1 = [];
                for loc = 1:n_cued_stimlocs % cued vs uncued

                    % shuffle orientation every 8 trials
                    shuffle_motdir = shuffle_concat([1:n_motdir],n_coh);
                    shuffle_motdir = shuffle_motdir' + repelem([0:n_motdir:(n_unique_cases-1)],n_motdir)';

                    % shuffle coherence every 3 trials
                    shuffle_c = shuffle_concat(1:n_coh,n_unique_cases/n_coh);
                    case_vec  = reshape(1:n_unique_cases,[],n_coh);

                    % Reshape number of unique image cases, such that we can
                    % allocate unique image nrs across combinations of coh and
                    % motdir
                    case_vec = case_vec';
                    case_vec = case_vec(:);
                    shuffled_c = case_vec(shuffle_c + repelem([0:n_coh:(n_unique_cases-1)],n_coh),:);

                    % SHUFFLE ORDER OF UNIQUE IMAGES
                    conds_shuffle0 = cond_table(shuffle_motdir,:); % <-- first shuffle based on unique nr of orientations
                    conds_shuffle0 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of coh, such that new trial order prioritized coh order

                    % get covert spatial attention cuing direction vector
                    [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, stimloc_cues, taskClass);

                    % Add cue vec column
                    conds_shuffle0.is_cued = cue_vec;

                    % Add catch trial vector
                    conds_shuffle0.is_catch = false(size(conds_shuffle0,1),1);

                    % Do some checks:
                    assert(isequal(length(conds_shuffle0.unique_im_nr),n_unique_cases)); % Check unique image nr
                    assert(isequal(sum(conds_shuffle0.stimloc==1),n_unique_cases/2));  % Check left stim loc
                    assert(isequal(sum(conds_shuffle0.stimloc==2),n_unique_cases/2));  % Check right stim loc
                    assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),repelem(params.stim.rdk.dots_direction,n_coh)')); % Check motdir
                    if ~ismember(taskClass{:}, {'ltm','img'})
                        assert(isequal(sort(conds_shuffle0.rdk_coherence,'ascend'),repelem(params.stim.rdk.dots_coherence,n_motdir)')); % Check coherence levels
                        assert(isequal(sort(reshape(conds_shuffle0.rdk_coherence,n_coh,[])),repmat(params.stim.rdk.dots_coherence',1,n_motdir))); % Check coherence order
                    else
                        assert(isequal(conds_shuffle0.rdk_coherence,repelem(params.stim.rdk.dots_coherence(3),n_motdir)')); % Check contrast
                    end


                    % Add catch trials
                    if sum(isnan(cue_vec)) == size(conds_shuffle0.stimloc,1)
                        cue_dir_for_catch_trial = [];
                    elseif sum(conds_shuffle0.stimloc==3)==size(conds_shuffle0.stimloc,1)
                        cue_dir_for_catch_trial = [];
                    elseif sum(conds_shuffle0.stimloc==1 & conds_shuffle0.is_cued==1) > sum(conds_shuffle0.stimloc==2 & conds_shuffle0.is_cued==1)  % if we have more lefts than rights
                        cue_dir_for_catch_trial = 2; % catch trial will cue right
                    elseif sum(conds_shuffle0.stimloc==2 & conds_shuffle0.is_cued==1) > sum(conds_shuffle0.stimloc==1 & conds_shuffle0.is_cued==1)   % if we have more rights than lefts
                        cue_dir_for_catch_trial = 1; % catch trial will cue right
                    elseif sum(conds_shuffle0.stimloc==2 & conds_shuffle0.is_cued==1) == sum(conds_shuffle0.stimloc==1 & conds_shuffle0.is_cued==1) 
                        if ~isempty(prior_cue_dir_for_catch_trial)
                            cue_dir_for_catch_trial = setdiff([1,2],prior_cue_dir_for_catch_trial); % do the opposite of previous trial
                        else
                            cue_dir_for_catch_trial = shuffle_concat([1,2],1); % Define catch trial stim loc cue dir
                            cue_dir_for_catch_trial = cue_dir_for_catch_trial(1);
                        end
                    end

                    if catch_trial_flag
                        % Check if we want to skip a catch trial
                        if ismember(taskClass{:}, {'ltm','img'})
                            n_catch_trials = catch_trial_occurance(rep);
                        else
                            n_catch_trials = 1;
                        end
                    else
                        n_catch_trials = 0;
                    end
                    
                    if  n_catch_trials > 0 
                        % Create catch trial
                        conds_shuffle1 = vcd_addCatchTrials(conds_shuffle0, n_catch_trials, cue_dir_for_catch_trial);
                        prior_cue_dir_for_catch_trial = cue_dir_for_catch_trial; % update cue dir
                    else
                        conds_shuffle1 = conds_shuffle0;
                    end
                    conds_master_single_rep = [conds_master_single_rep;conds_shuffle1];

                    clear conds_shuffle0 conds_shuffle1
                end

                % Merge unique im into trials and add fix cue thickening direction.
                conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep);
                clear conds_master_single_rep

                % Add WM change.
                if strcmp(taskClass{:},'wm')
                    n_deltas            = length(params.stim.rdk.delta_from_ref);
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    delta_vec           = NaN(size(conds_single_rep_merged.orient_dir));
                    delta_vec(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1,1) = shuffle_concat(repelem(params.stim.rdk.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left cued
                    delta_vec(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1),1) = shuffle_concat(repelem(params.stim.rdk.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left uncued
                    delta_vec(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2,2) = shuffle_concat(repelem(params.stim.rdk.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right cued
                    delta_vec(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2),2) = shuffle_concat(repelem(params.stim.rdk.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right uncued
                     
                    % Calculate absolute orientation of stim 2
                    orient_dir2 = conds_single_rep_merged.orient_dir + delta_vec;

                    % find unique WM test im nr
                    idx_wm = reshape(params.stim.rdk.unique_im_nrs_wm_test,4,[]);
                    [~,idx_l] = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.rdk.unique_im_nrs_core);
                    [~,idx_r] = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.rdk.unique_im_nrs_core);
                    assert(isequal(length(idx_l),length(idx_r)))
                    
                    stim2_im_nr = NaN(size(conds_single_rep_merged,1),2);
                    [~,shuffle_delta] = ismember(delta_vec,params.stim.rdk.delta_from_ref);
                    non_nan_trials = find(noncatch_trials_idx);
                    clear tmp;
                    for kk = 1:length(idx_l)
                        idx_ori_l = find(params.stim.rdk.dots_direction==conds_single_rep_merged.orient_dir(non_nan_trials(kk),1));
                        idx_ori_r = find(params.stim.rdk.dots_direction==conds_single_rep_merged.orient_dir(non_nan_trials(kk),2));
                        tmp(kk,:) = [idx_wm(shuffle_delta(kk,1),idx_ori_l); idx_wm(shuffle_delta(kk,2),idx_ori_r)]';
                    end
                    stim2_im_nr(noncatch_trials_idx,:)  = tmp;
                    % add to table
                    conds_single_rep_merged.stim2_delta       = delta_vec; 
                    conds_single_rep_merged.stim2_im_nr       = stim2_im_nr;
                    conds_single_rep_merged.stim2_orient_dir  = orient_dir2;
                    conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged,1),2); % should be boolean


                % Add IMG text prompt and quiz dots
                elseif strcmp(taskClass{:},'img')
                    n_quiz_images       = length(unique(params.stim.rdk.imagery_quiz_images));
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    shuffle_delta       = [shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images)); ...
                                            shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images))]';
                    
                    % create yes/no overlap vector 
                    delta_vec                        = NaN(size(conds_single_rep_merged,1),2);
                    delta_vec(noncatch_trials_idx,:) = shuffle_delta;
                    
                    % find unique IMG test im nr                    
                    [~,idx_l] = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.rdk.unique_im_nrs_specialcore);
                    [~,idx_r] = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.rdk.unique_im_nrs_specialcore);
                    idx_img   = reshape(params.stim.rdk.unique_im_nrs_img_test,length(params.stim.rdk.imagery_quiz_images),[]); % 10 yes + 10 no overlap test images per unique image
                    idx_img_yes = idx_img(params.stim.rdk.imagery_quiz_images==1,:);
                    idx_img_no  = idx_img(params.stim.rdk.imagery_quiz_images==2,:);
                    tmp = NaN(size(delta_vec(noncatch_trials_idx,:)));
                    for ll = unique(idx_l)'
                        tmp( (delta_vec(noncatch_trials_idx,1)==1 & idx_l==ll),1) = idx_img_no(randi(size(idx_img_no,1),[sum(delta_vec(noncatch_trials_idx,1)==1 & idx_l==ll),1]),ll);
                        tmp( (delta_vec(noncatch_trials_idx,1)==2 & idx_l==ll),1) = idx_img_yes(randi(size(idx_img_yes,1),[sum(delta_vec(noncatch_trials_idx,1)==2 & idx_l==ll),1]),ll);
                    end
                    for ll = unique(idx_r)'
                        tmp( (delta_vec(noncatch_trials_idx,2)==1 & idx_r==ll),2) = idx_img_no(randi(size(idx_img_no,1),[sum(delta_vec(noncatch_trials_idx,2)==1 & idx_r==ll),1]),ll);
                        tmp( (delta_vec(noncatch_trials_idx,2)==2 & idx_r==ll),2) = idx_img_yes(randi(size(idx_img_yes,1),[sum(delta_vec(noncatch_trials_idx,2)==2 & idx_r==ll),1]),ll);
                    end
                    
                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx,:)==2),idx_img_yes)));
                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx,:)==1),idx_img_no)));
                    
                    stim2_im_nr = NaN(size(conds_single_rep_merged.orient_dir,1),2);
                    stim2_im_nr(noncatch_trials_idx,:) = tmp;
                    
                    conds_single_rep_merged.stim2_delta    = delta_vec;
                    conds_single_rep_merged.stim2_im_nr    = stim2_im_nr;
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged,1),2); % should be boolean


                % Add ltm pair.
                elseif strcmp(taskClass{:},'ltm')
                    conds_single_rep_merged.stim2_im_nr   = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged,1),2); % should be boolean
                else
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_im_nr       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged,1),2); % should be boolean
                end
                % Keep track of repeat
                conds_single_rep_merged.repeat_nr = rep.*ones(size(conds_single_rep_merged,1),1);

                % Accummulate
                cond_master = [cond_master; conds_single_rep_merged];
            end

            % thickening direction doesn't have to match
            % between left and right, or do they??...
            if ~strcmp(taskClass{:},'fix')
                assert(isequal(sum(cond_master.is_cued==1),sum(cond_master.is_cued==2)))
            end


            
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        case 'dot'

            % Get stim manipulations
            n_dot_loc  = length(unique(cond_table.orient_dir));

            if strcmp(taskClass{:},'ltm') || strcmp(taskClass{:},'img')
                assert(isequal(n_dot_loc,length(params.stim.dot.unique_im_nrs_specialcore)));
            else
                assert(isequal(n_dot_loc,size(params.stim.dot.ang_deg,2)));
            end

            if strcmp(taskClass{:},'fix')
                stimloc_cues = NaN;  % {1:cued,0:uncued, NaN:neutral/no cue} or 3??
                n_cued_stimlocs = length(stimloc_cues);
            else
                stimloc_cues = [1,0];  % {1:cued,0:uncued, NaN:neutral/no cue}
                n_cued_stimlocs = length(stimloc_cues);
            end

            n_unique_cases = n_dot_loc;

            % Keep track of stim loc cue loc for catch trials.
            prior_cue_dir_for_catch_trial = [];

            for rep = 1:nr_reps

                conds_master_single_rep = [];

                cue_vec = []; im_cue_vec_loc1 = [];
                for loc = 1:n_cued_stimlocs % cued vs uncued

                    % shuffle dot angle every 8 trials
                    shuffle_loc = shuffle_concat([1:n_dot_loc],1);

                    % SHUFFLE TRIALS
                    conds_shuffle0 = cond_table(shuffle_loc,:); % <-- shuffle based on unique nr of orientations

                    % get covert spatial attention cuing direction vector
                    [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, stimloc_cues, taskClass);

                    % Add cue vec
                    conds_shuffle0.is_cued = cue_vec;

                    % Add catch trial vector
                    conds_shuffle0.is_catch = false(size(conds_shuffle0,1),1);

                    % Do some checks:
                    assert(isequal(length(conds_shuffle0.unique_im_nr),n_unique_cases)); % Check unique image nr
                    assert(isequal(sum(conds_shuffle0.stimloc==1),n_unique_cases/2));  % Check left stim loc
                    assert(isequal(sum(conds_shuffle0.stimloc==2),n_unique_cases/2));  % Check right stim loc
                    if ~ismember(taskClass{:}, {'ltm','img'})
                        assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),repelem(params.stim.dot.ang_deg,1)')); % Check dot angle
                    else
                        assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),params.stim.dot.ang_deg(ismember(params.stim.dot.unique_im_nrs_core,params.stim.dot.unique_im_nrs_specialcore))')); % Check contrast
                    end


                    % Add catch trials
                    if sum(isnan(cue_vec)) == size(conds_shuffle0.stimloc,1)
                        cue_dir_for_catch_trial = [];
                    elseif sum(conds_shuffle0.stimloc==3)==size(conds_shuffle0.stimloc,1)
                        cue_dir_for_catch_trial = [];
                    elseif sum(conds_shuffle0.stimloc==1 & conds_shuffle0.is_cued==1) > sum(conds_shuffle0.stimloc==2 & conds_shuffle0.is_cued==1)  % if we have more lefts than rights
                        cue_dir_for_catch_trial = 2; % catch trial will cue right
                    elseif sum(conds_shuffle0.stimloc==2 & conds_shuffle0.is_cued==1) > sum(conds_shuffle0.stimloc==1 & conds_shuffle0.is_cued==1)   % if we have more rights than lefts
                        cue_dir_for_catch_trial = 1; % catch trial will cue right
                    elseif sum(conds_shuffle0.stimloc==2 & conds_shuffle0.is_cued==1) == sum(conds_shuffle0.stimloc==1 & conds_shuffle0.is_cued==1) 
                        if ~isempty(prior_cue_dir_for_catch_trial)
                            cue_dir_for_catch_trial = setdiff([1,2],prior_cue_dir_for_catch_trial); % do the opposite of previous trial
                        else
                            cue_dir_for_catch_trial = shuffle_concat([1,2],1); % Define catch trial stim loc cue dir
                            cue_dir_for_catch_trial = cue_dir_for_catch_trial(1);
                        end
                    end

                    if catch_trial_flag
                        % Check if we want to skip a catch trial
                        if ismember(taskClass{:}, {'ltm','img'})
                            n_catch_trials = catch_trial_occurance(rep);
                        else
                            n_catch_trials = 1;
                        end
                    else
                        n_catch_trials = 0;
                    end
                    
                    if  n_catch_trials > 0
                        % create catch trial
                        conds_shuffle1 = vcd_addCatchTrials(conds_shuffle0, n_catch_trials, cue_dir_for_catch_trial);
                        prior_cue_dir_for_catch_trial = cue_dir_for_catch_trial;
                    else
                        conds_shuffle1 = conds_shuffle0;
                    end

                    % Vert cat tables
                    conds_master_single_rep = [conds_master_single_rep;conds_shuffle1];

                    clear conds_shuffle0 conds_shuffle1
                end

                % Merge trials and add fix cue thickening direction
                conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep);
                clear conds_master_single_rep

                % Add WM change.
                if strcmp(taskClass{:},'wm')
                    n_deltas      = length(params.stim.dot.delta_from_ref);
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    delta_vec     = NaN(size(conds_single_rep_merged.orient_dir));

                    while 1
                        % create shuffled delta vector
                        shuffle_delta       = NaN(size(conds_single_rep_merged.orient_dir));
                        delta_vec(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1,1) = shuffle_concat(repelem(params.stim.dot.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left cued
                        delta_vec(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1),1) = shuffle_concat(repelem(params.stim.dot.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left uncued
                        delta_vec(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2,2) = shuffle_concat(repelem(params.stim.dot.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right cued
                        delta_vec(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2),2) = shuffle_concat(repelem(params.stim.dot.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right uncued
                        
                        % calculate absolute angle of stim 2 (test stim)
                        orient_dir2   = conds_single_rep_merged.orient_dir + delta_vec;
                    
                        % Ensure that absolute angles of left and right test
                        % stimuli don't overlap (here we use a threshold of
                        % an angle that is more than 20 degrees difference between the two stimuli)
                        dot_too_close = diff(orient_dir2,[],2) < params.stim.dot.min_ang_distance_test_stim;
                        
                        if sum(dot_too_close)==0
                            break;
                        end
                    end
                    
                    % find unique WM test im nr
                    idx_wm      = reshape(params.stim.dot.unique_im_nrs_wm_test,4,[]);
                    [~,idx_l]   = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.dot.unique_im_nrs_core);
                    [~,idx_r]   = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.dot.unique_im_nrs_core);
                    assert(isequal(length(idx_l),length(idx_r)))

                    stim2_im_nr = NaN(size(conds_single_rep_merged,1),2);
                    [~,shuffle_delta] = ismember(delta_vec,params.stim.dot.delta_from_ref);
                    non_nan_trials = find(noncatch_trials_idx);
                    clear tmp;
                    for kk = 1:length(idx_l)
                        idx_ori_l = find(params.stim.dot.ang_deg==conds_single_rep_merged.orient_dir(non_nan_trials(kk),1));
                        idx_ori_r = find(params.stim.dot.ang_deg==conds_single_rep_merged.orient_dir(non_nan_trials(kk),2));
                        tmp(kk,:) = [idx_wm(shuffle_delta(kk,1),idx_ori_l); idx_wm(shuffle_delta(kk,2),idx_ori_r)]';
                    end
                    stim2_im_nr(noncatch_trials_idx,:)  = tmp;

                    % add conditions to table
                    conds_single_rep_merged.stim2_delta       = delta_vec;
                    conds_single_rep_merged.stim2_im_nr       = stim2_im_nr;
                    conds_single_rep_merged.stim2_orient_dir  = orient_dir2;
                    conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged,1),2); % should be boolean
                

                % Add IMG text prompt and quiz dots
                elseif strcmp(taskClass{:},'img')
                    n_quiz_images       = length(unique(params.stim.dot.imagery_quiz_images));
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    shuffle_delta       = [shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images)); ...
                                            shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images))]';
                    
                    % create yes/no overlap vector 
                    delta_vec                        = NaN(size(conds_single_rep_merged,1),2);
                    delta_vec(noncatch_trials_idx,:) = shuffle_delta;
                    
                    % find unique IMG test im nr                    
                    [~,idx_l] = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.dot.unique_im_nrs_specialcore);
                    [~,idx_r] = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.dot.unique_im_nrs_specialcore);
                    idx_img   = reshape(params.stim.dot.unique_im_nrs_img_test,length(params.stim.dot.imagery_quiz_images),[]); % 10 yes + 10 no overlap test images per unique image
                    idx_img_yes = idx_img(params.stim.dot.imagery_quiz_images==1,:);
                    idx_img_no  = idx_img(params.stim.dot.imagery_quiz_images==2,:);
                    tmp = NaN(size(delta_vec(noncatch_trials_idx,:)));
                    for ll = unique(idx_l)'
                        tmp( (delta_vec(noncatch_trials_idx,1)==1 & idx_l==ll),1) = idx_img_no(randi(size(idx_img_no,1),[sum(delta_vec(noncatch_trials_idx,1)==1 & idx_l==ll),1]),ll);
                        tmp( (delta_vec(noncatch_trials_idx,1)==2 & idx_l==ll),1) = idx_img_yes(randi(size(idx_img_yes,1),[sum(delta_vec(noncatch_trials_idx,1)==2 & idx_l==ll),1]),ll);
                    end
                    for ll = unique(idx_r)'
                        tmp( (delta_vec(noncatch_trials_idx,2)==1 & idx_r==ll),2) = idx_img_no(randi(size(idx_img_no,1),[sum(delta_vec(noncatch_trials_idx,2)==1 & idx_r==ll),1]),ll);
                        tmp( (delta_vec(noncatch_trials_idx,2)==2 & idx_r==ll),2) = idx_img_yes(randi(size(idx_img_yes,1),[sum(delta_vec(noncatch_trials_idx,2)==2 & idx_r==ll),1]),ll);
                    end

                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx,:)==2),idx_img_yes)));
                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx,:)==1),idx_img_no)));
                    
                    stim2_im_nr = NaN(size(conds_single_rep_merged.orient_dir,1),2);
                    stim2_im_nr(noncatch_trials_idx,:) = tmp;
                    
                    conds_single_rep_merged.stim2_delta       = delta_vec;
                    conds_single_rep_merged.stim2_im_nr       = stim2_im_nr;
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged,1),2); % should be boolean
                    
                % Add ltm pair.
                elseif strcmp(taskClass{:},'ltm')
                    conds_single_rep_merged.stim2_im_nr       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged,1),2); % should be boolean
                else
                    conds_single_rep_merged.stim2_im_nr       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged,1),2); % should be boolean
                end

                % Keep track of repeat
                conds_single_rep_merged.repeat_nr = rep.*ones(size(conds_single_rep_merged,1),1);

                % Accummulate
                cond_master = [cond_master; conds_single_rep_merged];
                clear conds_single_rep_merged2 conds_master_single_rep
            end
            % thickening direction doesn't have to match
            % between left and right, or do they??...
            if ~strcmp(taskClass{:},'fix')
                assert(isequal(sum(cond_master.is_cued==1),sum(cond_master.is_cued==2)))
            end


            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'obj'

            % Get stim manipulations
            n_super_cat = length(unique(cond_table.super_cat));
            assert(isequal(length(params.stim.obj.super_cat),n_super_cat));

            % UNIQUE COBJ array dims: 8 basic categories
            basic_cat_vec = []; n_basic_cat = [];
            for ni = 1:n_super_cat
                tmp = intersect(unique(cond_table.basic_cat_name,'stable'),unique(params.stim.obj.basic_cat{ni}, 'stable')','stable');
                n_basic_cat(ni) = length(tmp);

                basic_tmp = [];
                for nj = 1:length(tmp)
                    if ismember(taskClass{:},{'ltm','img'})
                        [~,idx] = intersect(unique(params.stim.obj.basic_cat{ni},'stable'),tmp(nj));
                        basic_tmp = cat(2,basic_tmp,idx);
                    else
                        basic_tmp = cat(1,basic_tmp,nj.*(arrayfun(@(x) strcmp(x, tmp(nj)),params.stim.obj.basic_cat{ni})));
                    end
                end
                basic_cat_vec = cat(2, basic_cat_vec, sum(basic_tmp,1));
            end
            assert(isequal(basic_cat_vec',cond_table.basic_cat))

            % Get super and sub category info
            super_cat_vec = []; sub_cat_vec = []; n_sub_cat = [];
            for ii = 1:length(n_basic_cat)
                if ismember(taskClass{:},{'ltm','img'})
                    [~,idx] = intersect(unique(params.stim.obj.sub_cat{ii}, 'stable'),unique(cond_table.sub_cat_name,'stable'),'stable');
                    n_sub_cat(ii) = length(idx);
                    super_cat_vec = cat(2, super_cat_vec, repelem(ii,n_sub_cat(ii)));
                    sub_cat_vec = cat(2, sub_cat_vec, idx');
                else
                    n_sub_cat(ii) = length(intersect(unique(cond_table.sub_cat_name,'stable'),unique(params.stim.obj.sub_cat{ii}, 'stable'))); %length(p.stim.obj.sub_cat{ii});
                    super_cat_vec = cat(2, super_cat_vec, repelem(ii,n_sub_cat(ii)));
                    sub_cat_vec = cat(2, sub_cat_vec, 1:n_sub_cat(ii));
                end
            end
            assert(isequal(sub_cat_vec',cond_table.sub_cat));

            if strcmp(taskClass{:},'fix')
                stimloc_cues = NaN;  % {1:cued,0:uncued, NaN:neutral/no cue} or 3??
                n_cued_stimlocs = length(stimloc_cues);
            else
                stimloc_cues = [1,0]; % {1:cued,0:uncued, NaN:neutral/no cue}
                n_cued_stimlocs = length(stimloc_cues);
            end

            % Keep track of stim loc cue loc for catch trials.
            prior_cue_dir_for_catch_trial = [];

            % Shuffle trials for each repeat
            for rep = 1:nr_reps

                conds_master_single_rep = [];

                % We prioritize super category level, but will shuffle unique object images within basic category
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
                assert(isequal(1:length(super_cat_vec),unique(basic_cat_shuffle_idx)));

                % Get cueing vector
                cue_vec = []; im_cue_vec_loc1 = [];
                for loc = 1:n_cued_stimlocs % cued vs uncued

                    % Make sure every 8 trials, subject sees images from at least 4 super
                    % classes
                    super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                    shuffled_super_cat    = super_cat_vec(super_cat_shuffle_idx);
                    for jj = 1:n_trials_per_block:(length(shuffled_super_cat))

                        while 1
                            if isempty(setdiff([1:n_super_cat],unique(shuffled_super_cat(jj:(jj+n_trials_per_block-1)),'legacy'))) || ...
                                    isequal(setdiff([1:n_super_cat],unique(shuffled_super_cat(jj:(jj+n_trials_per_block-1)),'legacy')),1)
                                break;
                            end
                            super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                            shuffled_super_cat = super_cat_vec(super_cat_shuffle_idx);

                        end
                    end
                    % Check if n_sub_cat still holds
                    assert(isequal(histcounts(shuffled_super_cat),n_sub_cat));

                    % NOW SHUFFLE!
                    conds_shuffle0 = cond_table(basic_cat_shuffle_idx,:); % <-- shuffle based on basic category
                    conds_shuffle1 = conds_shuffle0(super_cat_shuffle_idx,:); % <-- shuffle based on super category, last one to give it priority

                    % get covert spatial attention cuing direction vector
                    [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle1, stimloc_cues, taskClass);

                    % Add cue vec
                    conds_shuffle1.is_cued = cue_vec;

                    % Add catch trial vector
                    conds_shuffle1.is_catch = false(size(conds_shuffle1,1),1);

                    % Do some checks:
                    assert(isequal(length(conds_shuffle0.unique_im_nr),n_unique_cases)); % Check unique image nr
                    assert(isequal(sum(conds_shuffle1.stimloc==1),n_unique_cases/2));  % Check left stim loc
                    assert(isequal(sum(conds_shuffle1.stimloc==2),n_unique_cases/2));  % Check right stim loc
                    assert(isequal(sort(conds_shuffle1.super_cat_name),sort(repelem(params.stim.obj.super_cat,n_sub_cat)'))); % Check supercat
                    if ~ismember(taskClass{:}, {'ltm','img'})
                        assert(isequal(sort(conds_shuffle1.orient_dir,'ascend'),sort(params.stim.obj.facing_dir_deg,'ascend')')); % Check object facing dir
                    else
                        assert(isequal(sort(conds_shuffle1.orient_dir,'ascend'),sort(params.stim.obj.facing_dir_deg(ismember(params.stim.obj.unique_im_nrs_core,params.stim.obj.unique_im_nrs_specialcore))','ascend'))); % Check contrast
                    end


                     % Add catch trials
                    if sum(isnan(cue_vec)) == size(conds_shuffle1.stimloc,1)
                        cue_dir_for_catch_trial = [];
                    elseif sum(conds_shuffle1.stimloc==3)==size(conds_shuffle1.stimloc,1)
                        cue_dir_for_catch_trial = [];
                    elseif sum(conds_shuffle1.stimloc==1 & conds_shuffle1.is_cued==1) > sum(conds_shuffle1.stimloc==2 & conds_shuffle1.is_cued==1)  % if we have more lefts than rights
                        cue_dir_for_catch_trial = 2; % catch trial will cue right
                    elseif sum(conds_shuffle1.stimloc==2 & conds_shuffle1.is_cued==1) > sum(conds_shuffle1.stimloc==1 & conds_shuffle1.is_cued==1)   % if we have more rights than lefts
                        cue_dir_for_catch_trial = 1; % catch trial will cue right
                    elseif sum(conds_shuffle1.stimloc==2 & conds_shuffle1.is_cued==1) == sum(conds_shuffle1.stimloc==1 & conds_shuffle1.is_cued==1) 
                        if ~isempty(prior_cue_dir_for_catch_trial)
                            cue_dir_for_catch_trial = setdiff([1,2],prior_cue_dir_for_catch_trial); % do the opposite of previous trial
                        else
                            cue_dir_for_catch_trial = shuffle_concat([1,2],1); % Define catch trial stim loc cue dir
                            cue_dir_for_catch_trial = cue_dir_for_catch_trial(1);
                        end
                    end

                    if catch_trial_flag
                        % Check if we want to skip a catch trial
                        if ismember(taskClass{:}, {'ltm','img'})
                            n_catch_trials = catch_trial_occurance(rep);
                        else
                            n_catch_trials = 1;
                        end
                    else
                        n_catch_trials = 0;
                    end
                    
                    if  n_catch_trials > 0
                        % create catch trial
                        conds_shuffle2 = vcd_addCatchTrials(conds_shuffle1, n_catch_trials, cue_dir_for_catch_trial);
                        prior_cue_dir_for_catch_trial = cue_dir_for_catch_trial;
                    else
                        conds_shuffle2 = conds_shuffle1;
                    end

                    % Vert cat tables
                    conds_master_single_rep = [conds_master_single_rep;conds_shuffle2];

                    clear conds_shuffle0 conds_shuffle1 conds_shuffle2
                end

                % Merge trials and add fix cue thickening direction
                conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep);
                clear conds_master_single_rep;

                % Add WM change.
                if strcmp(taskClass{:},'wm')
                    % create shuffled delta vector
                    n_deltas            = length(params.stim.obj.delta_from_ref);
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    shuffle_delta       = NaN(size(conds_single_rep_merged.orient_dir));
                    delta_vec(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1,1) = shuffle_concat(repelem(params.stim.obj.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left cued
                    delta_vec(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1),1) = shuffle_concat(repelem(params.stim.obj.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left uncued
                    delta_vec(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2,2) = shuffle_concat(repelem(params.stim.obj.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right cued
                    delta_vec(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2),2) = shuffle_concat(repelem(params.stim.obj.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right uncued

                    facing_dir2   = conds_single_rep_merged.orient_dir + delta_vec;
                   
                    % find unique WM test im nr
                    idx_wm = reshape(params.stim.obj.unique_im_nrs_wm_test,4,[]);
                    [~,idx_l] = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.obj.unique_im_nrs_core);
                    [~,idx_r] = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.obj.unique_im_nrs_core);
                    assert(isequal(length(idx_l),length(idx_r)))
                    
                    stim2_im_nr = NaN(size(conds_single_rep_merged,1),2);
                    [~,shuffle_delta] = ismember(delta_vec,params.stim.obj.delta_from_ref);
                    non_nan_trials = find(noncatch_trials_idx);
                    clear tmp;
                    for kk = 1:length(idx_l)
                        idx_ori_l = find(params.stim.obj.facing_dir_deg==conds_single_rep_merged.orient_dir(non_nan_trials(kk),1));
                        idx_ori_r = find(params.stim.obj.facing_dir_deg==conds_single_rep_merged.orient_dir(non_nan_trials(kk),2));
                        tmp(kk,:) = [idx_wm(shuffle_delta(kk,1),idx_ori_l); idx_wm(shuffle_delta(kk,2),idx_ori_r)]';
                    end
                    stim2_im_nr(noncatch_trials_idx,:)       = tmp;
                    conds_single_rep_merged.stim2_delta      = delta_vec;
                    conds_single_rep_merged.stim2_im_nr      = stim2_im_nr; 
                    conds_single_rep_merged.stim2_orient_dir = facing_dir2; 
                    conds_single_rep_merged.is_lure          = false(size(conds_single_rep_merged,1),2); % should be boolean
                

                % Add IMG text prompt and quiz dots
                elseif strcmp(taskClass{:},'img')
                    n_quiz_images       = length(unique(params.stim.obj.imagery_quiz_images));
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    shuffle_delta       = [shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images)); ...
                                            shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images))]';
                    
                    % create yes/no overlap vector
                    delta_vec                        = NaN(size(conds_single_rep_merged,1),2);
                    delta_vec(noncatch_trials_idx,:) = shuffle_delta;
                    
                    % find unique IMG test im nr
                    [~,idx_l] = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.obj.unique_im_nrs_specialcore);
                    [~,idx_r] = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.obj.unique_im_nrs_specialcore);
                    idx_img   = reshape(params.stim.obj.unique_im_nrs_img_test,length(params.stim.obj.imagery_quiz_images),[]); % 10 yes + 10 no overlap test images per unique image
                    idx_img_yes = idx_img(params.stim.obj.imagery_quiz_images==1,:);
                    idx_img_no  = idx_img(params.stim.obj.imagery_quiz_images==2,:);
                    tmp = NaN(size(delta_vec(noncatch_trials_idx,:)));
                    for ll = unique(idx_l)'
                        tmp( (delta_vec(noncatch_trials_idx,1)==1 & idx_l==ll),1) = idx_img_no(randi(size(idx_img_no,1),[sum(delta_vec(noncatch_trials_idx,1)==1 & idx_l==ll),1]),ll);
                        tmp( (delta_vec(noncatch_trials_idx,1)==2 & idx_l==ll),1) = idx_img_yes(randi(size(idx_img_yes,1),[sum(delta_vec(noncatch_trials_idx,1)==2 & idx_l==ll),1]),ll);
                    end
                    for ll = unique(idx_r)'
                        tmp( (delta_vec(noncatch_trials_idx,2)==1 & idx_r==ll),2) = idx_img_no(randi(size(idx_img_no,1),[sum(delta_vec(noncatch_trials_idx,2)==1 & idx_r==ll),1]),ll);
                        tmp( (delta_vec(noncatch_trials_idx,2)==2 & idx_r==ll),2) = idx_img_yes(randi(size(idx_img_yes,1),[sum(delta_vec(noncatch_trials_idx,2)==2 & idx_r==ll),1]),ll);
                    end
                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx,:)==2),idx_img_yes)));
                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx,:)==1),idx_img_no)));

                    stim2_im_nr = NaN(size(conds_single_rep_merged.orient_dir,1),2);
                    stim2_im_nr(noncatch_trials_idx,:) = tmp;
                    
                    conds_single_rep_merged.stim2_delta    = delta_vec;
                    conds_single_rep_merged.stim2_im_nr    = stim2_im_nr;
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure          = false(size(conds_single_rep_merged,1),2); % should be boolean
                
                    % Add ltm pair.
                elseif strcmp(taskClass{:},'ltm')
                    conds_single_rep_merged.stim2_im_nr      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure          = false(size(conds_single_rep_merged,1),2); % should be boolean
                else
                    conds_single_rep_merged.stim2_im_nr      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure          = false(size(conds_single_rep_merged,1),2); % should be boolean
                end

                % Keep track of repeat
                conds_single_rep_merged.repeat_nr = rep.*ones(size(conds_single_rep_merged,1),1);

                % Accummulate
                cond_master = [cond_master; conds_single_rep_merged];
            end

            % thickening direction doesn't have to match
            % between left and right, or do they??...
            if ~strcmp(taskClass{:},'fix')
                assert(isequal(sum(cond_master.is_cued==1),sum(cond_master.is_cued==2)));
            end

            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'ns'

            % Get stim manipulations
            n_super_cat     = length(params.stim.ns.super_cat);
            assert(isequal(length(unique(cond_table.super_cat)),n_super_cat));

            % Get super and sub category info
            bsc_cat = unique(cond_table.basic_cat_name);
            sub_cat = unique(cond_table.sub_cat_name);

            % we want indoor and outdoor scenes, regardless of task
            assert(isequal(bsc_cat', unique(cat(2,params.stim.ns.basic_cat{:}))));

            if ~ismember(taskClass{:},{'ltm','img'})
                assert(isequal(sub_cat', unique(cat(2,params.stim.ns.sub_cat{:}))));
            else
                tmp_sub = cat(2,reshape(cat(1,params.stim.ns.sub_cat{1,1},params.stim.ns.sub_cat{1,2}),1,[]),...
                                reshape(cat(1,params.stim.ns.sub_cat{2,1},params.stim.ns.sub_cat{2,2}),1,[]),...
                                reshape(cat(1,params.stim.ns.sub_cat{3,1},params.stim.ns.sub_cat{3,2}),1,[]),...
                                reshape(cat(1,params.stim.ns.sub_cat{4,1},params.stim.ns.sub_cat{4,2}),1,[]),...
                                reshape(cat(1,params.stim.ns.sub_cat{5,1},params.stim.ns.sub_cat{5,2}),1,[]))';
                [~,idx] = intersect(params.stim.ns.unique_im_nrs_core,params.stim.ns.unique_im_nrs_specialcore);
                assert(isequal(sub_cat, unique(tmp_sub(idx))));
            end


            n_basic_cat = []; n_sub_cat = []; super_cat_vec = []; sub_cat_vec = [];
            if ismember(taskClass{:},{'ltm','img'})
                for ni = 1:n_super_cat
                    n_basic_cat(ni) = length(params.stim.ns.basic_cat{ni});
                    n_sub_cat(ni) = length(cond_table.basic_cat_name(cond_table.super_cat==ni));
                    super_cat_vec = cat(2, super_cat_vec, repelem(ni,n_sub_cat(ni)));
                    sub_cat_vec =  cat(2, sub_cat_vec,[1:n_sub_cat(ni)]);
                end
                assert(isequal(cond_table.sub_cat,sub_cat_vec'));
            else
                for ni = 1:n_super_cat
                    n_basic_cat(ni) = length(params.stim.ns.basic_cat{ni});
                    for nj = 1:n_basic_cat(ni)
                        n_sub_cat(ni,nj) = length(params.stim.ns.sub_cat{ni,nj});
                        super_cat_vec = cat(2, super_cat_vec, repelem(ni,n_sub_cat(ni,nj)));
                        sub_cat_vec   = cat(2, sub_cat_vec, [1:n_sub_cat(ni,nj)]);
                    end
                end
            end
            assert(isequal(cond_table.super_cat,super_cat_vec')); 

            for rep = 1:nr_reps

                conds_master_single_rep = [];

                % Make sure every block, subject sees images from 5 super
                % classes
                super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                shuffled_super_cat    = super_cat_vec(super_cat_shuffle_idx);

                good_blocks = 0;
                block_start_idx = 1:n_trials_per_block:(length(shuffled_super_cat));
                if n_trials_per_block < n_super_cat
                    allowed_mismatch0 = diff([n_trials_per_block,n_super_cat])*2;
                else
                    allowed_mismatch0 = 0;
                end

                while 1
                    super_cat_shuffle_idx = shuffle_concat(super_cat_shuffle_idx,1);
                    shuffled_super_cat = super_cat_vec(super_cat_shuffle_idx);

                    good_blocks = 0;

                    for jj = block_start_idx
                        idx = jj:(jj+n_trials_per_block-1);

                        if idx(end)>length(shuffled_super_cat)
                            idx(idx>length(shuffled_super_cat)) = [];
                            mismatch = diff([length(idx),n_super_cat])*2;
                        else
                            mismatch = allowed_mismatch0;
                        end
                        missing_conditions = setdiff([1:n_super_cat],unique(shuffled_super_cat(idx),'legacy'));

                        if isempty(missing_conditions) || length(missing_conditions) <= mismatch
                            good_blocks = good_blocks+1;
                        end
                    end
                    if good_blocks == length(block_start_idx)-1
                        break;
                    end
                end            

                % Check if n_sub_cat still holds
                if ismember(taskClass{:},{'ltm','img'})
                    assert(isequal(histcounts(shuffled_super_cat),n_sub_cat))
                else
                    assert(isequal(histcounts(shuffled_super_cat),sum(n_sub_cat,2)'))
                end


                % shuffle basic category (so shuffle every other image)
    %             basic_cat_vec = cond_table.basic_cat;
                basic_cat_shuffle_idx = [];
                for ii = 1:length(n_basic_cat)
                    % create [ 0 0 2 2 4 4 ] vector
                    curr_im = repelem(length(basic_cat_shuffle_idx): 2 : length(basic_cat_shuffle_idx)+2*n_basic_cat(ii),n_basic_cat(ii));
                    if ismember(taskClass{:},{'ltm','img'})
                        basic_cat_shuffle_idx = cat(2, basic_cat_shuffle_idx, curr_im + shuffle_concat([1:n_basic_cat(ii)],n_sub_cat(ii)));
                    else
                        % add shuffled [1,2,1,2,1,2] vector
                        basic_cat_shuffle_idx = cat(2, basic_cat_shuffle_idx, curr_im + shuffle_concat([1:n_basic_cat(ii)],n_sub_cat(ii,1)));
                    end
                end
                if ismember(taskClass{:},{'ltm','img'})
                    basic_cat_shuffle_idx = basic_cat_shuffle_idx(ismember(basic_cat_shuffle_idx,1:length(basic_cat_shuffle_idx)/2));
                end
                % NOW SHUFFLE!
                conds_shuffle0 = cond_table(basic_cat_shuffle_idx,:); % <-- shuffle based on basic category
                conds_master_single_rep = conds_shuffle0(super_cat_shuffle_idx,:); % <-- shuffle based on super category, last one to give it priority

                % Add catch trial vector
                conds_master_single_rep.is_catch = false(size(conds_master_single_rep,1),1);

                if catch_trial_flag
                    % Check if we want to skip a catch trial
                    if ismember(taskClass{:}, {'ltm','img'})
                        n_catch_trials = catch_trial_occurance(rep);
                    else
                        n_catch_trials = 1;
                    end
                else
                    n_catch_trials = 0;
                end
                
                if  n_catch_trials > 0
                    % Add catch trials
                    conds_master_single_rep2 = vcd_addCatchTrials(conds_master_single_rep, n_catch_trials, []);
                else
                    conds_master_single_rep2 = conds_master_single_rep;
                end

                % No need to merge unique im into trials because there is only one
                % stim.
                conds_master_single_rep2.unique_trial_nr = [1:size(conds_master_single_rep2,1)]';

                conds_master_single_rep2.stim_nr_left  = conds_master_single_rep2.unique_im_nr;
                conds_master_single_rep2.stim_nr_right = NaN(size(conds_master_single_rep2,1),1);

                conds_master_single_rep2.unique_im_nr = [];
                conds_master_single_rep2.stimloc = [];
                conds_master_single_rep2.stimloc_name = [];

                % Both sides will be thickened given single central image
                conds_master_single_rep2.is_cued = 3.*ones(size(conds_master_single_rep2,1),1);

                % fill second column with nans, so we can merge across stim
                % classes
                fnames = {'orient_dir','contrast','gbr_phase','rdk_coherence',...
                    'super_cat','basic_cat','sub_cat','affordance_cat',...
                    'super_cat_name','basic_cat_name','sub_cat_name','affordance_name','stim_class_name','is_special_core'};

                for fn = 1:length(fnames)
                    if size(conds_master_single_rep2.(fnames{fn}),2) == 1
                        if iscell(conds_master_single_rep2.(fnames{fn}))
                            conds_master_single_rep2.(fnames{fn}) = [conds_master_single_rep2.(fnames{fn}), num2cell(NaN(size(conds_master_single_rep2,1),1))];
                        else
                            conds_master_single_rep2.(fnames{fn}) = [conds_master_single_rep2.(fnames{fn}), NaN(size(conds_master_single_rep2,1),1)];
                        end
                    end
                end

                % reorder columns:
                newOrder = {'unique_trial_nr','stim_class', 'task_class', ...
                    'stim_nr_left','stim_nr_right','is_cued','is_catch', 'stim_class_name','task_class_name'};
                for ii = 1:length(newOrder)
                    newOrderIdx(ii) = find(strcmp(conds_master_single_rep2.Properties.VariableNames,newOrder{ii}));
                end
                newOrderIdx = [newOrderIdx, setdiff([1:size(conds_master_single_rep2,2)],newOrderIdx)];
                conds_master_single_rep2 = conds_master_single_rep2(:,newOrderIdx);


                % add WM change
                if strcmp(taskClass{:},'wm')
                    n_changes            = length(params.stim.ns.change_im);
                    noncatch_trials_idx  = ~conds_master_single_rep2.is_catch(:,1);

                    while 1
                        shuffle_delta    = shuffle_concat(params.stim.ns.change_im, ceil(size(conds_master_single_rep2.unique_trial_nr,1)/n_changes))';
                        delta_vec        = NaN(length(noncatch_trials_idx),2);
                        shuffle_delta    = shuffle_delta(1:sum(noncatch_trials_idx)); % crop length as we generated more deltas than needed (mod(30,4)==2)
                        delta_vec(noncatch_trials_idx,1) = shuffle_delta;
                        n = histcounts(shuffle_delta); n(n==0)=[];
                        % while we can't perfectly balance all delta across
                        % 30 unique core images, we at least try to have
                        % equal nr of "added" and "removed" changes..
                        if sum(n([1,2]))==sum(n([3,4]))
                            break; 
                        end
                    end

                    % find unique WM test im nr
                    [~,idx_c] = ismember(conds_master_single_rep2.stim_nr_left(noncatch_trials_idx),params.stim.ns.unique_im_nrs_core);
                    idx_wm = params.stim.ns.unique_im_nrs_wm_test(1:4:end);
                    stim2_im_nr = NaN(size(conds_master_single_rep2.orient_dir,1),2);
                    stim2_im_nr(noncatch_trials_idx,1)        = [idx_wm(idx_c) + (shuffle_delta(1:n_unique_cases)-1)']';
                    conds_master_single_rep2.stim2_delta      = delta_vec;
                    conds_master_single_rep2.stim2_im_nr      = stim2_im_nr;
                    conds_master_single_rep2.stim2_orient_dir = NaN(size(conds_master_single_rep2.orient_dir,1),2);
                    conds_master_single_rep2.is_lure          = false(size(conds_master_single_rep2,1),2); % should be boolean


                % take subset of IMG trials, add text prompt and quiz dots loc
                elseif strcmp(taskClass{:},'img')
                    n_quiz_images       = length(unique(params.stim.ns.imagery_quiz_images));
                    noncatch_trials_idx = ~(conds_master_single_rep2.is_catch);
                    shuffle_delta       = shuffle_concat(1:n_quiz_images, 1+(size(conds_master_single_rep2(noncatch_trials_idx,:),1)/n_quiz_images))';
                    shuffle_delta       = shuffle_delta(noncatch_trials_idx);
                    % create yes/no overlap vector
                    delta_vec                        = NaN(size(conds_master_single_rep2,1),2);
                    delta_vec(noncatch_trials_idx,1) = shuffle_delta;
                    
                    % find unique IMG test im nr
                    [~,idx_c] = ismember(conds_master_single_rep2.stim_nr_left(noncatch_trials_idx),params.stim.ns.unique_im_nrs_specialcore);
                    idx_img   = reshape(params.stim.ns.unique_im_nrs_img_test,length(params.stim.ns.imagery_quiz_images),[]); % 10 yes + 10 no overlap test images per unique image
                    idx_img_yes = idx_img(params.stim.ns.imagery_quiz_images==1,:);
                    idx_img_no  = idx_img(params.stim.ns.imagery_quiz_images==2,:);
                    
                    non_nan_trials = find(noncatch_trials_idx);
                    tmp = NaN(size(idx_c));
                    for ll = unique(idx_c)'
                        tmp( (delta_vec(noncatch_trials_idx,1)==1 & idx_c==non_nan_trials(ll)),1) = idx_img_no(randi(size(idx_img_no,1),[sum(delta_vec(noncatch_trials_idx,1)==1 & idx_c==non_nan_trials(ll)),1]),non_nan_trials(ll));
                        tmp( (delta_vec(noncatch_trials_idx,1)==2 & idx_c==non_nan_trials(ll)),1) = idx_img_yes(randi(size(idx_img_yes,1),[sum(delta_vec(noncatch_trials_idx,1)==2 & idx_c==non_nan_trials(ll)),1]),non_nan_trials(ll));
                    end
                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx)==2),idx_img_yes)));
                    assert(all(ismember(tmp(delta_vec(noncatch_trials_idx)==1),idx_img_no)));
                    
                    stim2_im_nr = NaN(size(conds_master_single_rep2.stim_nr_left,1),2);
                    stim2_im_nr(noncatch_trials_idx,1) = tmp;
                    
                    conds_master_single_rep2.stim2_delta    = delta_vec;
                    conds_master_single_rep2.stim2_im_nr    = stim2_im_nr;
                    conds_master_single_rep2.stim2_orient_dir  = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.is_lure          = false(size(conds_master_single_rep2,1),2); % should be boolean
                                    

                % Add ltm pair.
                elseif strcmp(taskClass{:},'ltm')
                    conds_master_single_rep2.stim2_im_nr      = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_orient_dir = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_delta      = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.is_lure          = false(size(conds_master_single_rep2,1),2); % should be boolean
                else
                    conds_master_single_rep2.stim2_im_nr   = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_orient_dir = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_delta      = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.is_lure          = false(size(conds_master_single_rep2,1),2); % should be boolean
                end

                % Keep track of repeat
                conds_master_single_rep2.repeat_nr = rep.*ones(size(conds_master_single_rep2,1),1);


                % Accummulate
                cond_master = [cond_master; conds_master_single_rep2];

                clear conds_master_single_rep conds_master_single_rep2

            end

    end
end

return