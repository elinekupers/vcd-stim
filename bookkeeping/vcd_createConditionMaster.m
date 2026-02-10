function cond_master = vcd_createConditionMaster(params, cond_table, env_type)
% VCD function to create condition master table.
%
%  cond_master = vcd_createConditionMaster(params, cond_table, n_trials_per_block)
%
% Purpose:
% This function creates a table with N rows (trials) by M columns (stimulus
% class conditions) for the requested stimulus class. The trials are
% shuffled per set of unique stimuli x cuing options, such that the subject 
% sees all K unique conditions every K trials. 
%
%
% INPUTS:
%  params               : (struct) stimulus and experimental session parameters
%  cond_table           : (table) table with core stimuli information for a
%                          particular stimulus class. Created by
%                          vcd_defineUniqueImages.m.
%  n_trials_per_block   : (int) number of trials per block (4 or 8)
%  env_type             : (char) character string defining if this is for 
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
%   1: {'unique_stim_nr'  } unique images (parafoveal stimulus patches: 1 trial occupies 2 rows)
%   2: {'stimloc'         } stimulus location relative to fixation dot: 1 = left, 2 = right, 3 = central
%   3: {'stimloc_name'    } same as column two but then in text 
%   4: {'orient_dir'      } orientation (gabor), motion direction (rdk), or angle (dot), or facing direction (obj) in deg
%   5: {'contrast'        } stimulus contrast (Michelson fraction)
%   6: {'gbr_phase'       } gabor stimulus phase  (NaN for non gabor stim)
%   7: {'rdk_coherence'   } rdk stimulus dot coherence (fraction of dots, NaN for non rdk stim)
%   8: {'super_cat'       } superordinate object category level (for obj and ns) 
%   9: {'basic_cat'       } basic object category level (for obj and ns)
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
%   19: {'is_catch'       } is this a catch trial where no stimulus was shown in trial (1) or not (0)?
%   20: {'is_objectcatch' } is this a objectcatch trial where a randomly 
%                           chosen non-core rotation was shown (1) or not 
%                           (0). Object catch trials only occur for OBJ stimuli.
%                           stimuli in PC-task crossing. Condition master
%                           only preallocates NaN at this stage for each
%                           trial because we determine the nr of object
%                           catch trials (20%) across the entire session.
%   21: {'unique_trial_nr'} unique trial nr
%   22: {'stim2_delta'    } for  WM task: what delta gabor tilt/
%                           rdk motion dir / dot angle / object rotation /
%                           scene was used for test stim in second stimulus array (after the delay) (see
%                           p.stim.(stim_class).delta_from_ref).
%                           for LTM task: is this A-B pair a match (1) or
%                           non-match (0).
%                           for IMG task: filled with NaN
%   23: {'stim2_im_nr'    } for  WM/IMG/LTM task: unique test image nr in stimulus array 2
%                           (after the delay). Defined for double epoch
%                           tasks only.
%   24: {'stim2_orient_dir'} for  WM/IMG/LTM task: absolute gabor tilt/rdk motion direction/dot
%                           angle/obj rotation of test stimulus (so after
%                           delta is applied to core stimulus feature)
%   25: {'is_lure'        } for LTM task: is stim2_im_nr a lure stimulus (1) or
%                           not (0). If lure stimulus, then stim2_delta must 
%                           be 0 (i.e., incorrect match). 
%                           For non LTM task classes, cols are filled with NaN.
%   26: {'repeat_nr'      } Keep track how many times has this unique image has
%                           been repeated thusfar in the experiment
%
% Written by Eline Kupers Feb 2025 @ UMN


%% Get variables from input table


if (nargin < 4 && ~exist('env_type','var')) || isempty(env_type)
    % Infer environment if we haven't env_type
    if ~exist('env_type','var') || isempty(env_type)
        if strcmp(params.disp.name,'7TAS_BOLDSCREEN32')
            env_type = 'MRI';
        elseif ismember(params.disp.name,{'CCNYU_VIEWPIXX3D','PPROOM_EIZOFLEXSCAN','KKOFFICE_AOCQ3277','EKHOME_ASUSVE247'})
            env_type = 'BEHAVIOR';
        end
    end
end
    
% which stimulus and task classes are we dealing with (cond_table is
% created by "vcd_defineUniqueImages.m"
stimClass = unique(cond_table.stim_class_name);
taskClass = unique(cond_table.task_class_name);
n_unique_cases = length(unique(cond_table.unique_im_nr));

% Check if there is only one stimulus class in provided 
assert(length(stimClass)==1)
assert(length(taskClass)==1)

% Prepare cond_master
cond_master = [];

% get nr of repetitions (depends on session type: mri or behavior)
[~,~,~,~, ~,~, ~, ~, ~, ~, unique_trial_repeats, ~,catch_trial_flag] = vcd_getSessionEnvironmentParams(params, env_type);
nr_reps = ceil(unique_trial_repeats(strcmp(stimClass, params.exp.stimclassnames),strcmp(taskClass, params.exp.taskclassnames)));

% only insert one catch trial every other block, given that special core is 50% or less stimuli per class
if catch_trial_flag && (strcmp(taskClass{:},'ltm') || strcmp(taskClass{:},'img'))
    catch_trial_occurance = shuffle_concat([0,1],ceil(nr_reps));
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

                    % SHUFFLE ORDER OF UNIQUE IMAGES
                    % we shuffle all conditions, without constraints.
                    shuffle_c      = randperm(size(cond_table,1),size(cond_table,1));
                    conds_shuffle0 = cond_table(shuffle_c,:);

                    % Add catch trial vector
                    conds_shuffle0.is_catch = zeros(size(conds_shuffle0,1),1);

                    % get covert spatial attention cuing direction vector
                    [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, stimloc_cues, taskClass);

                    % Add cue vec column
                    conds_shuffle0.is_cued = cue_vec;

                    % Preallocate space for is_objectcatch trial vector 
                    % (randomly chosen non-core rotation, only for OBJ)
                    conds_shuffle0.is_objectcatch = NaN(size(conds_shuffle0,1),1);
                    
                    % Do some checks on unique trial condition master table:
                    assert(isequal(length(conds_shuffle0.unique_im_nr),n_unique_cases)); % Check unique image nr
                    assert(isequal(sum(conds_shuffle0.stimloc==1),size(conds_shuffle0,1)/2));  % Check left stim loc
                    assert(isequal(sum(conds_shuffle0.stimloc==2),size(conds_shuffle0,1)/2));  % Check right stim loc
                    assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),repelem(params.stim.gabor.ori_deg,n_contrasts)')); % Check orientations
                    assert(isequal(sort(conds_shuffle0.gbr_phase,'ascend'),repelem(params.stim.gabor.ph_deg,n_unique_cases/length(params.stim.gabor.ph_deg))')); % Check phase
                    if ~ismember(taskClass{:}, {'ltm','img'})
                        assert(isequal(sort(conds_shuffle0.contrast,'ascend'),repelem(params.stim.gabor.contrast,n_ori_bins)')); % Check contrast levels
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
                conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(params, conds_single_rep);

                clear conds_master_single_rep

                % Add WM change
                if strcmp(taskClass{:},'wm')
                    
                    % create offset vector
                    n_deltas            = length(params.stim.gabor.delta_from_ref);
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    delta_vec           = NaN(size(conds_single_rep_merged.orient_dir));
                    delta_vec0          = NaN(size(conds_single_rep_merged.orient_dir(noncatch_trials_idx,:)));
                    orient_dir2         = NaN(size(conds_single_rep_merged.orient_dir));
                    
                    delta_vec0(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1,1)    = shuffle_concat(repelem(params.stim.gabor.delta_from_ref, n_unique_cases/2/n_deltas),1); % left cued
                    delta_vec0(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1),1) = shuffle_concat(repelem(params.stim.gabor.delta_from_ref, n_unique_cases/2/n_deltas),1); % left uncued
                    delta_vec0(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2,2)    = shuffle_concat(repelem(params.stim.gabor.delta_from_ref, n_unique_cases/2/n_deltas),1); % right cued
                    delta_vec0(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2),2) = shuffle_concat(repelem(params.stim.gabor.delta_from_ref, n_unique_cases/2/n_deltas),1); % right uncued
                    
                    delta_vec(noncatch_trials_idx,:) = delta_vec0;
                    % calculate absolute orientation of stim 2
                    orient_dir2(noncatch_trials_idx,:) = conds_single_rep_merged.orient_dir(noncatch_trials_idx,:) + delta_vec0; 
                    
                    % find unique WM test im nr
                    idx_wm      = reshape(params.stim.gabor.unique_im_nrs_wm_test,4,[]);
                    [~,idx_l]   = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.gabor.unique_im_nrs_core);
                    [~,idx_r]   = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.gabor.unique_im_nrs_core);
                    assert(isequal(length(idx_l),length(idx_r)))
                    
                    stim2_im_nr = NaN(size(conds_single_rep_merged,1),2);
                    [~,shuffle_delta] = ismember(delta_vec0,params.stim.gabor.delta_from_ref); clear delta_vec0;
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
                     conds_single_rep_merged.stim2_orient_dir  = mod(orient_dir2,360);
                     conds_single_rep_merged.is_lure           = NaN(size(conds_single_rep_merged,1),2);  
                
                elseif strcmp(taskClass{:},'img')
                    % Add IMG text prompt and quiz dots
                    n_quiz_images       = length(unique(params.stim.gabor.imagery_quiz_images)); % 2 types: overlap or not.
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    shuffle_delta       = [shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images)); ... % 1 = overlap, 2 = no overlap
                                            shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images))]'; % 50/50 chance of seeing 1 or 2
                    
                    % create yes/no overlap vector:
                    % 1=yes overlap, 2=no overlap (50% chance)
                    % col 1 = left stim, col 2 = right stim
                    % left/right stim is not yoked.
                    img_vec                        = NaN(size(conds_single_rep_merged,1),2);
                    img_vec(noncatch_trials_idx,:) = shuffle_delta;
                    idx_l                          = NaN(size(conds_single_rep_merged,1),2);
                    idx_r                          = NaN(size(conds_single_rep_merged,1),2);

                    % find unique IMG test stim nrs that correspond to the special core stim nr in stim1                  
                    [~,idx_l_tmp]   = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.gabor.unique_im_nrs_specialcore);
                    [~,idx_r_tmp]   = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.gabor.unique_im_nrs_specialcore);
                    idx_img     = params.stim.gabor.imagery_quiz_dot_stim_nr; % 8x20 matrix:8 special core stim x 20 (10 yes + 10 no overlap) test images per unique image 
                    idx_img_yes = idx_img(:,params.stim.gabor.imagery_quiz_images==1); % nr stim x 10 yes overlap test images per unique image
                    idx_img_no  = idx_img(:,params.stim.gabor.imagery_quiz_images==2); % nr stim x 10 no overlap test images per unique image
                    assert(all(ismember(unique(params.stim.gabor.unique_im_nrs_specialcore([idx_l_tmp;idx_r_tmp])), params.stim.gabor.imagery_quiz_dot_specialcore_stim_nr)))
                    assert(isequal(size(idx_img_yes),size(idx_img_no)));
                    
                    idx_l(noncatch_trials_idx==1) = idx_l_tmp; clear idx_l_tmp
                    idx_r(noncatch_trials_idx==1) = idx_r_tmp; clear idx_r_tmp
                    
                    % To avoid that the table becomes massive, we provide
                    % one number that tells us something about two quiz dots:
                    % we populate the stim2_orient column with the "orientation" 
                    % of the imagery quiz dots test images, which is the
                    % the orientation of the line connecting the two
                    % imagery quiz dots dots, relative to 0 deg = 12 o'clock. 
                    % if quiz dots "overlap" they have the same orientation
                    % as the gabor tilt.
                    ori_img_yes = params.stim.gabor.imagery_quiz_dot_orient_deg(:,params.stim.gabor.imagery_quiz_images==1,:); % 8 orientations x 10 yes overlap quiz images x 2 dots
                    ori_img_no  = params.stim.gabor.imagery_quiz_dot_orient_deg(:,params.stim.gabor.imagery_quiz_images==2,:); % 8 orientations x 10 no overlap quiz images x 2 dots
                    
                    assert(isequal(sum(circulardiff(ori_img_yes(:,:,1),ori_img_yes(:,:,2),180),'all'),0)); % dots fall on a line
                    assert(isequal(sum(circulardiff(ori_img_no(:,:,1),ori_img_no(:,:,2),180),'all'),0)); % dots fall on a line
                    
                    qz_im           = NaN(size(img_vec)); % note: with catch trials
                    ori_dir         = NaN(size(img_vec)); % note: with catch trials
                    yes_trial_idx_l = false(size(img_vec,1),1);
                    no_trial_idx_l  = false(size(img_vec,1),1);
                    yes_trial_idx_r = false(size(img_vec,1),1); 
                    no_trial_idx_r  = false(size(img_vec,1),1);
                    
                    % convert to 4 logical vectors, for l/r stim x yes/no overlap
                    yes_trial_idx_l(img_vec(:,1)==1)=true; % Quiz dots overlap (1=yes) or not (2=no) left
                    no_trial_idx_l(img_vec(:,1)==2)=true; % Quiz dots overlap (1=yes) or not (2=no) left
                    yes_trial_idx_r(img_vec(:,2)==1)=true; % Quiz dots overlap (1=yes) or not (2=no) right 
                    no_trial_idx_r(img_vec(:,2)==2)=true; % Quiz dots overlap (1=yes) or not (2=no) right
                    
                    % randomly sample from the 10 exemplar quiz dot images (WITHOUT replacement)
                    randomly_selected_yes_test_images_l  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_l), false);  % left stim yes 
                    randomly_selected_no_test_images_l   = randsample(size(idx_img_no,2),sum(no_trial_idx_l), false); % left stim no 
                    randomly_selected_yes_test_images_r  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_r), false); % right stim yes 
                    randomly_selected_no_test_images_r   = randsample(size(idx_img_no,2),sum(no_trial_idx_r), false); % right stim no
                    
                    yes_counter = 1; no_counter = 1;
                    for ll = 1:length(idx_l) % loop over unique images on the left
                        if conds_single_rep_merged.is_cued(ll)==1
                            
                            if yes_trial_idx_l(ll)==1
                                qz_im(ll,1)   = idx_img_yes(idx_l(ll),randomly_selected_yes_test_images_l(yes_counter));
                                ori_dir(ll,1) = NaN; % squeeze(ori_img_yes(idx_l(ll),randomly_selected_yes_test_images_l(yes_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left is second angle
                                yes_counter   = yes_counter+1;
                            elseif no_trial_idx_l(ll)==1
                                qz_im(ll,1)   = idx_img_no(idx_l(ll),randomly_selected_no_test_images_l(no_counter));
                                ori_dir(ll,1) = NaN; % squeeze(ori_img_no(idx_l(ll),randomly_selected_no_test_images_l(no_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left is second angle
                                no_counter    = no_counter+1;
                            elseif isnan(img_vec(ll,1))
                                qz_im(ll,1)   = NaN;
                                ori_dir(ll,1) = NaN;
                            else
                                error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                            end
                        else
                            qz_im(ll,1)   = NaN;
                            ori_dir(ll,1) = NaN;
                            % remove uncued right side from table
                            conds_single_rep_merged.stim_nr_right(ll) = NaN;
                            conds_single_rep_merged.stim_class_name(ll,2) = {NaN};
                            conds_single_rep_merged.orient_dir(ll,2) = NaN;
                            conds_single_rep_merged.contrast(ll,2) = NaN;
                            conds_single_rep_merged.gbr_phase(ll,2) = NaN;
                            conds_single_rep_merged.is_special_core(ll,2) = NaN;
                        end
                    end
                    
                    yes_counter = 1; no_counter = 1;
                    for ll = 1:length(idx_r) % loop over unique images on the right
                        if conds_single_rep_merged.is_cued(ll)==2
                            if yes_trial_idx_r(ll)==1
                                qz_im(ll,2)   = idx_img_yes(idx_r(ll),randomly_selected_yes_test_images_r(yes_counter));
                                ori_dir(ll,2) = NaN; % squeeze(ori_img_yes(idx_r(ll),randomly_selected_yes_test_images_r(yes_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note right is first angle
                                yes_counter   = yes_counter+1;
                            elseif no_trial_idx_r(ll)==1
                                qz_im(ll,2)   = idx_img_no(idx_r(ll),randomly_selected_no_test_images_r(no_counter));
                                ori_dir(ll,2) = NaN; % squeeze(ori_img_no(idx_r(ll),randomly_selected_no_test_images_r(no_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note right is first angle
                                no_counter = no_counter+1;
                            elseif isnan(img_vec(ll,1))
                                qz_im(ll,2)   = NaN;
                                ori_dir(ll,2) = NaN;
                            else
                                error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                            end
                        else
                            qz_im(ll,2)   = NaN;
                            ori_dir(ll,2) = NaN;
                            % remove uncued left side from table
                            conds_single_rep_merged.stim_nr_left(ll) = NaN;
                            conds_single_rep_merged.stim_class_name(ll,1) = {NaN};
                            conds_single_rep_merged.orient_dir(ll,1) = NaN;
                            conds_single_rep_merged.contrast(ll,1) = NaN;
                            conds_single_rep_merged.gbr_phase(ll,1) = NaN;
                            conds_single_rep_merged.is_special_core(ll,1) = NaN;
                        end
                    end
                    
                    % all aligned quiz dots should match the tilt of the gabor shown in stim1
                    %assert(all(isequal(ori_dir(yes_trial_idx_l,1),conds_single_rep_merged.orient_dir(yes_trial_idx_l,1)))); % left side
                    %assert(all(ismember(ori_dir(yes_trial_idx_r,2),conds_single_rep_merged.orient_dir(yes_trial_idx_r,2)))); % right side

                    assert(length(unique(qz_im(:)))==length(qz_im(:)))
                    assert(all(ismember(qz_im(yes_trial_idx_l,1),idx_img_yes))); % check if yes quiz dots image nr matches for left stim
                    assert(all(ismember(qz_im(no_trial_idx_l,1),idx_img_no))); % check if no quiz dots image nr matches for left stim
                    assert(all(ismember(qz_im(yes_trial_idx_r,2),idx_img_yes))); % check if yes quiz dots image nr matches for right stim
                    assert(all(ismember(qz_im(no_trial_idx_r,2),idx_img_no))); % check if no quiz dots image nr matches for right stim
                    
                    assert(isequal(sum(img_vec(:,1)==1),sum(img_vec(:,1)==2))); % right stim side: equal yes/no overlap?
                    assert(isequal(sum(img_vec(:,2)==1),sum(img_vec(:,2)==2))); % left stim side: equal yes/no overlap?
                    
                    % fill in vectors & condition table
                    conds_single_rep_merged.stim2_delta      = img_vec; % yes (1) /no (2) dots are aligned
                    conds_single_rep_merged.stim2_im_nr      = qz_im;   % unique img quiz dot stim number
                    conds_single_rep_merged.stim2_orient_dir = ori_dir; % orientation of dots
                    
                    % NO lures in IMG, set to NaN
                    conds_single_rep_merged.is_lure          = NaN(size(conds_single_rep_merged,1),2);
                    
                % Add LTM pairs
                elseif strcmp(taskClass{:}, 'ltm')
                    % preallocate space
                    conds_single_rep_merged.stim2_im_nr       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure           = NaN(size(conds_single_rep_merged,1),2);
                    
                    % define non-catch trials
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    
                    [stim2_nr, stim2_is_lure, stim2_match] = vcd_getLTMTestStim(params, conds_single_rep_merged);
                    assert(isequal(isnan(stim2_nr(:,1)),~noncatch_trials_idx))
                    assert(isequal(isnan(stim2_nr(:,2)),~noncatch_trials_idx))
                    
                    conds_single_rep_merged.stim2_im_nr = stim2_nr;
                    conds_single_rep_merged.is_lure     = stim2_is_lure;
                    conds_single_rep_merged.stim2_delta = stim2_match;
                    
                    % add orientation info
                    tmp1 = vcd('fullinfo', stim2_nr(noncatch_trials_idx,1));
                    tmp2 = vcd('fullinfo', stim2_nr(noncatch_trials_idx,2));
                    conds_single_rep_merged.stim2_orient_dir(noncatch_trials_idx,:) = cat(2, cat(1, tmp1.orient_dir), cat(1, tmp2.orient_dir));
                    clear tmp1 tmp2 
                    
                    % We want to test both A->B pairings and B->A pairings.
                    % Make a copy of conditions table
                    conds_single_rep_mergedB = conds_single_rep_merged;
                    
                    % flip order A->B to B->A 
                    conds_single_rep_mergedB.stim_nr_left  = conds_single_rep_merged.stim2_im_nr(:,1);
                    conds_single_rep_mergedB.stim_nr_right = conds_single_rep_merged.stim2_im_nr(:,2);
                    
                    % reset gabor phase, orient_dir (we'll add if needed)
                    conds_single_rep_mergedB.gbr_phase     = NaN(size(conds_single_rep_mergedB.gbr_phase));
                    conds_single_rep_mergedB.orient_dir    = NaN(size(conds_single_rep_mergedB.orient_dir));
                    conds_single_rep_mergedB.contrast      = NaN(size(conds_single_rep_mergedB.contrast));
                    
                     
                    % update stim2
                    conds_single_rep_mergedB.stim2_im_nr      = [conds_single_rep_merged.stim_nr_left, conds_single_rep_merged.stim_nr_right];
                    conds_single_rep_mergedB.stim2_orient_dir = conds_single_rep_merged.orient_dir;
                    
                    % update gabor phase info (if lure stim)
                    for ii = 1:size(conds_single_rep_mergedB,1)
                        % Get stim nr for left stim loc
                        tmp_nr_l = conds_single_rep_mergedB.stim_nr_left(ii);
                        
                        if isnan(tmp_nr_l) || (tmp_nr_l==0) % catch trials
                            % keep as "gabor"
                            assert(isequal(conds_single_rep_mergedB.is_catch(ii),1))
                            % set stim1 nr to 0000 and stim2 nr to NaN
                            conds_single_rep_mergedB.stim_nr_left(ii)  = 0;
                            conds_single_rep_mergedB.stim2_im_nr(ii,1) = NaN;
                        else 
                            % get stiminfo
                            tmp = vcd('fullinfo', tmp_nr_l);
                            assert(isequal(tmp.stim_loc,1)); % should be left
                            
                            % update stimclass name, orient_dir and contrast
                            conds_single_rep_mergedB.stim_class_name{ii,1} = vcd('stimtostimclassname', tmp_nr_l);
                            conds_single_rep_mergedB.orient_dir(ii,1)      = tmp.orient_dir;
                            conds_single_rep_mergedB.contrast(ii,1)        = tmp.contrast;
                            
                            if ismember(tmp_nr_l, params.stim.gabor.unique_im_nrs_core) % if gabor
                                assert(isequal(tmp.stim_class,1))
                                % add gabor stim info
                                conds_single_rep_mergedB.gbr_phase(ii,1)  = tmp.gbr_phase;
                            elseif ismember(tmp_nr_l, params.stim.rdk.unique_im_nrs_core) % if rdk
                                assert(isequal(tmp.stim_class,2)) % should be rdk
                                conds_single_rep_mergedB.rdk_coherence(ii,1) = tmp.rdk_coherence; % update motion coherence
                            elseif ismember(tmp_nr_l, params.stim.dot.unique_im_nrs_core) % if dot
                                assert(isequal(tmp.stim_class,3)) % should be dot
                            elseif ismember(tmp_nr_l, params.stim.obj.unique_im_nrs_core) % if obj
                                assert(isequal(tmp.stim_class,4)); 
                                % update onbject categories
                                conds_single_rep_mergedB.super_cat(ii,1)        = tmp.super_cat;
                                conds_single_rep_mergedB.basic_cat(ii,1)        = tmp.basic_cat;
                                conds_single_rep_mergedB.sub_cat(ii,1)          = tmp.sub_cat;
                                conds_single_rep_mergedB.affordance_cat(ii,1)   = tmp.affordance_cat;
                                conds_single_rep_mergedB.super_cat_name(ii,1)   = tmp.super_cat_name;
                                conds_single_rep_mergedB.basic_cat_name(ii,1)   = tmp.basic_cat_name;
                                conds_single_rep_mergedB.sub_cat_name(ii,1)     = tmp.sub_cat_name;
                                conds_single_rep_mergedB.affordance_name(ii,1)  = tmp.affordance_name;
                            end
                        end
                    end
                    
                    for ii = 1:size(conds_single_rep_mergedB,1)
                        % Get stim nr for right stim loc
                        tmp_nr_r = conds_single_rep_mergedB.stim_nr_right(ii);
                        
                        if isnan(tmp_nr_r) || (tmp_nr_r==0) % catch trials
                            % keep as "gabor"
                            % set stim1 nr to 0000 and stim2 nr to NaN
                            conds_single_rep_mergedB.stim_nr_right(ii)  = 0;
                            conds_single_rep_mergedB.stim2_im_nr(ii,2) = NaN;
                        else
                            % get stim info
                            tmp = vcd('fullinfo', tmp_nr_r);
                            assert(isequal(tmp.stim_loc,2)); % must be right stim loc
                            % update stim class name
                            conds_single_rep_mergedB.stim_class_name{ii,2} = vcd('stimtostimclassname', tmp_nr_r);
                            conds_single_rep_mergedB.orient_dir(ii,2)      = tmp.orient_dir; % update tilt/mot dir/facing dir info
                            conds_single_rep_mergedB.contrast(ii,2)        = tmp.contrast; % update contrast
                            if ismember(tmp_nr_r, params.stim.gabor.unique_im_nrs_core) % if gabor
                                assert(isequal(tmp.stim_class,1))
                                conds_single_rep_mergedB.gbr_phase(ii,2)  = tmp.gbr_phase;
                            elseif ismember(tmp_nr_r, params.stim.rdk.unique_im_nrs_core) % if rdk  
                                assert(isequal(tmp.stim_class,2)) % must be rdk
                                conds_single_rep_mergedB.rdk_coherence(ii,2) = tmp.rdk_coherence; % add rdk motion coherence
                            elseif ismember(tmp_nr_r, params.stim.dot.unique_im_nrs_core) % if dot
                                assert(isequal(tmp.stim_class,3)) % should be dot    
                            elseif ismember(tmp_nr_r, params.stim.obj.unique_im_nrs_core) % if obj
                                % add object category info
                                conds_single_rep_mergedB.super_cat(ii,2)        = tmp.super_cat;
                                conds_single_rep_mergedB.basic_cat(ii,2)        = tmp.basic_cat;
                                conds_single_rep_mergedB.sub_cat(ii,2)          = tmp.sub_cat;
                                conds_single_rep_mergedB.affordance_cat(ii,2)   = tmp.affordance_cat;
                                conds_single_rep_mergedB.super_cat_name(ii,2)   = tmp.super_cat_name;
                                conds_single_rep_mergedB.basic_cat_name(ii,2)   = tmp.basic_cat_name;
                                conds_single_rep_mergedB.sub_cat_name(ii,2)     = tmp.sub_cat_name;
                                conds_single_rep_mergedB.affordance_name(ii,2)  = tmp.affordance_name;
                            end
                        end
                    end

                    % update unique_trial_nr
                    conds_single_rep_mergedB.unique_trial_nr = max(conds_single_rep_mergedB.unique_trial_nr) + conds_single_rep_mergedB.unique_trial_nr;
                    
                    % check unique trial number (should continue counting)
                    assert(isequal(conds_single_rep_mergedB.unique_trial_nr, ...
                        [min(conds_single_rep_mergedB.unique_trial_nr):max(conds_single_rep_mergedB.unique_trial_nr)]'));
                    
                    % merge condition tables
                    conds_single_rep_merged = cat(1,conds_single_rep_merged,conds_single_rep_mergedB);
                    
                     % check unique trial number (should continue counting)
                    assert(isequal(sum(~noncatch_trials_idx)*2, sum(conds_single_rep_merged.is_catch)));
                    assert(isequal(conds_single_rep_merged.unique_trial_nr, ...
                                    [min(conds_single_rep_merged.unique_trial_nr):max(conds_single_rep_merged.unique_trial_nr)]'));
                    
                else
                    conds_single_rep_merged.stim2_delta      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_im_nr      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure          = NaN(size(conds_single_rep_merged,1),2); 
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
                    
                    % SHUFFLE ORDER OF UNIQUE IMAGES
                    % we shuffle all conditions, without constraints.
                    shuffle_c      = randperm(size(cond_table,1),size(cond_table,1));
                    conds_shuffle0 = cond_table(shuffle_c,:);
                    
                    % get covert spatial attention cuing direction vector
                    [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, stimloc_cues, taskClass);

                    % Add cue vec column
                    conds_shuffle0.is_cued = cue_vec;

                    % Add catch trial vector
                    conds_shuffle0.is_catch = zeros(size(conds_shuffle0,1),1);

                    % Preallocate space for is_objectcatch trial vector
                    % (randomly chosen non-core rotation, only for OBJ)
                    conds_shuffle0.is_objectcatch = NaN(size(conds_shuffle0,1),1);
                    
                    % Do some checks:
                    assert(isequal(length(conds_shuffle0.unique_im_nr),n_unique_cases)); % Check unique image nr
                    assert(isequal(sum(conds_shuffle0.stimloc==1),n_unique_cases/2));  % Check left stim loc
                    assert(isequal(sum(conds_shuffle0.stimloc==2),n_unique_cases/2));  % Check right stim loc
                    assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),repelem(params.stim.rdk.dots_direction,n_coh)')); % Check motdir
                    if ~ismember(taskClass{:}, {'ltm','img'})
                        assert(isequal(sort(conds_shuffle0.rdk_coherence,'ascend'),repelem(params.stim.rdk.dots_coherence,n_motdir)')); % Check coherence levels
                    else
                        assert(isequal(conds_shuffle0.rdk_coherence,repelem(params.stim.rdk.dots_coherence(3),n_motdir)')); % Check coherence
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
                conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(params, conds_master_single_rep);
                clear conds_master_single_rep

                % Add WM change.
                if strcmp(taskClass{:},'wm')
                    n_deltas            = length(params.stim.rdk.delta_from_ref);
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    delta_vec           = NaN(size(conds_single_rep_merged.orient_dir));
                    delta_vec0          = NaN(size(conds_single_rep_merged.orient_dir(noncatch_trials_idx,:)));
                    orient_dir2         = NaN(size(conds_single_rep_merged.orient_dir));
                    
                    delta_vec0(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1,1) = shuffle_concat(repelem(params.stim.rdk.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left cued
                    delta_vec0(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1),1) = shuffle_concat(repelem(params.stim.rdk.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left uncued
                    delta_vec0(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2,2) = shuffle_concat(repelem(params.stim.rdk.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right cued
                    delta_vec0(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2),2) = shuffle_concat(repelem(params.stim.rdk.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right uncued
                    delta_vec(noncatch_trials_idx,:) = delta_vec0; 
                    
                    % Calculate absolute orientation of stim 2
                    orient_dir2(noncatch_trials_idx,:)   = conds_single_rep_merged.orient_dir(noncatch_trials_idx,:) + delta_vec0;

                    % find unique WM test im nr
                    idx_wm = reshape(params.stim.rdk.unique_im_nrs_wm_test,4,[]);
                    [~,idx_l] = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.rdk.unique_im_nrs_core);
                    [~,idx_r] = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.rdk.unique_im_nrs_core);
                    assert(isequal(length(idx_l),length(idx_r)))
                    
                    stim2_im_nr = NaN(size(conds_single_rep_merged,1),2);
                    [~,shuffle_delta] = ismember(delta_vec0,params.stim.rdk.delta_from_ref); clear delta_vec0;
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
                    conds_single_rep_merged.stim2_orient_dir  = mod(orient_dir2,360); % apply circular wrap
                    conds_single_rep_merged.is_lure           = NaN(size(conds_single_rep_merged,1),2);


                % Add IMG text prompt and quiz dots
                elseif strcmp(taskClass{:},'img')
                    % Add IMG text prompt and quiz dots
                    n_quiz_images       = length(unique(params.stim.rdk.imagery_quiz_images)); % 2 types: overlap or not.
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    shuffle_delta       = [shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images)); ... % 1 = overlap, 2 = no overlap
                                            shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images))]'; % 50/50 chance of seeing 1 or 2
                    
                    % create yes/no overlap vector:
                    % 1=yes overlap, 2=no overlap (50% chance)
                    % col 1 = left stim, col 2 = right stim
                    % left/right stim is not yoked.
                    img_vec                        = NaN(size(conds_single_rep_merged,1),2);
                    img_vec(noncatch_trials_idx,:) = shuffle_delta;
                    idx_l                          = NaN(size(conds_single_rep_merged,1),2);
                    idx_r                          = NaN(size(conds_single_rep_merged,1),2);
                    
                    % find unique IMG test stim nrs that correspond to the special core stim nr in stim1                  
                    [~,idx_l_tmp]   = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.rdk.unique_im_nrs_specialcore);
                    [~,idx_r_tmp]   = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.rdk.unique_im_nrs_specialcore);
                    idx_img     = params.stim.rdk.imagery_quiz_dot_stim_nr; % 8x20 matrix:8 special core stim x 20 (10 yes + 10 no overlap) test images per unique image 
                    idx_img_yes = idx_img(:,params.stim.rdk.imagery_quiz_images==1); % nr stim x 10 yes overlap test images per unique image
                    idx_img_no  = idx_img(:,params.stim.rdk.imagery_quiz_images==2); % nr stim x 10 no overlap test images per unique image
                    assert(all(ismember(unique(params.stim.rdk.unique_im_nrs_specialcore([idx_l_tmp;idx_r_tmp])), params.stim.rdk.imagery_quiz_dot_specialcore_stim_nr)))
                    assert(isequal(size(idx_img_yes),size(idx_img_no)));
                    
                    idx_l(noncatch_trials_idx==1) = idx_l_tmp; clear idx_l_tmp
                    idx_r(noncatch_trials_idx==1) = idx_r_tmp; clear idx_r_tmp
                    
                    % To avoid that the table becomes massive, we provide
                    % one number that tells us something about two quiz dots:
                    % we populate the stim2_orient column with the "orientation" 
                    % of the imagery quiz dots test images, which is the
                    % the orientation of the line connecting the two
                    % imagery quiz dots dots, relative to 0 deg = 12 o'clock. 
                    % if quiz dots "overlap" they have the same orientation
                    % as the rdk motion direction.
                    ori_img_yes = params.stim.rdk.imagery_quiz_dot_orient_deg(:,params.stim.rdk.imagery_quiz_images==1,:); % 8 orientations x 10 yes overlap quiz images x 2 dots
                    ori_img_no  = params.stim.rdk.imagery_quiz_dot_orient_deg(:,params.stim.rdk.imagery_quiz_images==2,:); % 8 orientations x 10 no overlap quiz images x 2 dots
                    
                    assert(isequal(sum(circulardiff(ori_img_yes(:,:,1),ori_img_yes(:,:,2),180),'all'),0)); % dots fall on a line
                    assert(isequal(sum(circulardiff(ori_img_no(:,:,1),ori_img_no(:,:,2),180),'all'),0)); % dots fall on a line
                    
                    qz_im           = NaN(size(img_vec)); % note: with catch trials
                    ori_dir         = NaN(size(img_vec)); % note: with catch trials
                    yes_trial_idx_l = false(size(img_vec,1),1);
                    no_trial_idx_l  = false(size(img_vec,1),1);
                    yes_trial_idx_r = false(size(img_vec,1),1); 
                    no_trial_idx_r  = false(size(img_vec,1),1);
                    
                    % convert to 4 logical vectors, for l/r stim x yes/no overlap
                    yes_trial_idx_l(img_vec(:,1)==1)=true; % Quiz dots overlap (1=yes) or not (2=no) left
                    no_trial_idx_l(img_vec(:,1)==2)=true; % Quiz dots overlap (1=yes) or not (2=no) left
                    yes_trial_idx_r(img_vec(:,2)==1)=true; % Quiz dots overlap (1=yes) or not (2=no) right 
                    no_trial_idx_r(img_vec(:,2)==2)=true; % Quiz dots overlap (1=yes) or not (2=no) right
                    
                    % randomly sample from the 10 exemplar quiz dot images (WITHOUT replacement)
                    randomly_selected_yes_test_images_l  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_l), false);  % left stim yes 
                    randomly_selected_no_test_images_l   = randsample(size(idx_img_no,2),sum(no_trial_idx_l), false); % left stim no 
                    randomly_selected_yes_test_images_r  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_r), false); % right stim yes 
                    randomly_selected_no_test_images_r   = randsample(size(idx_img_no,2),sum(no_trial_idx_r), false); % right stim no
                    
                    yes_counter = 1; no_counter = 1;
                    for ll = 1:length(idx_l) % loop over unique images on the left
                        if conds_single_rep_merged.is_cued(ll)==1
                            if yes_trial_idx_l(ll)==1
                                qz_im(ll,1)   = idx_img_yes(idx_l(ll),randomly_selected_yes_test_images_l(yes_counter));
                                ori_dir(ll,1) = NaN; % squeeze(ori_img_yes(idx_l(ll),randomly_selected_yes_test_images_l(yes_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left is second angle
                                yes_counter   = yes_counter+1;
                            elseif no_trial_idx_l(ll)==1
                                qz_im(ll,1)   = idx_img_no(idx_l(ll),randomly_selected_no_test_images_l(no_counter));
                                ori_dir(ll,1) = NaN; % squeeze(ori_img_no(idx_l(ll),randomly_selected_no_test_images_l(no_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left is second angle
                                no_counter    = no_counter+1;
                            elseif isnan(img_vec(ll,1))
                                qz_im(ll,1)   = NaN;
                                ori_dir(ll,1) = NaN;
                            else
                                error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                            end
                        else
                            qz_im(ll,1)   = NaN;
                            ori_dir(ll,1) = NaN;
                            % remove uncued right side from table
                            conds_single_rep_merged.stim_nr_right(ll) = NaN;
                            conds_single_rep_merged.stim_class_name(ll,2) = {NaN};
                            conds_single_rep_merged.orient_dir(ll,2) = NaN;
                            conds_single_rep_merged.contrast(ll,2) = NaN;
                            conds_single_rep_merged.rdk_coherence(ll,2) = NaN;
                            conds_single_rep_merged.is_special_core(ll,2) = NaN;
                        end
                    end
                    
                    yes_counter = 1; no_counter = 1;
                    for ll = 1:length(idx_r) % loop over unique images on the right
                        if conds_single_rep_merged.is_cued(ll)==2
                            if yes_trial_idx_r(ll)==1
                                qz_im(ll,2)   = idx_img_yes(idx_r(ll),randomly_selected_yes_test_images_r(yes_counter));
                                ori_dir(ll,2) = NaN; % squeeze(ori_img_yes(idx_r(ll),randomly_selected_yes_test_images_r(yes_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note right is first angle
                                yes_counter   = yes_counter+1;
                            elseif no_trial_idx_r(ll)==1
                                qz_im(ll,2)   = idx_img_no(idx_r(ll),randomly_selected_no_test_images_r(no_counter));
                                ori_dir(ll,2) = NaN; % squeeze(ori_img_no(idx_r(ll),randomly_selected_no_test_images_r(no_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note right is first angle
                                no_counter = no_counter+1;
                            elseif isnan(img_vec(ll,1))
                                qz_im(ll,2)   = NaN;
                                ori_dir(ll,2) = NaN;
                            else
                                error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                            end
                        else
                            qz_im(ll,2)   = NaN;
                            ori_dir(ll,2) = NaN;
                            % remove uncued left side from table
                            conds_single_rep_merged.stim_nr_left(ll) = NaN;
                            conds_single_rep_merged.stim_class_name(ll,1) = {NaN};
                            conds_single_rep_merged.orient_dir(ll,1) = NaN;
                            conds_single_rep_merged.contrast(ll,1) = NaN;
                            conds_single_rep_merged.rdk_coherence(ll,1) = NaN;
                            conds_single_rep_merged.is_special_core(ll,1) = NaN;
                        end
                    end
                    
                    % all aligned quiz dots should match the tilt of the rdk motion direction shown in stim1
                    %assert(all(isequal(ori_dir(yes_trial_idx_l,1),conds_single_rep_merged.orient_dir(yes_trial_idx_l,1)))); % left side
                    %assert(all(ismember(ori_dir(yes_trial_idx_r,2),conds_single_rep_merged.orient_dir(yes_trial_idx_r,2)))); % right side

                    assert(length(unique(qz_im(:)))==length(qz_im(:)))
                    assert(all(ismember(qz_im(yes_trial_idx_l,1),idx_img_yes))); % check if yes quiz dots image nr matches for left stim
                    assert(all(ismember(qz_im(no_trial_idx_l,1),idx_img_no))); % check if no quiz dots image nr matches for left stim
                    assert(all(ismember(qz_im(yes_trial_idx_r,2),idx_img_yes))); % check if yes quiz dots image nr matches for right stim
                    assert(all(ismember(qz_im(no_trial_idx_r,2),idx_img_no))); % check if no quiz dots image nr matches for right stim
                    
                    assert(isequal(sum(img_vec(:,1)==1),sum(img_vec(:,1)==2))); % right stim side: equal yes/no overlap?
                    assert(isequal(sum(img_vec(:,2)==1),sum(img_vec(:,2)==2))); % left stim side: equal yes/no overlap?
                    
                    % fill in vectors & condition table
                    conds_single_rep_merged.stim2_delta      = img_vec; % yes (1) /no (2) dots are aligned
                    conds_single_rep_merged.stim2_im_nr      = qz_im;   % unique img quiz dot stim number
                    conds_single_rep_merged.stim2_orient_dir = ori_dir; % orientation of dots
                    
                    % NO lures in IMG, set to NaN
                    conds_single_rep_merged.is_lure          = NaN(size(conds_single_rep_merged,1),2);

                % Add ltm pair.
                elseif strcmp(taskClass{:},'ltm')
                    conds_single_rep_merged.stim2_im_nr       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure           = NaN(size(conds_single_rep_merged,1),2);

                    [stim2_nr, stim2_is_lure, stim2_match] = vcd_getLTMTestStim(params, conds_single_rep_merged);
                    
                    % check if catch trials match
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    assert(isequal(isnan(stim2_nr(:,1)),~noncatch_trials_idx))
                    assert(isequal(isnan(stim2_nr(:,2)),~noncatch_trials_idx))
                    
                    % Add paired stim, if its a matching pair and if it is
                    % a lure when not a match
                    conds_single_rep_merged.stim2_im_nr = stim2_nr;
                    conds_single_rep_merged.is_lure     = stim2_is_lure;
                    conds_single_rep_merged.stim2_delta = stim2_match;
                    
                    % add orientation info
                    tmp1 = vcd('fullinfo', stim2_nr(noncatch_trials_idx,1));
                    tmp2 = vcd('fullinfo', stim2_nr(noncatch_trials_idx,2));
                    conds_single_rep_merged.stim2_orient_dir(noncatch_trials_idx,:) = cat(2, cat(1, tmp1.orient_dir), cat(1, tmp2.orient_dir));
                    clear tmp1 tmp2
                    
                    % We want to test both A->B pairings and B->A pairings.
                    % Make a copy of conditions table
                    conds_single_rep_mergedB = conds_single_rep_merged;
                    
                     % flip order A->B to B->A 
                    conds_single_rep_mergedB.stim_nr_left  = conds_single_rep_merged.stim2_im_nr(:,1);
                    conds_single_rep_mergedB.stim_nr_right = conds_single_rep_merged.stim2_im_nr(:,2);
                    
                    % Reset orientation, add back later
                    conds_single_rep_mergedB.orient_dir    = NaN(size(conds_single_rep_merged.orient_dir,1),2);
                    conds_single_rep_mergedB.contrast    = NaN(size(conds_single_rep_merged.orient_dir,1),2);
                    
                    % update stim2
                    conds_single_rep_mergedB.stim2_im_nr      = [conds_single_rep_merged.stim_nr_left, conds_single_rep_merged.stim_nr_right];
                    conds_single_rep_mergedB.stim2_orient_dir = conds_single_rep_merged.orient_dir;
                    
                    % update stim info 
                    for ii = 1:size(conds_single_rep_mergedB,1)
                        % Get stim nr for left stim loc
                        tmp_nr_l = conds_single_rep_mergedB.stim_nr_left(ii);
                        
                        if isnan(tmp_nr_l) || (tmp_nr_l==0) % catch trials
                            % keep as "rdk"
                            assert(isequal(conds_single_rep_mergedB.is_catch(ii),1))
                            % set stim1 nr to 0000 and stim2 nr to NaN
                            conds_single_rep_mergedB.stim_nr_left(ii) = 0;
                            conds_single_rep_mergedB.stim2_im_nr(ii,1)   = NaN;
                        else
                            % get stim info
                            tmp = vcd('fullinfo', tmp_nr_l);
                            % Update stimclass, orientation info
                            conds_single_rep_mergedB.stim_class_name{ii,1} = vcd('stimtostimclassname', tmp_nr_l); % update stim class name
                            conds_single_rep_mergedB.orient_dir(ii,1)      = tmp.orient_dir;
                            conds_single_rep_mergedB.rdk_coherence(ii,1)   = NaN; % set rdk coh to NaN (we'll add back when needed)
                            conds_single_rep_mergedB.contrast(ii,1)        = tmp.contrast; % add contrast
                            assert(isequal(tmp.stim_loc,1)); % must be left
                            
                            if ismember(tmp_nr_l, params.stim.gabor.unique_im_nrs_core) % if gabor
                                assert(isequal(tmp.stim_class,1)) % must be gabor
                                conds_single_rep_mergedB.gbr_phase(ii,1) = tmp.gbr_phase; % add gabor phase
                            elseif ismember(tmp_nr_l, params.stim.rdk.unique_im_nrs_core) % if rdk
                                assert(isequal(tmp.stim_class,2)); % must be rdk
                                conds_single_rep_mergedB.rdk_coherence(ii,1) = tmp.rdk_coherence; % add rdk motion coherence
                            elseif ismember(tmp_nr_l, params.stim.dot.unique_im_nrs_core) % if dot
                                assert(isequal(tmp.stim_class,3)); % must be dot
                            elseif ismember(tmp_nr_l, params.stim.obj.unique_im_nrs_core) % if obj
                                assert(isequal(tmp.stim_class,4)); % must be object
                                % add category info
                                conds_single_rep_mergedB.super_cat(ii,1)        = tmp.super_cat;
                                conds_single_rep_mergedB.basic_cat(ii,1)        = tmp.basic_cat;
                                conds_single_rep_mergedB.sub_cat(ii,1)          = tmp.sub_cat;
                                conds_single_rep_mergedB.affordance_cat(ii,1)   = tmp.affordance_cat;
                                conds_single_rep_mergedB.super_cat_name(ii,1)   = tmp.super_cat_name;
                                conds_single_rep_mergedB.basic_cat_name(ii,1)   = tmp.basic_cat_name;
                                conds_single_rep_mergedB.sub_cat_name(ii,1)     = tmp.sub_cat_name;
                                conds_single_rep_mergedB.affordance_name(ii,1)  = tmp.affordance_name;
                            end
                        end
                    end
                    
                    % update stim info 
                    for ii = 1:size(conds_single_rep_mergedB,1)
                        % Get stim nr for right stim loc
                        tmp_nr_r = conds_single_rep_mergedB.stim_nr_right(ii);
                        
                        if isnan(tmp_nr_r) || (tmp_nr_r==0) % catch trials
                            % keep as "rdk"
                            % set stim1 nr to 0000 and stim2 nr to NaN
                            conds_single_rep_mergedB.stim_nr_right(ii) = 0;
                            conds_single_rep_mergedB.stim2_im_nr(ii,2) = NaN;
                        else
                            tmp = vcd('fullinfo', tmp_nr_r);
                            assert(isequal(tmp.stim_loc,2)); % must be right
                            % update stimclass name, orientation
                            conds_single_rep_mergedB.stim_class_name{ii,2} = vcd('stimtostimclassname', tmp_nr_r);
                            conds_single_rep_mergedB.orient_dir(ii,2)    = tmp.orient_dir;
                            conds_single_rep_mergedB.rdk_coherence(ii,2) = NaN; % set rdk coh to NaN
                             conds_single_rep_mergedB.contrast(ii,2)  = tmp.contrast; % add contrast
                            if ismember(tmp_nr_r, params.stim.gabor.unique_im_nrs_core)  % if gabor
                                assert(isequal(tmp.stim_class,1)) % must be gabor
                                conds_single_rep_mergedB.gbr_phase(ii,2) = tmp.gbr_phase; % add gabor phase
                            elseif ismember(tmp_nr_r, params.stim.rdk.unique_im_nrs_core) % if rdk
                                assert(isequal(tmp.stim_class,2)); % must be rdk
                                conds_single_rep_mergedB.rdk_coherence(ii,2) = tmp.rdk_coherence; % add rdk motion coherence back
                            elseif ismember(tmp_nr_r, params.stim.dot.unique_im_nrs_core) % if dot
                                assert(isequal(tmp.stim_class,3)); % must be dot
                            elseif ismember(tmp_nr_r, params.stim.obj.unique_im_nrs_core) % if obj
                                assert(isequal(tmp.stim_class,4)); % must be object
                                % add object category info
                                conds_single_rep_mergedB.super_cat(ii,2)        = tmp.super_cat;
                                conds_single_rep_mergedB.basic_cat(ii,2)        = tmp.basic_cat;
                                conds_single_rep_mergedB.sub_cat(ii,2)          = tmp.sub_cat;
                                conds_single_rep_mergedB.affordance_cat(ii,2)   = tmp.affordance_cat;
                                conds_single_rep_mergedB.super_cat_name(ii,2)   = tmp.super_cat_name;
                                conds_single_rep_mergedB.basic_cat_name(ii,2)   = tmp.basic_cat_name;
                                conds_single_rep_mergedB.sub_cat_name(ii,2)     = tmp.sub_cat_name;
                                conds_single_rep_mergedB.affordance_name(ii,2)  = tmp.affordance_name;
                            end
                        end
                    end
                                                
                    % update unique_trial_nr
                    conds_single_rep_mergedB.unique_trial_nr = max(conds_single_rep_mergedB.unique_trial_nr) + conds_single_rep_mergedB.unique_trial_nr;
                    
                    % check unique trial number (should continue counting)
                    assert(isequal(conds_single_rep_mergedB.unique_trial_nr, ...
                        [min(conds_single_rep_mergedB.unique_trial_nr):max(conds_single_rep_mergedB.unique_trial_nr)]'));
                    
                    % merge condition tables
                    conds_single_rep_merged = cat(1,conds_single_rep_merged,conds_single_rep_mergedB);
                    
                    % check unique trial number (should continue counting)
                    assert(isequal(sum(~noncatch_trials_idx)*2, sum(conds_single_rep_merged.is_catch)));
                    assert(isequal(conds_single_rep_merged.unique_trial_nr, ...
                        [min(conds_single_rep_merged.unique_trial_nr):max(conds_single_rep_merged.unique_trial_nr)]'));
                else
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_im_nr       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure           = NaN(size(conds_single_rep_merged,1),2); 
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

                    % we shuffle all conditions, without constraints.
                    shuffle_c      = randperm(size(cond_table,1),size(cond_table,1));
                    conds_shuffle0 = cond_table(shuffle_c,:);
                    
                    % get covert spatial attention cuing direction vector
                    [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, stimloc_cues, taskClass);

                    % Add cue vec
                    conds_shuffle0.is_cued = cue_vec;

                    % Add catch trial vector
                    conds_shuffle0.is_catch = zeros(size(conds_shuffle0,1),1);

                    % Preallocate space for is_objectcatch trial vector
                    % (randomly chosen non-core rotation, only for OBJ)
                    conds_shuffle0.is_objectcatch = NaN(size(conds_shuffle0,1),1);
                    
                    % Do some checks:
                    assert(isequal(length(conds_shuffle0.unique_im_nr),n_unique_cases)); % Check unique image nr
                    assert(isequal(sum(conds_shuffle0.stimloc==1),n_unique_cases/2));  % Check left stim loc
                    assert(isequal(sum(conds_shuffle0.stimloc==2),n_unique_cases/2));  % Check right stim loc
                    if ~ismember(taskClass{:}, {'ltm','img'})
                        assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),sort(params.stim.dot.ang_deg)')); % Check dot angle
                    else
                        assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),sort(params.stim.dot.ang_deg(ismember(params.stim.dot.unique_im_nrs_core,params.stim.dot.unique_im_nrs_specialcore)))')); % Check contrast
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
                conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(params, conds_master_single_rep);
                clear conds_master_single_rep

                % Add WM change.
                if strcmp(taskClass{:},'wm')
                    n_deltas      = length(params.stim.dot.delta_from_ref);
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    delta_vec     = NaN(size(conds_single_rep_merged.orient_dir));
                    orient_dir2   = NaN(size(conds_single_rep_merged.orient_dir));
                    while 1
                        % create shuffled delta vector
                        delta_vec0(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1,1) = shuffle_concat(repelem(params.stim.dot.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left cued
                        delta_vec0(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1),1) = shuffle_concat(repelem(params.stim.dot.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left uncued
                        delta_vec0(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2,2) = shuffle_concat(repelem(params.stim.dot.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right cued
                        delta_vec0(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2),2) = shuffle_concat(repelem(params.stim.dot.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right uncued
                        
                        % calculate absolute angle of stim 2 (test stim)
                        orient_dir2(noncatch_trials_idx,:)   = conds_single_rep_merged.orient_dir(noncatch_trials_idx,:) + delta_vec0;
                        
                        % apply circular wrap 
                        orient_dir2 = mod(orient_dir2, 360*ones(size(orient_dir2)));

                        % Ensure that absolute angles of left and right test stimuli don't overlap. 
                        % Any combination of two dot test stimuli that are
                        % less than 70 degrees apart are considered too close.
                        dot_too_close = abs(circulardiff(orient_dir2(:,1),orient_dir2(:,2),360)) <= params.stim.dot.min_ang_distance_test_stim;
                        
                        if sum(dot_too_close)==0
                            break;
                        end
                    end
                    
                    delta_vec(noncatch_trials_idx,:) = delta_vec0; 
                    
                    % find unique WM test im nr
                    idx_wm      = reshape(params.stim.dot.unique_im_nrs_wm_test,4,[]);
                    [~,idx_l]   = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.dot.unique_im_nrs_core);
                    [~,idx_r]   = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.dot.unique_im_nrs_core);
                    assert(isequal(length(idx_l),length(idx_r)))

                    stim2_im_nr = NaN(size(conds_single_rep_merged,1),2);
                    [~,shuffle_delta] = ismember(delta_vec0,params.stim.dot.delta_from_ref); clear delta_vec0;
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
                    conds_single_rep_merged.is_lure           = NaN(size(conds_single_rep_merged,1),2); 
                

                % Add IMG text prompt and quiz dots
                elseif strcmp(taskClass{:},'img')
                    % Add IMG text prompt and quiz dots
                    n_quiz_images       = length(unique(params.stim.dot.imagery_quiz_images)); % 2 types: straddle single dot or not.
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    shuffle_delta       = [shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images)); ... % 1 = straddle, 2 = both are left or right from single dot
                                            shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images))]'; % 50/50 chance of seeing 1 or 2
                    
                    % create yes/no overlap vector:
                    % 1=yes straddle, 2=no straddle (50% chance)
                    % col 1 = left stim, col 2 = right stim
                    % left/right stim is not yoked.
                    img_vec                        = NaN(size(conds_single_rep_merged,1),2);
                    img_vec(noncatch_trials_idx,:) = shuffle_delta;
                    idx_l                          = NaN(size(conds_single_rep_merged,1),2);
                    idx_r                          = NaN(size(conds_single_rep_merged,1),2);
                    
                    % find unique IMG test stim nrs that correspond to the special core stim nr in stim1                  
                    [~,idx_l_tmp]   = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.dot.unique_im_nrs_specialcore);
                    [~,idx_r_tmp]   = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.dot.unique_im_nrs_specialcore);
                    idx_img     = params.stim.rdk.imagery_quiz_dot_stim_nr; % 8x20 matrix:8 special core stim x 20 (10 yes + 10 no overlap) test images per unique image 
                    idx_img_yes = idx_img(:,params.stim.dot.imagery_quiz_images==1); % nr stim x 10 yes overlap test images per unique image
                    idx_img_no  = idx_img(:,params.stim.dot.imagery_quiz_images==2); % nr stim x 10 no overlap test images per unique image
                    assert(all(ismember(unique(params.stim.dot.unique_im_nrs_specialcore([idx_l_tmp;idx_r_tmp])), params.stim.dot.imagery_quiz_dot_specialcore_stim_nr)))
                    assert(isequal(size(idx_img_yes),size(idx_img_no)));
                    
                    idx_l(noncatch_trials_idx==1) = idx_l_tmp; clear idx_l_tmp
                    idx_r(noncatch_trials_idx==1) = idx_r_tmp; clear idx_r_tmp
                    
                    % To avoid that the table becomes massive, we provide
                    % one number that tells us something about two quiz dots:
                    % we populate the stim2_orient column with the orientation
                    % of the imagery quiz dots test images, which is the
                    % the orientation 1 imagery quiz dots (?) relative to 0 deg = 12 o'clock. 
                    % if quiz dots "straddle" they are left and right of the single dot location.
                    ori_img_yes = params.stim.dot.imagery_quiz_dot_orient_deg(:,params.stim.rdk.imagery_quiz_images==1,:); % 8 locations x 10 yes straddle quiz images x 2 dots
                    ori_img_no  = params.stim.dot.imagery_quiz_dot_orient_deg(:,params.stim.rdk.imagery_quiz_images==2,:); % 8 locations x 10 no straddle quiz images x 2 dots
                    
                    assert(all(reshape(ori_img_yes(1:4,:,:),1,[]) > 180)); % left quiz dots locations fall within 180 and 360 deg
                    assert(all(reshape(ori_img_no(5:end,:,:),1,[]) < 180)); % right quiz dots locations fall within 0 and 180 deg
                    
                    
                    yes_trial_idx_l = false(size(img_vec,1),1);
                    no_trial_idx_l  = false(size(img_vec,1),1);
                    yes_trial_idx_r = false(size(img_vec,1),1); 
                    no_trial_idx_r  = false(size(img_vec,1),1);
                    
                    % convert to 4 logical vectors, for l/r stim x yes/no overlap
                    yes_trial_idx_l(img_vec(:,1)==1)=true; % Quiz dots straddle (1=yes) or not (2=no) left
                    no_trial_idx_l(img_vec(:,1)==2)=true; % Quiz dots straddle (1=yes) or not (2=no) left
                    yes_trial_idx_r(img_vec(:,2)==1)=true; % Quiz dots straddle (1=yes) or not (2=no) right 
                    no_trial_idx_r(img_vec(:,2)==2)=true; % Quiz dots straddle (1=yes) or not (2=no) right
                    
                    while 1
                        dot_too_close = [];
                        qz_im           = NaN(size(img_vec)); % note: with catch trials
                        ori_dir         = NaN(size(img_vec)); % note: with catch trials
                        
                        % randomly sample from the 10 exemplar quiz dot images (WITHOUT replacement)
                        randomly_selected_yes_test_images_l  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_l), false);  % left stim yes
                        randomly_selected_no_test_images_l   = randsample(size(idx_img_no,2),sum(no_trial_idx_l), false); % left stim no
                        randomly_selected_yes_test_images_r  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_r), false); % right stim yes
                        randomly_selected_no_test_images_r   = randsample(size(idx_img_no,2),sum(no_trial_idx_r), false); % right stim no
                        
                        tmp_ori = NaN(size(ori_img_yes,1),2,2);
                        
                        yes_counter = 1; no_counter = 1;
                        for ll = 1:length(idx_l) % loop over unique images on the left
                            if conds_single_rep_merged.is_cued(ll)==1
                                if yes_trial_idx_l(ll)==1
                                    qz_im(ll,1)   = idx_img_yes(idx_l(ll),randomly_selected_yes_test_images_l(yes_counter));
                                    tmp_ori(ll,1,:) = squeeze(ori_img_yes(idx_l(ll),randomly_selected_yes_test_images_l(yes_counter),:)); % store both dot locations to check if quiz dots are too close.
                                    ori_dir(ll,1) = NaN; % rad2deg(circ_mean(squeeze(deg2rad(ori_img_yes(idx_l(ll),randomly_selected_yes_test_images_l(yes_counter),:))),[],1)); % average angle in deg of dot 1 and 2 (0 deg = 12 o'clock) -- note left col using first angle
                                    yes_counter   = yes_counter+1;
                                elseif no_trial_idx_l(ll)==1
                                    qz_im(ll,1)   = idx_img_no(idx_l(ll),randomly_selected_no_test_images_l(no_counter));
                                    tmp_ori(ll,1,:) = squeeze(ori_img_no(idx_l(ll),randomly_selected_no_test_images_l(no_counter),:)); % store both dot locations to check if quiz dots are too close.
                                    ori_dir(ll,1) = NaN; % rad2deg(circ_mean(squeeze(deg2rad(ori_img_no(idx_l(ll),randomly_selected_no_test_images_l(no_counter),:))),[],1)); % average angle in deg of dot 1 and 2 (0 deg = 12 o'clock) -- note left col, using first angle
                                    no_counter    = no_counter+1;
                                elseif isnan(img_vec(ll,1))
                                    qz_im(ll,1)   = NaN;
                                    ori_dir(ll,1) = NaN;
                                    tmp_ori(ll,1,:) = NaN(2,1);
                                else
                                    error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                                end
                            else
                                qz_im(ll,1)   = NaN;
                                ori_dir(ll,1) = NaN;
                                % remove uncued right side from table
                                conds_single_rep_merged.stim_nr_right(ll) = NaN;
                                conds_single_rep_merged.stim_class_name(ll,2) = {NaN};
                                conds_single_rep_merged.orient_dir(ll,2) = NaN;
                                conds_single_rep_merged.contrast(ll,2) = NaN;
                                conds_single_rep_merged.is_special_core(ll,2) = NaN;
                            end
                        end
                        
                        yes_counter = 1; no_counter = 1;
                        for ll = 1:length(idx_r) % loop over unique images on the right
                            if conds_single_rep_merged.is_cued(ll)==2
                                if yes_trial_idx_r(ll)==1
                                    qz_im(ll,2)   = idx_img_yes(idx_r(ll),randomly_selected_yes_test_images_r(yes_counter));
                                    tmp_ori(ll,2,:) = squeeze(ori_img_yes(idx_r(ll),randomly_selected_yes_test_images_r(yes_counter),:)); % store both dot locations to check if quiz dots are too close.
                                    ori_dir(ll,2) = NaN; % rad2deg(circ_mean(squeeze(deg2rad(ori_img_yes(idx_r(ll),randomly_selected_yes_test_images_r(yes_counter),:))),[],1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note right side, using first angle
                                    yes_counter   = yes_counter+1;
                                elseif no_trial_idx_r(ll)==1
                                    qz_im(ll,2)   = idx_img_no(idx_r(ll),randomly_selected_no_test_images_r(no_counter));
                                    tmp_ori(ll,2,:) = squeeze(ori_img_no(idx_r(ll),randomly_selected_no_test_images_r(no_counter),:)); % store both dot locations to check if quiz dots are too close.
                                    ori_dir(ll,2) = NaN; % rad2deg(circ_mean(squeeze(deg2rad(ori_img_no(idx_r(ll),randomly_selected_no_test_images_r(no_counter),:))),[],1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note right side, using first angle
                                    no_counter = no_counter+1;
                                elseif isnan(img_vec(ll,1))
                                    qz_im(ll,2)   = NaN;
                                    ori_dir(ll,2) = NaN;
                                    tmp_ori(ll,2,:) = NaN(2,1);
                                else
                                    error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                                end
                            else
                                qz_im(ll,2)   = NaN;
                                ori_dir(ll,2) = NaN;
                                % remove uncued left side from table
                                conds_single_rep_merged.stim_nr_left(ll) = NaN;
                                conds_single_rep_merged.stim_class_name(ll,1) = {NaN};
                                conds_single_rep_merged.orient_dir(ll,1) = NaN;
                                conds_single_rep_merged.contrast(ll,1) = NaN;
                                conds_single_rep_merged.is_special_core(ll,1) = NaN;
                            end
                        end

                   
                        % calculate absolute angle of stim 2 (test stim), apply circular wrap 
                        upper_left = max(squeeze(tmp_ori(:,1,:)),[],2); % quiz dot angle closests to upper vertical meridian.
                        lower_left = min(squeeze(tmp_ori(:,1,:)),[],2); % quiz dot angle closests to lower vertical meridian.
                        
                        lower_right = max(squeeze(tmp_ori(:,2,:)),[],2); % quiz dot angle closests to lower vertical meridian.
                        upper_right = min(squeeze(tmp_ori(:,2,:)),[],2); % quiz dot angle closests to upper vertical meridian.

                        abs_angle_upper = abs(circulardiff(upper_left,upper_right,360));
                        abs_angle_lower = abs(circulardiff(lower_left,lower_right,360));

                        % Ensure that absolute angles of left and right test stimuli don't overlap. 
                        % Any combination of two dot test stimuli that are
                        % less than 70 degrees apart are considered too close.
                        dot_too_close_upper = abs_angle_upper <= params.stim.dot.min_ang_distance_test_stim;
                        dot_too_close_lower = abs_angle_lower <= params.stim.dot.min_ang_distance_test_stim;
                        
                        dot_too_close = cat(2,dot_too_close_upper,dot_too_close_lower);
                        
                        if sum(dot_too_close, 'all')==0
                            break;
                        end 
                    end
%                     % all aligned quiz dots should match the angle shown in stim1
%                     assert(all(isequal(ori_dir(yes_trial_idx_l,1),conds_single_rep_merged.orient_dir(yes_trial_idx_l,1)))); % left side
%                     assert(all(ismember(ori_dir(yes_trial_idx_r,2),conds_single_rep_merged.orient_dir(yes_trial_idx_r,2)))); % right side

                    assert(length(unique(qz_im(:)))==length(qz_im(:)))
                    assert(all(ismember(qz_im(yes_trial_idx_l,1),idx_img_yes))); % check if yes quiz dots image nr matches for left stim
                    assert(all(ismember(qz_im(no_trial_idx_l,1),idx_img_no))); % check if no quiz dots image nr matches for left stim
                    assert(all(ismember(qz_im(yes_trial_idx_r,2),idx_img_yes))); % check if yes quiz dots image nr matches for right stim
                    assert(all(ismember(qz_im(no_trial_idx_r,2),idx_img_no))); % check if no quiz dots image nr matches for right stim
                    
                    assert(isequal(sum(img_vec(:,1)==1),sum(img_vec(:,1)==2))); % right stim side: equal yes/no overlap?
                    assert(isequal(sum(img_vec(:,2)==1),sum(img_vec(:,2)==2))); % left stim side: equal yes/no overlap?
                    
                    % fill in vectors & condition table
                    conds_single_rep_merged.stim2_delta      = img_vec; % yes (1) /no (2) dots are aligned
                    conds_single_rep_merged.stim2_im_nr      = qz_im;   % unique img quiz dot stim number
                    conds_single_rep_merged.stim2_orient_dir = ori_dir; % angle of dot locations
                    
                    % NO lures in IMG, set to NaN
                    conds_single_rep_merged.is_lure          = NaN(size(conds_single_rep_merged,1),2);

                    
                % Add ltm pair.
                elseif strcmp(taskClass{:},'ltm')
                    % preallocate space
                    conds_single_rep_merged.stim2_im_nr       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure           = NaN(size(conds_single_rep_merged,1),2); 
                    
                    [stim2_nr, stim2_is_lure, stim2_match] = vcd_getLTMTestStim(params, conds_single_rep_merged);

                    % check if catch trials match
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    assert(isequal(isnan(stim2_nr(:,1)),~noncatch_trials_idx))
                    assert(isequal(isnan(stim2_nr(:,2)),~noncatch_trials_idx))
                    
                    conds_single_rep_merged.stim2_im_nr = stim2_nr;
                    conds_single_rep_merged.is_lure     = stim2_is_lure;
                    conds_single_rep_merged.stim2_delta = stim2_match;
                    
                    % add orientation info to A->B trials
                    tmp1 = vcd('fullinfo', stim2_nr(noncatch_trials_idx,1));
                    tmp2 = vcd('fullinfo', stim2_nr(noncatch_trials_idx,2));
                    conds_single_rep_merged.stim2_orient_dir(noncatch_trials_idx,:) = cat(2, cat(1, tmp1.orient_dir), cat(1, tmp2.orient_dir));
                    clear tmp1 tmp2
                    
                    % We want to test both A->B pairings and B->A pairings.
                    % Make a copy of conditions table
                    conds_single_rep_mergedB = conds_single_rep_merged;
                    
                     % flip order A->B to B->A 
                    conds_single_rep_mergedB.stim_nr_left  = conds_single_rep_merged.stim2_im_nr(:,1);
                    conds_single_rep_mergedB.stim_nr_right = conds_single_rep_merged.stim2_im_nr(:,2);
                    
                    % update stim2
                    conds_single_rep_mergedB.stim2_im_nr      = [conds_single_rep_merged.stim_nr_left, conds_single_rep_merged.stim_nr_right];
                    conds_single_rep_mergedB.stim2_orient_dir = conds_single_rep_merged.orient_dir; % stim2 orient_dir B = stim1 orient_dir A
                    conds_single_rep_mergedB.orient_dir       = conds_single_rep_merged.stim2_orient_dir; % stim 1 orient_dir B = orient_dir A
                    
                    % update stimclass name
                    conds_single_rep_mergedB.stim_class_name(noncatch_trials_idx,1) = vcd('stimtostimclassname', conds_single_rep_mergedB.stim_nr_left(noncatch_trials_idx));
                    conds_single_rep_mergedB.stim_class_name(noncatch_trials_idx,2) = vcd('stimtostimclassname', conds_single_rep_mergedB.stim_nr_right(noncatch_trials_idx));
                    
                    % update stim info 
                    for ii = 1:size(conds_single_rep_mergedB,1)
                        % Get stim nr for left loc
                        tmp_nr_l = conds_single_rep_mergedB.stim_nr_left(ii);
                        if isnan(tmp_nr_l) || (tmp_nr_l==0) % catch trials
                            % keep as "dot"
                            assert(isequal(conds_single_rep_mergedB.is_catch(ii),1))
                            % set stim1 nr to 0000 and stim2 nr to NaN
                            conds_single_rep_mergedB.stim_nr_left(ii)  = 0;
                            conds_single_rep_mergedB.stim2_im_nr(ii,1) = NaN;
                        else % other stim classes
                            tmp = vcd('fullinfo', tmp_nr_l);
                            assert(~isnan(tmp.orient_dir));
                            % conds_single_rep_mergedB.orient_dir(ii,1) = tmp.orient_dir; % check orientation
                            assert(isequal(tmp.orient_dir,conds_single_rep_mergedB.orient_dir(ii,1))); 
                            assert(isequal(tmp.stim_loc, 1)); % must be left
                            if ismember(tmp_nr_l, params.stim.gabor.unique_im_nrs_core) % if gabor
                                conds_single_rep_mergedB.gbr_phase(ii,1) = tmp.gbr_phase; % add gabor phase
                                conds_single_rep_mergedB.contrast(ii,1)  = tmp.contrast; % add gabor contrast (0.8)
                                assert(isequal(tmp.stim_class, 1)); % must be gabor
                            elseif ismember(tmp_nr_l, params.stim.rdk.unique_im_nrs_core) % if rdk
                                conds_single_rep_mergedB.rdk_coherence(ii,1) = tmp.rdk_coherence; % add rdk motion coherence
                                assert(isequal(tmp.stim_class, 2)); % must be rdk
                            elseif ismember(tmp_nr_l, params.stim.dot.unique_im_nrs_core) % if dot
                                assert(~isnan(conds_single_rep_mergedB.orient_dir(ii,1)));
                                assert(isequal(tmp.stim_class, 3)); % must be dot    
                            elseif ismember(tmp_nr_l, params.stim.obj.unique_im_nrs_core) % if obj 
                                assert(isequal(tmp.stim_class, 4)); % must be obj
                                % add object category info
                                conds_single_rep_mergedB.super_cat(ii,1)        = tmp.super_cat;
                                conds_single_rep_mergedB.basic_cat(ii,1)        = tmp.basic_cat;
                                conds_single_rep_mergedB.sub_cat(ii,1)          = tmp.sub_cat;
                                conds_single_rep_mergedB.affordance_cat(ii,1)   = tmp.affordance_cat;
                                conds_single_rep_mergedB.super_cat_name(ii,1)   = tmp.super_cat_name;
                                conds_single_rep_mergedB.basic_cat_name(ii,1)   = tmp.basic_cat_name;
                                conds_single_rep_mergedB.sub_cat_name(ii,1)     = tmp.sub_cat_name;
                                conds_single_rep_mergedB.affordance_name(ii,1)  = tmp.affordance_name;
                            end
                        end
                    end
                    
                    % update stim info 
                    for ii = 1:size(conds_single_rep_mergedB,1)
                        % Get stim nr for right loc
                        tmp_nr_r = conds_single_rep_mergedB.stim_nr_right(ii); 
                        if isnan(tmp_nr_r) || (tmp_nr_r==0) % catch trials
                            % keep as "dot"
                            assert(isequal(conds_single_rep_mergedB.is_catch(ii),1))
                            % set stim1 nr to 0000 and stim2 nr to NaN
                            conds_single_rep_mergedB.stim_nr_right(ii)  = 0;
                            conds_single_rep_mergedB.stim2_im_nr(ii,2)  = NaN;
                        else % incorrect non-lure
                            tmp = vcd('fullinfo', tmp_nr_r);
                            assert(isequal(tmp.stim_loc, 2)); % must be right
                            assert(isequal(tmp.orient_dir,conds_single_rep_mergedB.orient_dir(ii,2))); % check orientation
                            % conds_single_rep_mergedB.orient_dir(ii,2) = tmp.orient_dir; % check orientation
                            if ismember(tmp_nr_r, params.stim.gabor.unique_im_nrs_core) % if gabor
                                conds_single_rep_mergedB.gbr_phase(ii,2) = tmp.gbr_phase; % add gabor phase
                                conds_single_rep_mergedB.contrast(ii,2)  = tmp.contrast; % add gabor contrast (0.8)
                                assert(isequal(tmp.stim_class, 1)); % must be gabor
                            elseif ismember(tmp_nr_r, params.stim.rdk.unique_im_nrs_core) % if rdk
                                conds_single_rep_mergedB.rdk_coherence(ii,2) = tmp.rdk_coherence; % add rdk motion coherence
                                assert(isequal(tmp.stim_class, 2)); % must be rdk
                            elseif ismember(tmp_nr_r, params.stim.dot.unique_im_nrs_core) % if dot
                                assert(isequal(tmp.stim_class, 3)); % must be dot
                                assert(~isnan(conds_single_rep_mergedB.orient_dir(ii,2)));
                            elseif ismember(tmp_nr_r, params.stim.obj.unique_im_nrs_core) % if obj
                                assert(isequal(tmp.stim_class, 4)); % must be obj
                                % add object category info
                                conds_single_rep_mergedB.super_cat(ii,2)        = tmp.super_cat;
                                conds_single_rep_mergedB.basic_cat(ii,2)        = tmp.basic_cat;
                                conds_single_rep_mergedB.sub_cat(ii,2)          = tmp.sub_cat;
                                conds_single_rep_mergedB.affordance_cat(ii,2)   = tmp.affordance_cat;
                                conds_single_rep_mergedB.super_cat_name(ii,2)   = tmp.super_cat_name;
                                conds_single_rep_mergedB.basic_cat_name(ii,2)   = tmp.basic_cat_name;
                                conds_single_rep_mergedB.sub_cat_name(ii,2)     = tmp.sub_cat_name;
                                conds_single_rep_mergedB.affordance_name(ii,2)  = tmp.affordance_name;
                            end
                        end
                    end
                    % update unique_trial_nr
                    conds_single_rep_mergedB.unique_trial_nr = max(conds_single_rep_mergedB.unique_trial_nr) + conds_single_rep_mergedB.unique_trial_nr;
                    % check unique trial number (should continue counting)
                    assert(isequal(conds_single_rep_mergedB.unique_trial_nr, ...
                        [min(conds_single_rep_mergedB.unique_trial_nr):max(conds_single_rep_mergedB.unique_trial_nr)]'));
                    
                    % merge condition tables
                    conds_single_rep_merged = cat(1,conds_single_rep_merged,conds_single_rep_mergedB);
                    
                    % check unique trial number (should continue counting)
                    assert(isequal(sum(~noncatch_trials_idx)*2, sum(conds_single_rep_merged.is_catch)));
                    assert(isequal(conds_single_rep_merged.unique_trial_nr, ...
                        [min(conds_single_rep_merged.unique_trial_nr):max(conds_single_rep_merged.unique_trial_nr)]'));
                else
                    conds_single_rep_merged.stim2_im_nr       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir  = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta       = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure           = NaN(size(conds_single_rep_merged,1),2);
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

                    % we shuffle all conditions, without constraints.
                    shuffle_c      = randperm(size(cond_table,1),size(cond_table,1));
                    conds_shuffle0 = cond_table(shuffle_c,:);
                    
                    % get covert spatial attention cuing direction vector
                    [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1,conds_shuffle0, stimloc_cues, taskClass);

                    % Add cue vec
                    conds_shuffle0.is_cued = cue_vec;

                    % Add general catch trial vector (no stim)
                    conds_shuffle0.is_catch = zeros(size(conds_shuffle0,1),1);
                    
                    % Preallocate space for is_objectcatch trial vector
                    % (randomly chosen non-core rotation, only for OBJ-PC)
                    conds_shuffle0.is_objectcatch = NaN(size(conds_shuffle0,1),1);

                    % Do some checks:
                    assert(isequal(length(conds_shuffle0.unique_im_nr),n_unique_cases)); % Check unique image nr
                    assert(isequal(sum(conds_shuffle0.stimloc==1),n_unique_cases/2));  % Check left stim loc
                    assert(isequal(sum(conds_shuffle0.stimloc==2),n_unique_cases/2));  % Check right stim loc
                    assert(isequal(sort(conds_shuffle0.super_cat_name),sort(repelem(params.stim.obj.super_cat,n_sub_cat)'))); % Check supercat
                    if ~ismember(taskClass{:}, {'ltm','img'})
                        assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),sort(params.stim.obj.facing_dir_deg,'ascend')')); % Check object facing dir
                    else
                        assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),sort(params.stim.obj.facing_dir_deg(ismember(params.stim.obj.unique_im_nrs_core,params.stim.obj.unique_im_nrs_specialcore))','ascend'))); % Check contrast
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

                    clear conds_shuffle0 conds_shuffle1 conds_shuffle2
                end

                % Merge trials and add fix cue thickening direction
                conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(params, conds_master_single_rep);
                clear conds_master_single_rep;

                % Add WM change.
                if strcmp(taskClass{:},'wm')
                    % create shuffled delta vector
                    n_deltas            = length(params.stim.obj.delta_from_ref);
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    facing_dir2         = NaN(size(conds_single_rep_merged.orient_dir));
                    delta_vec           = NaN(size(conds_single_rep_merged.orient_dir));
                    delta_vec0(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1,1) = shuffle_concat(repelem(params.stim.obj.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left cued
                    delta_vec0(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==1),1) = shuffle_concat(repelem(params.stim.obj.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % left uncued
                    delta_vec0(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2,2) = shuffle_concat(repelem(params.stim.obj.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right cued
                    delta_vec0(~(conds_single_rep_merged.is_cued(noncatch_trials_idx)==2),2) = shuffle_concat(repelem(params.stim.obj.delta_from_ref, (n_unique_cases/2/n_deltas)),1); % right uncued
                    
                    facing_dir2(noncatch_trials_idx,:) = conds_single_rep_merged.orient_dir(noncatch_trials_idx,:) + delta_vec0;
                    delta_vec(noncatch_trials_idx,:)   = delta_vec0; 
                    
                    % find unique WM test im nr
                    idx_wm = reshape(params.stim.obj.unique_im_nrs_wm_test,4,[]);
                    [~,idx_l] = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.obj.unique_im_nrs_core);
                    [~,idx_r] = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.obj.unique_im_nrs_core);
                    assert(isequal(length(idx_l),length(idx_r)))
                    
                    stim2_im_nr = NaN(size(conds_single_rep_merged,1),2);
                    [~,shuffle_delta] = ismember(delta_vec0,params.stim.obj.delta_from_ref); clear delta_vec0;
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
                    conds_single_rep_merged.is_lure          = NaN(size(conds_single_rep_merged,1),2);
                

                % Add IMG text prompt and quiz dots
                elseif strcmp(taskClass{:},'img')
                    % Add IMG text prompt and quiz dots
                    n_quiz_images       = length(unique(params.stim.obj.imagery_quiz_images)); % 2 types: overlap or not.
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    shuffle_delta       = [shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images)); ... % 1 = overlap, 2 = no overlap
                                            shuffle_concat(1:n_quiz_images, (size(conds_single_rep_merged(noncatch_trials_idx,:),1)/n_quiz_images))]'; % 50/50 chance of seeing 1 or 2
                    
                    % create yes/no overlap vector:
                    % 1=yes overlap, 2=no overlap (50% chance)
                    % col 1 = left stim, col 2 = right stim
                    % left/right stim is not yoked.
                    img_vec                        = NaN(size(conds_single_rep_merged,1),2);
                    img_vec(noncatch_trials_idx,:) = shuffle_delta;
                    idx_l                          = NaN(size(conds_single_rep_merged,1),2);
                    idx_r                          = NaN(size(conds_single_rep_merged,1),2);
                    
                    % find unique IMG test stim nrs that correspond to the special core stim nr in stim1                  
                    [~,idx_l_tmp]   = ismember(conds_single_rep_merged.stim_nr_left(noncatch_trials_idx),params.stim.obj.unique_im_nrs_specialcore);
                    [~,idx_r_tmp]   = ismember(conds_single_rep_merged.stim_nr_right(noncatch_trials_idx),params.stim.obj.unique_im_nrs_specialcore);
                    idx_img     = params.stim.obj.imagery_quiz_dot_stim_nr; % 8x20 matrix:8 special core stim x 20 (10 yes + 10 no overlap) test images per unique image 
                    idx_img_yes = idx_img(:,params.stim.obj.imagery_quiz_images==1); % nr stim x 10 yes overlap test images per unique image
                    idx_img_no  = idx_img(:,params.stim.obj.imagery_quiz_images==2); % nr stim x 10 no overlap test images per unique image
                    assert(all(ismember(unique(params.stim.obj.unique_im_nrs_specialcore([idx_l_tmp;idx_r_tmp])), params.stim.obj.imagery_quiz_dot_specialcore_stim_nr)))
                    assert(isequal(size(idx_img_yes),size(idx_img_no)));
                    
                    idx_l(noncatch_trials_idx==1) = idx_l_tmp; clear idx_l_tmp
                    idx_r(noncatch_trials_idx==1) = idx_r_tmp; clear idx_r_tmp
                    
                    % To avoid that the table becomes massive, we provide
                    % one number that tells us something about two quiz dots:
                    % we populate the stim2_orient column with the "orientation" 
                    % of the imagery quiz dots test images, which is the
                    % the orientation of the line connecting the two
                    % imagery quiz dots dots, relative to 0 deg = 12 o'clock. 
                    % when quiz dots overlap, their location is within the 
                    % outline of the object.
                    ori_img_yes = params.stim.obj.imagery_quiz_dot_orient_deg(:,params.stim.obj.imagery_quiz_images==1,:); % 8 locations x 10 yes overlap quiz images x 2 dots
                    ori_img_no  = params.stim.obj.imagery_quiz_dot_orient_deg(:,params.stim.obj.imagery_quiz_images==2,:); % 8 locations x 10 no overlap quiz images x 2 dots
                    
                    
                    qz_im           = NaN(size(img_vec)); % note: with catch trials
                    ori_dir         = NaN(size(img_vec)); % note: with catch trials
                    yes_trial_idx_l = false(size(img_vec,1),1);
                    no_trial_idx_l  = false(size(img_vec,1),1);
                    yes_trial_idx_r = false(size(img_vec,1),1); 
                    no_trial_idx_r  = false(size(img_vec,1),1);
                    
                    % convert to 4 logical vectors, for l/r stim x yes/no overlap
                    yes_trial_idx_l(img_vec(:,1)==1)=true; % Quiz dots overlap (1=yes) or not (2=no) left
                    no_trial_idx_l(img_vec(:,1)==2)=true; % Quiz dots overlap (1=yes) or not (2=no) left
                    yes_trial_idx_r(img_vec(:,2)==1)=true; % Quiz dots overlap (1=yes) or not (2=no) right 
                    no_trial_idx_r(img_vec(:,2)==2)=true; % Quiz dots overlap (1=yes) or not (2=no) right
                    
                    % randomly sample from the 10 exemplar quiz dot images (WITHOUT replacement)
                    randomly_selected_yes_test_images_l  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_l), false);  % left stim yes 
                    randomly_selected_no_test_images_l   = randsample(size(idx_img_no,2),sum(no_trial_idx_l), false); % left stim no 
                    randomly_selected_yes_test_images_r  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_r), false); % right stim yes 
                    randomly_selected_no_test_images_r   = randsample(size(idx_img_no,2),sum(no_trial_idx_r), false); % right stim no
                    
                    yes_counter = 1; no_counter = 1;
                    for ll = 1:length(idx_l) % loop over unique images on the left
                        if conds_single_rep_merged.is_cued(ll)==1
                            if yes_trial_idx_l(ll)==1
                                qz_im(ll,1)   = idx_img_yes(idx_l(ll),randomly_selected_yes_test_images_l(yes_counter));
                                ori_dir(ll,1) = NaN; % squeeze(ori_img_yes(idx_l(ll),randomly_selected_yes_test_images_l(yes_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left col using first angle
                                yes_counter   = yes_counter+1;
                            elseif no_trial_idx_l(ll)==1
                                qz_im(ll,1)   = idx_img_no(idx_l(ll),randomly_selected_no_test_images_l(no_counter));
                                ori_dir(ll,1) = NaN; % squeeze(ori_img_no(idx_l(ll),randomly_selected_no_test_images_l(no_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left col, using first angle
                                no_counter    = no_counter+1;
                            elseif isnan(img_vec(ll,1))
                                qz_im(ll,1)   = NaN;
                                ori_dir(ll,1) = NaN;
                            else
                                error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                            end
                        else
                            qz_im(ll,1)   = NaN;
                            ori_dir(ll,1) = NaN;
                            % remove uncued right side from table
                            conds_single_rep_merged.stim_nr_left(ll) = NaN;
                            conds_single_rep_merged.stim_class_name(ll,2) = {NaN};
                            conds_single_rep_merged.orient_dir(ll,2) = NaN;
                            conds_single_rep_merged.contrast(ll,2) = NaN;
                            conds_single_rep_merged.is_special_core(ll,2) = NaN;
                            conds_single_rep_merged.super_cat(ll,2) = NaN;
                            conds_single_rep_merged.basic_cat(ll,2) = NaN;
                            conds_single_rep_merged.sub_cat(ll,2) = NaN;
                            conds_single_rep_merged.affordance_cat(ll,2) = NaN;
                            conds_single_rep_merged.super_cat_name(ll,2) = {NaN};
                            conds_single_rep_merged.basic_cat_name(ll,2) = {NaN};
                            conds_single_rep_merged.sub_cat_name(ll,2) = {NaN};
                            conds_single_rep_merged.affordance_name(ll,2) = {NaN};
                        end
                    end
                    
                    yes_counter = 1; no_counter = 1;
                    for ll = 1:length(idx_r) % loop over unique images on the right
                        if conds_single_rep_merged.is_cued(ll)==2
                            if yes_trial_idx_r(ll)==1
                                qz_im(ll,2)   = idx_img_yes(idx_r(ll),randomly_selected_yes_test_images_r(yes_counter));
                                ori_dir(ll,2) = NaN; % squeeze(ori_img_yes(idx_r(ll),randomly_selected_yes_test_images_r(yes_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note right side, using first angle
                                yes_counter   = yes_counter+1;
                            elseif no_trial_idx_r(ll)==1
                                qz_im(ll,2)   = idx_img_no(idx_r(ll),randomly_selected_no_test_images_r(no_counter));
                                ori_dir(ll,2) = NaN; % squeeze(ori_img_no(idx_r(ll),randomly_selected_no_test_images_r(no_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note right side, using first angle
                                no_counter = no_counter+1;
                            elseif isnan(img_vec(ll,1))
                                qz_im(ll,2)   = NaN;
                                ori_dir(ll,2) = NaN;
                            else
                                error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                            end
                        else
                            qz_im(ll,2)   = NaN;
                            ori_dir(ll,2) = NaN;
                            % remove uncued right side from table
                            conds_single_rep_merged.stim_nr_left(ll) = NaN;
                            conds_single_rep_merged.stim_class_name(ll,1) = {NaN};
                            conds_single_rep_merged.orient_dir(ll,1) = NaN;
                            conds_single_rep_merged.contrast(ll,1) = NaN;
                            conds_single_rep_merged.is_special_core(ll,1) = NaN;
                            conds_single_rep_merged.super_cat(ll,1) = NaN;
                            conds_single_rep_merged.basic_cat(ll,1) = NaN;
                            conds_single_rep_merged.sub_cat(ll,1) = NaN;
                            conds_single_rep_merged.affordance_cat(ll,1) = NaN;
                            conds_single_rep_merged.super_cat_name(ll,1) = {NaN};
                            conds_single_rep_merged.basic_cat_name(ll,1) = {NaN};
                            conds_single_rep_merged.sub_cat_name(ll,1) = {NaN};
                            conds_single_rep_merged.affordance_name(ll,1) = {NaN};
                        end
                    end
                    

                    assert(length(unique(qz_im(:)))==length(qz_im(:)))
                    assert(all(ismember(qz_im(yes_trial_idx_l,1),idx_img_yes))); % check if yes quiz dots image nr matches for left stim
                    assert(all(ismember(qz_im(no_trial_idx_l,1),idx_img_no))); % check if no quiz dots image nr matches for left stim
                    assert(all(ismember(qz_im(yes_trial_idx_r,2),idx_img_yes))); % check if yes quiz dots image nr matches for right stim
                    assert(all(ismember(qz_im(no_trial_idx_r,2),idx_img_no))); % check if no quiz dots image nr matches for right stim
                    
                    assert(isequal(sum(img_vec(:,1)==1),sum(img_vec(:,1)==2))); % right stim side: equal yes/no overlap?
                    assert(isequal(sum(img_vec(:,2)==1),sum(img_vec(:,2)==2))); % left stim side: equal yes/no overlap?
                    
                    % fill in vectors & condition table
                    conds_single_rep_merged.stim2_delta      = img_vec; % yes (1) /no (2) dots are aligned
                    conds_single_rep_merged.stim2_im_nr      = qz_im;   % unique img quiz dot stim number
                    conds_single_rep_merged.stim2_orient_dir = ori_dir; % angle of line connecting two quiz dot locations
                    
                    % NO lures in IMG, set to NaN
                    conds_single_rep_merged.is_lure          = NaN(size(conds_single_rep_merged,1),2);
                
                    % Add ltm pair.
                elseif strcmp(taskClass{:},'ltm')
                    % preallocate space
                    conds_single_rep_merged.stim2_im_nr      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure          = NaN(size(conds_single_rep_merged,1),2); 
                    
                    [stim2_nr, stim2_is_lure, stim2_match] = vcd_getLTMTestStim(params, conds_single_rep_merged);

                    % check if catch trials match
                    noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                    assert(isequal(isnan(stim2_nr(:,1)),~noncatch_trials_idx))
                    assert(isequal(isnan(stim2_nr(:,2)),~noncatch_trials_idx))
                    
                    conds_single_rep_merged.stim2_im_nr = stim2_nr;
                    conds_single_rep_merged.is_lure     = stim2_is_lure;
                    conds_single_rep_merged.stim2_delta = stim2_match;
                    
                    % add orientation info
                    tmp1 = vcd('fullinfo', stim2_nr(noncatch_trials_idx,1));
                    tmp2 = vcd('fullinfo', stim2_nr(noncatch_trials_idx,2));
                    conds_single_rep_merged.stim2_orient_dir(noncatch_trials_idx,:) = cat(2, cat(1, tmp1.orient_dir), cat(1, tmp2.orient_dir));
                    clear tmp1 tmp2
                    
                    % We want to test both A->B pairings and B->A pairings.
                    % Make a copy of conditions table
                    conds_single_rep_mergedB = conds_single_rep_merged;
                    
                     % flip order A->B to B->A 
                    conds_single_rep_mergedB.stim_nr_left  = conds_single_rep_merged.stim2_im_nr(:,1);
                    conds_single_rep_mergedB.stim_nr_right = conds_single_rep_merged.stim2_im_nr(:,2);
                    conds_single_rep_mergedB.orient_dir    = NaN(size(conds_single_rep_merged.orient_dir,1),2);
                    conds_single_rep_mergedB.contrast      = NaN(size(conds_single_rep_merged.orient_dir,1),2);
                    
                    % update stim2
                    conds_single_rep_mergedB.stim2_im_nr      = [conds_single_rep_merged.stim_nr_left, conds_single_rep_merged.stim_nr_right];
                    conds_single_rep_mergedB.stim2_orient_dir = conds_single_rep_merged.orient_dir;
                    
                    % update stimclass name
                    conds_single_rep_mergedB.stim_class_name(noncatch_trials_idx,1) = vcd('stimtostimclassname', conds_single_rep_mergedB.stim_nr_left(noncatch_trials_idx));
                    conds_single_rep_mergedB.stim_class_name(noncatch_trials_idx,2) = vcd('stimtostimclassname', conds_single_rep_mergedB.stim_nr_right(noncatch_trials_idx));
                    
                    % update stim info (we only need to do this for the
                    % table with B->A trials, now that "B" test stim will
                    % become "A" reference stim and we do not store this
                    % info fotr "B" stim in table with A->B trials.
                    for ii = 1:size(conds_single_rep_mergedB,1)
                        % Get stim nr for left stim loc
                        tmp_nr_l = conds_single_rep_mergedB.stim_nr_left(ii);
                        
                        % set category info to nan (we'll add it if required)
                        conds_single_rep_mergedB.super_cat(ii,1)        = NaN;
                        conds_single_rep_mergedB.basic_cat(ii,1)        = NaN;
                        conds_single_rep_mergedB.sub_cat(ii,1)          = NaN;
                        conds_single_rep_mergedB.affordance_cat(ii,1)   = NaN;
                        conds_single_rep_mergedB.super_cat_name{ii,1}   = NaN;
                        conds_single_rep_mergedB.basic_cat_name{ii,1}   = NaN;
                        conds_single_rep_mergedB.sub_cat_name{ii,1}     = NaN;
                        conds_single_rep_mergedB.affordance_name{ii,1}  = NaN;
                        
                        if isnan(tmp_nr_l) || (tmp_nr_l==0) % catch trials
                            % keep as "obj"
                            assert(isequal(conds_single_rep_mergedB.is_catch(ii),1))
                            % set stim1 nr to 0000 and stim2 nr to NaN
                            conds_single_rep_mergedB.stim_nr_left(ii)  = 0;
                            conds_single_rep_mergedB.stim2_im_nr(ii,1) = NaN;

                        else % classic stim
                            tmp = vcd('fullinfo', tmp_nr_l); % Get stim info
                            conds_single_rep_mergedB.orient_dir(ii,1) = tmp.orient_dir; % update gabor tilt/rdk motion direction/dot angle
                            conds_single_rep_mergedB.contrast(ii,1)   = tmp.contrast; % add contrast
                            assert(isequal(tmp.stim_loc,1)); % must be left
                            if ismember(tmp_nr_l, params.stim.gabor.unique_im_nrs_core) % if gabor
                                conds_single_rep_mergedB.gbr_phase(ii,1) = tmp.gbr_phase; % add gabor phase
                                assert(isequal(tmp.stim_class,1)); % must be gabor
                            elseif ismember(tmp_nr_l, params.stim.rdk.unique_im_nrs_core) % if rdk
                                conds_single_rep_mergedB.rdk_coherence(ii,1) = tmp.rdk_coherence; % add motion coherence
                                assert(isequal(tmp.stim_class,2)); % must be rdk
                            elseif ismember(tmp_nr_l, params.stim.dot.unique_im_nrs_core) % if dot
                                assert(isequal(tmp.stim_class,3)); % must be dot
                            elseif ismember(tmp_nr_l, params.stim.obj.unique_im_nrs_core) % if obj
                                assert(isequal(tmp.stim_class,4)); % must be obj
                                % update category info
                                conds_single_rep_mergedB.super_cat(ii,1)        = tmp.super_cat;
                                conds_single_rep_mergedB.basic_cat(ii,1)        = tmp.basic_cat;
                                conds_single_rep_mergedB.sub_cat(ii,1)          = tmp.sub_cat;
                                conds_single_rep_mergedB.affordance_cat(ii,1)   = tmp.affordance_cat;
                                conds_single_rep_mergedB.super_cat_name(ii,1)   = tmp.super_cat_name;
                                conds_single_rep_mergedB.basic_cat_name(ii,1)   = tmp.basic_cat_name;
                                conds_single_rep_mergedB.sub_cat_name(ii,1)     = tmp.sub_cat_name;
                                conds_single_rep_mergedB.affordance_name(ii,1)  = tmp.affordance_name;
                            end
                        end
                    end
                        
                    % update stim info
                    for ii = 1:size(conds_single_rep_mergedB,1)
                        % Get stim nr for right stim loc
                        tmp_nr_r = conds_single_rep_mergedB.stim_nr_right(ii);
                        
                        % set category info to nan (we'll add it if required)
                        conds_single_rep_mergedB.super_cat(ii,2)        = NaN;
                        conds_single_rep_mergedB.basic_cat(ii,2)        = NaN;
                        conds_single_rep_mergedB.sub_cat(ii,2)          = NaN;
                        conds_single_rep_mergedB.affordance_cat(ii,2)   = NaN;
                        conds_single_rep_mergedB.super_cat_name{ii,2}   = NaN;
                        conds_single_rep_mergedB.basic_cat_name{ii,2}   = NaN;
                        conds_single_rep_mergedB.sub_cat_name{ii,2}     = NaN;
                        conds_single_rep_mergedB.affordance_name{ii,2}  = NaN;
                        
                        if isnan(tmp_nr_r) || (tmp_nr_r==0) % catch trials
                            % keep as "obj"
                            assert(isequal(conds_single_rep_mergedB.is_catch(ii),1))
                            % set stim1 nr to 0000 and stim2 nr to NaN
                            conds_single_rep_mergedB.stim_nr_right(ii)  = 0;
                            conds_single_rep_mergedB.stim2_im_nr(ii,2) = NaN;
                        else % other stim class
                            tmp = vcd('fullinfo', tmp_nr_r); % Get stim info
                            conds_single_rep_mergedB.orient_dir(ii,2) = tmp.orient_dir; % add gabor tilt/rdk mot dir/dot angle/obj facing direction
                            conds_single_rep_mergedB.contrast(ii,2)   = tmp.contrast; % add contrast info
                            assert(isequal(tmp.stim_loc,2)); % must be right
                            if ismember(tmp_nr_r, params.stim.gabor.unique_im_nrs_core) % if gabor
                                    conds_single_rep_mergedB.gbr_phase(ii,2) = tmp.gbr_phase; % add phase info
                                    assert(isequal(tmp.stim_class,1)); % must be gabor
                            elseif ismember(tmp_nr_r, params.stim.rdk.unique_im_nrs_core) % if rdk
                                    conds_single_rep_mergedB.rdk_coherence(ii,2) = tmp.rdk_coherence; % add coherence info
                                    assert(isequal(tmp.stim_class,2)); % must be rdk
                            elseif ismember(tmp_nr_r, params.stim.dot.unique_im_nrs_core) % if dot
                                assert(isequal(tmp.stim_class,3)); % must be dot
                            elseif ismember(tmp_nr_r, params.stim.obj.unique_im_nrs_core) % if obj
                                assert(isequal(tmp.stim_class,4)); % must be obj
                                conds_single_rep_mergedB.super_cat(ii,2)        = tmp.super_cat;
                                conds_single_rep_mergedB.basic_cat(ii,2)        = tmp.basic_cat;
                                conds_single_rep_mergedB.sub_cat(ii,2)          = tmp.sub_cat;
                                conds_single_rep_mergedB.affordance_cat(ii,2)   = tmp.affordance_cat;
                                conds_single_rep_mergedB.super_cat_name(ii,2)   = tmp.super_cat_name;
                                conds_single_rep_mergedB.basic_cat_name(ii,2)   = tmp.basic_cat_name;
                                conds_single_rep_mergedB.sub_cat_name(ii,2)     = tmp.sub_cat_name;
                                conds_single_rep_mergedB.affordance_name(ii,2)  = tmp.affordance_name;
                            end
                        end
                    end
                    
                    % update unique_trial_nr
                    conds_single_rep_mergedB.unique_trial_nr = max(conds_single_rep_mergedB.unique_trial_nr) + conds_single_rep_mergedB.unique_trial_nr;
                    
                    % check unique trial number (should continue counting)
                    assert(isequal(conds_single_rep_mergedB.unique_trial_nr, ...
                        [min(conds_single_rep_mergedB.unique_trial_nr):max(conds_single_rep_mergedB.unique_trial_nr)]'));

                    % merge condition tables
                    conds_single_rep_merged = cat(1,conds_single_rep_merged,conds_single_rep_mergedB);
                    
                    % check unique trial number (should continue counting)
                    assert(isequal(sum(~noncatch_trials_idx)*2, sum(conds_single_rep_merged.is_catch)));
                    assert(isequal(conds_single_rep_merged.unique_trial_nr, ...
                        [min(conds_single_rep_merged.unique_trial_nr):max(conds_single_rep_merged.unique_trial_nr)]'));

                else
                    conds_single_rep_merged.stim2_im_nr      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_delta      = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.stim2_orient_dir = NaN(size(conds_single_rep_merged,1),2);
                    conds_single_rep_merged.is_lure          = NaN(size(conds_single_rep_merged,1),2); 
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
            for ni = 1:n_super_cat
                n_basic_cat(ni) = length(params.stim.ns.basic_cat{ni});
                for nj = 1:n_basic_cat(ni)
                    n_sub_cat(ni,nj) = length(params.stim.ns.sub_cat{ni,nj});
                    super_cat_vec = cat(2, super_cat_vec, repelem(ni,n_sub_cat(ni,nj)));
                end
                sub_cat_vec = cat(2, sub_cat_vec, repelem(1:n_sub_cat(ni,nj),n_basic_cat(ni)));
            end
            if ismember(taskClass{:},{'ltm','img'})
                super_cat_vec = super_cat_vec(ismember(params.stim.ns.unique_im_nrs_core,params.stim.ns.unique_im_nrs_specialcore));
                sub_cat_vec = sub_cat_vec(ismember(params.stim.ns.unique_im_nrs_core,params.stim.ns.unique_im_nrs_specialcore));
            end
            assert(isequal(cond_table.super_cat,super_cat_vec')); 
            assert(isequal(cond_table.sub_cat,sub_cat_vec'));
            
            for rep = 1:nr_reps

                conds_master_single_rep = [];

                % we shuffle all conditions, without constraints, regardless of superordinate category class
                shuffle_c               = randperm(size(cond_table,1),size(cond_table,1));
                conds_master_single_rep = cond_table(shuffle_c,:);

                % Add catch trial vector
                conds_master_single_rep.is_catch = zeros(size(conds_master_single_rep,1),1);
                
                % Preallocate space for is_objectcatch trial vector
                % (randomly chosen non-core rotation, only for OBJ-PC)
                conds_master_single_rep.is_objectcatch = NaN(size(conds_master_single_rep,1),1);
                
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

                % No need to merge unique im into trials because there is only one stim.
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
                    noncatch_trials_idx  = conds_master_single_rep2.is_catch(:,1)==0;

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
                    [~,stim2_im_delta_idx] = ismember(shuffle_delta(1:n_unique_cases),params.stim.ns.change_im);
                    stim2_im_nr(noncatch_trials_idx,1)        = idx_wm(idx_c)' + stim2_im_delta_idx-1;
                    conds_master_single_rep2.stim2_delta      = delta_vec;
                    conds_master_single_rep2.stim2_im_nr      = stim2_im_nr;
                    conds_master_single_rep2.stim2_orient_dir = NaN(size(conds_master_single_rep2.orient_dir,1),2);
                    conds_master_single_rep2.is_lure          = NaN(size(conds_master_single_rep2,1),2);

                % take subset of IMG trials, add text prompt and quiz dots loc
                elseif strcmp(taskClass{:},'img')
                    n_quiz_images       = length(unique(params.stim.ns.imagery_quiz_images));
                    noncatch_trials_idx = conds_master_single_rep2.is_catch==0;
                    shuffle_delta       = shuffle_concat(1:n_quiz_images, (size(conds_master_single_rep2(noncatch_trials_idx,:),1)/n_quiz_images))';
                    
                    % create yes/no overlap vector:
                    % 1=yes overlap, 2=no overlap (50% chance)
                    % col 1 = left stim, col 2 = right stim
                    % left/right stim is not yoked.
                    img_vec                        = NaN(size(conds_master_single_rep2,1),1);
                    img_vec(noncatch_trials_idx,:) = shuffle_delta;
                    
                    % find unique IMG test stim nrs that correspond to the special core stim nr in stim1
                    idx_c     = NaN(size(conds_master_single_rep2,1),1); % <--- 14 trials x 1
                    [~,idx_c_tmp] = ismember(conds_master_single_rep2.stim_nr_left(noncatch_trials_idx),params.stim.ns.unique_im_nrs_specialcore); % central image only
                    idx_img   = reshape(params.stim.ns.unique_im_nrs_img_test,length(params.stim.ns.imagery_quiz_images),[])'; % 20x15 matrix: 10 yes + 10 no overlap test images for each of the 15 special core images
                    idx_img_yes = idx_img(:,params.stim.ns.imagery_quiz_images==1); % Quiz dots overlap (1=yes) or not (2=no)
                    idx_img_no  = idx_img(:,params.stim.ns.imagery_quiz_images==2); % Quiz dots overlap (1=yes) or not (2=no)
                    
                    assert(isequal(size(idx_img_yes),size(idx_img_no)))
                    assert(all(ismember(unique(params.stim.ns.unique_im_nrs_specialcore(idx_c_tmp)), params.stim.ns.imagery_quiz_dot_specialcore_stim_nr)))
                    
                    idx_c(noncatch_trials_idx==1) = idx_c_tmp; clear idx_c_tmp
                    
                    % To avoid that the table becomes massive, we provide
                    % one number that tells us something about two quiz dots:
                    % we populate the stim2_orient column with the "orientation" 
                    % of the imagery quiz dots test images, which is the
                    % the orientation of the line connecting the two
                    % imagery quiz dots dots, relative to 0 deg = 12 o'clock. 
                    % when quiz dots overlap, their location is within the 
                    % outline of the object in the scene.
                    ori_img_yes = params.stim.ns.imagery_quiz_dot_orient_deg(:,params.stim.ns.imagery_quiz_images==1,1,:); % 8 locations x 10 yes overlap quiz images x 2 dots
                    ori_img_no  = params.stim.ns.imagery_quiz_dot_orient_deg(:,params.stim.ns.imagery_quiz_images==2,1,:); % 8 locations x 10 no overlap quiz images x 2 dots
                    
                    
                    qz_im           = NaN(size(img_vec)); % note: with catch trials
                    ori_dir         = NaN(size(img_vec)); % note: with catch trials
                    yes_trial_idx_c = false(size(img_vec,1),1);
                    no_trial_idx_c  = false(size(img_vec,1),1);
                    
                    
                    % convert to 4 logical vectors, for l/r stim x yes/no overlap
                    yes_trial_idx_c(img_vec(:,1)==1)=true; % Quiz dots overlap (1=yes) or not (2=no) 
                    no_trial_idx_c(img_vec(:,1)==2)=true; % Quiz dots overlap (1=yes) or not (2=no) 
                    
                    % randomly sample from the 10 exemplar quiz dot images (WITHOUT replacement)
                    randomly_selected_yes_test_images_c  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_c), false);  % left stim yes 
                    randomly_selected_no_test_images_c   = randsample(size(idx_img_no,2),sum(no_trial_idx_c), false); % left stim no 
                    
                    yes_counter = 1; no_counter = 1;
                    for ll = 1:length(idx_c) % loop over unique images on the left
                        if yes_trial_idx_c(ll)==1
                            qz_im(ll)   = idx_img_yes(idx_c(ll),randomly_selected_yes_test_images_c(yes_counter));
                            ori_dir(ll) = NaN; % squeeze(ori_img_yes(idx_c(ll),randomly_selected_yes_test_images_c(yes_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left col using first angle
                            yes_counter   = yes_counter+1;
                        elseif no_trial_idx_c(ll)==1
                            qz_im(ll)   = idx_img_no(idx_c(ll),randomly_selected_no_test_images_c(no_counter));
                            ori_dir(ll) = NaN; % squeeze(ori_img_no(idx_c(ll),randomly_selected_no_test_images_c(no_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left col, using first angle
                            no_counter    = no_counter+1;
                        elseif isnan(img_vec(ll))
                            qz_im(ll)   = NaN;
                            ori_dir(ll) = NaN;
                        else
                            error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                        end
                    end
                    
                    assert(length(unique(qz_im(:)))==length(qz_im(:)))
                    assert(all(ismember(qz_im(yes_trial_idx_c,1),idx_img_yes))); % check if yes quiz dots image nr matches for central stim
                    assert(all(ismember(qz_im(no_trial_idx_c,1),idx_img_no))); % check if no quiz dots image nr matches for central stim                    
                    assert(isequal(sum(img_vec(:,1)==1),sum(img_vec(:,1)==2))); % equal yes/no overlap?
                    
                    % fill in vectors & condition table
                    conds_master_single_rep2.stim2_im_nr      = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_orient_dir = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_delta      = NaN(size(conds_master_single_rep2,1),2);
                    
                    conds_master_single_rep2.stim2_delta(:,1)      = img_vec; % yes (1) /no (2) dots are aligned
                    conds_master_single_rep2.stim2_im_nr(:,1)      = qz_im;   % unique img quiz dot stim number
                    conds_master_single_rep2.stim2_orient_dir(:,1) = ori_dir; % angle of line connecting two quiz dot locations
                    
                    % NO lures in IMG, set to NaN
                    conds_master_single_rep2.is_lure          = NaN(size(conds_master_single_rep2,1),2);          

                
                elseif strcmp(taskClass{:},'img')
                    % take subset of IMG trials, add text prompt and quiz dots loc
                    n_quiz_images       = length(unique(params.stim.ns.imagery_quiz_images));
                    noncatch_trials_idx = conds_master_single_rep2.is_catch==0;
                    shuffle_delta       = shuffle_concat(1:n_quiz_images, (size(conds_master_single_rep2(noncatch_trials_idx,:),1)/n_quiz_images))';
                    
                    % create yes/no overlap vector:
                    % 1=yes overlap, 2=no overlap (50% chance)
                    % col 1 = left stim, col 2 = right stim
                    % left/right stim is not yoked.
                    img_vec                        = NaN(size(conds_master_single_rep2,1),1);
                    img_vec(noncatch_trials_idx,:) = shuffle_delta;
                    
                    % find unique IMG test stim nrs that correspond to the special core stim nr in stim1
                    idx_c     = NaN(size(conds_master_single_rep2,1),1); % <--- EK start here! (what should dim be here? 10 or 14?)
                    [~,idx_c_tmp] = ismember(conds_master_single_rep2.stim_nr_left(noncatch_trials_idx),params.stim.ns.unique_im_nrs_specialcore); % central image only
                    idx_img   = reshape(params.stim.ns.unique_im_nrs_img_test,length(params.stim.ns.imagery_quiz_images),[])'; % 20x15 matrix: 10 yes + 10 no overlap test images for each of the 15 special core images
                    idx_img_yes = idx_img(:,params.stim.ns.imagery_quiz_images==1); % Quiz dots overlap (1=yes) or not (2=no)
                    idx_img_no  = idx_img(:,params.stim.ns.imagery_quiz_images==2); % Quiz dots overlap (1=yes) or not (2=no)
                    
                    assert(isequal(size(idx_img_yes),size(idx_img_no)))
                    assert(all(ismember(unique(params.stim.ns.unique_im_nrs_specialcore(idx_c_tmp)), params.stim.ns.imagery_quiz_dot_specialcore_stim_nr)))
                    
                    idx_c(noncatch_trials_idx==1) = idx_c_tmp; clear idx_c_tmp
                    
                    % To avoid that the table becomes massive, we provide
                    % one number that tells us something about two quiz dots:
                    % we populate the stim2_orient column with the "orientation" 
                    % of the imagery quiz dots test images, which is the
                    % the orientation of the line connecting the two
                    % imagery quiz dots dots, relative to 0 deg = 12 o'clock. 
                    % when quiz dots overlap, their location is within the 
                    % outline of the object in the scene.
                    ori_img_yes = params.stim.ns.imagery_quiz_dot_orient_deg(:,params.stim.ns.imagery_quiz_images==1,1,:); % 8 locations x 10 yes overlap quiz images x 2 dots
                    ori_img_no  = params.stim.ns.imagery_quiz_dot_orient_deg(:,params.stim.ns.imagery_quiz_images==2,1,:); % 8 locations x 10 no overlap quiz images x 2 dots
                    
                    
                    qz_im           = NaN(size(img_vec)); % note: with catch trials
                    ori_dir         = NaN(size(img_vec)); % note: with catch trials
                    yes_trial_idx_c = false(size(img_vec,1),1);
                    no_trial_idx_c  = false(size(img_vec,1),1);
                    
                    
                    % convert to 4 logical vectors, for l/r stim x yes/no overlap
                    yes_trial_idx_c(img_vec(:,1)==1)=true; % Quiz dots overlap (1=yes) or not (2=no) 
                    no_trial_idx_c(img_vec(:,1)==2)=true; % Quiz dots overlap (1=yes) or not (2=no) 
                    
                    % randomly sample from the 10 exemplar quiz dot images (WITHOUT replacement)
                    randomly_selected_yes_test_images_c  = randsample(size(idx_img_yes,2),sum(yes_trial_idx_c), false);  % left stim yes 
                    randomly_selected_no_test_images_c   = randsample(size(idx_img_no,2),sum(no_trial_idx_c), false); % left stim no 
                    
                    yes_counter = 1; no_counter = 1;
                    for ll = 1:length(idx_c) % loop over unique images on the left
                        if yes_trial_idx_c(ll)==1
                            qz_im(ll)   = idx_img_yes(idx_c(ll),randomly_selected_yes_test_images_c(yes_counter));
                            ori_dir(ll) = squeeze(ori_img_yes(idx_c(ll),randomly_selected_yes_test_images_c(yes_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left col using first angle
                            yes_counter   = yes_counter+1;
                        elseif no_trial_idx_c(ll)==1
                            qz_im(ll)   = idx_img_no(idx_c(ll),randomly_selected_no_test_images_c(no_counter));
                            ori_dir(ll) = squeeze(ori_img_no(idx_c(ll),randomly_selected_no_test_images_c(no_counter),1)); % angle in deg for dot 1 and 2 (0 deg = 12 o'clock) -- note left col, using first angle
                            no_counter    = no_counter+1;
                        elseif isnan(img_vec(ll))
                            qz_im(ll)   = NaN;
                            ori_dir(ll) = NaN;
                        else
                            error('[%s]: Imagery test image must be "yes" or "no".',mfilename);
                        end
                    end
                    
                    assert(length(unique(qz_im(:)))==length(qz_im(:)))
                    assert(all(ismember(qz_im(yes_trial_idx_c,1),idx_img_yes))); % check if yes quiz dots image nr matches for central stim
                    assert(all(ismember(qz_im(no_trial_idx_c,1),idx_img_no))); % check if no quiz dots image nr matches for central stim                    
                    assert(isequal(sum(img_vec(:,1)==1),sum(img_vec(:,1)==2))); % equal yes/no overlap?
                    
                    % fill in vectors & condition table
                    conds_master_single_rep2.stim2_im_nr      = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_orient_dir = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_delta      = NaN(size(conds_master_single_rep2,1),2);
                    
                    conds_master_single_rep2.stim2_delta(:,1)      = img_vec; % yes (1) /no (2) dots are aligned
                    conds_master_single_rep2.stim2_im_nr(:,1)      = qz_im;   % unique img quiz dot stim number
                    conds_master_single_rep2.stim2_orient_dir(:,1) = ori_dir; % angle of line connecting two quiz dot locations
                    
                    % NO lures in IMG, set to NaN
                    conds_master_single_rep2.is_lure          = NaN(size(conds_master_single_rep2,1),2);          

                    
                elseif strcmp(taskClass{:},'ltm') % Add ltm pair.
                    
                    % check if catch trials match
                    noncatch_trials_idx = conds_master_single_rep2.is_catch==0;

                    % preallocate space
                    conds_master_single_rep2.stim2_im_nr      = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_orient_dir = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_delta      = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.is_lure     	  = NaN(size(conds_master_single_rep2,1),2);
                    
                    % get paired LTM stim
                    [stim2_nr, stim2_is_lure, stim2_match] = vcd_getLTMTestStim(params, conds_master_single_rep2);
                    assert(isequal(isnan(stim2_nr(:,1)),~noncatch_trials_idx))
                    
                    % insert pairs, if they are lures and correct/incorrect
                    % pairings (stim2_match)
                    conds_master_single_rep2.stim2_im_nr = stim2_nr;
                    conds_master_single_rep2.is_lure     = stim2_is_lure;
                    conds_master_single_rep2.stim2_delta = stim2_match;
                    
                    % We want to test both A->B pairings and B->A pairings.
                    % Make a copy of conditions table
                    conds_single_rep_mergedB = conds_master_single_rep2;
                    
                    % flip order A->B to B->A
                    conds_single_rep_mergedB.stim_nr_left       = conds_master_single_rep2.stim2_im_nr(:,1);
                    conds_single_rep_mergedB.stim2_im_nr(:,1)   = conds_master_single_rep2.stim_nr_left;
                    
                    % central stim, so right stim loc = nan
                    assert(all(isnan(conds_single_rep_mergedB.stim2_im_nr(:,2))))
                    assert(all(isnan(conds_single_rep_mergedB.stim_nr_right)))
                    
                    % update scene category
                    for ii = 1:size(conds_single_rep_mergedB,1)
                        tmp_stim_c = conds_single_rep_mergedB.stim_nr_left(ii);
                        if isnan(tmp_stim_c) || (tmp_stim_c==0) % if catch trial
                            conds_single_rep_mergedB.stim_nr_left(ii)  = 0;
                            conds_single_rep_mergedB.stim2_im_nr(ii,1) = NaN;
                            % set category info to NaN
                            conds_single_rep_mergedB.super_cat(ii,1)        = NaN;
                            conds_single_rep_mergedB.basic_cat(ii,1)        = NaN;
                            conds_single_rep_mergedB.sub_cat(ii,1)          = NaN;
                            conds_single_rep_mergedB.affordance_cat(ii,1)   = NaN;
                            conds_single_rep_mergedB.super_cat_name{ii,1}   = NaN;
                            conds_single_rep_mergedB.basic_cat_name{ii,1}   = NaN;
                            conds_single_rep_mergedB.sub_cat_name{ii,1}     = NaN;
                            conds_single_rep_mergedB.affordance_name{ii,1}  = NaN;
                        else
                            tmp = vcd('fullinfo', tmp_stim_c); % Get stim info
                            assert(isequal(tmp.stim_class,5)); % must be ns
                            assert(isequal(tmp.stim_loc,3)); % must be center
                            % update category info
                            conds_single_rep_mergedB.super_cat(ii,1)        = tmp.super_cat;
                            conds_single_rep_mergedB.basic_cat(ii,1)        = tmp.basic_cat;
                            conds_single_rep_mergedB.sub_cat(ii,1)          = tmp.sub_cat;
                            conds_single_rep_mergedB.affordance_cat(ii,1)   = tmp.affordance_cat;
                            conds_single_rep_mergedB.super_cat_name(ii,1)   = tmp.super_cat_name;
                            conds_single_rep_mergedB.basic_cat_name(ii,1)   = tmp.basic_cat_name;
                            conds_single_rep_mergedB.sub_cat_name(ii,1)     = tmp.sub_cat_name;
                            conds_single_rep_mergedB.affordance_name(ii,1)  = tmp.affordance_name;
                        end
                    end
                    
                    % update unique_trial_nr
                    conds_single_rep_mergedB.unique_trial_nr = max(conds_single_rep_mergedB.unique_trial_nr) + conds_single_rep_mergedB.unique_trial_nr;
                    
                    % check unique trial number (should continue counting)
                    assert(isequal(conds_single_rep_mergedB.unique_trial_nr, ...
                        [min(conds_single_rep_mergedB.unique_trial_nr):max(conds_single_rep_mergedB.unique_trial_nr)]'));
                    
                    % merge condition tables
                    conds_master_single_rep2 = cat(1,conds_master_single_rep2,conds_single_rep_mergedB);
                    
                    % check unique trial number (should continue counting)
                    assert(isequal(sum(~noncatch_trials_idx)*2, sum(conds_master_single_rep2.is_catch)));
                    assert(isequal(conds_master_single_rep2.unique_trial_nr, ...
                        [min(conds_master_single_rep2.unique_trial_nr):max(conds_master_single_rep2.unique_trial_nr)]'));
                    
                else % fill with nans
                    conds_master_single_rep2.stim2_im_nr      = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_orient_dir = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.stim2_delta      = NaN(size(conds_master_single_rep2,1),2);
                    conds_master_single_rep2.is_lure          = NaN(size(conds_master_single_rep2,1),2);
                end
                
                % Keep track of repeats of unique trials
                conds_master_single_rep2.repeat_nr = rep.*ones(size(conds_master_single_rep2,1),1);

                % Accummulate
                cond_master = [cond_master; conds_master_single_rep2];

                clear conds_master_single_rep conds_master_single_rep2

            end

    end
end

return