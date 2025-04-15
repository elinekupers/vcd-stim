function cond_master = vcd_createConditionMaster(p, cond_table, n_trials_per_block)
% VCD function to create condition master table.
%
%  cond_master = vcd_createConditionMaster(p, cond_table, n_trials_per_block)
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
%  p                    : (struct) stimulus and experimental session parameters
%  unique_im            : (matrix) unique stimulus images x M conditions
%  n_unique_cases       : (int) number of unique cases. NOTE: this
%                           number is NOT the same as size(unique_im,1), as
%                           unique cases only counts the fully-crossed
%                           stim features of interest.
%  n_trials_per_block   : (int) number of trials per block (4 or 8)
%  stimClass            : (str) stimulus super class
%  taskClass            : (str) task super class
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
%   1: {'unique_im_nr'   } unique images (parafoveal stimulus patches: 1 trial occupies 2 rows)
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
%   18: {'is_cued'         } cue status: 0 = uncued, 1 = cued
%   19: {'unique_trial_nr'} unique trial nr
%   20: {'thickening_dir' } thickening direction of fixation dot rim (for spatial attention cue) 1 = left, 2 = right, 3 = both/neutral
%   21: {'stim2_delta'    } for double epoch tasks, what predefined delta between stim 1 and stim 2 did we choose from p.stim.(xx).delta_ref
%   22: {'stim2'          } consequence of 21, so updated stim feature of stim2 (e.g., 80 deg orientation for a stim1: 95 - 15 deg)
%   23: {'ltm_stim_pair'  } for LTM task: each unique image nr is associated with
%                           another stimulus
%   24: {'is_lure'         } for LTM task: XX of the trials we show a lure stimulus
%   25: {'repeat_nr'      } Keep track how many times has this unique image has
%                           been repeated thusfar in the experiment
%
% Written by Eline Kupers Feb 2025 @ UMN


%% Get variables from input table
stimClass = unique(cond_table.stim_class_name);
taskClass = unique(cond_table.task_class_name);
n_unique_cases = length(unique(cond_table.unique_im_nr));

assert(length(stimClass)==1)
assert(length(taskClass)==1)

% Prepare cond_master
cond_master = [];

switch stimClass{:}
    
    case 'gabor'
        
        % Get gabor stim manipulations
        n_contrasts    = length(p.stim.gabor.contrast);
        assert(isequal(length(unique(cond_table.contrast)),n_contrasts));
        
        n_ori_bins     = length(p.stim.gabor.ori_deg);
        assert(isequal(length(unique(cond_table.orient_dir)),n_ori_bins));
        
        if strcmp(taskClass,'fix')
            stimloc_cues   = NaN; % {1:cued,0:uncued, NaN:no cue} or 3??
            n_cued_stimlocs = length(stimloc_cues);
        else
            stimloc_cues   = [1,0]; %{1:cued,0:uncued, NaN:no cue, 3:both sides cued?}
            n_cued_stimlocs = length(stimloc_cues);
        end
        

        % loop over repeats of unique trials
        for rep = 1:p.exp.n_unique_trial_repeats
            
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
                assert(isequal(sort(conds_shuffle0.unique_im_nr,'ascend'),[1:n_unique_cases]')); % Check unique image nr
                assert(isequal(sum(conds_shuffle0.stimloc==1),size(conds_shuffle0,1)/2));  % Check left stim loc
                assert(isequal(sum(conds_shuffle0.stimloc==2),size(conds_shuffle0,1)/2));  % Check right stim loc
                assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),repelem(p.stim.gabor.ori_deg,n_contrasts)')); % Check orientations
                assert(isequal(sort(conds_shuffle0.contrast,'ascend'),repelem(p.stim.gabor.contrast,n_ori_bins)')); % Check contrast
                assert(isequal(sort(conds_shuffle0.gbr_phase,'ascend'),repelem(p.stim.gabor.ph_deg,n_unique_cases/length(p.stim.gabor.ph_deg))')); % Check phase
                assert(isequal(sort(reshape(conds_shuffle0.contrast,n_contrasts,[])),repmat(p.stim.gabor.contrast',1,n_ori_bins))); % Check contrast order
                
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
                
                conds_shuffle1 = vcd_addCatchTrials(conds_shuffle0, 1, cue_dir_for_catch_trial);
                prior_cue_dir_for_catch_trial = cue_dir_for_catch_trial;

                % Vertcat single repeat to create this master table
                conds_single_rep = [conds_single_rep;conds_shuffle1];
                
                clear conds_shuffle0 conds_shuffle1
            end
            
            % Merge trials and add fix cue thickening direction
            conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_single_rep);
            
            clear conds_master_single_rep
            
            % Add WM change
            if strcmp(taskClass{:},'wm')
                n_deltas      = length(p.stim.gabor.delta_from_ref);
                shuffle_delta = shuffle_concat(1:length(p.stim.gabor.delta_from_ref), (n_unique_cases/(n_deltas/2)));
                noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                delta_vec     = NaN(size(conds_single_rep_merged.orient_dir));
                delta_vec(noncatch_trials_idx,:) = reshape(p.stim.gabor.delta_from_ref(shuffle_delta)',[],2);
                orient2_vec = conds_single_rep_merged.orient_dir + delta_vec;
                conds_single_rep_merged.stim2_delta = delta_vec;
                conds_single_rep_merged.stim2       = num2cell(orient2_vec);
            else
                delta_vec   = NaN(size(conds_single_rep_merged.orient_dir));
                orient2_vec = NaN(size(conds_single_rep_merged.orient_dir));
                conds_single_rep_merged.stim2_delta = delta_vec;
                conds_single_rep_merged.stim2       = num2cell(orient2_vec);
            end
            
            % Add IMG text prompt and quiz dots
            if strcmp(taskClass{:},'img')
                n_quiz_images    = length(unique(p.exp.trial.img.quiz_images));
                quiz_img_type_shuffle = [shuffle_concat(1:n_quiz_images, size(conds_single_rep_merged,1)/n_quiz_images); ...
                                         shuffle_concat(1:n_quiz_images, size(conds_single_rep_merged,1)/n_quiz_images)]';
                img_quiz_im_vec  = quiz_img_type_shuffle;
                img_quiz_txt_vec = NaN(size(conds_single_rep_merged,1),1);
                conds_single_rep_merged.img_txt_prompt = num2cell(img_quiz_txt_vec);
                conds_single_rep_merged.stim2          = num2cell(img_quiz_im_vec);
            else
                img_quiz_im_vec   = NaN(size(conds_single_rep_merged,1),2);
                img_quiz_txt_vec  = NaN(size(conds_single_rep_merged,1),1);
                conds_single_rep_merged.stim2          = num2cell(img_quiz_im_vec);
                conds_single_rep_merged.img_txt_prompt = num2cell(img_quiz_txt_vec);
            end
 
            %% TODO Preallocate space for LTM pairing (we'll do this later)
            pair_vec = NaN(size(conds_single_rep_merged.stim2));
            lure_vec = false(size(pair_vec)); % should be boolean
            
            conds_single_rep_merged.ltm_stim_pair = pair_vec;
            conds_single_rep_merged.is_lure        = lure_vec;
            
            % Keep track of the times a unique condition has been repeated
            rep_vec = rep.*ones(size(conds_single_rep_merged,1),1);
            
            conds_single_rep_merged.repeat_nr = rep_vec;
            unique_trial_nr = conds_single_rep_merged.unique_trial_nr + ((rep-1)*max(conds_single_rep_merged.unique_trial_nr));
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

    case 'rdk'
        
        % Get stim manipulations
        n_coh       = length(p.stim.rdk.dots_coherence); % levels of dot coherence
        n_motdir    = length(p.stim.rdk.dots_direction); % number of motion direction bins per hemifield
        
        if strcmp(taskClass,'fix')
            stimloc_cues   = NaN; % {1:cued,0:uncued, NaN:neutral/no cue} or 3??
            n_cued_stimlocs = length(stimloc_cues);
        else
            stimloc_cues = [1,0]; % {1:cued,0:uncued, NaN:neutral/no cue}
            n_cued_stimlocs = length(stimloc_cues);
        end
        
        % Keep track of stim loc cue loc for catch trials.
        prior_cue_dir_for_catch_trial = [];
        
        for rep = 1:p.exp.n_unique_trial_repeats
            
            conds_master_single_rep = [];
            
            cue_vec = []; im_cue_vec_loc1 = [];
            for loc = 1:n_cued_stimlocs % cued vs uncued
                
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
                assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),repelem(p.stim.rdk.dots_direction,n_coh)')); % Check motdir
                assert(isequal(sort(conds_shuffle0.rdk_coherence,'ascend'),repelem(p.stim.rdk.dots_coherence,n_motdir)')); % Check coherence
                assert(isequal(sort(reshape(conds_shuffle0.rdk_coherence,n_coh,[])),repmat(p.stim.rdk.dots_coherence',1,n_motdir))); % Check coh order
                
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
                
                % Create catch trial
                conds_shuffle1 = vcd_addCatchTrials(conds_shuffle0, 1, cue_dir_for_catch_trial);
                prior_cue_dir_for_catch_trial = cue_dir_for_catch_trial;
                
                conds_master_single_rep = [conds_master_single_rep;conds_shuffle1];

                clear conds_shuffle0 conds_shuffle1
            end

            % Merge unique im into trials and add fix cue thickening direction.
            conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep);
            clear conds_master_single_rep

            % Add WM change.
            if strcmp(taskClass,'wm')
                n_deltas      = length(p.stim.rdk.delta_from_ref);
                shuffle_delta = shuffle_concat(1:length(p.stim.rdk.delta_from_ref), (n_unique_cases/(n_deltas/2)));
                noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                delta_vec     = NaN(size(conds_single_rep_merged.orient_dir));
                delta_vec(noncatch_trials_idx,:) = reshape(p.stim.rdk.delta_from_ref(shuffle_delta)',[],2);
                motdir2       = conds_single_rep_merged.orient_dir + delta_vec;
                conds_single_rep_merged.stim2_delta = delta_vec;
                conds_single_rep_merged.stim2 = num2cell(motdir2);
            else
                delta_vec     = NaN(size(conds_single_rep_merged.orient_dir));
                motdir2       = NaN(size(conds_single_rep_merged.orient_dir));
                conds_single_rep_merged.stim2_delta = delta_vec;
                conds_single_rep_merged.stim2 = num2cell(motdir2);
            end

            
            % Add IMG text prompt and quiz dots
            if strcmp(taskClass{:},'img')
                n_quiz_images    = length(unique(p.exp.trial.img.quiz_images));
                quiz_img_type_shuffle = [shuffle_concat(1:n_quiz_images, size(conds_single_rep_merged,1)/n_quiz_images); ...
                    shuffle_concat(1:n_quiz_images, size(conds_single_rep_merged,1)/n_quiz_images)]';
                img_quiz_im_vec  = quiz_img_type_shuffle;
                img_quiz_txt_vec  = NaN(size(conds_single_rep_merged,1),1);
                conds_single_rep_merged.img_txt_prompt = num2cell(img_quiz_txt_vec);
                conds_single_rep_merged.stim2          = num2cell(img_quiz_im_vec);
            else
                img_quiz_im_vec   = NaN(size(conds_single_rep_merged,1),2);
                img_quiz_txt_vec  = NaN(size(conds_single_rep_merged,1),1);
                conds_single_rep_merged.stim2          = num2cell(img_quiz_im_vec);
                conds_single_rep_merged.img_txt_prompt = num2cell(img_quiz_txt_vec);
            end

            
            % Add ltm pair.
            %                         if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
            % %                             pair_match = [];
            % %                             delta_vec = shuffle_concat(p.stim.rdk.ltm_pairs, (n_unique_cases/(n_deltas/2)));
            % %                             pair_vec = conds_master_single_rep3(:,7) + delta_vec';
            %                         else
            conds_single_rep_merged.ltm_stim_pair = NaN(size(conds_single_rep_merged.stim2_delta));
            conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged.ltm_stim_pair)); % should be boolean
            
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
        
        
    case 'dot'
        
        % Get stim manipulations
        n_dot_loc  = size(p.stim.dot.loc_deg,2);
        
        if strcmp(taskClass,'fix')
            stimloc_cues = NaN;  % {1:cued,0:uncued, NaN:neutral/no cue} or 3??
            n_cued_stimlocs = length(stimloc_cues);
        else
            stimloc_cues = [1,0];  % {1:cued,0:uncued, NaN:neutral/no cue}
            n_cued_stimlocs = length(stimloc_cues);
        end
        
        n_unique_cases = n_dot_loc;
        
        % Keep track of stim loc cue loc for catch trials.
        prior_cue_dir_for_catch_trial = [];
        
        for rep = 1:p.exp.n_unique_trial_repeats
            
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
                assert(isequal(sort(conds_shuffle0.orient_dir,'ascend'),repelem(p.stim.dot.ang_deg,1)')); % Check dot angle
                
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
                
                conds_shuffle1 = vcd_addCatchTrials(conds_shuffle0, 1, cue_dir_for_catch_trial);
                prior_cue_dir_for_catch_trial = cue_dir_for_catch_trial;

                % Vert cat tables
                conds_master_single_rep = [conds_master_single_rep;conds_shuffle1];
                
                clear conds_shuffle0 conds_shuffle1
            end
            
            % Merge trials and add fix cue thickening direction
            conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep);
            clear conds_master_single_rep
            
            % Add WM change.
            if strcmp(taskClass{:},'wm')
                n_deltas     = length(p.stim.dot.delta_from_ref);
                shuffle_delta = shuffle_concat(1:length(p.stim.dot.delta_from_ref), (n_unique_cases/(n_deltas/2)));
                noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                delta_vec     = NaN(size(conds_single_rep_merged.orient_dir));
                delta_vec(noncatch_trials_idx,:) = reshape(p.stim.dot.delta_from_ref(shuffle_delta)',[],2);
                loc_deg2       = conds_single_rep_merged.orient_dir + delta_vec;
            else
                delta_vec   = NaN(size(conds_single_rep_merged.orient_dir));
                loc_deg2    = NaN(size(conds_single_rep_merged.orient_dir));
            end
            conds_single_rep_merged.stim2_delta = delta_vec;
            conds_single_rep_merged.stim2       = num2cell(loc_deg2);

            % Add IMG text prompt and quiz dots
            if strcmp(taskClass{:},'img')
               n_quiz_images    = length(unique(p.exp.trial.img.quiz_images));
                quiz_img_type_shuffle = [shuffle_concat(1:n_quiz_images, size(conds_single_rep_merged,1)/n_quiz_images); ...
                    shuffle_concat(1:n_quiz_images, size(conds_single_rep_merged,1)/n_quiz_images)]';
                img_quiz_im_vec   = quiz_img_type_shuffle;
                img_quiz_txt_vec  = NaN(size(conds_single_rep_merged,1),1);
                conds_single_rep_merged.img_txt_prompt = num2cell(img_quiz_txt_vec);
                conds_single_rep_merged.stim2          = num2cell(img_quiz_im_vec);
            else
                img_quiz_im_vec   = NaN(size(conds_single_rep_merged,1),2);
                img_quiz_txt_vec  = NaN(size(conds_single_rep_merged,1),1);
                conds_single_rep_merged.stim2          = num2cell(img_quiz_im_vec);
                conds_single_rep_merged.img_txt_prompt = num2cell(img_quiz_txt_vec);
            end

            
            % Add ltm pair.
            %                         if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
            % %                             pair_match = [];
            % %                             delta_vec = shuffle_concat(p.stim.obj.delta_from_ref, (n_unique_cases/(n_deltas/2)));
            % %                             pair_vec = conds_master_single_rep3(:,7) + delta_vec';
            %                         else
            conds_single_rep_merged.ltm_stim_pair = NaN(size(conds_single_rep_merged.stim2_delta));
            conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged.ltm_stim_pair)); % should be boolean
            %                         end
            
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
        
        
    case 'obj'
        
        % Get stim manipulations
        n_super_cat = length(p.stim.obj.super_cat);
        stimloc_vec = [1,2]; %{'left','right'}
        
        % UNIQUE COBJ array dims: 8 basic categories
        basic_cat_vec = [];
        for ni = 1:n_super_cat
            tmp = unique(p.stim.obj.basic_cat{ni}, 'stable');
            n_basic_cat(ni) = length(tmp);
            
            basic_tmp = [];
            for nj = 1:length(tmp)
                basic_tmp = cat(1,basic_tmp,nj.*(arrayfun(@(x) strcmp(x, tmp(nj)),p.stim.obj.basic_cat{ni})));
            end
            basic_cat_vec = cat(2, basic_cat_vec, sum(basic_tmp,1));
        end
        
        % Get super and sub category info
        super_cat_vec = []; sub_cat_vec = [];
        for ii = 1:length(n_basic_cat)
            n_sub_cat(ii) = length(p.stim.obj.sub_cat{ii});
            super_cat_vec = cat(2, super_cat_vec, repelem(ii,n_sub_cat(ii)));
            sub_cat_vec = cat(2, sub_cat_vec, 1:n_sub_cat(ii));
        end
        
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
        for rep = 1:p.exp.n_unique_trial_repeats

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
            assert(isequal(1:length(super_cat_vec),unique(basic_cat_shuffle_idx)))
            
            % Get cueing vector
            cue_vec = []; im_cue_vec_loc1 = [];
            for loc = 1:n_cued_stimlocs % cued vs uncued
                
                % Make sure every 8 trials, subject sees images from 4 super
                % classes
                super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                shuffled_super_cat    = super_cat_vec(super_cat_shuffle_idx);
                for jj = 1:n_trials_per_block:(length(shuffled_super_cat))
                    
                    while ~isempty(setdiff([1:n_super_cat],unique(shuffled_super_cat(jj:(jj+n_trials_per_block-1)),'legacy')))
                        super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                        shuffled_super_cat = super_cat_vec(super_cat_shuffle_idx);
                        
                    end
                end
                % Check if n_sub_cat still holds
                assert(isequal(histcounts(shuffled_super_cat),n_sub_cat))
                
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
                assert(isequal(sort(conds_shuffle1.super_cat_name),sort(repelem(p.stim.obj.super_cat,n_sub_cat)'))); % Check supercat
                
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
                
                conds_shuffle2 = vcd_addCatchTrials(conds_shuffle1, 1, cue_dir_for_catch_trial);
                prior_cue_dir_for_catch_trial = cue_dir_for_catch_trial;

                % Vert cat tables
                conds_master_single_rep = [conds_master_single_rep;conds_shuffle2];
                
                clear conds_shuffle0 conds_shuffle1 conds_shuffle2
            end
            
            % Merge trials and add fix cue thickening direction
            conds_single_rep_merged = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep);
            clear conds_master_single_rep;
            
            % Add WM change.
            if strcmp(taskClass{:},'wm')
                n_deltas      = length(p.stim.obj.delta_from_ref);
                shuffle_delta = shuffle_concat(1:length(p.stim.obj.delta_from_ref), (n_unique_cases/(n_deltas/2)));
                noncatch_trials_idx = ~isnan(conds_single_rep_merged.orient_dir(:,1));
                delta_vec     = NaN(size(conds_single_rep_merged.orient_dir));
                delta_vec(noncatch_trials_idx,:) = reshape(p.stim.obj.delta_from_ref(shuffle_delta)',[],2);
                facing_dir2   = conds_single_rep_merged.orient_dir + delta_vec;
                conds_single_rep_merged.stim2_delta = delta_vec;
                conds_single_rep_merged.stim2       = num2cell(facing_dir2);
            else
                delta_vec     = NaN(size(conds_single_rep_merged.orient_dir));
                facing_dir2   = NaN(size(conds_single_rep_merged.orient_dir));
                conds_single_rep_merged.stim2_delta = delta_vec;
                conds_single_rep_merged.stim2       = num2cell(facing_dir2);
            end

            
            % Add IMG text prompt and quiz dots
            if strcmp(taskClass{:},'img')
                n_quiz_images    = length(unique(p.exp.trial.img.quiz_images));
                quiz_img_type_shuffle = [shuffle_concat(1:n_quiz_images, size(conds_single_rep_merged,1)/n_quiz_images); ...
                    shuffle_concat(1:n_quiz_images, size(conds_single_rep_merged,1)/n_quiz_images)]';
                img_quiz_im_vec  = quiz_img_type_shuffle;
                img_quiz_txt_vec  = NaN(size(conds_single_rep_merged,1),1);
                conds_single_rep_merged.img_txt_prompt = num2cell(img_quiz_txt_vec);
                conds_single_rep_merged.stim2          = num2cell(img_quiz_im_vec);
            else
                img_quiz_im_vec   = NaN(size(conds_single_rep_merged,1),2);
                img_quiz_txt_vec  = NaN(size(conds_single_rep_merged,1),1);
                conds_single_rep_merged.stim2          = num2cell(img_quiz_im_vec);
                conds_single_rep_merged.img_txt_prompt = num2cell(img_quiz_txt_vec);
            end
            
            % Add ltm pair.
            %                         if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
            % %                             pair_match = [];
            % %                             delta_vec = shuffle_concat(p.stim.obj.delta_from_ref, (n_unique_cases/(n_deltas/2)));
            % %                             pair_vec = conds_master_single_rep3(:,7) + delta_vec';
            %                         else
            conds_single_rep_merged.ltm_stim_pair = NaN(size(conds_single_rep_merged.stim2_delta));
            conds_single_rep_merged.is_lure = false(size(conds_single_rep_merged.ltm_stim_pair));
            
            % Keep track of repeat
            conds_single_rep_merged.repeat_nr = rep.*ones(size(conds_single_rep_merged,1),1);
            
            % Accummulate
            cond_master = [cond_master; conds_single_rep_merged];
        end
        
        % thickening direction doesn't have to match
        % between left and right, or do they??...
        if ~strcmp(taskClass,'fix')
            assert(isequal(sum(cond_master.is_cued==1),sum(cond_master.is_cued==2)))
        end
        
    case 'ns'
        
        % Get stim manipulations
        stimloc_cues    = 3; %{1:'cued',0:'uncued',3:'bothcued'}; % should we make bothcued=2?
        n_cued_stimlocs = length(stimloc_cues);
        n_super_cat     = length(p.stim.ns.super_cat);
        
       
        % Get super and sub category info
        n_basic_cat = [];   n_sub_cat = [];
        super_cat_vec = []; sub_cat_vec = [];
        for ni = 1:n_super_cat
            n_basic_cat(ni) = length(p.stim.ns.basic_cat{ni});
            
            for nj = 1:n_basic_cat(ni)

                n_sub_cat(ni,nj) = length(p.stim.ns.sub_cat{ni,nj});
                
                super_cat_vec = cat(2, super_cat_vec, repelem(ni,n_sub_cat(ni,nj)));
%                 basic_cat_vec0 = cat(2, basic_cat_vec0, repelem(nj, n_sub_cat(ni,nj)));
                sub_cat_vec   = cat(2, sub_cat_vec, [1:n_sub_cat(ni,nj)]);

            end
        end
        

        
        for rep = 1:p.exp.n_unique_trial_repeats
            
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
            assert(isequal(histcounts(shuffled_super_cat),sum(n_sub_cat,2)'))
            
            % shuffle basic category (so shuffle every other image)
%             basic_cat_vec = cond_table.basic_cat;
            basic_cat_shuffle_idx = [];
            for ii = 1:length(n_basic_cat)
                % create [ 0 0 2 2 4 4 ] vector
                curr_im = repelem(length(basic_cat_shuffle_idx): 2 : length(basic_cat_shuffle_idx)+2*n_basic_cat(ii),n_basic_cat(ii));
                % add shuffled [1,2,1,2,1,2] vector
                basic_cat_shuffle_idx = cat(2, basic_cat_shuffle_idx, curr_im + shuffle_concat([1:n_basic_cat(ii)],n_sub_cat(ii,1)));
            end

            
            % NOW SHUFFLE!
            conds_shuffle0 = cond_table(basic_cat_shuffle_idx,:); % <-- shuffle based on basic category
            conds_master_single_rep = conds_shuffle0(super_cat_shuffle_idx,:); % <-- shuffle based on super category, last one to give it priority
            
            % Add catch trial vector
            conds_master_single_rep.is_catch = false(size(conds_master_single_rep,1),1);

            % Add catch trials
            conds_master_single_rep2 = vcd_addCatchTrials(conds_master_single_rep, 1, []);
            
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
            fnames = {'orient_dir','contrast','gbr_phase','rdk_coherence','super_cat','basic_cat','sub_cat','super_cat_name','basic_cat_name','sub_cat_name','stim_class_name','is_in_img'};
            
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
                n_changes = length(p.stim.ns.change_im);
                change_blindness_vec = shuffle_concat(1:length(p.stim.ns.change_im), ceil(size(conds_master_single_rep2.unique_trial_nr,1)/n_changes))';
                noncatch_trials_idx  = ~conds_master_single_rep2.is_catch(:,1);
                delta_vec            = NaN(length(noncatch_trials_idx),2);
                change_blindness_vec = change_blindness_vec(1:sum(noncatch_trials_idx)); % check length 
                delta_vec(noncatch_trials_idx,1) = change_blindness_vec;
                cblind_im_name       = num2cell(NaN(length(noncatch_trials_idx),2)); % add column of nans to match with l/r stim trials
                cblind_im_name(noncatch_trials_idx,1) = p.stim.ns.change_im(change_blindness_vec)';
                conds_master_single_rep2.stim2_delta = delta_vec;
                conds_master_single_rep2.stim2 = num2cell(cblind_im_name);
            else
                delta_vec            = NaN(size(conds_master_single_rep2,1),2);
                cblind_im_name       = delta_vec;
                conds_master_single_rep2.stim2_delta = delta_vec;
                conds_master_single_rep2.stim2 = num2cell(cblind_im_name);
            end

            
            % take subset of IMG trials, add text prompt and quiz dots loc
            if strcmp(taskClass{:},'img')
                n_quiz_images    = length(unique(p.exp.trial.img.quiz_images));
                quiz_img_type_shuffle = shuffle_concat(1:n_quiz_images, ceil(size(conds_master_single_rep2,1)/n_quiz_images))'; ...
                quiz_img_type_shuffle = quiz_img_type_shuffle(1:size(conds_master_single_rep2,1));
                quiz_img_type_shuffle = [quiz_img_type_shuffle, NaN(size(conds_master_single_rep2,1),1)];
                img_quiz_im_vec  = quiz_img_type_shuffle;
                img_quiz_txt_vec  = NaN(size(conds_master_single_rep2,1),1);
                conds_master_single_rep2.img_txt_prompt = num2cell(img_quiz_txt_vec);
                conds_master_single_rep2.stim2          = num2cell(img_quiz_im_vec);
            else
                conds_master_single_rep2.stim2           = num2cell(NaN(size(conds_master_single_rep2,1),2));
                conds_master_single_rep2.img_txt_prompt  = num2cell(NaN(size(conds_master_single_rep2,1),1));
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
            conds_master_single_rep2.ltm_stim_pair = NaN(size(conds_master_single_rep2.stim2_delta));
            conds_master_single_rep2.is_lure = false(size(conds_master_single_rep2.ltm_stim_pair));
            
            % Keep track of repeat
            conds_master_single_rep2.repeat_nr = rep.*ones(size(conds_master_single_rep2,1),1);
            
            
            % Accummulate
            cond_master = [cond_master; conds_master_single_rep2];
            
            clear conds_master_single_rep conds_master_single_rep2
            
        end
        
end

return