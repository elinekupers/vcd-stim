function cond_master = vcd_shuffleUniqueImages(p, unique_im, stimClass, use_fix_flag)


cond_master = [];

switch stimClass
    
    case 'gabor'
        
        % Get gabor stim manipulations
        n_contrasts    = length(p.stim.gabor.contrast);
        n_ori_bins     = length(p.stim.gabor.ori_deg);
        loc_stim       = [1,2]; %{'L','R'};
        n_stimloc_cues = length(loc_stim);
        n_unique_cases = n_contrasts*n_ori_bins;
        
        for rep = 1:p.exp.n_unique_trial_repeats
            
            conds_single_rep = [];
            
            % If this is a fixation task
            if use_fix_flag
                
                % shuffle orientation every 8 trials
                shuffle_ori = shuffle_concat([1:n_ori_bins],n_contrasts*2);
                shuffle_ori = shuffle_ori' + repmat(repelem([0:n_ori_bins:(n_unique_cases-1)],n_ori_bins),1,2)';
                
                % shuffle contrasts every 3 trials
                shuffle_c = shuffle_concat(1:n_contrasts,(n_unique_cases/n_contrasts)*2);
                
                % Reshape number of unique image cases, such that we can
                % allocate unique image nrs the combinations of contrast and orientation
                case_vec  = reshape(1:size(unique_im,1),[],3);
                case_vec  = case_vec';
                case_vec  = case_vec(:);
                
                %
                shuffled_c = case_vec(shuffle_c + repmat(repelem([0:n_contrasts:(n_unique_cases-1)],n_contrasts),1,2),:);
                
                % SHUFFLE UNIQUE IMAGES!!
                conds_shuffle0 = unique_im(shuffle_ori,:); % <-- first shuffle based on unique nr of orientations
                conds_shuffle0 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of contrasts, such that new trial order prioritized contrast order
                
                % we don't have spatial cues in fixation condition
                cue_vec = NaN(size(unique_im,1),1);
                
                % Horz cat cue vec to shuffled condition vector
                conds_shuffle1 = [conds_shuffle0, cue_vec];
                
                % Vertcat single repeat to create this master table
                conds_single_rep = [conds_single_rep;conds_shuffle1];
                
            else % non-fixation tasks
                
                for loc = 1:n_stimloc_cues % cued vs uncued
                    
                    % shuffle orientation every 8 trials
                    shuffle_ori = shuffle_concat([1:n_ori_bins],n_contrasts);
                    shuffle_ori = shuffle_ori' + repelem([0:n_ori_bins:(n_unique_cases-1)],n_ori_bins)';
                    assert(isequal(unique(shuffle_ori),[1:(n_ori_bins*n_contrasts)]'))
                    
                    % shuffle contrasts every 3 trials
                    shuffle_c = shuffle_concat(1:n_contrasts,n_unique_cases/n_contrasts);
                    case_vec = reshape(1:n_unique_cases,[],3);
                    case_vec = case_vec';
                    case_vec = case_vec(:);
                    
                    shuffled_c = case_vec(shuffle_c + repelem([0:n_contrasts:(n_unique_cases-1)],n_contrasts),:);
                    assert(isequal(unique(shuffle_c + repelem([0:n_contrasts:(n_unique_cases-1)],n_contrasts)),[1:n_unique_cases]))
                    
                    % TRICKY STUFF: WE PRIORITISE SHUFFLE ORDER OF UNIQUE IMAGES
                    conds_shuffle0 = unique_im(shuffle_ori,:); % <-- first shuffle based on unique nr of orientations
                    conds_shuffle0 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of contrasts, such that new trial order prioritized contrast order
                    
                    % Add spatial cueing vector
                    % MORE TRICKY STUFF: we need to ensure each unique image is
                    % presented the same number of times while being cued
                    % or uncued. To keep cueing balanced, we randomly cue
                    % half of the unique images on the left (and thus create
                    % a set of uncued unique images on the right), and then
                    % ensure those uncued images on the right still get
                    % cued in a second repeat.
                    
                    if loc == 1 % left stim
                        % generate random cueing vector
                        cue_vec = shuffle_concat(stimloc_cues,n_unique_cases/n_stimloc_cues)';
                        
                        % Keep a copy for next round (i.e., stim loc 2)
                        cue_vec1 = cue_vec;
                        im_cue_vec1 = conds_shuffle0(find(cue_vec1),1);
                        
                    elseif loc == 2 % right stim
                        % We create an anti-cue vector of the loc1 cue_vec,
                        % such that we cue those unique images that haven't
                        % been cued left.
                        
                        % preallocate space
                        cue_vec_anti = zeros(size(cue_vec));
                        
                        % Get copy, these unique images were cued previously
                        prev_cue_vec = im_cue_vec1;
                        
                        % These unique images still need to be cued:
                        im_to_be_cued = setdiff([1:n_unique_cases],prev_cue_vec);
                        [~,im_to_be_cued_i] = intersect(conds_shuffle0(:,1),im_to_be_cued','stable');
                        
                        % check if this is true
                        assert(isequal(sort(conds_shuffle0(im_to_be_cued_i)),im_to_be_cued'))
                        
                        cue_vec_anti(im_to_be_cued_i) = 1;
                        cue_vec = cue_vec_anti;
                    end
                    
                    % Horz cat cue vec to shuffled condition vector
                    conds_shuffle1 = [conds_shuffle0, cue_vec];
                    
                    % Vertcat single repeat to create this master table
                    conds_single_rep = [conds_single_rep;conds_shuffle1];
                    
                end
            end
            
            
            %% reorganize column order
            conds_master_single_rep2 = NaN(size(conds_single_rep));
            
            conds_master_single_rep2(:,1) = conds_single_rep(:,1); % unique im nr
            conds_master_single_rep2(:,2) = conds_single_rep(:,2); % stim loc
            conds_master_single_rep2(:,3) = conds_single_rep(:,6); % cue stat
            conds_master_single_rep2(:,4) = conds_single_rep(:,3); % orient
            conds_master_single_rep2(:,5) = conds_single_rep(:,4); % contrast
            conds_master_single_rep2(:,6) = conds_single_rep(:,5); % phase
            
            % Merge trials and add fix cue thickening direction
            conds_master_single_rep3 = vcd_mergeUniqueImageIntoTrials(conds_master_single_rep2,use_fix_flag);
            
            % Add WM change
            if strcmp(tasks(task_crossings(curr_task)).name,'wm')
                n_deltas = length(p.stim.gabor.delta_from_ref);
                delta_vec = shuffle_concat(p.stim.gabor.delta_from_ref, (n_unique_cases/(n_deltas/2)));
                orient2 = conds_master_single_rep3(:,6) + delta_vec';
            else
                orient2 = NaN(size(conds_master_single_rep3,1),1);
            end
            
            %% %%% TODO: Add ltm pair!!!
            %                         if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
            % %                             pair_match = [];
            % %                             delta_vec = shuffle_concat(p.stim.rdk.ltm_pairs, (n_unique_cases/(n_deltas/2)));
            % %                             pair_vec = conds_master_single_rep3(:,7) + delta_vec';
            %                         else
            pair_vec = NaN(size(conds_master_single_rep3,1),1);
            lure_vec = pair_vec; % should be boolean
            %                         end
            %%
            
            % Keep track of the times a unique condition has been repeated
            rep_vec = rep.*ones(size(conds_master_single_rep3,1),1);
            
            conds_master_single_rep3 = [conds_master_single_rep3, orient2, pair_vec, lure_vec, rep_vec];
            
            % Accummulate
            cond_master = [cond_master; conds_master_single_rep3];
            
            % thickening direction doesn't have to match
            % between left and right, from my
            % understanding...
            % assert(isequal(sum(cond_master(:,2)==1),sum(cond_master(:,2)==2)))
            
            % Clean up
            clear  conds_master_single_rep2 conds_master_single_rep3
            
        end
        
        
    case 'rdk'   
        
        % Get stim manipulations
        n_coh       = length(p.stim.rdk.dots_coherence);
        n_motdir    = length(p.stim.rdk.dots_direction); % number of motion direction bins per hemifield
        loc_stim    = [1,2];%{'L','R'};
        n_stim_loc  = length(loc_stim);
        
        % gabors array dims: 8 orientations x 4 phases x 3 contrasts x 4 deltas
        n_unique_cases = n_coh*n_motdir;

        for rep = 1:n_repeats
            
            conds_master_single_rep = [];
            if use_fix_flag
                
                % shuffle orientation every 8 trials
                shuffle_motdir = shuffle_concat([1:n_motdir],n_coh*2);
                shuffle_motdir = shuffle_motdir' + repmat(repelem([0:n_motdir:(n_unique_cases-1)],n_motdir),1,2)';
                
                % shuffle coherence every 3 trials
                shuffle_c = shuffle_concat(1:n_coh,(n_unique_cases/n_coh)*2);
                case_vec = reshape(1:size(unique_im,1),[],3);
                case_vec = case_vec';
                case_vec = case_vec(:);
                shuffled_c = case_vec(shuffle_c + repmat(repelem([0:n_coh:(n_unique_cases-1)],n_coh),1,2),:);
                
                % TRICKY STUFF: WE PRIORITISE SHUFFLE ORDER OF UNIQUE IMAGES              
                conds_shuffle0 = unique_im(shuffle_motdir,:); % <-- first shuffle based on unique nr of motion directions
                conds_shuffle0 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of coherence levels, such that new trial order prioritizes coherence order
                
                % we don't have spatial cues in fixation condition
                cue_vec = NaN(size(unique_im,1),1);
                
                % Horz cat cue vec to shuffled condition vector
                conds_shuffle1 = [conds_shuffle0, cue_vec];
                
                % Vertcat single repeat to create this master table
                conds_master_single_rep = [conds_master_single_rep;conds_shuffle1];
                
            else
                
                for loc = 1:n_stimloc_cues % cued vs uncued
                    
                    % shuffle orientation every 8 trials
                    shuffle_motdir = shuffle_concat([1:n_motdir],n_coh);
                    shuffle_motdir = shuffle_motdir' + repelem([0:n_motdir:(n_unique_cases-1)],n_motdir)';
                    
                    % shuffle coherence every 3 trials
                    shuffle_c = shuffle_concat(1:n_coh,n_unique_cases/n_coh);
                    case_vec = reshape(1:n_unique_cases,[],3);
                    case_vec = case_vec';
                    case_vec = case_vec(:);
                    shuffled_c = case_vec(shuffle_c + repelem([0:n_coh:(n_unique_cases-1)],n_coh),:);
                    
                    conds_shuffle0 = unique_im(shuffle_motdir,:); % <-- first shuffle based on unique nr of orientations
                    conds_shuffle0 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of coh, such that new trial order prioritized coh order
                    
                    if loc == 1
                        cue_vec = shuffle_concat(stimloc_cues,n_unique_cases/n_stimloc_cues)';
                        
                        % store for next round
                        cue_vec1 = cue_vec;
                        im_cue_vec1 = conds_shuffle0(find(cue_vec1),1);
                        
                    elseif loc == 2
                        cue_vec_anti = zeros(size(cue_vec));
                        % These unique images were cued previously
                        prev_cue_vec = im_cue_vec1;
                        % These unique images still need to be
                        % cued to keep cueing balanced
                        im_to_be_cued = setdiff([1:n_unique_cases],prev_cue_vec);
                        [~,im_to_be_cued_i] = intersect(conds_shuffle0(:,1),im_to_be_cued','stable');
                        assert(isequal(sort(conds_shuffle0(im_to_be_cued_i)),im_to_be_cued'))
                        cue_vec_anti(im_to_be_cued_i) = 1;
                        cue_vec = cue_vec_anti;
                    end
                    
                    conds_shuffle1 = [conds_shuffle0, cue_vec];
                    
                    conds_master_single_rep = [conds_master_single_rep;conds_shuffle1];
                end
                
                clear conds_shuffle0 conds_shuffle1 conds_shuffle2
            end
            
            % Reorganize column order
            conds_master_single_rep2 = NaN(size(conds_master_single_rep));
            
            conds_master_single_rep2(:,1) = conds_master_single_rep(:,1); % unique im nr
            conds_master_single_rep2(:,2) = conds_master_single_rep(:,2); % stim loc
            conds_master_single_rep2(:,3) = conds_master_single_rep(:,5); % cue stat
            conds_master_single_rep2(:,4) = conds_master_single_rep(:,3); % mot dir
            conds_master_single_rep2(:,5) = conds_master_single_rep(:,4); % coh
            
            % Merge unique im into trials and add fix cue thickening direction.
            conds_master_single_rep3 = create_condition_master_trials(conds_master_single_rep2,use_fix_flag);
            
            % Add WM change.
            if strcmp(tasks(task_crossings(curr_task)).name,'wm')
                n_deltas = length(p.stim.rdk.delta_from_ref);
                delta_vec = shuffle_concat(p.stim.rdk.delta_from_ref, (n_unique_cases/(n_deltas/2)));
                motdir2 = conds_master_single_rep3(:,6) + delta_vec';
            else
                motdir2 = NaN(size(conds_master_single_rep3,1),1);
            end
            % Add ltm pair.
            %                         if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
            % %                             pair_match = [];
            % %                             delta_vec = shuffle_concat(p.stim.rdk.ltm_pairs, (n_unique_cases/(n_deltas/2)));
            % %                             pair_vec = conds_master_single_rep3(:,7) + delta_vec';
            %                         else
            pair_vec = NaN(size(conds_master_single_rep3,1),1);
            lure_vec = pair_vec; % should be boolean
            %                         end
            
            % Keep track of repeat
            rep_vec = rep.*ones(size(conds_master_single_rep3,1),1);
            conds_master_single_rep3 = [conds_master_single_rep3, motdir2, pair_vec, lure_vec, rep_vec];
            
            % Accummulate
            cond_master = [cond_master; conds_master_single_rep3];
            
            % Clean up
            clear  conds_master_single_rep2 conds_master_single_rep3
        end
end

return