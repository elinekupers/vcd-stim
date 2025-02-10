function all_trials = vcd_makeTrials(p)


% We need to define all possible combinations
% Depending on the session, run, and subject, we need to define the
% combinations we want to show, converted into miniblocks, where trials are
% sampling the manipulations in a pseudo-random way.

n_repeats = p.exp.n_unique_trial_repeats;

%%
% In each run, we have manipulations that we prioritize to fully sample,
% otherwise it is difficult to compare conditions (e.g., we want to sample
% all contrast levels within the run).
priority_stim_manip = struct('name',{},'priority',{},'other',{});

priority_stim_manip(1).name     = {'gabor'};
priority_stim_manip(1).priority = {'contrast','delta_ref'}; % Priority manipulation
priority_stim_manip(1).other    = {'ori_bin'};                 % Other manipulations
priority_stim_manip(2).name     = {'rdk'};
priority_stim_manip(2).priority = {'coherence','delta_ref'};
priority_stim_manip(2).other    = {'ori_bin'};
priority_stim_manip(3).name     = {'dot'};
priority_stim_manip(3).priority = {'ori_bin'};
priority_stim_manip(3).other    = {};
priority_stim_manip(4).name     = {'cobj'};
priority_stim_manip(4).priority = {'super_cat','basic_cat'};
priority_stim_manip(4).other    = {'sub_cat'};
priority_stim_manip(5).name     = {'ns'};
priority_stim_manip(5).priority = {'super_cat','basic_cat'};
priority_stim_manip(5).other    = {'sub_cat'}; % ERK: do we sample/index these or does it not matter since stim are selected to be balanced?

%% Loop over stimulus super classes
all_trials = struct();

for stimClass_idx = 1:length(p.exp.stimClassLabels)
    
    stimClass_name = p.exp.stimClassLabels{stimClass_idx};
    
    switch stimClass_name
        
        %% GABOR
        case 'gabor'
            
            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.exp.stimClassLabels,x), {stimClass_name}, 'UniformOutput', false));
            task_crossings = find(p.exp.crossings(bsc_idx,:));
            
            % Get stim manipulations
            n_contrasts = length(p.stim.gabor.contrast);
            n_ori_bins  = length(p.stim.gabor.ori_deg);                     % number of bins per hemifield
            n_phases    = length(p.stim.gabor.ph_deg);
            loc_stim    = [1,2]; %{'L','R'};
            n_stim_loc  = length(loc_stim);
            
            tasks = [];
            
            for curr_task = 1:length(task_crossings)
                
                tasks(task_crossings(curr_task)).name = p.exp.taskClassLabels{task_crossings(curr_task)};
                use_fix_flag = strcmp(tasks(task_crossings(curr_task)).name,'fix');
                
                if use_fix_flag
                    stimloc_cues = 3;%{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                    tasks(task_crossings(curr_task)).has_loc_cue = false;
                else
                    tasks(task_crossings(curr_task)).has_loc_cue =  p.exp.trial.stim_LR_loc_cue(stimClass_idx);
                    stimloc_cues = [1,0]; %{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                end
                
                % gabors array dims: 8 orientations x 4 phases x 3 contrasts x 4 deltas
                n_unique_cases = n_contrasts*n_ori_bins;
                
                % Preallocate space
                trial = ...
                    struct('trial_nr',[],'unique_im_nr',[],'contrast',[],'phase',[],'orient',[],...
                    'ref_delta',[], 'stim_loc', [],'cue_status',[], 'covert_att_loc_cue', []);
                
                if p.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.exp.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.exp.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end
                
                if use_fix_flag
                    % Get vectors for each stimulus manipulation in one repeat of the unique cases
                    ori_vec     = [repmat(p.stim.gabor.ori_deg, 1, n_unique_cases/n_ori_bins), ...
                                    repmat(p.stim.gabor.ori_deg, 1, n_unique_cases/n_ori_bins)];
                    ph_vec      = [repmat(p.stim.gabor.ph_deg, 1, n_unique_cases/n_phases), ...
                                    repmat(p.stim.gabor.ph_deg, 1, n_unique_cases/n_phases)];
                    stimloc_vec = [repmat(loc_stim, 1, n_unique_cases/n_stim_loc), ...
                                    repmat(fliplr(loc_stim), 1, n_unique_cases/n_stim_loc)];
                    con_vec     = [repelem(p.stim.gabor.contrast, n_unique_cases/n_contrasts), ...
                                    repelem(p.stim.gabor.contrast, n_unique_cases/n_contrasts)]; 
                else
                    % Get vectors for each stimulus manipulation in one repeat of the unique cases
                    ori_vec = repmat(p.stim.gabor.ori_deg, 1, n_unique_cases/n_ori_bins);
                    ph_vec  = repmat(p.stim.gabor.ph_deg, 1, n_unique_cases/n_phases);
                    stimloc_vec = repmat(loc_stim, 1, n_unique_cases/n_stim_loc);
                    con_vec = repelem(p.stim.gabor.contrast,n_unique_cases/n_contrasts);
                end
                
                % give each unique image a nr, define it's properties
                unique_im = cat(1, 1:length(ori_vec), stimloc_vec, ori_vec, con_vec, ph_vec)';
                
                % Start shuffling the conditions and add to master matrix
                cond_master = [];
                for rep = 1:n_repeats
                    
                    conds_master_single_rep = []; 
                    
                    if use_fix_flag
                        
                        % shuffle orientation every 8 trials
                        shuffle_ori = shuffle_concat([1:n_ori_bins],n_contrasts*2);
                        shuffle_ori = shuffle_ori' + repmat(repelem([0:n_ori_bins:(n_unique_cases-1)],n_ori_bins),1,2)';
                        
                        % shuffle contrasts every 3 trials
                        shuffle_c = shuffle_concat(1:n_contrasts,(n_unique_cases/n_contrasts)*2);
                        case_vec = reshape(1:size(unique_im,1),[],3);
                        case_vec = case_vec';
                        case_vec = case_vec(:);
                        shuffled_c = case_vec(shuffle_c + repmat(repelem([0:n_contrasts:(n_unique_cases-1)],n_contrasts),1,2),:);
                        
                        conds_shuffle0 = unique_im(shuffle_ori,:); % <-- first shuffle based on unique nr of orientations
                        conds_shuffle1 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of contrasts, such that new trial order prioritized contrast order
 
                        
                        cue_vec = NaN(size(unique_im,1),1);
                        
                        conds_shuffle2 = [conds_shuffle1, cue_vec];
                        
                        conds_master_single_rep = [conds_master_single_rep;conds_shuffle2];
                        
                    else
                        for loc = 1:n_stimloc_cues % cued vs uncued
                            % shuffle orientation every 8 trials
                            shuffle_ori = shuffle_concat([1:n_ori_bins],n_contrasts);
                            shuffle_ori = shuffle_ori' + repelem([0:n_ori_bins:(n_unique_cases-1)],n_ori_bins)';
                            
                            % shuffle contrasts every 3 trials
                            shuffle_c = shuffle_concat(1:n_contrasts,n_unique_cases/n_contrasts);
                            case_vec = reshape(1:n_unique_cases,[],3);
                            case_vec = case_vec';
                            case_vec = case_vec(:);
                            shuffled_c = case_vec(shuffle_c + repelem([0:n_contrasts:(n_unique_cases-1)],n_contrasts),:);
                            
                            conds_shuffle0 = unique_im(shuffle_ori,:); % <-- first shuffle based on unique nr of orientations
                            conds_shuffle1 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of contrasts, such that new trial order prioritized contrast order
                            
                            if loc == 1
                                cue_vec = shuffle_concat(stimloc_cues,n_unique_cases/n_stimloc_cues)';
                                
                            elseif loc == 2
                                cue_vec_anti = zeros(size(cue_vec));
                                prev_cue_vec = conds_shuffle1(find(cue_vec),1);
                                
                                conds_to_be_cued = setdiff([1:n_unique_cases],prev_cue_vec);
                                conds_to_be_cued_i = intersect(conds_shuffle1(:,1),conds_to_be_cued','stable');
                                cue_vec_anti(conds_to_be_cued_i) = 1;
                                cue_vec = cue_vec_anti;
                            end
                        end
                        
                        conds_shuffle2 = [conds_shuffle1, cue_vec];
                        
                        conds_master_single_rep = [conds_master_single_rep;conds_shuffle2];
                    end
                    
                    clear conds_shuffle0 conds_shuffle1 conds_shuffle2
                    
                    % reorganize column order
                    conds_master_single_rep2 = NaN(size(conds_master_single_rep));
                    
                    conds_master_single_rep2(:,1) = conds_master_single_rep(:,1); % unique im nr
                    conds_master_single_rep2(:,2) = conds_master_single_rep(:,2); % stim loc
                    conds_master_single_rep2(:,3) = conds_master_single_rep(:,6); % cue stat
                    conds_master_single_rep2(:,4) = conds_master_single_rep(:,3); % orient
                    conds_master_single_rep2(:,5) = conds_master_single_rep(:,4); % contrast
                    conds_master_single_rep2(:,6) = conds_master_single_rep(:,5); % phase
                    
                    % Merge trials and add fix cue thickening direction
                    conds_master_single_rep3 = create_condition_master_trials(conds_master_single_rep2,use_fix_flag);
                    
                    % Add WM change
                    if strcmp(tasks(task_crossings(curr_task)).name,'wm')
                        n_deltas = length(p.stim.gabor.delta_from_ref);
                        delta_vec = shuffle_concat(p.stim.gabor.delta_from_ref, (n_unique_cases/n_deltas));
                        orient2 = conds_master_single_rep3(:,6) + delta_vec';
                    else
                        orient2 = NaN(size(conds_master_single_rep3,1),1);
                    end
                    
                    % Keep track of repeat
                    rep_vec = rep.*ones(size(conds_master_single_rep3,1),1);
                    conds_master_single_rep3 = [conds_master_single_rep3, orient2, rep_vec];
                    
                    % Accummulate
                    cond_master = [cond_master; conds_master_single_rep3];
                    
                    % Clean up
                    clear  conds_master_single_rep2 conds_master_single_rep3
                    
                end
                
                % conds_master_reordered :   n/2 trials x 9 matrix with the following columns
                %                               1: trial nr (every trial occupies two rows)
                %                               2: thickening_dir (1=left, 2=right, 3=both)
                %                               3: unique im nr
                %                               4: stim loc (1=left, 2=right, 3=central)
                %                               5: cue status (0=uncued, 1 =cued)
                %                        OTHER columns contain stim info, e.g. for Gabors:
                %                               6: angle
                %                               7: contrast
                %                               8: phase
                %                               9: delta ref
                %                               10: repeat nr
                
                %% Divide trials into miniblocks
                % if nr of unique cases doesn't fit in the number of trials
                % per block, then shave of 2 trials.. hopefully that works.
                lingering_trials = mod(size(cond_master,1)/2,n_trials_per_block);
                if lingering_trials ~= 0
                    n_trials_per_block2 = n_trials_per_block - 2;
                else
                    n_trials_per_block2 = n_trials_per_block;
                end
                miniblock_trials = reshape(1:2:size(cond_master,1),n_trials_per_block2,[]);
                
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        for stimloc = [1,2]
                            trial(local_trial_nr,stimloc).trial_nr           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),1);
                            trial(local_trial_nr,stimloc).covert_att_loc_cue = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),2);
                            trial(local_trial_nr,stimloc).unique_im_nr       = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),3);
                            trial(local_trial_nr,stimloc).stim_loc           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),4);
                            trial(local_trial_nr,stimloc).cue_status         = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),5);
                            
                            trial(local_trial_nr,stimloc).orient             = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),6);
                            trial(local_trial_nr,stimloc).contrast           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),7);
                            trial(local_trial_nr,stimloc).phase              = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),8);
                            trial(local_trial_nr,stimloc).ref_delta          = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),9);
                            trial(local_trial_nr,stimloc).repeat_nr          = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),10);

                        end
                        
                    end
                    
                    tasks(task_crossings(curr_task)).miniblock(mm).trial = trial;
                end
            end
            clear cond_master
            all_trials.stim(stimClass_idx).name  = stimClass_name;
            all_trials.stim(stimClass_idx).tasks = tasks;
            
            %% RDK
        case 'rdk'
            
            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.exp.stimClassLabels,x), {'rdk'}, 'UniformOutput', false));
            task_crossings = find(p.exp.crossings(bsc_idx,:));
            
            % Get stim manipulations
            n_coh       = length(p.stim.rdk.dots_coherence);
            n_motdir    = length(p.stim.rdk.dots_direction);                 % number of motion direction bins per hemifield
            loc_stim    = [1,2];%{'L','R'};
            n_stim_loc  = length(loc_stim);
            
            % gabors array dims: 8 orientations x 4 phases x 3 contrasts x 4 deltas
            n_unique_cases = n_coh*n_motdir;
            
            tasks = [];
            
            for curr_task = 1:length(task_crossings)
                
                tasks(task_crossings(curr_task)).name = p.exp.taskClassLabels{task_crossings(curr_task)};
                
                use_fix_flag = strcmp(tasks(task_crossings(curr_task)).name,'fix');
                
                if use_fix_flag
                    stimloc_cues = 3;%{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                    tasks(task_crossings(curr_task)).has_loc_cue = false;
                else
                    tasks(task_crossings(curr_task)).has_loc_cue =  p.exp.trial.stim_LR_loc_cue(stimClass_idx);
                    stimloc_cues = [1,0]; %{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                end
                
                
                % Preallocate space
                trial = ...
                    struct('trial_nr',[],'unique_im_nr',[],'cue_status',[],'coh',[],'motdir',[],...
                    'ref_delta',[],'covert_att_loc_cue', []);
                
                if p.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.exp.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.exp.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end
                
                if use_fix_flag
                    % Get vectors for each stimulus manipulation in one repeat of the unique cases
                    motdir_vec = [repmat(p.stim.rdk.dots_direction, 1, n_unique_cases/n_motdir), ...
                                    repmat(p.stim.rdk.dots_direction, 1, n_unique_cases/n_motdir)];
                    stimloc_vec = [repmat(loc_stim, 1, n_unique_cases/n_stim_loc), ...
                                    repmat(fliplr(loc_stim), 1, n_unique_cases/n_stim_loc)];
                    coh_vec     = [repelem(p.stim.rdk.dots_coherence,n_unique_cases/n_coh), ...
                                    repelem(p.stim.rdk.dots_coherence,n_unique_cases/n_coh)];
                else
                    % Get vectors for each stimulus manipulation in one repeat of the unique cases
                    motdir_vec = repmat(p.stim.rdk.dots_direction, 1, n_unique_cases/n_motdir);
                    stimloc_vec = repmat(loc_stim, 1, n_unique_cases/n_stim_loc);
                    coh_vec = repelem(p.stim.rdk.dots_coherence,n_unique_cases/n_coh);
                end
                
                % give each unique image a nr, define it's properties
                unique_im = cat(1, 1:length(motdir_vec), stimloc_vec, motdir_vec, coh_vec)';
                
                % Start shuffling the conditions and add to master matrix
                cond_master = [];
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
                        
                        conds_shuffle0 = unique_im(shuffle_motdir,:); % <-- first shuffle based on unique nr of orientations
                        conds_shuffle1 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of coherence levels, such that new trial order prioritizes coherence order
 
                        cue_vec = NaN(size(unique_im,1),1);
                        
                        conds_shuffle2 = [conds_shuffle1, cue_vec];
                        
                        conds_master_single_rep = [conds_master_single_rep;conds_shuffle2];
                        
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
                            conds_shuffle1 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of coh, such that new trial order prioritized coh order

                            if loc == 1
                                cue_vec = shuffle_concat(stimloc_cues,n_unique_cases/n_stimloc_cues)';
                                
                            elseif loc == 2
                                cue_vec_anti = zeros(size(cue_vec));
                                prev_cue_vec = conds_shuffle1(find(cue_vec),1);
                                
                                conds_to_be_cued = setdiff([1:n_unique_cases],prev_cue_vec);
                                conds_to_be_cued_i = intersect(conds_shuffle1(:,1),conds_to_be_cued','stable');
                                cue_vec_anti(conds_to_be_cued_i) = 1;
                                cue_vec = cue_vec_anti;
                            end
                            
                            conds_shuffle2 = [conds_shuffle1, cue_vec];

                            conds_master_single_rep = [conds_master_single_rep;conds_shuffle2];
                        end
                    end                    
                    clear conds_shuffle0 conds_shuffle1 conds_shuffle2
                    
                                        
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
                    
                    % Keep track of repeat
                    rep_vec = rep.*ones(size(conds_master_single_rep3,1),1);
                    conds_master_single_rep3 = [conds_master_single_rep3, motdir2, rep_vec];
                    
                    % Accummulate
                    cond_master = [cond_master; conds_master_single_rep3];
                    
                    % Clean up
                    clear  conds_master_single_rep2 conds_master_single_rep3
                end
                
                %% Divide trials into miniblocks
                % if nr of unique cases doesn't fit in the number of trials
                % per block, then shave of 2 trials.. hopefully that works.
                lingering_trials = mod(size(cond_master,1)/2,n_trials_per_block);
                if lingering_trials ~= 0
                    n_trials_per_block2 = n_trials_per_block - 2;
                else
                    n_trials_per_block2 = n_trials_per_block;
                end
                miniblock_trials = reshape(1:2:size(cond_master,1),n_trials_per_block2,[]);
                
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        for stimloc = [1,2]
                            trial(local_trial_nr,stimloc).trial_nr           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),1);
                            trial(local_trial_nr,stimloc).covert_att_loc_cue = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),2);
                            trial(local_trial_nr,stimloc).unique_im_nr       = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),3);
                            trial(local_trial_nr,stimloc).stim_loc           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),4);
                            trial(local_trial_nr,stimloc).cue_status         = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),5);
                            
                            trial(local_trial_nr,stimloc).motdir         	 = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),6);
                            trial(local_trial_nr,stimloc).coh                = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),7);
                            trial(local_trial_nr,stimloc).ref_delta          = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),8);
                            trial(local_trial_nr,stimloc).repeat_nr          = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),9);

                        end
                        
                    end
                    
                    tasks(task_crossings(curr_task)).miniblock(mm).trial = trial;
                end
            end
            clear cond_master
            all_trials.stim(stimClass_idx).name  = stimClass_name;
            all_trials.stim(stimClass_idx).tasks = tasks;
            
            %% SIMPLE DOTS
        case 'dot'
            
            % Get stim manipulations
            n_dot_loc      = size(p.stim.dot.loc_deg,2);
            loc_stim = [1,2];%{'L','R'};
            n_stim_loc = length(loc_stim);
            
            n_unique_cases = n_dot_loc;
            
            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.exp.stimClassLabels,x), {'dot'}, 'UniformOutput', false));
            task_crossings = find(p.exp.crossings(bsc_idx,:));
  
            tasks = [];
            
            for curr_task = 1:length(task_crossings)
                
                tasks(task_crossings(curr_task)).name = p.exp.taskClassLabels{task_crossings(curr_task)};

                use_fix_flag = strcmp(tasks(task_crossings(curr_task)).name,'fix');
                
                if use_fix_flag
                    stimloc_cues = 3;%{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                    tasks(task_crossings(curr_task)).has_loc_cue = false;
                else
                    tasks(task_crossings(curr_task)).has_loc_cue =  p.exp.trial.stim_LR_loc_cue(stimClass_idx);
                    stimloc_cues = [1,0]; %{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                end
                
                % Preallocate space
                trial = ...
                    struct('loc_deg',[],'cue_status',[],'unique_im_nr',[],...
                    'trial_nr',[],'stim_loc',[],'ref_delta',[],'covert_att_loc_cue', []);
                
                if p.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.exp.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.exp.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end
            
                if use_fix_flag % double stim locations for fix task
                    % Get vectors for each stimulus manipulation in one repeat of the unique cases
                    dot_loc_vec = [repmat(p.stim.dot.loc_deg, 1, n_unique_cases/n_dot_loc), ...
                        repmat(p.stim.dot.loc_deg, 1, n_unique_cases/n_dot_loc)];
                    stimloc_vec = [repmat(loc_stim, 1, n_unique_cases/n_stim_loc), ...
                        repmat(fliplr(loc_stim), 1, n_unique_cases/n_stim_loc)];
                else
                    % Get vectors for each stimulus manipulation in one repeat of the unique cases
                    dot_loc_vec = repmat(p.stim.dot.loc_deg, 1, n_unique_cases/n_dot_loc);
                    stimloc_vec = repmat(loc_stim, 1, n_unique_cases/n_stim_loc);
                end
                
                % give each unique image a nr, define it's properties
                unique_im = cat(1, 1:length(dot_loc_vec), stimloc_vec, dot_loc_vec)';

                % Start shuffling the conditions and add to master matrix
                cond_master = [];
                
                for rep = 1:n_repeats
                    
                    conds_master_single_rep = [];
                    
                    if use_fix_flag 
                     
                        % shuffle dot angle every 8 trials
                        shuffle_loc = shuffle_concat([1:n_dot_loc],2);
                        
                        conds_shuffle1 = unique_im(shuffle_loc,:); % <-- shuffle based on unique nr of orientations
                        
                        cue_vec = NaN(size(unique_im,1),1);
                        
                        conds_shuffle2 = [conds_shuffle1, cue_vec];
                        
                        conds_master_single_rep = [conds_master_single_rep;conds_shuffle2];
                        
                    else
                        
                        for loc = 1:n_stimloc_cues % cued vs uncued
                            
                            % shuffle dot angle every 8 trials
                            shuffle_loc = shuffle_concat([1:n_dot_loc],1);
                            
                            conds_shuffle1 = unique_im(shuffle_loc,:); % <-- shuffle based on unique nr of orientations
                            
                            if loc == 1
                                cue_vec = shuffle_concat(stimloc_cues,n_unique_cases/n_stimloc_cues)';
                                
                            elseif loc == 2
                                cue_vec_anti = zeros(size(cue_vec));
                                prev_cue_vec = conds_shuffle1(find(cue_vec),1);
                                
                                conds_to_be_cued = setdiff([1:n_unique_cases],prev_cue_vec);
                                conds_to_be_cued_i = intersect(conds_shuffle1(:,1),conds_to_be_cued','stable');
                                cue_vec_anti(conds_to_be_cued_i) = 1;
                                cue_vec = cue_vec_anti;
                            end
                            
                            conds_shuffle2 = [conds_shuffle1, cue_vec];
                            
                            conds_master_single_rep = [conds_master_single_rep;conds_shuffle2];
                        end
                    end
                    clear conds_shuffle0 conds_shuffle1 conds_shuffle2
                    
                    % Reorganize columns
                    conds_master_single_rep2 = NaN(size(conds_master_single_rep));
                    
                    conds_master_single_rep2(:,1) = conds_master_single_rep(:,1); % unique im nr
                    conds_master_single_rep2(:,2) = conds_master_single_rep(:,2); % stim loc
                    conds_master_single_rep2(:,3) = conds_master_single_rep(:,4); % cue stat
                    conds_master_single_rep2(:,4) = conds_master_single_rep(:,3); % dot loc
                    
                    % Merge unique im into trials, and add fix cue thickening direction.                     
                    conds_master_single_rep3 = create_condition_master_trials(conds_master_single_rep2,use_fix_flag);
                    
                    % Add WM change.
                    if strcmp(tasks(task_crossings(curr_task)).name,'wm')
                        n_deltas = length(p.stim.dot.delta_from_ref);
                        delta_vec = shuffle_concat(p.stim.dot.delta_from_ref, (n_unique_cases/(n_deltas/2)));
                        angle2 = conds_master_single_rep3(:,6) + delta_vec';
                    else
                        angle2 = NaN(size(conds_master_single_rep3,1),1);
                    end
                    
                    % Keep track of repeat
                    rep_vec = rep.*ones(size(conds_master_single_rep3,1),1);
                    conds_master_single_rep3 = [conds_master_single_rep3, angle2, rep_vec];
                    
                    % Accummulate
                    cond_master = [cond_master; conds_master_single_rep3];
                    
                    % Clean up
                    clear  conds_master_single_rep2 conds_master_single_rep3
                end
                
                %% Divide trials into miniblocks
                % if nr of unique cases doesn't fit in the number of trials
                % per block, then shave of 2 trials.. hopefully that works.
                lingering_trials = mod(size(cond_master,1)/2,n_trials_per_block);
                if lingering_trials ~= 0
                    n_trials_per_block2 = n_trials_per_block - 2;
                else
                    n_trials_per_block2 = n_trials_per_block;
                end
                miniblock_trials = reshape(1:2:size(cond_master,1),n_trials_per_block2,[]);
                
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        for stimloc = [1,2]
                            trial(local_trial_nr,stimloc).trial_nr           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),1);
                            trial(local_trial_nr,stimloc).covert_att_loc_cue = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),2);
                            trial(local_trial_nr,stimloc).unique_im_nr       = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),3);
                            trial(local_trial_nr,stimloc).stim_loc           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),4);
                            trial(local_trial_nr,stimloc).cue_status         = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),5);
                            
                            trial(local_trial_nr,stimloc).loc_deg            = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),6);
                            trial(local_trial_nr,stimloc).ref_delta          = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),7);
                            trial(local_trial_nr,stimloc).repeat_nr          = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),8);
                        end
                        
                    end
                    
                    tasks(task_crossings(curr_task)).miniblock(mm).trial = trial;
                end
            end
            clear cond_master
            
            all_trials.stim(stimClass_idx).name  = stimClass_name;
            all_trials.stim(stimClass_idx).tasks = tasks;
            
            
            %% COMPLEX OBJECTS
        case 'cobj'
            
            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.exp.stimClassLabels,x), {'cobj'}, 'UniformOutput', false));
            task_crossings = find(p.exp.crossings(bsc_idx,:));

            % Get stim manipulations
           
            n_super_cat = length(p.stim.cobj.super_cat);
            loc_stim = [1,2];%{'L','R'};
            n_stim_loc = length(loc_stim);

            tasks = [];
            
            for curr_task = 1:length(task_crossings)
                
                tasks(task_crossings(curr_task)).name = p.exp.taskClassLabels{task_crossings(curr_task)};
                use_fix_flag = strcmp(tasks(task_crossings(curr_task)).name,'fix');

                if use_fix_flag
                    stimloc_cues = 3;%{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                    tasks(task_crossings(curr_task)).has_loc_cue = false;
                else
                    tasks(task_crossings(curr_task)).has_loc_cue =  p.exp.trial.stim_LR_loc_cue(stimClass_idx);
                    stimloc_cues = [1,0]; %{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                end
                
                % Preallocate space
                trial = ...
                    struct('loc_deg',[],'cue_status',[],'unique_im_nr',[],...
                    'trial_nr',[],'stim_loc',[], 'super_cat',[],'super_cat_name',[],...
                    'basic_cat',[],'basic_cat_name',[],...
                    'sub_cat',[],'sub_cat_name',[],...
                    'ref_delta',[],'covert_att_loc_cue', []);
                
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
                
                n_unique_cases = length(basic_cat_vec);
                
                % Check type of trial
                if p.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.exp.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.exp.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end
                
                if use_fix_flag % double stim locations for fix task
                    super_cat_vec = []; sub_cat_vec = [];
                    for ii = 1:length(n_basic_cat)
                        n_sub_cat(ii) = length(p.stim.cobj.sub_cat{ii});
                        super_cat_vec = cat(2, super_cat_vec, repelem(ii,n_sub_cat(ii)));
                        sub_cat_vec = cat(2, sub_cat_vec, 1:n_sub_cat(ii));
                    end
                    super_cat_vec = repmat(super_cat_vec,1,2);
                    sub_cat_vec   = repmat(sub_cat_vec,1,2);
                    basic_cat_vec = repmat(basic_cat_vec,1,2);
                    stimloc_vec   = [repmat(loc_stim, 1, n_unique_cases/n_stim_loc), ...
                        repmat(fliplr(loc_stim), 1, n_unique_cases/n_stim_loc)];
                    facing_dir_vec = repmat(p.stim.cobj.facing_dir_deg,1,2);
                    
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
                    
                    
                else
                    
                    % Get vectors for each stimulus manipulation in one repeat of the unique cases
                    super_cat_vec = []; sub_cat_vec = [];
                    for ii = 1:length(n_basic_cat)
                        n_sub_cat(ii) = length(p.stim.cobj.sub_cat{ii});
                        super_cat_vec = cat(2, super_cat_vec, repelem(ii,n_sub_cat(ii)));
                        sub_cat_vec = cat(2, sub_cat_vec, 1:n_sub_cat(ii));
                    end
                    stimloc_vec = repmat(loc_stim, 1, n_unique_cases/n_stim_loc);
                    facing_dir_vec = p.stim.cobj.facing_dir_deg;
                    
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

                end
                
                % give each unique image a nr, define it's properties
                unique_im = cat(1, 1:length(stimloc_vec), stimloc_vec, super_cat_vec, basic_cat_vec, sub_cat_vec, facing_dir_vec)';
                
                
                % Start shuffling the conditions and add to master matrix
                cond_master = [];
                
                for rep = 1:n_repeats
                    
                    conds_master_single_rep = [];
                    if use_fix_flag
                        
                        cue_vec = NaN(size(unique_im,1),1);
                        
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
                        assert(isequal(histcounts(shuffled_super_cat),n_sub_cat*2))
                        
                        % NOW SHUFFLE!
                        conds_shuffle1 = unique_im(basic_cat_shuffle_idx,:); % <-- shuffle based on basic category
                        conds_shuffle2 = conds_shuffle1(super_cat_shuffle_idx,:); % <-- shuffle based on super category, last one to give it priority
                        
                        conds_shuffle3 = [conds_shuffle2, cue_vec];
                        
                        conds_master_single_rep = [conds_master_single_rep;conds_shuffle3];
                    else
                        
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
                            conds_shuffle1 = unique_im(basic_cat_shuffle_idx,:); % <-- shuffle based on basic category
                            conds_shuffle2 = conds_shuffle1(super_cat_shuffle_idx,:); % <-- shuffle based on super category, last one to give it priority
                            
                            % shuffle cued vs uncued status
                            if loc == 1
                                cue_vec = shuffle_concat(stimloc_cues,n_unique_cases/n_stimloc_cues)';
                                
                            elseif loc == 2
                                cue_vec_anti = zeros(size(cue_vec));
                                prev_cue_vec = conds_shuffle2(find(cue_vec),1);
                                
                                conds_to_be_cued = setdiff([1:n_unique_cases],prev_cue_vec);
                                conds_to_be_cued_i = intersect(conds_shuffle2(:,1),conds_to_be_cued','stable');
                                cue_vec_anti(conds_to_be_cued_i) = 1;
                                cue_vec = cue_vec_anti;
                            end
                            
                            conds_shuffle3 = [conds_shuffle2, cue_vec];
                            
                            conds_master_single_rep = [conds_master_single_rep;conds_shuffle3];
                        end
                    end
                    % save some memory
                    clear conds_shuffle0 conds_shuffle1 conds_shuffle2 conds_shuffle3
                    
                    
                    % Reorganize columns
                    conds_master_single_rep2 = NaN(size(conds_master_single_rep));
                    
                    conds_master_single_rep2(:,1) = conds_master_single_rep(:,1); % unique im nr
                    conds_master_single_rep2(:,2) = conds_master_single_rep(:,2); % stim loc
                    conds_master_single_rep2(:,3) = conds_master_single_rep(:,7); % cue stat
                    conds_master_single_rep2(:,4) = conds_master_single_rep(:,3); % super category
                    conds_master_single_rep2(:,5) = conds_master_single_rep(:,4); % basic category
                    conds_master_single_rep2(:,6) = conds_master_single_rep(:,5); % sub category
                    conds_master_single_rep2(:,7) = conds_master_single_rep(:,6); % facing direction
                    
                    % Merge unique im into trials, and add fix cue thickening direction.                     
                    use_fix_flag = strcmp(tasks(task_crossings(curr_task)).name,'fix');
                    conds_master_single_rep3 = create_condition_master_trials(conds_master_single_rep2,use_fix_flag);
                    
                    % Add WM change.
                    if strcmp(tasks(task_crossings(curr_task)).name,'wm')
                        n_deltas = length(p.stim.cobj.delta_from_ref);
                        delta_vec = shuffle_concat(p.stim.cobj.delta_from_ref, (n_unique_cases/(n_deltas/2)));
                        angle2 = conds_master_single_rep3(:,7) + delta_vec';
                    else
                        angle2 = NaN(size(conds_master_single_rep3,1),1);
                    end
                    
                    % Keep track of repeat
                    rep_vec = rep.*ones(size(conds_master_single_rep3,1),1);
                    conds_master_single_rep3 = [conds_master_single_rep3, angle2, rep_vec];
                    
                    % Accummulate
                    cond_master = [cond_master; conds_master_single_rep3];
                    
                    % Clean up
                    clear  conds_master_single_rep2 conds_master_single_rep3
                    
                end
                
                %% Divide trials into miniblocks
                % if nr of unique cases doesn't fit in the number of trials
                % per block, then shave of 2 trials.. hopefully that works.
                lingering_trials = mod(size(cond_master,1)/2,n_trials_per_block);
                if lingering_trials ~= 0
                    n_trials_per_block2 = n_trials_per_block - 2;
                else
                    n_trials_per_block2 = n_trials_per_block;
                end
                miniblock_trials = reshape(1:2:size(cond_master,1),n_trials_per_block2,[]);
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        for stimloc = [1,2]
                            trial(local_trial_nr,stimloc).trial_nr           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),1);
                            trial(local_trial_nr,stimloc).covert_att_loc_cue = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),2);
                            trial(local_trial_nr,stimloc).unique_im_nr       = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),3);
                            trial(local_trial_nr,stimloc).stim_loc           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),4);
                            trial(local_trial_nr,stimloc).cue_status         = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),5);
                            
                            trial(local_trial_nr,stimloc).super_cat         = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),6);
                            trial(local_trial_nr,stimloc).basic_cat         = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),7);
                            trial(local_trial_nr,stimloc).sub_cat           = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),8);
                            trial(local_trial_nr,stimloc).super_cat_name    = p.stim.cobj.super_cat(trial(local_trial_nr,stimloc).super_cat);
                            trial(local_trial_nr,stimloc).basic_cat_name    = p.stim.cobj.basic_cat{trial(local_trial_nr,stimloc).super_cat}(trial(local_trial_nr,stimloc).basic_cat);
                            trial(local_trial_nr,stimloc).sub_cat_name      = p.stim.cobj.sub_cat{trial(local_trial_nr,stimloc).super_cat}(trial(local_trial_nr,stimloc).sub_cat);
                            
                            trial(local_trial_nr,stimloc).facing_dir        = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),9);
                            trial(local_trial_nr,stimloc).ref_delta         = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),10);
                            trial(local_trial_nr,stimloc).repeat_nr         = cond_master(miniblock_trials(local_trial_nr,mm)+(stimloc-1),11);

                        end
                        
                    end
                    
                    tasks(task_crossings(curr_task)).miniblock(mm).trial = trial;
                end
            end
            clear cond_master
            
            all_trials.stim(stimClass_idx).name  = stimClass_name;
            all_trials.stim(stimClass_idx).tasks = tasks;
            
            %% NATURALISTIC SCENES
            
        case 'ns'
            
            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.exp.stimClassLabels,x), {'ns'}, 'UniformOutput', false));
            task_crossings = find(p.exp.crossings(bsc_idx,:));
            
            % Get stim manipulations
            n_basic_cat = []; n_sub_cat = [];
            n_super_cat = length(p.stim.ns.super_cat);
            
            for ni = 1:n_super_cat
                n_basic_cat(ni) = length(p.stim.ns.basic_cat{ni});
                for nj = 1:n_basic_cat(ni)
                    n_sub_cat(ni,nj) = length(p.stim.ns.sub_cat{ni,nj});
                end
            end
            
            n_unique_cases = sum(n_sub_cat(:));
            
            loc_stim = 1; %{'central'};
            
            % Get vectors for each stimulus manipulation in one repeat of the unique cases
            super_cat_vec = []; basic_cat_vec = []; sub_cat_vec = [];
            basic_cat_shuffle_idx = [];
            for ii = 1:length(n_basic_cat)
                super_cat_vec = cat(2, super_cat_vec, repelem(ii,sum(n_sub_cat(ii,:))));
                basic_cat_vec = cat(2, basic_cat_vec, repmat(1:n_basic_cat(ii),1,n_sub_cat(ii,1)));
                sub_cat_vec = cat(2, sub_cat_vec, [1:n_sub_cat(ii,1), 1:n_sub_cat(ii,2)]);
                curr_im = repelem(length(basic_cat_shuffle_idx): 2 : length(basic_cat_shuffle_idx)+2*n_basic_cat(ii),n_basic_cat(ii));
                basic_cat_shuffle_idx = cat(2, basic_cat_shuffle_idx, curr_im + shuffle_concat([1:n_basic_cat(ii)],n_sub_cat(ii,1)));
            end
            stimloc_vec = repmat(loc_stim, 1, n_unique_cases);
            
            % give each unique image a nr, define it's properties
            unique_im = cat(1, 1:n_unique_cases, stimloc_vec, super_cat_vec, basic_cat_vec, sub_cat_vec)';
            
            tasks = struct('name',[],'miniblock',[],'has_loc_cue',[]);
            
            for curr_task = 1:length(task_crossings)
                
                tasks(task_crossings(curr_task)).name = p.exp.taskClassLabels{task_crossings(curr_task)};
                use_fix_flag = strcmp(tasks(task_crossings(curr_task)).name,'fix');
                tasks(task_crossings(curr_task)).has_loc_cue =  use_fix_flag;
                
                stimloc_cues = 3;%{'cued','uncued'};
                n_stimloc_cues = length(stimloc_cues);
                tasks(task_crossings(curr_task)).has_loc_cue = false;
                
                % Preallocate space
                trial = ...
                    struct(...
                    'stim_loc',[],'cue_status',[],'unique_im_nr',[],...
                    'trial_nr',[],'super_cat',[],'super_cat_name',[],...
                    'basic_cat',[],'basic_cat_name',[],...
                    'sub_cat',[],'sub_cat_name',[],...
                    'change_im',[],'covert_att_loc_cue',[]);
                
                if p.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.exp.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.exp.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end
                
                % Start shuffling the conditions and add to master matrix
                cond_master = [];
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
                    conds_shuffle1 = unique_im(basic_cat_shuffle_idx,:); % <-- shuffle based on basic category
                    conds_shuffle2 = conds_shuffle1(super_cat_shuffle_idx,:); % <-- shuffle based on super category, last one to give it priority
                    
                    % shuffle cued vs uncued status
                    cue_vec = NaN(n_unique_cases,1);
                    cond_shuffle3 = [conds_shuffle2, cue_vec];
                    conds_master_single_rep = [conds_master_single_rep; cond_shuffle3];
                    
                    clear conds_shuffle0 conds_shuffle1 conds_shuffle2 conds_shuffle3
                    
                    % No need to merge unique im into trials because there is only one
                    % stim.
                    trial_vec = [1:size(conds_master_single_rep,1)]';
                    thickening_dir = repmat(stimloc_cues,1,size(conds_master_single_rep,1)/n_stimloc_cues)';
                    
                    % Reorganize columns
                    conds_master_single_rep2 = NaN(size(conds_master_single_rep,1),8);
                    conds_master_single_rep2(:,1) = trial_vec; % trial nr
                    conds_master_single_rep2(:,2) = conds_master_single_rep(:,6); % cue stat
                    conds_master_single_rep2(:,3) = conds_master_single_rep(:,1); % unique im nr
                    conds_master_single_rep2(:,4) = conds_master_single_rep(:,2); % stim loc
                    conds_master_single_rep2(:,5) = thickening_dir; % cue stat
                    
                    conds_master_single_rep2(:,6) = conds_master_single_rep(:,3); % super category
                    conds_master_single_rep2(:,7) = conds_master_single_rep(:,4); % basic category
                    conds_master_single_rep2(:,8) = conds_master_single_rep(:,5); % sub category
                    
                    conds_master_single_rep3 = conds_master_single_rep2;
                    
                    % add WM change
                    if strcmp(tasks(task_crossings(curr_task)).name,'wm')
                        n_changes = length(p.stim.ns.change_im);
                        change_blindness_vec = shuffle_concat(1:length(p.stim.ns.change_im), (size(trial_vec,1)/n_changes))';
                    else
                        change_blindness_vec = NaN(size(conds_master_single_rep3,1),1);
                    end
                                        
                    % Keep track of repeat
                    rep_vec = rep.*ones(size(conds_master_single_rep3,1),1);
                    conds_master_single_rep3 = [conds_master_single_rep3, change_blindness_vec, rep_vec]; % change blindness: 1 = easy, 2 = hard
                    
                    % Accummulate
                    cond_master = [cond_master; conds_master_single_rep3];
                    
                    % Clean up
                    clear  conds_master_single_rep2 conds_master_single_rep3
                end
                
                
                %% Divide trials into miniblocks
                %                 lingering_trials = mod(size(conds_master_reordered,1),n_trials_per_block);
                %                 if (n_trials_per_block == 8) && (lingering_trials == 6)
                %                     n_trials_per_block2 = n_trials_per_block - 2;
                %                 elseif (n_trials_per_block == 4) && (lingering_trials == 2)
                %                     n_trials_per_block2 = n_trials_per_block + 1;
                %                 else
                n_trials_per_block2 = n_trials_per_block;
                %                 end
                miniblock_trials = reshape(1:size(cond_master,1),n_trials_per_block,[]);
                
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        trial(1,local_trial_nr).trial_nr           = cond_master(miniblock_trials(local_trial_nr,mm),1);
                        trial(1,local_trial_nr).covert_att_loc_cue = cond_master(miniblock_trials(local_trial_nr,mm),2);
                        trial(1,local_trial_nr).unique_im_nr       = cond_master(miniblock_trials(local_trial_nr,mm),3);
                        trial(1,local_trial_nr).stim_loc           = cond_master(miniblock_trials(local_trial_nr,mm),4);
                        trial(1,local_trial_nr).cue_status         = cond_master(miniblock_trials(local_trial_nr,mm),5);
                        
                        trial(1,local_trial_nr).super_cat               = cond_master(miniblock_trials(local_trial_nr,mm),6);
                        trial(1,local_trial_nr).basic_cat               = cond_master(miniblock_trials(local_trial_nr,mm),7);
                        trial(1,local_trial_nr).sub_cat                 = cond_master(miniblock_trials(local_trial_nr,mm),8);
                        trial(1,local_trial_nr).super_cat_name          = p.stim.ns.super_cat(trial(local_trial_nr).super_cat);
                        trial(1,local_trial_nr).basic_cat_name          = p.stim.ns.basic_cat{trial(local_trial_nr).super_cat}(trial(local_trial_nr).basic_cat);
                        trial(1,local_trial_nr).sub_cat_name            = p.stim.ns.sub_cat{trial(local_trial_nr).super_cat,trial(local_trial_nr).basic_cat}(trial(local_trial_nr).sub_cat);
                        trial(1,local_trial_nr).change_blindness        = cond_master(miniblock_trials(local_trial_nr,mm),9);
                        trial(1,local_trial_nr).repeat_nr               = cond_master(miniblock_trials(local_trial_nr,mm),10);

                    end
                    
                    tasks(task_crossings(curr_task)).miniblock(mm).trial = trial;
                    
                end
            end
            
            clear cond_master
            
            all_trials.stim(stimClass_idx).name = stimClass_name;
            all_trials.stim(stimClass_idx).tasks = tasks;
    end
end

if p.store_params
    save(fullfile(vcd_rootPath,'workspaces','info',sprintf('trials_%s.mat',datestr(now,30))),'all_trials','priority_stim_manip')
end

return








