function all_trials = vcd_makeTrials(p)


% We need to define all possible combinations
% Depending on the session, run, and subject, we need to define the
% combinations we want to show, converted into miniblocks, where trials are
% sampling the manipulations in a pseudo-random way.

p.task = vcd_getTaskParams;

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

for ll = 1:length(p.task.stimClassLabels)
    
    switch p.task.stimClassLabels{ll}
        
        %% GABOR
        case 'gabor'
            
            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.task.stimClassLabels,x), {'gabor'}, 'UniformOutput', false));
            task_crossings = find(p.task.crossings(bsc_idx,:));
            
            % Get stim manipulations
            n_contrasts = length(p.stim.gabor.contrast);
            n_ori_bins  = length(p.stim.gabor.ori_deg);                     % number of bins per hemifield
            n_phases    = length(p.stim.gabor.ph_deg);
            loc_stim = [1,2];%{'L','R'};
            n_stim_loc = length(loc_stim);
            
            tasks = [];
            
            for curr_task = 1:length(task_crossings)
                
                tasks(curr_task).name = p.task.taskClassLabels{task_crossings(curr_task)};
                use_fix_flag = strcmp(tasks(curr_task).name,'FIX');

                if use_fix_flag
                    stimloc_cues = 3;%{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                    tasks(curr_task).has_loc_cue = false;
                else
                    tasks(curr_task).has_loc_cue =  p.task.trial.stim_LR_loc_cue(ll);
                    stimloc_cues = [1,0]; %{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                end
                
                % gabors array dims: 8 orientations x 4 phases x 3 contrasts x 4 deltas
                n_unique_cases = n_contrasts*n_ori_bins;
                
                % Preallocate space
                trial = ...
                    struct('trial_nr',[],'unique_im_nr',[],'contrast',[],'phase',[],'orient',[],...
                    'ref_delta',[], 'stim_loc', [],'cue_status',[], 'covert_att_loc_cue', []);
                
                if p.task.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.task.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.task.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end
                
                % Get vectors for each stimulus manipulation in one repeat of the unique cases 
                ori_vec = repmat(p.stim.gabor.ori_deg, 1, n_unique_cases/n_ori_bins);
                ph_vec  = repmat(p.stim.gabor.ph_deg, 1, n_unique_cases/n_phases);
                stimloc_vec = repmat(loc_stim, 1, n_unique_cases/n_stim_loc);
                con_vec = repelem(p.stim.gabor.contrast,n_unique_cases/n_contrasts);
                
                % give each unique image a nr, define it's properties
                unique_im = cat(1, 1:length(ori_vec), stimloc_vec, ori_vec, con_vec, ph_vec)';
                
                % Start shuffling the conditions and add to master matrix
                conds_master = [];
                
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
                    
                    if use_fix_flag
                       cue_vec = NaN(n_unique_cases,1); 
                    else
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

                    conds_master = [conds_master;conds_shuffle2];
                end
            
                clear conds_shuffle0 conds_shuffle1 conds_shuffle2

                
                conds_master2 = NaN(size(conds_master));
                
                conds_master2(:,1) = conds_master(:,1); % unique im nr
                conds_master2(:,2) = conds_master(:,2); % stim loc
                conds_master2(:,3) = conds_master(:,6); % cue stat
                conds_master2(:,4) = conds_master(:,3); % orient
                conds_master2(:,5) = conds_master(:,4); % contrast
                conds_master2(:,6) = conds_master(:,5); % phase

                conds_master_reordered = create_condition_master_trials(conds_master2,use_fix_flag);
                
                if strcmp(tasks(curr_task).name,'WM')
                    n_deltas = length(p.stim.gabor.delta_from_ref);
                    delta_vec = shuffle_concat(p.stim.gabor.delta_from_ref, 2*(n_unique_cases/n_deltas));
                    orient2 = conds_master_reordered(:,6) + delta_vec';
                else
                    orient2 = NaN(size(conds_master_reordered,1),1);
                end
                conds_master_reordered = [conds_master_reordered, orient2];
                
                                
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
               
                %% Divide trials into miniblocks
                % if nr of unique cases doesn't fit in the number of trials
                % per block, then shave of 2 trials.. hopefully that works.
                lingering_trials = mod(size(conds_master_reordered,1)/2,n_trials_per_block);
                if lingering_trials ~= 0
                    n_trials_per_block2 = n_trials_per_block - 2;
                else
                    n_trials_per_block2 = n_trials_per_block;
                end
                miniblock_trials = reshape(1:2:size(conds_master_reordered,1),n_trials_per_block2,[]);
                
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        for stimloc = [1,2]
                            trial(local_trial_nr,stimloc).trial_nr           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),1);
                            trial(local_trial_nr,stimloc).covert_att_loc_cue = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),2);
                            trial(local_trial_nr,stimloc).unique_im_nr       = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),3);
                            trial(local_trial_nr,stimloc).stim_loc           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),4);
                            trial(local_trial_nr,stimloc).cue_status         = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),5);

                            trial(local_trial_nr,stimloc).orient             = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),6);
                            trial(local_trial_nr,stimloc).contrast           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),7);
                            trial(local_trial_nr,stimloc).phase              = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),8);
                            trial(local_trial_nr,stimloc).ref_delta          = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),9);
                        end
                        
                    end
                    
                    tasks(curr_task).miniblock(mm).trial = trial;
                end
            end         

            all_trials.(p.task.stimClassLabels{ll}).tasks = tasks;
        
        %% RDK
        case 'rdk'
            
            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.task.stimClassLabels,x), {'rdk'}, 'UniformOutput', false));
            task_crossings = find(p.task.crossings(bsc_idx,:));
            
            % Get stim manipulations
            n_coh       = length(p.stim.rdk.dots_coherence);
            n_motdir    = length(p.stim.rdk.dots_direction);                 % number of motion direction bins per hemifield
            loc_stim    = [1,2];%{'L','R'};
            n_stim_loc  = length(loc_stim);
            
            % gabors array dims: 8 orientations x 4 phases x 3 contrasts x 4 deltas
            n_unique_cases = n_coh*n_motdir;

            tasks = [];
            
            for curr_task = 1:length(task_crossings)
                
                tasks(curr_task).name = p.task.taskClassLabels{task_crossings(curr_task)};
                
                use_fix_flag = strcmp(tasks(curr_task).name,'FIX');

                if use_fix_flag
                    stimloc_cues = 3;%{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                    tasks(curr_task).has_loc_cue = false;
                else
                    tasks(curr_task).has_loc_cue =  p.task.trial.stim_LR_loc_cue(ll);
                    stimloc_cues = [1,0]; %{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                end
                
                
                % Preallocate space
                trial = ...
                    struct('trial_nr',[],'unique_im_nr',[],'cue_status',[],'coh',[],'motdir',[],...
                            'ref_delta',[],'covert_att_loc_cue', []);

                if p.task.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.task.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.task.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end
                
                % Get vectors for each stimulus manipulation in one repeat of the unique cases 
                motdir_vec = repmat(p.stim.rdk.dots_direction, 1, n_unique_cases/n_motdir);
                stimloc_vec = repmat(loc_stim, 1, n_unique_cases/n_stim_loc);
                coh_vec = repelem(p.stim.rdk.dots_coherence,n_unique_cases/n_coh);
                
                % give each unique image a nr, define it's properties
                unique_im = cat(1, 1:length(motdir_vec), stimloc_vec, motdir_vec, coh_vec)';
                
                % Start shuffling the conditions and add to master matrix
                conds_master = [];
                
                for loc = 1:n_stimloc_cues % cued vs uncued

                    % shuffle orientation every 8 trials
                    shuffle_motdir = shuffle_concat([1:n_motdir],n_coh);
                    shuffle_motdir = shuffle_motdir' + repelem([0:n_motdir:(n_unique_cases-1)],n_motdir)';

                    % shuffle contrasts every 3 trials
                    shuffle_c = shuffle_concat(1:n_coh,n_unique_cases/n_coh);
                    case_vec = reshape(1:n_unique_cases,[],3);
                    case_vec = case_vec';
                    case_vec = case_vec(:);
                    shuffled_c = case_vec(shuffle_c + repelem([0:n_coh:(n_unique_cases-1)],n_coh),:);

                    conds_shuffle0 = unique_im(shuffle_motdir,:); % <-- first shuffle based on unique nr of orientations
                    conds_shuffle1 = conds_shuffle0(shuffled_c,:); % <-- then shuffle based on unique nr of coh, such that new trial order prioritized coh order
                    if use_fix_flag
                       cue_vec = NaN(n_unique_cases,1); 
                    else
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

                    conds_master = [conds_master;conds_shuffle2];
                end
            
                clear conds_shuffle0 conds_shuffle1 conds_shuffle2

                
                conds_master2 = NaN(size(conds_master));
                
                conds_master2(:,1) = conds_master(:,1); % unique im nr
                conds_master2(:,2) = conds_master(:,2); % stim loc
                conds_master2(:,3) = conds_master(:,5); % cue stat
                conds_master2(:,4) = conds_master(:,3); % mot dir
                conds_master2(:,5) = conds_master(:,4); % coh

                use_fix_flag = strcmp(tasks(curr_task).name,'FIX');
                conds_master_reordered = create_condition_master_trials(conds_master2,use_fix_flag);
                
                if strcmp(tasks(curr_task).name,'WM')
                    n_deltas = length(p.stim.rdk.delta_from_ref);
                    delta_vec = shuffle_concat(p.stim.rdk.delta_from_ref, (n_unique_cases/n_deltas)*2);
                    motdir2 = conds_master_reordered(:,6) + delta_vec';
                else
                    motdir2 = NaN(size(conds_master_reordered,1),1);
                end
                conds_master_reordered = [conds_master_reordered, motdir2];
                
                
                
                 %% Divide trials into miniblocks
                % if nr of unique cases doesn't fit in the number of trials
                % per block, then shave of 2 trials.. hopefully that works.
                lingering_trials = mod(size(conds_master_reordered,1)/2,n_trials_per_block);
                if lingering_trials ~= 0
                    n_trials_per_block2 = n_trials_per_block - 2;
                else
                    n_trials_per_block2 = n_trials_per_block;
                end
                miniblock_trials = reshape(1:2:size(conds_master_reordered,1),n_trials_per_block2,[]);
                
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        for stimloc = [1,2]
                            trial(local_trial_nr,stimloc).trial_nr           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),1);
                            trial(local_trial_nr,stimloc).covert_att_loc_cue = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),2);
                            trial(local_trial_nr,stimloc).unique_im_nr       = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),3);
                            trial(local_trial_nr,stimloc).stim_loc           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),4);
                            trial(local_trial_nr,stimloc).cue_status         = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),5);

                            trial(local_trial_nr,stimloc).motdir_bin         = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),6);
                            trial(local_trial_nr,stimloc).coh                = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),7);
                            trial(local_trial_nr,stimloc).ref_delta          = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),8);
                        end
                        
                    end
                    
                    tasks(curr_task).miniblock(mm).trial = trial;
                end
            end         
            
            all_trials.(p.task.stimClassLabels{ll}).tasks = tasks;
            
            %% SIMPLE DOTS
        case 'dot'
            
            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.task.stimClassLabels,x), {'dot'}, 'UniformOutput', false));
            task_crossings = find(p.task.crossings(bsc_idx,:));
            
            % Get stim manipulations
            n_dot_loc      = size(p.stim.dot.loc_deg,2);
            loc_stim = [1,2];%{'L','R'};
            n_stim_loc = length(loc_stim);
            
            
            n_unique_cases = n_dot_loc;
            
            % Get vectors for each stimulus manipulation in one repeat of the unique cases 
            dot_loc_vec = repmat(p.stim.dot.loc_deg, 1, n_unique_cases/n_dot_loc);
            stimloc_vec = repmat(loc_stim, 1, n_unique_cases/n_stim_loc);

            % give each unique image a nr, define it's properties
            unique_im = cat(1, 1:length(dot_loc_vec), stimloc_vec, dot_loc_vec)';

            tasks = [];
            
            for curr_task = 1:length(task_crossings)
                
                tasks(curr_task).name = p.task.taskClassLabels{task_crossings(curr_task)};
                use_fix_flag = strcmp(tasks(curr_task).name,'FIX');

                if use_fix_flag
                    stimloc_cues = 3;%{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                    tasks(curr_task).has_loc_cue = false;
                else
                    tasks(curr_task).has_loc_cue =  p.task.trial.stim_LR_loc_cue(ll);
                    stimloc_cues = [1,0]; %{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                end
                
                % Preallocate space
                trial = ...
                    struct('loc_deg',[],'cue_status',[],'unique_im_nr',[],...
                    'trial_nr',[],'stim_loc',[],'ref_delta',[],'covert_att_loc_cue', []);
                
                if p.task.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.task.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.task.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end
                
                % Start shuffling the conditions and add to master matrix
                conds_master = [];
                
                for loc = 1:n_stimloc_cues % cued vs uncued

                    % shuffle orientation every 8 trials
                    shuffle_loc = shuffle_concat([1:n_dot_loc],1);

                    conds_shuffle1 = unique_im(shuffle_loc,:); % <-- shuffle based on unique nr of orientations
                    if use_fix_flag
                       cue_vec = NaN(n_unique_cases,1); 
                    else
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

                    conds_master = [conds_master;conds_shuffle2];
                end
            
                clear conds_shuffle0 conds_shuffle1 conds_shuffle2

                
                conds_master2 = NaN(size(conds_master));
                
                conds_master2(:,1) = conds_master(:,1); % unique im nr
                conds_master2(:,2) = conds_master(:,2); % stim loc
                conds_master2(:,3) = conds_master(:,4); % cue stat
                conds_master2(:,4) = conds_master(:,3); % dot loc
 
                use_fix_flag = strcmp(tasks(curr_task).name,'FIX');
                conds_master_reordered = create_condition_master_trials(conds_master2,use_fix_flag);
                
                if strcmp(tasks(curr_task).name,'WM')
                    n_deltas = length(p.stim.dot.delta_from_ref);
                    delta_vec = shuffle_concat(p.stim.dot.delta_from_ref, (n_unique_cases/n_deltas)*2);
                    angle2 = conds_master_reordered(:,6) + delta_vec';
                else
                    angle2 = NaN(size(conds_master_reordered,1),1);
                end
                conds_master_reordered = [conds_master_reordered, angle2];
                
                
                
                 %% Divide trials into miniblocks
                % if nr of unique cases doesn't fit in the number of trials
                % per block, then shave of 2 trials.. hopefully that works.
                lingering_trials = mod(size(conds_master_reordered,1)/2,n_trials_per_block);
                if lingering_trials ~= 0
                    n_trials_per_block2 = n_trials_per_block - 2;
                else
                    n_trials_per_block2 = n_trials_per_block;
                end
                miniblock_trials = reshape(1:2:size(conds_master_reordered,1),n_trials_per_block2,[]);
                    
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        for stimloc = [1,2]
                            trial(local_trial_nr,stimloc).trial_nr           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),1);
                            trial(local_trial_nr,stimloc).covert_att_loc_cue = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),2);
                            trial(local_trial_nr,stimloc).unique_im_nr       = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),3);
                            trial(local_trial_nr,stimloc).stim_loc           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),4);
                            trial(local_trial_nr,stimloc).cue_status         = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),5);

                            trial(local_trial_nr,stimloc).loc_deg            = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),6);
                            trial(local_trial_nr,stimloc).ref_delta          = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),7);
                        end
                        
                    end
                    
                    tasks(curr_task).miniblock(mm).trial = trial;
                end
            end         
            
            all_trials.(p.task.stimClassLabels{ll}).tasks = tasks;
            
            
            %% COMPLEX OBJECTS
        case 'cobj'

            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.task.stimClassLabels,x), {'cobj'}, 'UniformOutput', false));
            task_crossings = find(p.task.crossings(bsc_idx,:));

            % Get stim manipulations
            basic_cat_vec = []; basic_cat_shuffle_idx = [];
            n_super_cat = length(p.stim.cobj.super_cat);
            loc_stim = [1,2];%{'L','R'};
            n_stim_loc = length(loc_stim);
                
            for ni = 1:n_super_cat
                tmp = unique(p.stim.cobj.basic_cat{ni}, 'stable');
                n_basic_cat(ni) = length(tmp);

                basic_tmp = [];

                for nj = 1:length(tmp)
                    basic_tmp = cat(1,basic_tmp,nj.*(arrayfun(@(x) strcmp(x, tmp(nj)),p.stim.cobj.basic_cat{ni})));
                end
                basic_cat_vec = cat(2, basic_cat_vec, sum(basic_tmp,1));
                basic_cat_shuffle_idx = cat(2, basic_cat_shuffle_idx,  numel(basic_cat_shuffle_idx) + shuffle_concat(1:size(basic_tmp,2),1));
            end
            assert(isequal(1:n_unique_cases,unique(basic_cat_shuffle_idx)))

            % Get vectors for each stimulus manipulation in one repeat of the unique cases
            super_cat_vec = []; sub_cat_vec = [];
            for ii = 1:length(n_basic_cat)
                n_sub_cat(ii) = length(p.stim.cobj.sub_cat{ii});
                super_cat_vec = cat(2, super_cat_vec, repelem(ii,n_sub_cat(ii)));
                sub_cat_vec = cat(2, sub_cat_vec, 1:n_sub_cat(ii));
            end
            stimloc_vec = repmat(loc_stim, 1, n_unique_cases/n_stim_loc);

            n_unique_cases = length(super_cat_vec);
            
            % give each unique image a nr, define it's properties
            unique_im = cat(1, 1:n_unique_cases, stimloc_vec, super_cat_vec, basic_cat_vec, sub_cat_vec, p.stim.cobj.facing_dir_deg)';

            tasks = [];
            
            for curr_task = 1:length(task_crossings)
                
                tasks(curr_task).name = p.task.taskClassLabels{task_crossings(curr_task)};
                use_fix_flag = strcmp(tasks(curr_task).name,'FIX');

                if use_fix_flag
                    stimloc_cues = 3;%{'cued','uncued'};
                    n_stimloc_cues = length(stimloc_cues);
                    tasks(curr_task).has_loc_cue = false;
                else
                    tasks(curr_task).has_loc_cue =  p.task.trial.stim_LR_loc_cue(ll);
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
                
                if p.task.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.task.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.task.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end

                 
                % Start shuffling the conditions and add to master matrix
                conds_master = [];
                
                for loc = 1:n_stimloc_cues % cued vs uncued

                    % Make sure every 8 trials, subject sees images from 4 super
                    % classes
                    super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                    shuffled_super_cat = super_cat_vec(super_cat_shuffle_idx);
                    for jj = 1:n_trials_per_block:(length(shuffled_super_cat))
                   
%                         [C,IA,IC] = unique(shuffle_super_cat(jj:(jj+7)),'legacy');

                        while ~isempty(setdiff([1:n_super_cat],unique(shuffled_super_cat(jj:(jj+n_trials_per_block-1)),'legacy')))
                            super_cat_shuffle_idx = randperm(length(super_cat_vec),length(super_cat_vec));
                            shuffled_super_cat = super_cat_vec(super_cat_shuffle_idx);

                            
                            
                            % attempt to swap indices  
                            
%                             missing_cat = setdiff([1:n_super_cat],C);
%                             missing_cat_loc = find(shuffle_super_cat == missing_cat);
%                             [N, edges] = histcounts(IC);
%                             edges = edges(1:4)+0.5;
%                             [~,fi] = max(N);
%                             
%                             swap_int = find(IC==edges(fi));
%                             IC(swap_int(1)) =  missing_cat;
%                             
%                             shuffle_super_cat(missing_cat_loc(1))=edges(fi);
% 
%                             
%                             shuffle_super_cat(jj:(jj+7)) = IC;
                        end
                    end
                    % Check if n_sub_cat still holds
                    assert(isequal(histcounts(shuffled_super_cat),n_sub_cat))
                    
                    % NOW SHUFFLE!
                    conds_shuffle1 = unique_im(basic_cat_shuffle_idx,:); % <-- shuffle based on basic category
                    conds_shuffle2 = conds_shuffle1(super_cat_shuffle_idx,:); % <-- shuffle based on super category, last one to give it priority

                    % shuffle cued vs uncued status
                    if use_fix_flag
                       cue_vec = NaN(n_unique_cases,1); 
                    else
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
                    end
                    conds_shuffle3 = [conds_shuffle2, cue_vec];

                    conds_master = [conds_master;conds_shuffle3];
                end
            
                clear conds_shuffle0 conds_shuffle1 conds_shuffle2 conds_shuffle3

                
                conds_master2 = NaN(size(conds_master));
                
                conds_master2(:,1) = conds_master(:,1); % unique im nr
                conds_master2(:,2) = conds_master(:,2); % stim loc
                conds_master2(:,3) = conds_master(:,7); % cue stat
                conds_master2(:,4) = conds_master(:,3); % super category
                conds_master2(:,5) = conds_master(:,4); % basic category
                conds_master2(:,6) = conds_master(:,5); % sub category
                conds_master2(:,7) = conds_master(:,6); % facing direction

                use_fix_flag = strcmp(tasks(curr_task).name,'FIX');
                conds_master_reordered = create_condition_master_trials(conds_master2,use_fix_flag);
                
                if strcmp(tasks(curr_task).name,'WM')
                    n_deltas = length(p.stim.cobj.delta_from_ref);
                    delta_vec = shuffle_concat(p.stim.cobj.delta_from_ref, (n_unique_cases/n_deltas)*2);
                    angle2 = conds_master_reordered(:,7) + delta_vec';
                else
                    angle2 = NaN(size(conds_master_reordered,1),1);
                end
                conds_master_reordered = [conds_master_reordered, angle2];
                
                
                
                 %% Divide trials into miniblocks
                % if nr of unique cases doesn't fit in the number of trials
                % per block, then shave of 2 trials.. hopefully that works.
                lingering_trials = mod(size(conds_master_reordered,1)/2,n_trials_per_block);
                if lingering_trials ~= 0
                    n_trials_per_block2 = n_trials_per_block - 2;
                else
                    n_trials_per_block2 = n_trials_per_block;
                end
                miniblock_trials = reshape(1:2:size(conds_master_reordered,1),n_trials_per_block2,[]);
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        for stimloc = [1,2]
                            trial(local_trial_nr,stimloc).trial_nr           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),1);
                            trial(local_trial_nr,stimloc).covert_att_loc_cue = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),2);
                            trial(local_trial_nr,stimloc).unique_im_nr       = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),3);
                            trial(local_trial_nr,stimloc).stim_loc           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),4);
                            trial(local_trial_nr,stimloc).cue_status         = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),5);

                            trial(local_trial_nr,stimloc).super_cat         = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),6);
                            trial(local_trial_nr,stimloc).basic_cat         = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),7);
                            trial(local_trial_nr,stimloc).sub_cat           = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),8);
                            trial(local_trial_nr,stimloc).facing_dir        = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),9);
                            trial(local_trial_nr,stimloc).ref_delta         = conds_master_reordered(miniblock_trials(local_trial_nr,mm)+(stimloc-1),10);
                        end
                        
                    end
                    
                    tasks(curr_task).miniblock(mm).trial = trial;
                end
            end
            
            all_trials.(p.task.stimClassLabels{ll}).tasks = tasks;
            
            %% NATURALISTIC SCENES
            
        case 'ns'
            
            % Get task crossings
            bsc_idx = cell2mat(cellfun(@(x) strcmp(p.task.stimClassLabels,x), {'ns'}, 'UniformOutput', false));
            task_crossings = find(p.task.crossings(bsc_idx,:));
            
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

            tasks = [];

            for curr_task = 1:length(task_crossings)
                
                tasks(curr_task).name = p.task.taskClassLabels{task_crossings(curr_task)};
                use_fix_flag = strcmp(tasks(curr_task).name,'FIX');
                tasks(curr_task).has_loc_cue =  use_fix_flag;

                stimloc_cues = 3;%{'cued','uncued'};
                n_stimloc_cues = length(stimloc_cues);
                tasks(curr_task).has_loc_cue = false;
                
                % Preallocate space
                trial = ...
                    struct(...
                    'stim_loc',[],'cue_status',[],'unique_im_nr',[],...
                    'trial_nr',[],'super_cat',[],'super_cat_name',[],...
                    'basic_cat',[],'basic_cat_name',[],...
                    'sub_cat',[],'sub_cat_name',[],...
                    'change_im',[],'covert_att_loc_cue',[]);
                
                if p.task.trial.single_epoch_tasks(task_crossings(curr_task))
                    n_trials_per_block    = p.task.miniblock.n_trials_single_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                else
                    n_trials_per_block    = p.task.miniblock.n_trials_double_epoch;
                    n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
                end

                % Start shuffling the conditions and add to master matrix
                conds_master = [];
                
                for iter = [1:4]
                
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
                    conds_master = [conds_master; cond_shuffle3];

                    clear conds_shuffle0 conds_shuffle1 conds_shuffle2 conds_shuffle3
                end
                
                trial_vec = [1:size(conds_master,1)]';
                thickening_dir = repmat(stimloc_cues,1,size(conds_master,1)/n_stimloc_cues)';
                
                conds_master2 = NaN(size(conds_master,1),8);
                conds_master2(:,1) = trial_vec; % trial nr
                conds_master2(:,2) = conds_master(:,6); % cue stat
                conds_master2(:,3) = conds_master(:,1); % unique im nr
                conds_master2(:,4) = conds_master(:,2); % stim loc
                conds_master2(:,5) = thickening_dir; % cue stat
                
                conds_master2(:,6) = conds_master(:,3); % super category
                conds_master2(:,7) = conds_master(:,4); % basic category
                conds_master2(:,8) = conds_master(:,5); % sub category
                
                conds_master_reordered = conds_master2;

                
                if strcmp(tasks(curr_task).name,'WM')
                    n_changes = length(p.stim.ns.change_im);
                    change_blindness_vec = shuffle_concat(1:length(p.stim.ns.change_im), (size(trial_vec,1)/n_changes))';
                else
                    change_blindness_vec = NaN(size(conds_master_reordered,1),1);
                end
                
                conds_master_reordered = [conds_master_reordered, change_blindness_vec]; % 1 = easy, 2 = hard

                
                
                 %% Divide trials into miniblocks
                lingering_trials = mod(size(conds_master_reordered,1)/2,n_trials_per_block);
                if lingering_trials ~= 0
                    n_trials_per_block2 = n_trials_per_block - 2;
                else
                    n_trials_per_block2 = n_trials_per_block;
                end
                miniblock_trials = reshape(1:size(conds_master_reordered,1),n_trials_per_block2,[]);
                    
                for mm = 1:size(miniblock_trials,2)
                    for local_trial_nr = 1:size(miniblock_trials,1)
                        trial(local_trial_nr).trial_nr           = conds_master_reordered(miniblock_trials(local_trial_nr,mm),1);
                        trial(local_trial_nr).covert_att_loc_cue = conds_master_reordered(miniblock_trials(local_trial_nr,mm),2);
                        trial(local_trial_nr).unique_im_nr       = conds_master_reordered(miniblock_trials(local_trial_nr,mm),3);
                        trial(local_trial_nr).stim_loc           = conds_master_reordered(miniblock_trials(local_trial_nr,mm),4);
                        trial(local_trial_nr).cue_status         = conds_master_reordered(miniblock_trials(local_trial_nr,mm),5);
                        
                        trial(local_trial_nr).super_cat         = conds_master_reordered(miniblock_trials(local_trial_nr,mm),6);
                        trial(local_trial_nr).basic_cat         = conds_master_reordered(miniblock_trials(local_trial_nr,mm),7);
                        trial(local_trial_nr).sub_cat           = conds_master_reordered(miniblock_trials(local_trial_nr,mm),8);
                        trial(local_trial_nr).change_blindness  = conds_master_reordered(miniblock_trials(local_trial_nr,mm),9);
                    end
                    
                    tasks(curr_task).miniblock(mm).trial = trial;
                end
            end
            
            
            all_trials.(p.task.stimClassLabels{ll}).tasks = tasks;
            
    end
end

if p.store_params
    save(fullfile(vcd_rootPath,'workspaces','trials.mat'),'all_trials','priority_stim_manip')
end

return








