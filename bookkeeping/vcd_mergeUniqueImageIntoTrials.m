function conds_master_reordered_merged = vcd_mergeUniqueImageIntoTrials(conds_master)
% VCD function to merge unique images into stimulus epochs where there is a
% left and right parafoveal stimulus, and add cuing direction column to
% condition master table.
% 
%   conds_master_reordered = vcd_mergeUniqueImageIntoTrials(conds_master,fix_task_flag)
% 
% INPUTS:
% conds_master  :       (table) with columns for stim-trial properties and 
%                       rows with n trials. 
% fix_task_flag :       logical flag to deal with FIX condition (thickening
%                       on both sides)
%
% OUTPUTS:
% conds_master_reordered :   n/2 trials x 9 matrix with the following columns
%                               1: trial nr (every trial occupies two rows)
%                               2: thickening_dir (1=left, 2=right, 3=both)
%                               3: unique im nr
%                               4: stim loc (1=left, 2=right)
%                               5: cue status (0 = uncued, 1 = cued)
%                               6: direction of central cue (1 = left, 2 = right).
%                        OTHER columns contain stim info, e.g. for Gabors:
%                               7: angle
%                               8: contrast
%                               9: phase

%% Determine variables from input table 
% Find stim locations (left or right)
left_stim  = find(conds_master.stimloc==1);
right_stim = find(conds_master.stimloc==2);

% Ensure we have equal left and right stim
assert(length(left_stim)==length(right_stim))

% Preallocate space: Pair left / right stimuli for each trial
trial_vec_i = NaN(size(conds_master,1)/2,2);

% Add the condition master to struct
trial_vec = repelem(1:size(trial_vec_i,1),2); % or     trial_vec = repelem(1:(n_unique_cases/2),2); ????

nr_cueing_conds = 4; %left/right x cued/uncued

% Define thickening of central cue
if ~strcmp(unique(conds_master.task_class_name),{'fix'})
    
    thickening_dir = zeros(length(trial_vec),1);
    th_count = 1;
    
    left_cued    = find((conds_master.iscued == 1) & (conds_master.stimloc == 1));
    left_uncued  = find((conds_master.iscued == 0) & (conds_master.stimloc == 1));
    right_cued   = find((conds_master.iscued == 1) & (conds_master.stimloc == 2));
    right_uncued = find((conds_master.iscued == 0) & (conds_master.stimloc == 2));
    
    cond_idx = shuffle_concat([1:nr_cueing_conds],size(trial_vec_i,1)/nr_cueing_conds);  
    
    
    for ii = 1:size(trial_vec_i,1)

        if cond_idx(ii) == 1 % left cued;
            leading_cue_status  = 1; 
            leading_stim_loc    = 1;
        elseif cond_idx(ii) == 2 % left uncued;
            leading_cue_status  = 0;
            leading_stim_loc    = 1;
        elseif cond_idx(ii) == 3 % right cued;
            leading_cue_status  = 1; 
            leading_stim_loc    = 2;
        elseif cond_idx(ii) == 4 % right uncued;
            leading_cue_status  = 0; 
            leading_stim_loc    = 2;
        end
            
        if leading_cue_status == 1 && leading_stim_loc == 1 % left cued
            thickening_dir(th_count)   = 1; % left thickening
            thickening_dir(th_count+1) = 1; % left thickening
            
            trial_vec_i(ii,1) = left_cued(1);
            left_cued(1) = [];
            trial_vec_i(ii,2) = right_uncued(1);
            right_uncued(1) = [];
            
        elseif leading_cue_status == 1 && leading_stim_loc == 2 % right cued
            thickening_dir(th_count) = 2; % right thickening
            thickening_dir(th_count+1) = 2; % right thickening
            
            trial_vec_i(ii,1) = left_uncued(1);
            left_uncued(1) = [];
            trial_vec_i(ii,2) = right_cued(1);
            right_cued(1) = [];
           
        elseif leading_cue_status == 0 && leading_stim_loc == 1 % left uncued
            thickening_dir(th_count) = 2; % right thickening
            thickening_dir(th_count+1) = 2; % right thickening
            
            trial_vec_i(ii,1) = left_uncued(1);
            left_uncued(1) = [];
            trial_vec_i(ii,2) = right_cued(1);
            right_cued(1) = [];
            
        elseif leading_cue_status == 0 && leading_stim_loc == 2 % right uncued
            thickening_dir(th_count) = 1; % left thickening
            thickening_dir(th_count+1) = 1; % left thickening
            
            trial_vec_i(ii,1) = left_cued(1);
            left_cued(1) = [];
            trial_vec_i(ii,2) = right_uncued(1);
            right_uncued(1) = [];
            
        elseif leading_cue_status == 0 || isnan(leading_cue_status)
            thickening_dir(th_count) = 3;  % both sides thickening
            thickening_dir(th_count+1) = 3;  % both sides thickening
        end
        th_count = th_count +2;
    end
    
    trial_vec_i = trial_vec_i';
    trial_vec_i = trial_vec_i(:);
    


elseif strcmp(unique(conds_master.task_class_name),{'fix'})

    thickening_dir = 3.*ones(size(trial_vec))'; % we thicken at both sides

    for ii = 1:length(conds_master.stimloc)
        stim_loc = conds_master.stimloc(ii);
        stim_to_allocate = find(isnan(trial_vec_i(:,stim_loc)));
        trial_vec_i(stim_to_allocate(1), stim_loc) = ii;
    end
    
    trial_vec_i = trial_vec_i';
    trial_vec_i = trial_vec_i(:);

end
  
if isfield(conds_master, 'unique_trial_nr')
    trial_vec = trial_vec' + max(conds_master.trial_vec);
else
    trial_vec = trial_vec';
end
conds_master_reordered = conds_master(trial_vec_i,:);
conds_master_reordered.unique_trial_nr = trial_vec;
conds_master_reordered.thickening_dir = thickening_dir;


% copy condition order table headers and scrub content
sz = [sum(conds_master_reordered.stimloc==3) + sum(sum(conds_master_reordered.stimloc==[1,2],2))/2, size(conds_master_reordered,2)];
conds_master_reordered_merged = vcd_preallocateNaNTable(sz(1), sz(2), conds_master_reordered(1,:), 2);
varTypes = varfun(@class,conds_master_reordered(1,:),'OutputFormat','cell');


% Loop over columns
for vt = 1:sz(2)
    % Reset trial nr (what row are we allocating)
    trial_nr = 1;
    
    % Merge left/right for each trial
    for tt = 1:sz(1)

        if conds_master_reordered.stimloc == 3
            if strcmp(varTypes(vt),'double')
                conds_master_reordered_merged.(conds_master_reordered.Properties.VariableNames{vt})(tt,:) = ...
                    [table2array(conds_master_reordered(trial_nr,vt)), NaN];
            elseif strcmp(varTypes(vt),'cell')
                conds_master_reordered_merged.(conds_master_reordered.Properties.VariableNames{vt})(tt,:) = ...
                    [{table2array(conds_master_reordered(trial_nr,vt))}, {NaN}];
            end
            trial_nr = trial_nr + 1;
            
        else
            if strcmp(varTypes(vt),'double')
                conds_master_reordered_merged.(conds_master_reordered.Properties.VariableNames{vt})(tt,:) = ...
                    [table2array(conds_master_reordered(trial_nr,vt)), table2array(conds_master_reordered(trial_nr+1,vt))];
            elseif strcmp(varTypes(vt),'cell')
                conds_master_reordered_merged.(conds_master_reordered.Properties.VariableNames{vt})(tt,:) = ...
                    [table2cell(conds_master_reordered(trial_nr,vt)), table2cell(conds_master_reordered(trial_nr+1,vt))];
            end
            
            trial_nr = trial_nr + 2;
        end
    end
end


% Randomize catch trial loc
sz = size(conds_master_reordered_merged);
catch_trials = find(conds_master_reordered_merged.iscatch(:,1)==true);
if ~isempty(catch_trials)
    
    while 1
        shuffle_trial_idx = randi(size(conds_master_reordered_merged,1),length(catch_trials));
        if shuffle_trial_idx ~= catch_trials
            break;
        end
    end
    
    for cc = 1:length(shuffle_trial_idx) % 11
        tmp = conds_master_reordered_merged;
            
        curr_catch_trial   = tmp(catch_trials(cc),:);  % 13 (will become 6)
        
        if shuffle_trial_idx(cc) < catch_trials(cc)
            conds_shuffle1_pre1  = tmp(1:(shuffle_trial_idx(cc)-1),:);  % 1-5 
            conds_shuffle1_post1 = tmp(shuffle_trial_idx(cc):(catch_trials(cc)-1),:); % 6-12
            conds_shuffle1_post2 = tmp((catch_trials(cc)+1):end,:); %post catch trial -end

           conds_master_reordered_merged = cat(1,conds_shuffle1_pre1,curr_catch_trial,conds_shuffle1_post1,conds_shuffle1_post2);

        elseif shuffle_trial_idx(cc) > catch_trials(cc)
            conds_shuffle1_pre1   = tmp(1:(catch_trials(cc)-1),:);  % 1-12
            conds_shuffle1_pre2   = tmp((catch_trials(cc)+1):shuffle_trial_idx(cc),:);
            conds_shuffle1_post1  = tmp((shuffle_trial_idx(cc)+1):end,:); %14-end
            conds_master_reordered_merged = cat(1,conds_shuffle1_pre1,conds_shuffle1_pre2,curr_catch_trial,conds_shuffle1_post1);
        end
        
        
    end
    assert(isequal(size(conds_master_reordered_merged),sz))
end


