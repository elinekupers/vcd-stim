function conds_master_reordered_merged = vcd_mergeUniqueImageIntoTrials(conds_master)
% VCD function to merge unique images into stimulus trials where there is a
% left and right parafoveal stimulus or single central stimulus.
%
%   conds_master_reordered = vcd_mergeUniqueImageIntoTrials(conds_master,fix_task_flag)
%
% This function does the following:
% Step 1: match left/right stimuli based on condition
% Step 2: add thickening direction of fixation circle into a column
% Step 3: recreate condition master table with merged unique images, where
% each row is now a single trial.
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
trial_vec = repelem(1:size(trial_vec_i,1),2); 

nr_cueing_conds = 4; %left/right x cued/uncued

% %%%% STEP 1 & 2 %%%%
% Define thickening of central cue
if ~strcmp(unique(conds_master.task_class_name),{'fix'})
    
    thickening_dir = zeros(length(trial_vec),1);
    th_count = 1;
    
    left_cued    = find((conds_master.is_cued == 1) & (conds_master.stimloc == 1));
    left_uncued  = find((conds_master.is_cued == 0) & (conds_master.stimloc == 1));
    right_cued   = find((conds_master.is_cued == 1) & (conds_master.stimloc == 2));
    right_uncued = find((conds_master.is_cued == 0) & (conds_master.stimloc == 2));
    
    assert(isequal(length(left_cued),length(right_uncued)))
    assert(isequal(length(left_uncued),length(right_cued)))
    
    if ismember(unique(conds_master.task_class_name),{'img','ltm'})  

        % double the nr of trials for imagery, but not catch trials, given that we only selected a subset
        left_cued    = [shuffle_concat(left_cued(ismember(left_cued,find(conds_master.is_special_core)))', 2)'; left_cued(ismember(left_cued,find(conds_master.is_catch)))];
        left_uncued  = [shuffle_concat(left_uncued(ismember(left_uncued,find(conds_master.is_special_core)))', 2)'; left_uncued(ismember(left_uncued,find(conds_master.is_catch)))];
        right_cued   = [shuffle_concat(right_cued(ismember(right_cued,find(conds_master.is_special_core)))', 2)'; right_cued(ismember(right_cued,find(conds_master.is_catch)))];
        right_uncued = [shuffle_concat(right_uncued(ismember(right_uncued,find(conds_master.is_special_core)))', 2)';right_uncued(ismember(right_uncued,find(conds_master.is_catch)))];

        % update trial vec and trial vec idx 
        trial_vec_i = NaN(length(left_cued)*2,2);
        trial_vec   = repelem(1:size(trial_vec_i,1),2); 

    end
    
    cuing_conditions = [ones(1,length(left_cued)),  2.*ones(1,length(left_uncued)), ...
                     3.*ones(1,length(right_cued)), 4.*ones(1,length(right_uncued))];
    cond_idx         = shuffle_concat([1:nr_cueing_conds],size(trial_vec_i,1)/nr_cueing_conds);
    
    if length(cond_idx) < length(cuing_conditions)
        cond_idx = cat(2,cond_idx, shuffle_concat(1:diff([length(cond_idx),size(trial_vec_i,1)]),1));
    end
    
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
            thickening_dir(th_count)   = 3;  % both sides thickening
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
conds_master_reordered.unique_trial_nr = trial_vec(1:size(conds_master_reordered,1));
conds_master_reordered.thickening_dir  = thickening_dir(1:size(conds_master_reordered,1));

% check that stim loc is left/right/left/right,etc
assert(isequal(conds_master_reordered.stimloc,repmat([1;2],size(conds_master_reordered.stimloc,1)/2,1)))

% %%%% STEP 3 %%%%

% get table size (rows x cols), and width inside a column
sz = [sum(conds_master_reordered.stimloc==3) + sum(sum(conds_master_reordered.stimloc==[1,2],2))/2, size(conds_master_reordered,2)];

% we assume most columns with stimulus are double width, except for the ones cherry-picked below:
col_width = 2.*ones(1,size(conds_master_reordered,2));
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'stim_class')) = 1;
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'unique_im_nr')) = 1;
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'stimloc')) = 1;
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'task_class')) = 1;
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'is_catch')) = 1;
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'is_cued')) = 1;
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'is_objectcatch')) = 1;
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'task_class_name')) = 1;
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'unique_trial_nr')) = 1;
col_width(strcmp(conds_master_reordered.Properties.VariableNames,'thickening_dir')) = 1;

% copy table headers and create an new table with the same header but scrubbed from its content
conds_master_reordered_merged = vcd_preallocateNaNTable(sz(1), sz(2), conds_master_reordered(1,:), col_width);

% Remove unnecesary columns:
conds_master_reordered_merged.stimloc_name = [];
conds_master_reordered_merged.is_cued = [];

% rename two columns: unique_im_nr->stim_nr_left and stimloc->stim_nr_right
conds_master_reordered_merged.Properties.VariableNames{strcmp(conds_master_reordered_merged.Properties.VariableNames,'unique_im_nr')} = 'stim_nr_left';
conds_master_reordered_merged.Properties.VariableNames{strcmp(conds_master_reordered_merged.Properties.VariableNames,'stimloc')} = 'stim_nr_right';
conds_master_reordered_merged.Properties.VariableNames{strcmp(conds_master_reordered_merged.Properties.VariableNames,'thickening_dir')} = 'is_cued';

varTypes = varfun(@class,conds_master_reordered(1,:),'OutputFormat','cell');


% Loop over columns
for vt = 1:sz(2)
    
    colName = conds_master_reordered.Properties.VariableNames{vt};
    if strcmp(colName,'unique_im_nr')
        colName2 = 'stim_nr_left';
    elseif strcmp(colName,'stimloc')
        colName2 = 'stim_nr_right';
    else
        colName2 = colName;
    end
    
    if ~ismember(colName2,{'stim_nr_right','is_cued','stimloc_name','thickening_dir'})
        
        % Reset trial nr (what row are we allocating)
        trial_nr = 1;
        
        % Merge left/right for each trial
        for tt = 1:sz(1)
            
            if conds_master_reordered.stimloc == 3
                if strcmp(colName2,'stim_nr_left')
                    conds_master_reordered_merged.(colName2)(tt) = conds_master_reordered.unique_im_nr(trial_nr,vt);
                elseif strcmp(colName2,'stim_nr_right')
                    conds_master_reordered_merged.(colName2)(tt) = NaN;
                elseif strcmp(colName2,'is_cued')
                    conds_master_reordered_merged.(colName2)(tt) = conds_master_reordered.thickening_dir(trial_nr,vt);
                else
                    if size(conds_master_reordered_merged.(colName2),2)==2 % double width columns
                        if strcmp(varTypes(vt),'double')
                            conds_master_reordered_merged.(colName2)(tt,:) = ...
                                [table2array(conds_master_reordered(trial_nr,vt)), NaN];
                        elseif strcmp(varTypes(vt),'cell')
                            conds_master_reordered_merged.(colName2)(tt,:) = ...
                                [{table2array(conds_master_reordered(trial_nr,vt))}, {NaN}];
                        end
                    elseif size(conds_master_reordered_merged.(colName2),2)==1 % single width columns
                        if strcmp(varTypes(vt),'double')
                            conds_master_reordered_merged.(colName2)(tt) = ...
                                table2array(conds_master_reordered(trial_nr,vt));
                        elseif strcmp(varTypes(vt),'cell')
                            conds_master_reordered_merged.(colName2)(tt) = ...
                                {table2array(conds_master_reordered(trial_nr,vt))};
                        end
                    end
                end
                trial_nr = trial_nr + 1;
            else
                if strcmp(colName2,'stim_nr_left')
                    assert(isequal(conds_master_reordered.stimloc(trial_nr),1));
                    conds_master_reordered_merged.stim_nr_left(tt) = conds_master_reordered.unique_im_nr(trial_nr);
                    
                    assert(isequal(conds_master_reordered.stimloc(trial_nr+1),2));
                    conds_master_reordered_merged.stim_nr_right(tt) = conds_master_reordered.unique_im_nr(trial_nr+1);
                    
                    assert(isequal(conds_master_reordered.thickening_dir(trial_nr),conds_master_reordered.thickening_dir(trial_nr+1)))
                    conds_master_reordered_merged.is_cued(tt) = conds_master_reordered.thickening_dir(trial_nr);
                else
                    if size(conds_master_reordered_merged.(colName2),2)==2 % double width columns
                        if strcmp(varTypes(vt),'double')
                            conds_master_reordered_merged.(colName2)(tt,:) = ...
                                [table2array(conds_master_reordered(trial_nr,vt)), table2array(conds_master_reordered(trial_nr+1,vt))];
                        elseif strcmp(varTypes(vt),'cell')
                            conds_master_reordered_merged.(colName2)(tt,:) = ...
                                [table2cell(conds_master_reordered(trial_nr,vt)), table2cell(conds_master_reordered(trial_nr+1,vt))];
                        end
                    elseif size(conds_master_reordered_merged.(colName2),2)==1 % single width columns
                        if strcmp(varTypes(vt),'double')
                            assert(isequal(table2array(conds_master_reordered(trial_nr,vt)), table2array(conds_master_reordered(trial_nr+1,vt))));
                            conds_master_reordered_merged.(colName2)(tt) = table2array(conds_master_reordered(trial_nr,vt));
                        elseif strcmp(varTypes(vt),'cell')
                            assert(isequal(table2cell(conds_master_reordered(trial_nr,vt)), table2cell(conds_master_reordered(trial_nr+1,vt))));
                            conds_master_reordered_merged.(colName2)(tt) = table2cell(conds_master_reordered(trial_nr,vt));
                        end
                    end
                end
                trial_nr = trial_nr + 2;
            end
        end
    end
end




% Randomize catch trial loc
sz = size(conds_master_reordered_merged);
catch_trials = find(conds_master_reordered_merged.is_catch(:,1)==true);
if ~isempty(catch_trials)
    
    while 1
        shuffle_trial_idx = randi(size(conds_master_reordered_merged,1),length(catch_trials)); % get a new trial number
        if abs(diff(shuffle_trial_idx)) ~= 1
            break;
        end
        if shuffle_trial_idx ~= catch_trials % we don't want to place it back in the same position
            break;
        end
    end
    
    for cc = 1:length(shuffle_trial_idx) % 11
        tmp = conds_master_reordered_merged;
        % grab catch trial
        curr_catch_trial   = tmp(catch_trials(cc),:);  %  13 will become new number, e.g., 6
        
        if shuffle_trial_idx(cc) < catch_trials(cc)
            conds_shuffle1_pre1  = tmp(1:(shuffle_trial_idx(cc)-1),:);  % e.g., 1-5
            conds_shuffle1_post1 = tmp(shuffle_trial_idx(cc):(catch_trials(cc)-1),:); % e.g., 6-12
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

% update trial nr, move more general columns to the left of the table
assert(isequal(sort(conds_master_reordered_merged.unique_trial_nr), [1:size(conds_master_reordered_merged.unique_trial_nr,1)]'))
conds_master_reordered_merged.unique_trial_nr = sort(conds_master_reordered_merged.unique_trial_nr);
[~,~,idx0] = intersect({'unique_trial_nr','stim_class','task_class','stim_nr_left','stim_nr_right','is_cued','is_catch','stim_class_name','task_class_name'}, ...
    conds_master_reordered_merged.Properties.VariableNames, 'stable');
newOrder = [idx0', setdiff([1:length(conds_master_reordered_merged.Properties.VariableNames)],idx0)];

tmp = conds_master_reordered_merged(:,newOrder);
conds_master_reordered_merged = tmp;

