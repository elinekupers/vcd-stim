function conds_master_reordered = create_condition_master_trials(conds_master,fix_task_flag)

% INPUTS:
% conds_master  :       col 1 has unique stim nr (int)
%                       col 2 has stim loc info (1=left,2=right)
%                       col 3 has stim cue status (0=uncued, 1=cued).
%                       col 4:end has stim feature info
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

if nargin < 2 || ~exist('fix_task_flag','var') || isempty(fix_task_flag)
    fix_task_flag = false;
end

left_stim = find(conds_master(:,2)==1);
right_stim = find(conds_master(:,2)==2);

assert(length(left_stim)==length(right_stim))

% pair left / right stimuli for each trial
trial_vec_i = NaN(size(conds_master,1)/2,2);

% add the condition master to struct
n_unique_cases = numel(unique(conds_master(:,1)));
trial_vec = repelem(1:size(trial_vec_i,1),2); % or     trial_vec = repelem(1:(n_unique_cases/2),2); ????


% Define thickening of central cue
if ~fix_task_flag
    
    thickening_dir = zeros(length(trial_vec),1);
    th_count = 1;
    
    left_cued    = find((conds_master(:,3) == 1) & (conds_master(:,2) == 1));
    left_uncued  = find((conds_master(:,3) == 0) & (conds_master(:,2) == 1));
    right_cued   = find((conds_master(:,3) == 1) & (conds_master(:,2) == 2));
    right_uncued = find((conds_master(:,3) == 0) & (conds_master(:,2) == 2));
    
    cond_idx = shuffle_concat([1:4],size(trial_vec_i,1)/4);  
    
    
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
    


elseif fix_task_flag

    thickening_dir = 3.*ones(size(trial_vec))'; % we thicken at both sides

    for ii = 1:length(conds_master(:,3))
        stim_loc = conds_master(ii,2);
        stim_to_allocate = find(isnan(trial_vec_i(:,stim_loc)));
        trial_vec_i(stim_to_allocate(1), stim_loc) = ii;
    end
    
    trial_vec_i = trial_vec_i';
    trial_vec_i = trial_vec_i(:);

end
  

conds_master_reordered = [trial_vec', thickening_dir, conds_master(trial_vec_i,:)];

end


%     
%     for tt = 1:2:length(trial_vec)
%         curr_trial = trial_vec(tt);
%         
%         if curr_trial == 1 % if first trial
%             leading_im = 1; % just pick the first one from the list of unique images
%         else
%             [leading_im,leading_im_i] = min([left_stim(1),right_stim(1)]); % or the first one that comes next
%         end
%         
%         if conds_master(leading_im,2) == 1 % if leading unique stim happens on the left
%             trial_vec_i(curr_trial,1) = leading_im; % set it left
%             assert(leading_im==left_stim(1))
%             left_stim(1) = [];
%         elseif conds_master(leading_im,2) == 2 % if leading unique  stim happens on the right
%             trial_vec_i(curr_trial,2) = leading_im; % set right
%             assert(leading_im==right_stim(1))
%             right_stim(1) = [];
%         end
%         
%         trial_im = trial_vec_i(curr_trial,:);
%         defined_side = ~isnan(trial_vec_i(curr_trial,:));
%         side_to_fill = isnan(trial_vec_i(curr_trial,:));
%         
%         if conds_master(trial_im(defined_side),2) == 1 % if stim is on the left
%             
%             if conds_master(trial_im(defined_side),3) == 1 % and if stim is cued
%                 thickening_dir(tt) = 1; % left
%                 
%                 counter = 1;
%                 while 1
%                     
%                     % then we want a right stim, that is uncued
%                     tmp = right_stim(counter);
%                     
%                     if conds_master(tmp,3) == 0
%                         trial_vec_i(curr_trial,find(side_to_fill)) = tmp;
%                         right_stim(counter) = [];
%                         break;
%                     else
%                         counter = counter+1;
%                     end
%                 end
%                 
%             elseif conds_master(trial_im(defined_side),3) == 0 %  if stim is uncued
%                 thickening_dir(tt) = 2; % If stim is uncued, the opposite rim side thickens (here right)
%                 
%                 counter = 1;
%                 while 1
%                     
%                     % then we want a right stim, that is cued
%                     tmp = right_stim(counter);
%                     
%                     if conds_master(tmp,3) == 1
%                         trial_vec_i(curr_trial,find(side_to_fill)) = tmp;
%                         right_stim(counter) = [];
%                         break;
%                     else
%                         counter = counter+1;
%                     end
%                 end
%                 
%             end
%             
%         elseif conds_master(trial_im(defined_side),2) == 2 % if stim is on the right
%             
%             if conds_master(trial_im(defined_side),3) == 1 % and if stim is cued
%                 thickening_dir(tt) = 2; % right
%                 
%                 counter = 1;
%                 while 1
%                     
%                     % then we want a left stim, that is uncued
%                     tmp = left_stim(counter);
%                     
%                     if conds_master(tmp,3) == 0 % if considered stim is uncued
%                         trial_vec_i(curr_trial,find(side_to_fill)) = tmp; % then we pair it
%                         left_stim(counter) = []; % and remove it from the list
%                         break;
%                     else
%                         counter = counter+1; % otherwise we go to the next
%                     end
%                 end
%                 
%             elseif conds_master(trial_im(defined_side),3) == 0 %  if stim is uncued
%                 thickening_dir(tt) = 1; % left
%                 
%                 counter = 1;
%                 while 1
%                     
%                     % then we want a right stim, that is cued
%                     tmp = left_stim(counter);
%                     
%                     if conds_master(tmp,3) == 1 % if  considered stim is cued
%                         trial_vec_i(curr_trial,find(side_to_fill)) = tmp;
%                         left_stim(counter) = [];
%                         break;
%                     else % if not cued, let's look at the next stim
%                         counter = counter+1;
%                     end
%                 end
%                 
%             end
%             
%         end
%         thickening_dir(tt+1) = thickening_dir(tt); % copy for completeness
%     end
%     clear counter
%     
% elseif fix_task_flag
%     thickening_dir = 3.*ones(size(trial_vec))'; % we thicken at both sides
% 
%     for tt = 1:2:length(trial_vec)
%         counter = 1;
%         curr_trial = trial_vec(tt);
%         
%         if curr_trial == 1 % if first trial
%             leading_im = 1; % just pick the first one from the list of unique images
%         else
%             [leading_im,leading_im_i] = min([left_stim(1),right_stim(1)]); % or the first one that comes next
%         end
%         
%         if conds_master(leading_im,2) == 1 % if leading unique stim happens on the left
%             trial_vec_i(curr_trial,1) = leading_im; % set it left
%             assert(leading_im==left_stim(1))
%             left_stim(1) = [];
%             
%         elseif conds_master(leading_im,2) == 2 % if leading unique  stim happens on the right
%             trial_vec_i(curr_trial,2) = leading_im; % set right
%             assert(leading_im==right_stim(1))
%             right_stim(1) = [];
%             
%         end
%         
%         trial_im = trial_vec_i(curr_trial,:);
%         defined_side = ~isnan(trial_vec_i(curr_trial,:));
%         side_to_fill = isnan(trial_vec_i(curr_trial,:));
%         
%         if conds_master(trial_im(defined_side),2) == 1 % if stim is on the left
%             
%             % then we want a right stim, that is uncued
%             trial_vec_i(curr_trial,find(side_to_fill))  = right_stim(counter);
%             right_stim(counter) = [];
%             
%         elseif conds_master(trial_im(defined_side),2) == 2 % if stim is on the right
%             
%             trial_vec_i(curr_trial,find(side_to_fill))  = left_stim(counter);
%             left_stim(counter) = []; % and remove it from the list
%             
%         end
%     end
%     clear counter
%     
% end
% 
% % Check if we used all stimuli
% assert(isempty(left_stim));
% assert(isempty(right_stim));




