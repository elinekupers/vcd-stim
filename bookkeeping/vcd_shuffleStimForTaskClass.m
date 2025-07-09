function master_table_out = vcd_shuffleStimForTaskClass(params, task_class_name_to_shuffle, master_table, nr_of_trials_per_block, session_type)
%% VCD function to combine and shuffle stimulus classes for SCC/LTM trials
%   master_table_out = vcd_shuffleStimForTaskClass(params, task_class_name_to_shuffle, master_table, nr_of_trials_per_block)
%
% Purpose: SCC and LTM crossings requires a mixture of stimulus classes:
% SCC task asks subjects to judge the stimulus class of the shown
% image. LTM task asks subjects to recall the learned pair associate.
% 
% Given that the master table is currently sorted per stimulus class, we
% need to reorganize the stimuli in scc/ltm blocks to have a (balanced)
% mixture of images from different classes. This function will do so, while
% making sure that we comply to even/uneven left/right trial order, as well
% as showing each stimulus class at least once per block. 
% NOTE: fully-crossed and balanced stimulus classes is not possible,
% because we have different number of unique images per stimulus class.
%
% INPUTS:
% * params                      : (struct) stimulus and experimental design
%                                   parameters. We need:
%     * params.stim.unique_im_nrs for each stimulus class to do some checks,
%     * params.exp.stimclassnames to get a list of stimulus classes.
% * task_class_name_to_shuffle  : (str) name of task class we want to
%                                   shuffle the stimuli (either 'scc' or
%                                   'ltm')
% * master_table                : (table) giant table with stimulus and trial definitions
% * nr_of_trials_per_block      : (int) nr of trials for each scc or ltm block.
%                                   Expecting either 8 or 4 trials/block.
% * session_type                : (str) is this a behavioral or mri
%                                   session? choose from 'MRI' or
%                                   'BEHAVIOR'
%
% OUTPUTS:
% * master_table_out            : (table) giant table, where scc blocks 
%                                  contain a mixture of gabor, rdk, dot,and
%                                  obj stimuli.
%
%% Check inputs
% we only want to shuffle scc and ltm task crossings
assert(any(strcmp(task_class_name_to_shuffle,{'scc','ltm'})));

% scc task crossings have 8 trials/block, ltm task crossings have 4 trials/block
assert(isequal(nr_of_trials_per_block, choose(strcmp(task_class_name_to_shuffle,'scc'), params.exp.block.n_trials_single_epoch,params.exp.block.n_trials_double_epoch)));


%% Get nr of cueing directions
cued_stim_loc = unique(master_table.is_cued(strcmp(master_table.task_class_name,{task_class_name_to_shuffle})));
    
% Get all the rows with corresponsing stim class name
all_task_rows = strcmp(master_table.task_class_name,{task_class_name_to_shuffle});

% get all indices for left and right stimuli, separate for
% cued/uncued/ncued conditions
for kk = 1:length(cued_stim_loc) % (1:left, 2:right, 3:neutral)

    % get task trials (filter out catch trials)
    m0_idx          = master_table.is_catch==0 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk);
    m0_sub          = find(m0_idx); % keep track of row nr
    
    if cued_stim_loc(kk) == 3 % scenes
        m0_center_img_nr  = master_table.stim_nr_left(m0_idx);
        
        % shape into an 8 trials x m blocks matrix, such that the
        % columns contain different stimulus classes
        n_leftovers = mod(size(m0_center_img_nr,1),nr_of_trials_per_block);
        if n_leftovers > 0
            leftovers_img_nr = m0_center_img_nr((end-(n_leftovers-1)):end);
            leftovers_sub    = m0_sub((end-(n_leftovers-1)):end);
            tmp_img_nr_center_rz  = reshape(m0_center_img_nr(1:(size(m0_center_img_nr,1)-n_leftovers)),[],nr_of_trials_per_block);
            tmp_sub_rz            = reshape(m0_sub(1:(size(m0_sub,1)-n_leftovers)),[],nr_of_trials_per_block);
        else
            tmp_img_nr_center_rz  = reshape(m0_center_img_nr,[],nr_of_trials_per_block);
            tmp_sub_rz            = reshape(m0_sub,[],nr_of_trials_per_block);            
        end
        % shuffle the order across columns such that gabors, rdks, dots, obj's get mixed
        shuffled_ind_center  = nr_of_trials_per_block*[0:(size(tmp_img_nr_center_rz,1)-1)]' ...
            + reshape(shuffle_concat([1:nr_of_trials_per_block],size(tmp_img_nr_center_rz,1)),nr_of_trials_per_block,[])';

        % transpose as matlab indexes over cols first
        tmp_img_nr_center1     = tmp_img_nr_center_rz';
        shuffled_ind_center1   = shuffled_ind_center';
        tmp_sub_rz1            = tmp_sub_rz';

        % shuffle trials based on shuffled indices
        tmp_img_nr_center2   = tmp_img_nr_center1(shuffled_ind_center1);
        tmp_sub_rz2_center2  = tmp_sub_rz1(shuffled_ind_center1);
        
        assert(isequal(sort(tmp_img_nr_center2,1),sort(tmp_img_nr_center1,1))) % sorted shuffled and sorted unshuffled left should be the same
        assert(~isequal(tmp_img_nr_center2,tmp_img_nr_center1)) % shuffled and unsort should not be the same
        
        img_nr_shuffle_center  = tmp_img_nr_center2(:); % image nrs
        sub_shuffle_center     = tmp_sub_rz2_center2(:); % row nrs
        
        % add leftovers at the end
        if n_leftovers > 0
            img_nr_shuffle_center = cat(1,img_nr_shuffle_center,leftovers_img_nr);
            sub_shuffle_center = cat(1,sub_shuffle_center,leftovers_sub);
        end
        
        % record catch trials
        if sum(master_table.is_catch==1) > 1
            m0_center_catch  = master_table.stim_nr_left(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
            sub_center_catch = find(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
        end
        
    else % gabor/rdk/dot/obj

        m0_left_img_nr  = master_table.stim_nr_left(m0_idx);
        m0_right_img_nr = master_table.stim_nr_right(m0_idx);
        
%         m0_left_catch   = master_table.stim_nr_left(catch_sub);
%         m0_right_catch  = master_table.stim_nr_right(catch_sub);
%         n_catch_l       = length(m0_left_catch);
%         n_catch_r       = length(m0_right_catch);
%         
        % try to shape nr of trials into an 8 trials x m blocks matrix, 
        % such that the columns contain different stimulus classes. If we
        % have a nr of trials that can't be divided by 8, we call them
        % "leftovers" and string them along
        n_leftovers_l = mod(size(m0_left_img_nr,1),nr_of_trials_per_block);
        n_leftovers_r = mod(size(m0_right_img_nr,1),nr_of_trials_per_block);
        
        if (n_leftovers_l > 0) && (n_leftovers_r > 0)
            leftovers_img_nr_l = m0_left_img_nr(end-(n_leftovers_l-1));
            leftovers_sub_l    = m0_sub(end-(n_leftovers_l-1));
            leftovers_img_nr_r = m0_right_img_nr(end-(n_leftovers_r-1));
            leftovers_sub_r    = m0_sub(end-(n_leftovers_r-1));
            tmp_img_nr_left_rz    = reshape(m0_left_img_nr(1:(size(m0_left_img_nr,1)-n_leftovers_l)),[],nr_of_trials_per_block);
            tmp_img_nr_right_rz   = reshape(m0_right_img_nr(1:(size(m0_right_img_nr,1)-n_leftovers_r)),[],nr_of_trials_per_block);
        else
            tmp_img_nr_left_rz    = reshape(m0_left_img_nr,[],nr_of_trials_per_block);
            tmp_img_nr_right_rz   = reshape(m0_right_img_nr,[],nr_of_trials_per_block);
        end
        
        % keep track of row nr 
        tmp_sub_rz            = reshape(m0_sub,[],nr_of_trials_per_block);
   
        % set a limit to the number of trial shuffles to avoid infinite loop
        fprintf('[%s]: Shuffle trials and try to get equally distributed nr of stimulus classes across trials\n',mfilename);
        attempts = 0;
        
        while 1
            attempts = attempts+1;
            
            if attempts > 10000
                error('[%s]: Can''t find a solution even though we tried more than 10,000 times! Aborting!', mfilename)
            end
            % shuffle the order across columns such that gabors, rdks, dots, obj's get mixed
            shuffled_ind_left  = nr_of_trials_per_block*[0:(size(tmp_img_nr_left_rz,1)-1)]' ...
                + reshape(shuffle_concat([1:nr_of_trials_per_block],size(tmp_img_nr_left_rz,1)),nr_of_trials_per_block,[])';
            
            shuffled_ind_right = nr_of_trials_per_block*[0:(size(tmp_img_nr_right_rz,1)-1)]' ...
                + reshape(shuffle_concat([1:nr_of_trials_per_block],size(tmp_img_nr_right_rz,1)),nr_of_trials_per_block,[])';
            
            assert(sum(sum(ismember(tmp_img_nr_left_rz,tmp_img_nr_right_rz)))==0);
            
            % transpose as matlab indexes over cols first
            tmp_img_nr_left1       = tmp_img_nr_left_rz';
            tmp_img_nr_right1      = tmp_img_nr_right_rz';
            shuffled_ind_left1     = shuffled_ind_left';
            shuffled_ind_right1    = shuffled_ind_right';
            assert(sum(sum(ismember(tmp_img_nr_left1,tmp_img_nr_right1)))==0);
            
            % keep track of row nr
            tmp_sub_rz1 = tmp_sub_rz';
            
            % shuffle trials based on shuffled indices
            tmp_img_nr_left2   = tmp_img_nr_left1(shuffled_ind_left1);
            tmp_img_nr_right2  = tmp_img_nr_right1(shuffled_ind_right1);
            assert(sum(sum(ismember(tmp_img_nr_left2,tmp_img_nr_right2)))==0);
            
            % keep track of row nr for left and right
            tmp_sub_rz2_l      = tmp_sub_rz1(shuffled_ind_left1);
            tmp_sub_rz2_r      = tmp_sub_rz1(shuffled_ind_right1);
            
            assert(isequal(sort(tmp_img_nr_left2,1),sort(tmp_img_nr_left1,1))) % sorted shuffled and sorted unshuffled left should be the same
            assert(isequal(sort(tmp_img_nr_right2,1),sort(tmp_img_nr_right1,1))) % sorted shuffled and sorted unshuffled right should be the same
            assert(~isequal(tmp_img_nr_left2,tmp_img_nr_left1)) % shuffled and unsort should not be the same
            assert(~isequal(tmp_img_nr_right2,tmp_img_nr_right1)) % shuffled and unsort should not be the same
            
            trial_order = [tmp_sub_rz2_l(:),tmp_sub_rz2_r(:)];

            % check if stimclasses are balanced 
            [~,stimclss_idx] = ismember(master_table.stim_class_name(trial_order),params.exp.stimclassnames);
            chance_of_stimclass = histcounts(stimclss_idx)./sum(histcounts(stimclss_idx));
            dt = ceil(max(1./chance_of_stimclass));
            n_start = 1:dt:size(trial_order,1);
            sample_ok = [];
            for nn = 1:length(n_start)
                if n_start(nn)+(dt-1) > size(trial_order,1)
                    curr_sample = stimclss_idx(n_start(nn):end,:);
                else
                    curr_sample = stimclss_idx(n_start(nn):(n_start(nn)+(dt-1)),:);
                end
                unique_stimclass = unique(curr_sample(:));
                if length(unique_stimclass)==length(unique(stimclss_idx(:)))
                    sample_ok(nn) = true;
                else
                    sample_ok(nn) = false;
                    break;
                end
            end
            if sum(sample_ok) == length(n_start)
                break;
            end
        end

        img_nr_shuffle_left(kk,:)  = tmp_img_nr_left2(:); % this will fail if we don't have equal nr of trials per stimloc/cueing status
        img_nr_shuffle_right(kk,:) = tmp_img_nr_right2(:); % this will fail if we don't have equal nr of trials per stimloc/cueing status
        assert(sum(sum(ismember(img_nr_shuffle_left,img_nr_shuffle_right)))==0);
        
        sub_shuffle_left(kk,:)  = tmp_sub_rz2_l(:);
        sub_shuffle_right(kk,:) = tmp_sub_rz2_r(:);
        
        % add leftovers at the end if there are any
        if (n_leftovers_l > 0) && (n_leftovers_r > 0)
            img_nr_shuffle_left(kk,:) = cat(1,img_nr_shuffle_left(kk,:),leftovers_img_nr_l);
            img_nr_shuffle_right(kk,:) = cat(1,img_nr_shuffle_right(kk,:),leftovers_img_nr_r);
            sub_shuffle_left(kk,:) = cat(1,sub_shuffle_left(kk,:),leftovers_sub_l);
            sub_shuffle_right(kk,:) = cat(1,sub_shuffle_right(kk,:),leftovers_sub_r);
        end

        % record catch trials
        if sum(master_table.is_catch==1) > 1
            img_nr_left_catch(kk,:)  = master_table.stim_nr_left(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
            img_nr_right_catch(kk,:) = master_table.stim_nr_right(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
            sub_shuffle_left_catch(kk,:)   = find(img_nr_left_catch(kk,:));
            sub_shuffle_right_catch(kk,:)  = find(img_nr_right_catch(kk,:));
        end

    end
end


% GOAL: 
% Versions A and B represent left cued/right uncued and right cued/left
% uncued. We combine these conditions first, then ensure left and right cued
% trials are mixed within a block, and then randomly assign the
% combined trials to the master table.

% sub_shuffle_left is indexing the table rows (not the actual image nrs)
% Note: sub_shuffle_left dim 1 order: left cued, right cued
%       sub_shuffle_right dim 1 order: left cued, right cued
if all(ismember(cued_stim_loc,[1,2]))
    combined_trial_shuffleA = [sub_shuffle_left(1,:); sub_shuffle_right(1,:)]'; % shuffled left cued, shuffled right uncued (N/2 trials x 2 loc (l/r))
    combined_trial_shuffleB = [sub_shuffle_left(2,:); sub_shuffle_right(2,:)]'; % shuffled left uncued, shuffled right cued (N/2 trials x loc (l/r))
    combined_trial_shuffleAB = cat(1, combined_trial_shuffleA, combined_trial_shuffleB); % N trials x 2 loc (l/r) (first half is cued left, second half is cued right)
    
    % Do the same for image number and get corresponding button responses
    combined_img_nr_shuffleA = [img_nr_shuffle_left(1,:);img_nr_shuffle_right(1,:)]'; % N/2 trials x 2 loc (l/r)
    combined_img_nr_shuffleB = [img_nr_shuffle_left(2,:);img_nr_shuffle_right(2,:)]'; % N/2 trials x 2 loc (l/r)
    combined_img_nr_shuffleAB = cat(1, combined_img_nr_shuffleA, combined_img_nr_shuffleB); % N trials x 2 loc (l/r) (first half is cued left, second half is cued right)
    
    % We want the number of cued left and cued right to be balanced within a block. 
    % right now we have combined_trial_shuffleAB:
        % * first column [1:N/2] = left cued, [(N/2)+1 : N] = left uncued
        % * second column [1:N/2] = right uncued, [(N/2)+1 : N] = right cued
        
    stimclass_left_cued  = vcd('stimtostimclassnumber',combined_img_nr_shuffleAB(1:floor(size(combined_img_nr_shuffleAB,1)/2),1));
    stimclass_right_cued = vcd('stimtostimclassnumber',combined_img_nr_shuffleAB((ceil(size(combined_img_nr_shuffleAB,1)/2)+1):end,2));
    correct_response = [stimclass_left_cued,stimclass_right_cued];
    cued_side        = [ones(size(stimclass_left_cued)),2*ones(size(stimclass_right_cued))];
    
    % for LTM: check that column lengths is the same as sum of repeats * special core
    % stimuli
    if strcmp(session_type, 'MRI')
        unique_trial_repeats = params.exp.n_unique_trial_repeats_mri;
    else
        unique_trial_repeats = params.exp.n_unique_trial_repeats_behavior;
    end
    
    if strcmp(task_class_name_to_shuffle,{'ltm'})
        assert(isequal(size(combined_trial_shuffleAB,1), ...
            unique(unique_trial_repeats(1:4,strcmp(task_class_name_to_shuffle,params.exp.taskclassnames)))  * ... 23 repeats
            (2*(length(params.stim.all_specialcore_im_nrs)-length(params.stim.ns.unique_im_nrs_specialcore))))); % 8*4 (=32) * 2 (double nr of repeats)
    elseif strcmp(task_class_name_to_shuffle,{'scc'})
        assert(isequal(size(combined_trial_shuffleAB,1), ...
            sum(unique_trial_repeats(1:4,strcmp(task_class_name_to_shuffle,params.exp.taskclassnames))'  .* ... nr of repeats
            [length(params.stim.gabor.unique_im_nrs_core),length(params.stim.rdk.unique_im_nrs_core),... nr of unique core images
            length(params.stim.dot.unique_im_nrs_core),length(params.stim.obj.unique_im_nrs_core) ]))); %
    else
        error('[%s]: Unique nr of stimulus classes doesn''t match expected number', mfilename)
    end
    
    stim_class_nrs = unique(correct_response);
    abs_stim_class_chance = chance_of_stimclass*params.exp.block.n_trials_single_epoch;
    
    max_attempts = 100000;
    attempt = 0;
    while 1
        attempt = attempt + 1;
        % Shuffle across all trials, using the following constraints:
        % - we want 50/50 left/right cued stimuli
        % - we want at least one of each stimulus class per block
        % - we want somewhat equal distribution of button responses.. (we
        % allow for a discrepancy of 2 or less across stim classes)
        shuffle_vec = [shuffle_concat(1:size(combined_trial_shuffleAB,1)/2,1); ... [1:(N/2) shuffled; (N/2):end shuffled]
                    (size(combined_trial_shuffleAB,1)/2) + shuffle_concat(1:size(combined_trial_shuffleAB,1)/2,1)]; % N trials (1-8 indices)
        shuffle_vec = shuffle_vec(:);
        cued_side_shuffled = cued_side(shuffle_vec'); 
        assert(isequal(cued_side_shuffled,repmat([1,2],1,length(cued_side_shuffled)/2))) % should be [1,2,1,2,1,2, etc]
        shuffled_responses = correct_response(shuffle_vec);
        block_start = 1:(params.exp.block.n_trials_single_epoch):length(shuffled_responses);
        
        % reset ocounters
        block_ok = zeros(1,length(block_start)); counter = 1; scc_ok = false;
        for cc = block_start
            n1 = histcounts(shuffled_responses(cc:(cc+(params.exp.block.n_trials_single_epoch-1))),[stim_class_nrs, stim_class_nrs(end)+1]); % check stim class within a block
            if all(n1>=1) % if we have at least one of each stimulus class per block, we are happy
                block_ok(counter) = 1;
                counter = counter + 1;
            elseif any(n1<1) 
                block_ok(counter) = 0;
                break;
            end
        end
        
        if ~isempty(block_ok) && sum(block_ok)==length(block_start)
            % if within block check completed, check is we have the
            % expected distribution of stimulus classes across all blocks
            % (given the unequal distribution of stimulus classes, i.e., we
            % have more gabors and rdks than dots and objects).
            n2 = histcounts(shuffled_responses,[stim_class_nrs, stim_class_nrs(end)+1]); % check stim class across all blocks.
            if isequal(n2,abs_stim_class_chance*(length(block_start)))
                scc_ok = true;
            else
                scc_ok = false;
            end
            
            % we also check if any two dot stimuli are considered too close
            % or not:
            updated_stim_im_nr = combined_img_nr_shuffleAB(shuffle_vec,:);
            updated_stim_im_idx = combined_trial_shuffleAB(shuffle_vec,:);
            left_stimclass_updated_stim_nr = vcd('stimtostimclassnumber',updated_stim_im_nr(:,1));
            right_stimclass_updated_stim_nr = vcd('stimtostimclassnumber',updated_stim_im_nr(:,2));
            stimclass_updated_stim_nr = [left_stimclass_updated_stim_nr; right_stimclass_updated_stim_nr]';
            dot_idx = find(sum(stimclass_updated_stim_nr==3,2)==2);
            diff_orient = abs(circulardiff(master_table.orient_dir(updated_stim_im_idx(dot_idx,1),1),master_table.orient_dir(updated_stim_im_idx(dot_idx,2),2),360));
            % if difference between dots is less than we allow, we reset
            % ssc_ok flag to false. If angles are ok, then we leave the
            % state of the scc_flag as set above.
            if any(diff_orient<=params.stim.dot.min_ang_distance_dot)
                scc_ok = false; 
            end
        end
        
        if scc_ok
            break;
        end
        
        if attempt > max_attempts
            error('\n[%s]: Can''t reach a shuffle that works with current constraints!',mfilename)
        end
        
        if mod(attempt,100)
            fprintf('.');
        end
    end
    fprintf('\n');
    combined_trial_shuffleAB3 = combined_trial_shuffleAB(shuffle_vec,:);
    combined_img_nr_shuffleAB3 = combined_img_nr_shuffleAB(shuffle_vec,:);
end

if strcmp(task_class_name_to_shuffle,{'ltm'})
    if exist('sub_shuffle_center','var')
        assert(isequal(size(sub_shuffle_center,1), (length(params.stim.ns.unique_im_nrs_specialcore)*unique_trial_repeats(5,strcmp(task_class_name_to_shuffle,params.exp.taskclassnames)))))
        % add column of nans to match two column structure for left/right
        % matrix
        sub_shuffle_center2 = cat(2,sub_shuffle_center,NaN(size(sub_shuffle_center)));
        
        combined_trial_shuffleABC  = NaN(size(combined_trial_shuffleAB3,1) + size(sub_shuffle_center2,1),2);
        combined_cueloc_shuffleABC = NaN(size(combined_trial_shuffleAB3,1) + size(sub_shuffle_center2,1),2);
        combined_im_nr_shuffleABC  = NaN(size(combined_trial_shuffleAB3,1) + size(sub_shuffle_center2,1),2);
        insert_scenes_here   = datasample(1:size(combined_trial_shuffleABC,1),size(sub_shuffle_center2,1), 'Replace',false);
        insert_classics_here = setdiff(1:size(combined_trial_shuffleABC,1),insert_scenes_here);
        combined_trial_shuffleABC(insert_scenes_here,:)   = sub_shuffle_center2;
        combined_trial_shuffleABC(insert_classics_here,:) = combined_trial_shuffleAB3;
        
        combined_cueloc_shuffleABC(insert_scenes_here,:)   = 3;
        combined_cueloc_shuffleABC(insert_classics_here,:) = cued_side_shuffled;
        
        combined_im_nr_shuffleABC(insert_scenes_here,:)   = img_nr_shuffle_center;
        combined_im_nr_shuffleABC(insert_classics_here,:) = combined_img_nr_shuffleAB3;
    else
        combined_trial_shuffleABC = combined_trial_shuffleAB3;
        combined_cueloc_shuffleABC = cued_side_shuffled;
        combined_im_nr_shuffleABC = combined_img_nr_shuffleAB3;
    end

    % shuffle catch trials
    if sum(master_table.is_catch==1) > 1
        shuffle_vec_catch = shuffle_concat(sub_center_catch(1,:),1);
        catch_table       = catch_table.stim_class_name(shuffle_vec_catch(:,1),1);
    end
else
    combined_trial_shuffleABC  = combined_trial_shuffleAB3;
    combined_cueloc_shuffleABC = cued_side_shuffled;
    combined_im_nr_shuffleABC  = combined_img_nr_shuffleAB3;
    
    % Shuffle catch trials
    if sum(master_table.is_catch==1) > 1
        shuffle_vec_catch = [shuffle_concat(sub_shuffle_left_catch(1,:),1); shuffle_concat(sub_shuffle_right_catch(2,:),1)];
        catch_table1 = master_table(shuffle_vec_catch(:,1),1);
        catch_table2 = master_table(shuffle_vec_catch(:,2),2);
    end
end
  

%% APPLY THE SHUFFLE!
shuffled_master_table = master_table([],:);
for xx = 1:size(master_table(1,:),2); col_widths(xx) = size(table2array(master_table(1,xx)),2); end
colNames = master_table.Properties.VariableNames;

double_width_cols = find(col_widths==2);

for ii = 1:size(combined_trial_shuffleABC,1)
    
    % Get left/center stimulus of this trial
    new_trial_l = master_table(combined_trial_shuffleABC(ii,1),:);
    shuffled_master_table(ii,:) = new_trial_l;
    
    % update block_nr, trial_nr, stim_class_nr, repeat nr according to new image/trial order
    shuffled_master_table.unique_trial_nr(ii) = ii;
    shuffled_master_table.stim_class(ii) = 99;
    shuffled_master_table.repeat_nr(ii)  = 99;
    shuffled_master_table.stim_class_unique_block_nr(ii) = ceil(ii/nr_of_trials_per_block);
    shuffled_master_table.trial_nr(ii) = mod(ii-1,nr_of_trials_per_block)+1;
    
    % check cuing & stim nr
    assert(isequal(new_trial_l.is_cued, combined_cueloc_shuffleABC(ii)));
    assert(isequal(new_trial_l.stim_nr_left, combined_im_nr_shuffleABC(ii,1)));
    
    % Get right stimulus of this trial if there is one
    if ~isnan(combined_trial_shuffleABC(ii,2))
        new_trial_r = master_table(combined_trial_shuffleABC(ii,2),:);
        assert(isequal(new_trial_l.is_cued,  new_trial_r.is_cued));
        assert(isequal(new_trial_l.is_catch, new_trial_r.is_catch));
        assert(isnan(new_trial_l.is_objectcatch)); assert(isnan(new_trial_r.is_objectcatch));
        assert(isnan(new_trial_l.cd_start)); assert(isnan(new_trial_r.cd_start));
        assert(isequal(new_trial_r.stim_nr_right, combined_im_nr_shuffleABC(ii,2)));
        shuffled_master_table.stim_nr_right(ii) = new_trial_r.stim_nr_right;
    else
        shuffled_master_table.stim_nr_right(ii) = NaN;
    end
    
    % Combine double columns
    for jj = 1:length(double_width_cols)
        colName = colNames{double_width_cols(jj)};
        shuffled_master_table.(colName)(ii,:) = [new_trial_l.(colName)(1), new_trial_r.(colName)(2)];
    end
end

% Check that we didn't loose any unique image nrs
for sc = 1:length(params.exp.stimclassnames)
    stimClass = params.exp.stimclassnames{sc};
    assert(all(ismember(shuffled_master_table.stim_nr_left(strcmp(shuffled_master_table.stim_class_name(:,1),stimClass)),params.stim.(stimClass).unique_im_nrs_core)))
    if sc < 5
        assert(all(ismember(shuffled_master_table.stim_nr_right(strcmp(shuffled_master_table.stim_class_name(:,2),stimClass)),params.stim.(stimClass).unique_im_nrs_core)))
    end
end
% Check that left and right image nrs are not overlapping
assert(sum(shuffled_master_table.stim_nr_left==shuffled_master_table.stim_nr_right)==0)

% Check that we have equal nr of left/right cuing conditions
assert(sum(shuffled_master_table.stim_nr_left==shuffled_master_table.stim_nr_right)==0)


%% CALCULATE REPEATS 

% Compute how often each new combination of stimuli x cue status is repeated
[uniquerow, ~, rowidx] = unique([shuffled_master_table.stim_nr_left, shuffled_master_table.stim_nr_right, shuffled_master_table.is_cued], 'rows'); 
nr_of_ccurrences = accumarray(rowidx, 1);

repeat_nr = zeros(1,length(nr_of_ccurrences));
for nrow = 1:size(shuffled_master_table,1)
    
    curr_comb = [shuffled_master_table.stim_nr_left(nrow), shuffled_master_table.stim_nr_right(nrow), shuffled_master_table.is_cued(nrow)];
    if isnan(curr_comb(2))
        [~,comb_nr] = intersect(uniquerow(:,[1,3]),curr_comb([1,3]),'rows');
    else
        [~,comb_nr] = intersect(uniquerow,curr_comb,'rows');
    end
    repeat_nr(comb_nr) = repeat_nr(comb_nr)+1;
    shuffled_master_table.repeat_nr(nrow) =  repeat_nr(comb_nr);
end

% copy master table without current task rows
master_table2     = master_table(~all_task_rows, :);

if sum(master_table.is_catch==1) > 1
    % Add catch trials back into shuffled_master_table
    shuffled_master_table2 = cat(1,shuffled_master_table,catch_table);
    shuffle_w_catch = randperm((size(shuffled_master_table,1)+1):(size(shuffled_master_table,1)+1)+length(shuffled_master_table2),length(catch_table));
    shuffled_master_table3 = shuffled_master_table2(shuffle_w_catch,:);
    shuffled_master_table = shuffled_master_table3;
end

% add shuffled SCC table, which combines all indiv stim class scc rows
master_table2 = cat(1,master_table2,shuffled_master_table);

% Give the new master table to the output var
master_table_out = master_table2;

return