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
    
% Get all the rows with SCC and the corresponsing stim classes
all_task_rows = strcmp(master_table.task_class_name,{task_class_name_to_shuffle});

% get all indices for left and right stimuli, separate for
% cued/uncued/ncued conditions
for kk = 1:length(cued_stim_loc) % left, right, neutral

    % get task trials (filter out catch trials)
    m0_idx          = ~master_table.is_catch & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk);
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
        if sum(master_table.is_catch) > 1
            m0_center_catch  = master_table.stim_nr_left(master_table.is_catch & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
            sub_center_catch = find(master_table.is_catch & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
        end
        
    else % gabor/rdk/dot/obj

        m0_left_img_nr  = master_table.stim_nr_left(m0_idx);
        m0_right_img_nr = master_table.stim_nr_right(m0_idx);
        
%         m0_left_catch   = master_table.stim_nr_left(catch_sub);
%         m0_right_catch  = master_table.stim_nr_right(catch_sub);
%         n_catch_l       = length(m0_left_catch);
%         n_catch_r       = length(m0_right_catch);
%         
        % shape into an 8 trials x m blocks matrix, such that the
        % columns contain different stimulus classes
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
        if sum(master_table.is_catch) > 1
            img_nr_left_catch(kk,:)  = master_table.stim_nr_left(master_table.is_catch & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
            img_nr_right_catch(kk,:) = master_table.stim_nr_right(master_table.is_catch & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
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

combined_trial_shuffleA = [sub_shuffle_left(1,:); sub_shuffle_right(1,:)]'; % shuffled left cued, shuffled right uncued (trials x loc (l/r))
combined_trial_shuffleB = [sub_shuffle_left(2,:); sub_shuffle_right(2,:)]'; % shuffled left uncued, shuffled right cued (trials x loc (l/r))

combined_trial_shuffleAB0 = cat(3, combined_trial_shuffleA, combined_trial_shuffleB); % 320 trials x 2 loc (l/r) x 2 cuedir (first left, then right)
combined_trial_shuffleAB1 = permute(combined_trial_shuffleAB0,[3,1,2]);  % 320 trials x 2 loc (l/r) x 2 cuedir --> 2 cuedir x 320 trials x 2 loc (l/r)

% We want the number of cued left and cued right to be balanced within a
% block. Right now: combined_scc_shuffleAB1(1,1,1) -> left cued, left stim
% combined_scc_shuffleAB1(2,1,1) -> right cued, left stim
% combined_scc_shuffleAB1(1,1,2) -> left cued, right stim
% combined_scc_shuffleAB1(2,1,2) -> right cued, right stim. 
% We want left/right cues to happen randomly across trials (rows), but keep
% indices for left and right stim separately (columns).
% To achieve this, we will first reshape combined_scc_shuffleAB1 along the
% first two dimensions, we will get a left cued trial in the first row, and
% a right cued trial in the second row.
combined_trial_shuffleAB2 = reshape(combined_trial_shuffleAB1, [], size(combined_trial_shuffleAB1,3)); % 2 cuedir x 320 trials x 2 loc (l/r) --> 640 trials x 2 loc (l/r)

% for LTM: check that column lengths is the same as sum of repeats * special core
% stimuli
if strcmp(session_type, 'MRI')
    unique_trial_repeats = params.exp.n_unique_trial_repeats_mri;
else
    unique_trial_repeats = params.exp.n_unique_trial_repeats_behavior;
end
    
if strcmp(task_class_name_to_shuffle,{'ltm'})
    assert(isequal(size(combined_trial_shuffleAB2,1), ...
    unique(unique_trial_repeats(1:4,strcmp(task_class_name_to_shuffle,params.exp.taskclassnames)))  * ... 23 repeats
    (2*(length(params.stim.all_specialcore_im_nrs)-length(params.stim.ns.unique_im_nrs_specialcore))))); % 8*4 (=32) * 2 (double nr of repeats)
else
    assert(isequal(size(combined_trial_shuffleAB2,1), ...
    sum(unique_trial_repeats(1:4,strcmp(task_class_name_to_shuffle,params.exp.taskclassnames))'  .* ... nr of repeats
    [length(params.stim.gabor.unique_im_nrs_core),length(params.stim.rdk.unique_im_nrs_core),... nr of unique core images
    length(params.stim.dot.unique_im_nrs_core),length(params.stim.obj.unique_im_nrs_core) ]))); % 
end
% shuffle order of uncued/cued
shuffle_vec = shuffle_concat(1:nr_of_trials_per_block,size(combined_trial_shuffleAB2,1)/nr_of_trials_per_block); % 640 trials (1-8 indices)
shuffle_vec = reshape(shuffle_vec, nr_of_trials_per_block,[]); % 8 trials x 80 blocks
shuffle_vec = shuffle_vec + [[0:(size(shuffle_vec,2)-1)].*nr_of_trials_per_block]; % 8 trials x 80 blocks (continuous counting)
shuffle_vec = shuffle_vec(:);
combined_trial_shuffleAB3 = combined_trial_shuffleAB2(shuffle_vec,:);

if exist('sub_shuffle_center','var')
    assert(isequal(size(sub_shuffle_center,1), (length(params.stim.ns.unique_im_nrs_specialcore)*unique_trial_repeats(5,strcmp(task_class_name_to_shuffle,params.exp.taskclassnames)))))
    % add column of nans to match two column structure for left/right
    % matrix
    sub_shuffle_center2 = cat(2,sub_shuffle_center,NaN(size(sub_shuffle_center)));
    
    if strcmp(task_class_name_to_shuffle,{'ltm'})
        combined_trial_shuffleABC = NaN(size(combined_trial_shuffleAB3,1)+size(sub_shuffle_center2,1),2);        
        
        insert_scenes_here = datasample(1:size(combined_trial_shuffleABC,1),size(sub_shuffle_center2,1), 'Replace',false);
        insert_classics_here = setdiff(1:size(combined_trial_shuffleABC,1),insert_scenes_here);
        combined_trial_shuffleABC(insert_scenes_here,:) = sub_shuffle_center2;
        combined_trial_shuffleABC(insert_classics_here,:) = combined_trial_shuffleAB3;
    end
    
    if sum(master_table.is_catch) > 1
        shuffle_vec_catch = shuffle_concat(sub_center_catch(1,:),1);
        catch_table = catch_table.stim_class_name(shuffle_vec_catch(:,1),1);
    end
else
    combined_trial_shuffleABC = combined_trial_shuffleAB3;
    
    if sum(master_table.is_catch) > 1
        shuffle_vec_catch = [shuffle_concat(sub_shuffle_left_catch(1,:),1); shuffle_concat(sub_shuffle_right_catch(2,:),1)];
        % Shuffle catch trials
        catch_table = master_table(shuffle_vec_catch(:,1),1);
        catch_table = master_table(shuffle_vec_catch(:,2),2);
    end
end
  

% APPLY THE SHUFFLE!
shuffled_master_table = table();
for xx = 1:size(master_table(1,:),2); col_widths(xx) = size(table2array(master_table(1,xx)),2); end
colNames = master_table.Properties.VariableNames;

double_width_cols = find(col_widths==2);

for ii = 1:size(combined_trial_shuffleABC,1)
    
    new_trial = master_table(combined_trial_shuffleABC(ii,1),:);
    
    new_trial.unique_trial_nr = ii;
    new_trial.stim_class = 99;
    new_trial.repeat_nr  = 99;
        
    % update block_nr,  and repeat nr according to
    % new image/trial order
    new_trial.stim_class_unique_block_nr = ceil(ii/nr_of_trials_per_block);
    new_trial.trial_nr = mod(ii-1,nr_of_trials_per_block)+1;

    new_trial.stim_nr_left = master_table.stim_nr_left(combined_trial_shuffleABC(ii,1));
    if ~isnan(combined_trial_shuffleABC(ii,2))
        new_trial.stim_nr_right = master_table.stim_nr_right(combined_trial_shuffleABC(ii,2));
        assert(isequal(master_table.is_cued(combined_trial_shuffleABC(ii,1)),master_table.is_cued(combined_trial_shuffleABC(ii,2))))
        assert(isequal(master_table.is_catch(combined_trial_shuffleABC(ii,1)),master_table.is_catch(combined_trial_shuffleABC(ii,2))))
    else
        new_trial.stim_nr_right = NaN;
    end

    for jj = 1:length(double_width_cols)
        colName = colNames{double_width_cols(jj)};
        if isnan(combined_trial_shuffleABC(ii,2))
            new_trial.(colName) = [master_table.(colName)(combined_trial_shuffleABC(ii,1),1), master_table.(colName)(combined_trial_shuffleABC(ii,1),2)];
        else
            new_trial.(colName) = [master_table.(colName)(combined_trial_shuffleABC(ii,1),1), master_table.(colName)(combined_trial_shuffleABC(ii,2),2)];
        end
    end
    shuffled_master_table  = cat(1,shuffled_master_table, new_trial);
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

if sum(master_table.is_catch) > 1
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