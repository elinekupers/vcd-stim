function master_table_out = vcd_shuffleStimForTaskClass(params, task_class_name_to_shuffle, master_table, nr_of_trials_per_block)
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
%     * params.exp.stimClassLabels to get a list of stimulus classes.
% * task_class_name_to_shuffle  : (str) name of task class we want to
%                                   shuffle the stimuli (either 'scc' or
%                                   'ltm')
% * master_table                : (table) giant table with stimulus and trial definitions
% * nr_of_trials_per_block      : (int) nr of trials for each scc or ltm block.
%                                   Expecting either 8 or 4 trials/block.
%
% OUTPUTS:
% * master_table_out            : (table) giant table, where scc blocks 
%                                  contain a mixture of gabor, rdk, dot,and
%                                  obj stimuli.
%
%% Check inputs
% we only want to shuffle scc and ltm task crossings
assert(ismember(strcmp(task_class_name_to_shuffle,{'scc','ltm'})));

% scc task crossings have 8 trials/block, ltm task crossings have 4 trials/block
assert(isequal(nr_of_trials_per_block, choose(strcmp(task_class_name_to_shuffle,'scc'), params.exp.block.n_trials_single_epoch,params.exp.block.n_trials_double_epoch)));


%% Get scc trials for left-cued, left-uncued, right-cued, right-uncued
if isequalwithequalnans(master_table.stim_nr_right(strcmp(master_table.task_class_name,{task_class_name_to_shuffle})), ...
        length(master_table.stim_nr_right(strcmp(master_table.task_class_name,{task_class_name_to_shuffle}))))
    scc_stim_loc = 3;
else
    scc_stim_loc = [1;2];
end
scc_cue_stat = unique(master_table.is_cued(strcmp(master_table.task_class_name,{task_class_name_to_shuffle})));

% Get all the rows with SCC and the corresponsing stim classes
all_task_rows = strcmp(master_table.task_class_name,{task_class_name_to_shuffle});

% get all indices for left and right stimuli, separate for cued/uncued
% conditions
for kk = 1:length(scc_cue_stat) % left and right

    % get scc trials (filter out catch trials)
    m0_idx          = ~master_table.is_catch & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == scc_cue_stat(kk);
    m0_sub          = find(m0_idx); % keep track of row nr
    m0_left_img_nr  = master_table.stim_nr_left(m0_idx);
    m0_right_img_nr = master_table.stim_nr_right(m0_idx);
    m0_img_class    = master_table.stim_class(m0_idx);


    % shape into an 8 trials x m blocks matrix, such that the
    % columns contain different stimulus classes
    tmp_img_nr_left_rz    = reshape(m0_left_img_nr,[],nr_of_trials_per_block);
    tmp_img_nr_right_rz   = reshape(m0_right_img_nr,[],nr_of_trials_per_block);
    tmp_sub_rz            = reshape(m0_sub,[],nr_of_trials_per_block);

    % shuffle the order across columns such that gabors, rdks, dots, obj's get mixed
    shuffled_ind_left  = nr_of_trials_per_block*[0:(size(tmp_img_nr_left_rz,1)-1)]' ...
                        + reshape(shuffle_concat([1:nr_of_trials_per_block],size(tmp_img_nr_left_rz,1)),nr_of_trials_per_block,[])';
    shuffled_ind_right = nr_of_trials_per_block*[0:(size(tmp_img_nr_right_rz,1)-1)]' ...
                        + reshape(shuffle_concat([1:nr_of_trials_per_block],size(tmp_img_nr_right_rz,1)),nr_of_trials_per_block,[])';
    
    % transpose as matlab indexes over cols first
    tmp_img_nr_left1       = tmp_img_nr_left_rz';
    tmp_img_nr_right1      = tmp_img_nr_right_rz';
    shuffled_ind_left1     = shuffled_ind_left';
    shuffled_ind_right1    = shuffled_ind_right';
    tmp_sub_rz1 = tmp_sub_rz';
    
    
    % shuffle trials based on shuffled indices
    tmp_img_nr_left2   = tmp_img_nr_left1(shuffled_ind_left1);
    tmp_img_nr_right2  = tmp_img_nr_right1(shuffled_ind_right1);
    tmp_sub_rz2_r      = tmp_sub_rz1(shuffled_ind_left1);
    tmp_sub_rz2_l      = tmp_sub_rz1(shuffled_ind_right1);

    assert(isequal(sort(tmp_img_nr_left2,1),sort(tmp_img_nr_left1,1))) % sorted shuffled and sorted unshuffled left should be the same
    assert(isequal(sort(tmp_img_nr_right2,1),sort(tmp_img_nr_right1,1))) % sorted shuffled and sorted unshuffled right should be the same
    assert(~isequal(tmp_img_nr_left2,tmp_img_nr_left1)) % shuffled and unsort should not be the same
    assert(~isequal(tmp_img_nr_right2,tmp_img_nr_right1)) % shuffled and unsort should not be the same

    img_nr_shuffle_left(kk,:)  = tmp_img_nr_left2(:); % this will fail if we don't have equal nr of trials per stimloc/cueing status
    img_nr_shuffle_right(kk,:) = tmp_img_nr_right2(:); % this will fail if we don't have equal nr of trials per stimloc/cueing status
    
    sub_shuffle_left(kk,:)  = tmp_sub_rz2_r(:); 
    sub_shuffle_right(kk,:) = tmp_sub_rz2_l(:); 
end


% GOAL: 
% Versions A and B represent left cued/right uncued and right cued/left
% uncued. We combine these conditions first, then ensure left and right cued
% trials are mixed within a block, and then randomly assign the
% combined trials to the master table.

% Note: sub_shuffle_left dim 1 order: left cued, right cued
%       sub_shuffle_right dim 1 order: left cued, right cued
combined_scc_shuffleA = [sub_shuffle_left(1,:); sub_shuffle_right(1,:)]'; % shuffled left cued, shuffled right uncued (trials x loc (l/r))
combined_scc_shuffleB = [sub_shuffle_left(2,:); sub_shuffle_right(2,:)]'; % shuffled left uncued, shuffled right cued (trials x loc (l/r))

combined_scc_shuffleAB0 = cat(3, combined_scc_shuffleA, combined_scc_shuffleB); % 320 trials x 2 loc (l/r) x 2 cuedir (first left, then right)
combined_scc_shuffleAB1 = permute(combined_scc_shuffleAB0,[3,1,2]);  % 320 trials x 2 loc (l/r) x 2 cuedir --> 2 cuedir x 320 trials x 2 loc (l/r)

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
combined_scc_shuffleAB2 = reshape(combined_scc_shuffleAB1, [], size(combined_scc_shuffleAB1,3)); % 2 cuedir x 320 trials x 2 loc (l/r) --> 640 trials x 2 loc (l/r)

% shuffle order of uncued/cued (but keep 
shuffle_vec = shuffle_concat(1:nr_of_trials_per_block,size(combined_scc_shuffleAB2,1)/nr_of_trials_per_block); % 640 trials (1-8 indices)
shuffle_vec = reshape(shuffle_vec, nr_of_trials_per_block,[]); % 8 trials x 80 blocks
shuffle_vec = shuffle_vec + [[0:(size(shuffle_vec,2)-1)].*nr_of_trials_per_block]; % 8 trials x 80 blocks (continuous counting)
shuffle_vec = shuffle_vec(:);
combined_scc_shuffleAB3 = combined_scc_shuffleAB2(shuffle_vec,:);

% APPLY THE SHUFFLE!
shuffled_master_table = table();
for xx = 1:size(master_table(1,:),2); col_widths(xx) = size(table2array(master_table(1,xx)),2); end
colNames = master_table.Properties.VariableNames;

double_width_cols = find(col_widths==2);

for ii = 1:size(combined_scc_shuffleAB3,1)
    
    new_trial = master_table(combined_scc_shuffleAB3(ii,1),:);
    
    new_trial.unique_trial_nr = ii;
    new_trial.stim_class = 99;
    new_trial.repeat_nr  = 99;
        
    % update block_nr, block_local_trial_nr and repeat nr according to
    % new image/trial order
    new_trial.stim_class_unique_block_nr = ceil(ii/nr_of_trials_per_block);
    new_trial.block_local_trial_nr = mod(ii-1,nr_of_trials_per_block)+1;

    new_trial.stim_nr_left = master_table.stim_nr_left(combined_scc_shuffleAB3(ii,1));
    new_trial.stim_nr_right = master_table.stim_nr_right(combined_scc_shuffleAB3(ii,2));
    
    assert(isequal(master_table.is_cued(combined_scc_shuffleAB3(ii,1)),master_table.is_cued(combined_scc_shuffleAB3(ii,2))))
    assert(isequal(master_table.is_catch(combined_scc_shuffleAB3(ii,1)),master_table.is_catch(combined_scc_shuffleAB3(ii,2))))

    
    for jj = 1:length(double_width_cols)
        colName = colNames{double_width_cols(jj)};
        
        new_trial.(colName) = [master_table.(colName)(combined_scc_shuffleAB3(ii,1),1), master_table.(colName)(combined_scc_shuffleAB3(ii,2),2)];
    end
    shuffled_master_table  = cat(1,shuffled_master_table, new_trial);
end

for sc = 1:length(params.exp.stimClassLabels)
    stimClass = params.exp.stimClassLabels{sc};
    assert(all(ismember(shuffled_master_table.stim_nr_left(strcmp(shuffled_master_table.stim_class_name(:,1),stimClass)),params.stim.(stimClass).unique_im_nrs)))
    assert(all(ismember(shuffled_master_table.stim_nr_right(strcmp(shuffled_master_table.stim_class_name(:,2),stimClass)),params.stim.(stimClass).unique_im_nrs)))
end

% Compute how often each new combination of stimuli x cue status is repeated
[uniquerow, ~, rowidx] = unique([shuffled_master_table.stim_nr_left, shuffled_master_table.stim_nr_right, shuffled_master_table.is_cued], 'rows'); 
noccurrences = accumarray(rowidx, 1);

repeat_nr = zeros(1,length(noccurrences));
for nrow = 1:size(shuffled_master_table,1)
    
    curr_comb = [shuffled_master_table.stim_nr_left(nrow), shuffled_master_table.stim_nr_right(nrow), shuffled_master_table.is_cued(nrow)];
    
    [~,comb_nr] = intersect(uniquerow,curr_comb,'rows');
    repeat_nr(comb_nr) = repeat_nr(comb_nr)+1;
    shuffled_master_table.repeat_nr(nrow) =  repeat_nr(comb_nr);
end

% copy master table without current task rows
master_table2 = master_table(~all_task_rows,:);

% add shuffled SCC table, which combines all indiv stim class scc rows
master_table2 = cat(1,master_table2,shuffled_master_table);

% Give the new master table to the output var
master_table_out = master_table2;

return