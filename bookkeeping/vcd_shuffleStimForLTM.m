function master_table_out = vcd_shuffleStimForLTM(master_table, nr_of_trials_per_block)
%% VCD function to combine and shuffle stimulus classes for SCC trials
%   master_table_out = vcd_shuffleStimForLTM(master_table, nr_of_trials_per_block)
%
% Purpose: LTM blocks ask subjects to judge the whether the test image 
% is the same as associated image pair learned the day before.  
% These pairs are mixed across stimulus classes, therefore we want to 
% (balanced) mixture of images from different classes for every LTM
% crossing. This function will do so, while making sure that we comply to
% even/uneven left/right trial order, as well as  showing each stimulus
% class at least once per block.
%
% INPUTS:
% * master_table                : (table) giant table with stimulus and trial definitions
% * nr_of_trials_per_block      : (int) nr of trials for each scc block
%
% OUTPUTS:
% * master_table_out            : (table) giant table, where scc blocks 
%                                  contain a mixture of gabor,rdk,dot, and
%                                  obj stimuli.
%
%
%% Get scc trials for left-cued, left-uncued, right-cued, right-uncued
ltm_sub_shuffle = []; 
ltm_idx = []; ltm_sub = [];

ltm_stim_loc = unique(master_table.stimloc(strcmp(master_table.task_class_name,{'ltm'})));
ltm_cue_stat = unique(master_table.iscued(strcmp(master_table.task_class_name,{'ltm'})));

% Get all the rows with SCC and the corresponsing stim classes
all_ltm_rows = strcmp(master_table.task_class_name,{'ltm'});
ltm_stim_rows = master_table.stim_class_name(all_ltm_rows);

count = 1;
for jj = 1:length(ltm_stim_loc) % 3 stim locations (l/r/c)
    for kk = 1:length(ltm_cue_stat) % cued and uncued (l/r/c)
        % get all indices for stim loc and cue status
        ltm_idx(count,:)    = (strcmp(master_table.task_class_name,{'ltm'}) & master_table.stimloc==ltm_stim_loc(jj) & master_table.iscued == ltm_cue_stat(kk));
        ltm_sub(count,:)    = find(ltm_idx(count,:));
        
        
        
        % shape into an 8 trials x m blocks matrix, such that the
        % columns contain different stimulus classes
        tmp_sub    = reshape(ltm_sub(count,:),[],nr_of_trials_per_block);
        
        % shuffle the order across columns such that gabors, rdks, dots, obj's get mixed
        ind        = nr_of_trials_per_block*[0:(size(tmp_sub,1)-1)]' ...
            + reshape(shuffle_concat([1:nr_of_trials_per_block],size(tmp_sub,1)),nr_of_trials_per_block,[])';
        
        % transpose as matlab indexes over cols first
        tmp_sub1   = tmp_sub';
        tmp_sub1   = tmp_sub1(ind');
        
        assert(isequal(sort(tmp_sub1,1),tmp_sub')) % sorted shuffled and unsort should be the same
        assert(~isequal(tmp_sub1,tmp_sub)) % shuffled and unsort should not be the same
        
        ltm_sub_shuffle(count,:) = tmp_sub1(:); % this will fail if we don't have equal nr of trials per stimloc/cueing status
        count = count + 1;
    end
end

% scc_cue_stat dim 1: is left uncued, left cued, right uncued, right cued.
% versions A and B represent left cued/right uncued and right cued/left
% uncued. We combine these conditions first, and then randomly assign the
% combined trials to the master table.

combined_scc_shuffleA = [ltm_sub_shuffle(1,:); ltm_sub_shuffle(4,:)]'; % left uncued, right uncued
combined_scc_shuffleB = [ltm_sub_shuffle(2,:); ltm_sub_shuffle(3,:)]'; % left cued, right uncued

combined_scc_shuffleAB0 = cat(3, combined_scc_shuffleA, combined_scc_shuffleB);
combined_scc_shuffleAB1 = permute(combined_scc_shuffleAB0,[3,1,2]); 
combined_scc_shuffleAB2 = reshape(combined_scc_shuffleAB1, [], size(combined_scc_shuffleAB1,3));

% shuffle order of uncued/cued
shuffle_vec = shuffle_concat(1:nr_of_trials_per_block,size(combined_scc_shuffleAB2,1)/nr_of_trials_per_block);
shuffle_vec = reshape(shuffle_vec, nr_of_trials_per_block,[]);
shuffle_vec = shuffle_vec + [[0:(size(shuffle_vec,2)-1)].*nr_of_trials_per_block];
shuffle_vec = shuffle_vec(:);
combined_scc_shuffleAB3 = combined_scc_shuffleAB2(shuffle_vec,:);

% APPLY THE SHUFFLE!
shuffled_SCC_table = table();
for ii = 1:size(combined_scc_shuffleAB3,1)
    new_trial = [master_table(combined_scc_shuffleAB3(ii,1),:); master_table(combined_scc_shuffleAB3(ii,2),:)];
    shuffled_SCC_table  = cat(1,shuffled_SCC_table, new_trial);
end
assert(isequal(shuffled_SCC_table.stimloc, repmat([1,2]',size(shuffled_SCC_table,1)/2,1)))

% for debugging: keep just n case
% old_scc_block_nr = shuffled_SCC_table.block_nr;
% old_scc_local_trial_nr = shuffled_SCC_table.block_local_trial_nr;
% old_scc_repeat_nr = shuffled_SCC_table.repeat_nr;
% old_scc_stimloc_nr = shuffled_SCC_table.stimloc;

% update block_nr, block_local_trial_nr and repeat nr according to
% new image/trial order
shuffled_SCC_table.stim_class_unique_block_nr = repelem(1:(size(shuffled_SCC_table,1)/(nr_of_trials_per_block*2)),(nr_of_trials_per_block*2))';
shuffled_SCC_table.block_local_trial_nr = repmat(repelem(1:nr_of_trials_per_block,2)',size(shuffled_SCC_table,1)/(nr_of_trials_per_block*2),1);

for ss = unique(shuffled_SCC_table.stim_class)'
    for ii = unique(shuffled_SCC_table.unique_im_nr(shuffled_SCC_table.stim_class==ss))'
        for cc = unique(shuffled_SCC_table.iscued)'
            selected_im = (shuffled_SCC_table.unique_im_nr==ii & shuffled_SCC_table.stim_class==ss & shuffled_SCC_table.iscued==cc);
            shuffled_SCC_table.repeat_nr(selected_im,:) = sort(shuffled_SCC_table.repeat_nr(selected_im,:));
        end
    end
end

% copy master table without current SCC rows
master_table2 = master_table(~all_scc_rows,:);

% add shuffled SCC table, which combines all indiv stim class scc rows
master_table2 = cat(1,master_table2,shuffled_SCC_table);

% Give the new master table to the output var
master_table_out = master_table2;

return