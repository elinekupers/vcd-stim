function [] = vcd_pairUniqueImagesLTM(p, unique_im)




if strcmp(tasks(task_crossings(curr_task)).name,'ltm')
    paired_stim = [];
    delta_vec = shuffle_concat(p.stim.gabor.ltm_pairs, (n_unique_cases/(n_deltas/2)));
    pair_vec = conds_master_single_rep3(:,7);
    
else
    
    
    
    % LTM stim-stim pairing will be probabilistic, with these tweaks:
    exp_session.session.ltm.prob_new_pairing  = 0.2;    % chance that LTM stim A will be match to stim C (instead of stim B), in a given session
    exp_session.session.ltm.prob_pair_order_flip = 0.2; % chance that LTM stim A -> B will flip to B -> A in a given session