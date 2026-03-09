function condition_master = vcd_createLTMIMGTrainingTrials(condition_master, params)

% only for deep demo sessions
assert(params.is_demo && ~params.is_wide)

nr_sessions = unique(condition_master.session_nr);

for ses = 1:length(nr_sessions)
    
    % Get LTM/IMG trials
    practice_trials_LTM = condition_master.session_type==2 & condition_master.session_nr == nr_sessions(ses) & condition_master.task_class==6;
    practice_trials_IMG = condition_master.session_type==2 & condition_master.session_nr == nr_sessions(ses) & condition_master.task_class==7;
    
    if sum(practice_trials_LTM)>1
        
        condition_master0 = condition_master(practice_trials_LTM,:);
        
        % loop over trials to correct the answer/remove uncued side
        for tt = 1:size(condition_master0,1)
            
            cued_side = condition_master0.is_cued(tt);
            
            if cued_side == 1
                if condition_master0.correct_response(tt) == 1
                    % ensure stim2 is correct
                    corr_stim = params.stim.all_ltm_pairs(ismember(params.stim.all_ltm_pairs(:,1),condition_master0.stim_nr_left(tt)),2);
                    assert(isequal(corr_stim,condition_master0.stim2_im_nr(tt,cued_side)));
                    assert(isequal(0,condition_master0.is_lure(tt,cued_side)));
                    assert(isequal(1,condition_master0.stim2_delta(tt,cued_side)));
                else
                    % make stim2 correct
                    corr_stim = params.stim.all_ltm_pairs(ismember(params.stim.all_ltm_pairs(:,1),condition_master0.stim_nr_left(tt)),2);
                    assert(~isequal(corr_stim,condition_master0.stim2_im_nr(tt,cued_side)));
                    
                    condition_master0.stim2_im_nr(tt,cued_side) = corr_stim;
                    condition_master0.correct_response(tt) = 1;
                    condition_master0.is_lure(tt,cued_side) = 0;
                    condition_master0.stim2_delta(tt,cued_side) = 1;
                    tmp = vcd('fullinfo',corr_stim);
                    condition_master0.stim2_orient_dir(tt,cued_side) = tmp.orient_dir;
                end
                
                % remove uncued right side from table
                condition_master0.stim_nr_right(tt) = NaN;
                condition_master0.stim_class_name(tt,2) = {NaN};
                condition_master0.condition_nr(tt,2) = NaN;
                condition_master0.condition_name(tt,2) = {NaN};
                condition_master0.orient_dir(tt,2) = NaN;
                condition_master0.contrast(tt,2) = NaN;
                condition_master0.gbr_phase(tt,2) = NaN;
                condition_master0.is_special_core(tt,2) = NaN;
                condition_master0.stim_class_name(tt,2) = {NaN};
                condition_master0.orient_dir(tt,2) = NaN;
                condition_master0.rdk_coherence(tt,2) = NaN;
                condition_master0.super_cat(tt,2) = NaN;
                condition_master0.basic_cat(tt,2) = NaN;
                condition_master0.sub_cat(tt,2) = NaN;
                condition_master0.affordance_cat(tt,2) = NaN;
                condition_master0.super_cat_name(tt,2) = {NaN};
                condition_master0.basic_cat_name(tt,2) = {NaN};
                condition_master0.sub_cat_name(tt,2) = {NaN};
                condition_master0.affordance_name(tt,2) = {NaN};
                condition_master0.stim2_im_nr(tt,2) = NaN;
                condition_master0.stim2_delta(tt,2) = NaN;
                condition_master0.stim2_orient_dir(tt,2) = NaN;
                condition_master0.is_lure(tt,2)     = NaN;
                
            elseif cued_side == 2
                
                if condition_master0.correct_response(tt) == 1
                    % ensure stim2 is correct
                    corr_stim = params.stim.all_ltm_pairs(ismember(params.stim.all_ltm_pairs(:,1),condition_master0.stim_nr_right(tt)),2);
                    assert(isequal(corr_stim,condition_master0.stim2_im_nr(tt,cued_side)));
                    assert(isequal(0,condition_master0.is_lure(tt,cued_side)));
                    assert(isequal(1,condition_master0.stim2_delta(tt,cued_side)));
                else
                    % make stim2 correct
                    corr_stim = params.stim.all_ltm_pairs(ismember(params.stim.all_ltm_pairs(:,1),condition_master0.stim_nr_right(tt)),2);
                    assert(~isequal(corr_stim,condition_master0.stim2_im_nr(tt,cued_side)));
                    
                    condition_master0.stim2_im_nr(tt,cued_side) = corr_stim;
                    condition_master0.correct_response(tt) = 1;
                    condition_master0.is_lure(tt,cued_side) = 0;
                    condition_master0.stim2_delta(tt,cued_side) = 1;
                    tmp = vcd('fullinfo',corr_stim);
                    condition_master0.stim2_orient_dir(tt,cued_side) = tmp.orient_dir;
                end
                
                % remove uncued left side from table
                condition_master0.stim_nr_left(tt) = NaN;
                condition_master0.stim_class_name(tt,1) = {NaN};
                condition_master0.condition_nr(tt,1) = NaN;
                condition_master0.condition_name(tt,1) = {NaN};
                condition_master0.orient_dir(tt,1) = NaN;
                condition_master0.contrast(tt,1) = NaN;
                condition_master0.gbr_phase(tt,1) = NaN;
                condition_master0.is_special_core(tt,1) = NaN;
                condition_master0.stim_class_name(tt,1) = {NaN};
                condition_master0.orient_dir(tt,1) = NaN;
                condition_master0.rdk_coherence(tt,1) = NaN;
                condition_master0.super_cat(tt,1) = NaN;
                condition_master0.basic_cat(tt,1) = NaN;
                condition_master0.sub_cat(tt,1) = NaN;
                condition_master0.affordance_cat(tt,1) = NaN;
                condition_master0.super_cat_name(tt,1) = {NaN};
                condition_master0.basic_cat_name(tt,1) = {NaN};
                condition_master0.sub_cat_name(tt,1) = {NaN};
                condition_master0.affordance_name(tt,1) = {NaN};
                condition_master0.stim2_im_nr(tt,1) = NaN;
                condition_master0.stim2_delta(tt,1) = NaN;
                condition_master0.stim2_orient_dir(tt,1) = NaN;
                condition_master0.is_lure(tt,1)     = NaN;
            
            
            elseif cued_side == 3
                if condition_master0.correct_response(tt) == 1
                    % ensure stim2 is correct
                    corr_stim = params.stim.all_ltm_pairs(ismember(params.stim.all_ltm_pairs(:,1),condition_master0.stim_nr_left(tt)),2);
                    assert(isequal(corr_stim,condition_master0.stim2_im_nr(tt,1)));
                    assert(isequal(0,condition_master0.is_lure(tt,1)));
                    assert(isequal(1,condition_master0.stim2_delta(tt,1)));
                else
                    % make stim2 correct
                    corr_stim = params.stim.all_ltm_pairs(ismember(condition_master0.stim_nr_left(tt),params.stim.all_ltm_pairs(:,1)),2);
                    assert(~isequal(corr_stim,condition_master0.stim2_im_nr(tt,1)));
                    
                    condition_master0.stim2_im_nr(tt,1) = corr_stim;
                    condition_master0.correct_response(tt)  = 1;
                    condition_master0.stim2_delta(tt,1)     = 1;
                    condition_master0.is_lure(tt,1)         = 0;
                end
                
                % remove uncued right side from table
                condition_master0.stim2_im_nr(tt,2) = NaN;
                condition_master0.stim2_delta(tt,2) = NaN;
                condition_master0.is_lure(tt,2)     = NaN;
                condition_master0.stim2_orient_dir(tt,2) = NaN;
            else
                % do nothing
            end
        end
        
        condition_master(practice_trials_LTM,:) = condition_master0;
        clear condition_master0
    end
    
    
    if sum(practice_trials_IMG)>1
        % Get IMG trials
        condition_master0 = condition_master(practice_trials_IMG,:);
        
        % loop over trials to insert special core stim into stim2 slot
        for tt = 1:size(condition_master0,1)

            % Insert special core image in stim2_im_nr
            cued_side = condition_master0.is_cued(tt);
            if ismember(cued_side,[1,3])
                condition_master0.stim2_im_nr(tt,1)  = condition_master0.stim_nr_left(tt);
                condition_master0.stim2_orient_dir(tt,1) = condition_master0.orient_dir(tt,1);
            elseif cued_side == 2
                condition_master0.stim2_im_nr(tt,2) = condition_master0.stim_nr_right(tt);
                condition_master0.stim2_orient_dir(tt,1) = condition_master0.orient_dir(tt,2);
            end
        end
        
        condition_master(practice_trials_IMG,:) = condition_master0;
        clear condition_master0
    end

end


return