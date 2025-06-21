function condition_master = vcd_updateGlobalCounters(params, condition_master, env_type)

% add global session, run, block, and trial nr
global_run_nrB     = 0;
global_block_nrB   = 0;
global_trial_nrB   = 0;

for ii = 1:size(condition_master,1)
    if ii == 1 
        global_run_nr     = 1;
        global_block_nr   = 1;
        global_trial_nr   = 1;
    
        global_run_counter   = global_run_nr;
        global_block_counter = global_block_nr;
        global_trial_counter = global_trial_nr;
        
        % MRI session type 1 uses regular counter. BEHAVIOR has no A/B versions, always use regular counter
    elseif (strcmp(env_type,'MRI') && condition_master.session_type(ii) == 1) || strcmp(env_type,'BEHAVIOR') 
        if condition_master.run_nr(ii-1) ~= condition_master.run_nr(ii)
            global_run_nr = global_run_nr + 1;
        end
        if condition_master.block_nr(ii-1) ~= condition_master.block_nr(ii)
            global_block_nr = global_block_nr + 1;
        end
        if condition_master.trial_nr(ii-1) ~= condition_master.trial_nr(ii)
            global_trial_nr = global_trial_nr + 1;
        end
        
        % Continue counting for counterB when it is not session 1 and 27
        if condition_master.session_nr(ii) > 1 & condition_master.session_type(ii) == 1 
            global_run_nrB     = global_run_nr;
            global_block_nrB   = global_block_nr;
            global_trial_nrB   = global_trial_nr;
        end
        
        global_run_counter   = global_run_nr;
        global_block_counter = global_block_nr;
        global_trial_counter = global_trial_nr;
        
        % session type 2 uses counterB for session 1 and 27
    elseif strcmp(env_type,'MRI') && condition_master.session_type(ii) == 2 && ismember(condition_master.session_nr(ii),[1,length(params.exp.session.n_mri_sessions)])
        if condition_master.run_nr(ii-1) ~= condition_master.run_nr(ii)
            global_run_nrB = global_run_nrB + 1;
        end
        if condition_master.block_nr(ii-1) ~= condition_master.block_nr(ii)
            global_block_nrB = global_block_nrB + 1;
        end
        if condition_master.trial_nr(ii-1) ~= condition_master.trial_nr(ii)
            global_trial_nrB = global_trial_nrB + 1;
        end
        global_run_counter = global_run_nrB;
        global_block_counter = global_block_nrB;
        global_trial_counter = global_trial_nrB;
        
    end
    
    condition_master.global_run_nr(ii)   = global_run_counter;
    condition_master.global_block_nr(ii) = global_block_counter;
    condition_master.global_trial_nr(ii) = global_trial_counter;
end  