function condition_master = vcd_updateGlobalCounters(params, condition_master, env_type)

% If we don't have is_wide field in parameters, assume it is not a wide MRI
% session.
if ~isfield(params, 'is_wide') || isempty(params.is_wide)
    params.is_wide = false;
end

if strcmp(env_type,'MRI')
    % Find which sessions have A and B versions:
    sessionB    = unique(condition_master.session_nr(condition_master.session_type == 2));
    nr_sessions = unique(condition_master.session_nr);
    assert(length(sessionB)==1); 
    
    % check if they correspond to the session we expect to have a B version)
    if params.is_wide
        assert(isequal(nr_sessions, params.exp.session.n_wide_sessions)) % since we separate the generation of wide and deep condition masters, we only expect one session
        assert(all(~isnan(params.exp.session.mri.wide.session_types)));
    else
        assert(isequal(max(nr_sessions), find(~isnan(params.exp.session.mri.deep.session_types(:,2)))));
        assert(all(~isnan(params.exp.session.mri.wide.session_types(end,:))));
    end
end

% Reset counters
global_run_nr      = 1;
global_block_nr    = 1;
global_trial_nr    = 1;

global_run_nrB     = 1;
global_block_nrB   = 1;
global_trial_nrB   = 1;

for st = unique(condition_master.session_type)'
    
    % if we have more than 1 session type, then update global counters separately.
    session_idx = find(condition_master.session_type==st);
    
    % add global session, run, block, and trial nr
    for ii = 1:length(session_idx)
        
        % MRI session type 1 uses regular counter. BEHAVIOR has no A/B versions, always use regular counter
        if strcmp(env_type,'BEHAVIOR') || (strcmp(env_type,'MRI') && st == 1)
            
            % We are dealing with all the deep MRI sessions, keep a record
            % of the penultimate session's last global numbers: we want to
            % use those as a starting point when we update global counters
            % for the B version of this session.
            if strcmp(env_type,'MRI') && params.is_wide ~=1 && max(condition_master.session_nr) > 1 && ...
                    isequal(condition_master.session_nr(session_idx(ii)), max(condition_master.session_nr))
                global_run_nrB   = global_run_nr;
                global_block_nrB = global_block_nr;
                global_trial_nrB = global_trial_nr;
            end
            
            if ii > 1
                if condition_master.run_nr(session_idx(ii-1)) ~= condition_master.run_nr(session_idx(ii))
                    global_run_nr   = global_run_nr + 1;
                end
                if condition_master.block_nr(session_idx(ii-1)) ~= condition_master.block_nr(session_idx(ii))
                    global_block_nr = global_block_nr + 1;
                end
                if condition_master.trial_nr(session_idx(ii-1)) ~= condition_master.trial_nr(session_idx(ii))
                    global_trial_nr = global_trial_nr + 1;
                end
            end
            
            % update the condition master columns
            condition_master.global_run_nr(session_idx(ii))   = global_run_nr;
            condition_master.global_block_nr(session_idx(ii)) = global_block_nr;
            condition_master.global_trial_nr(session_idx(ii)) = global_trial_nr;
            
        elseif strcmp(env_type,'MRI') && st == 2
            
            % if this is the B version of the wide session or last deep session
            if ismember(condition_master.session_nr(session_idx(ii)),sessionB)
                assert(condition_master.session_type(session_idx(ii)) == 2)
                if params.is_wide
                    assert(isequal(sessionB,1));
                else
                    assert(isequal(sessionB,max(params.exp.session.n_deep_sessions)));
                end
                
                if ii > 1
                    if condition_master.run_nr(session_idx(ii-1)) ~= condition_master.run_nr(session_idx(ii))
                        global_run_nrB   = global_run_nrB + 1;
                    end
                    if condition_master.block_nr(session_idx(ii-1)) ~= condition_master.block_nr(session_idx(ii))
                        global_block_nrB = global_block_nrB + 1;
                    end
                    if condition_master.trial_nr(session_idx(ii-1)) ~= condition_master.trial_nr(session_idx(ii))
                        global_trial_nrB = global_trial_nrB + 1;
                    end
                end
                % update the condition master columns
                condition_master.global_run_nr(session_idx(ii))   = global_run_nrB;
                condition_master.global_block_nr(session_idx(ii)) = global_block_nrB;
                condition_master.global_trial_nr(session_idx(ii)) = global_trial_nrB;
                
            end
        end
    end
end