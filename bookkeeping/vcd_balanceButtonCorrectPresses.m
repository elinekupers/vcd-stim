function condition_master = vcd_balanceButtonCorrectPresses(params, condition_master, env_type)
% VCD bookkeeping function to check if button presses are balanced to the
% extent possible, and shuffle trial order if needed.
%
%   condition_master = vcd_balanceButtonCorrectPresses(condition_master)
%
%
% check blocks per session type
[all_sessions,session_types] = vcd_getSessionEnvironmentParams(params, env_type);

if any(ismember(condition_master.Properties.VariableNames,'session_nr'))
    session_nrs = unique(condition_master.session_nr);
    session_flag = true;
else
    % create temporary session/session type/block nr vector
    session_flag = false;
    session_nrs   = 1:size(all_sessions,3);
    t = [condition_master.stim_class, condition_master.task_class];
    allocated_trials = zeros(size(condition_master,1),1);
end


for ses = session_nrs
    
    if session_flag
        ses_idx       = condition_master.session_nr==ses;
        session_types = unique(condition_master.session_type(ses_idx));
    end
    
    for st = session_types
        if ~isnan(session_types(ses,st))
            
            if session_flag
                % index rows if we have this information in the condition
                % master
                st_idx = (condition_master.session_nr==ses & condition_master.session_type==st);
            else
                stimtask_counter = all_sessions(:,:,ses,st);
                
                for sc = 1:size(stimtask_counter,1)
                    for tc = 1:size(stimtask_counter,2)
                        st_idx    = false(size(condition_master,1),1);
                        nr_blocks = stimtask_counter(sc,tc);
                        nr_trials = params.exp.nr_trials_per_block(sc,tc)*nr_blocks;
                        
                        if tc == 3
                            [~,ti] = intersect(t,[99,tc]);
                        else
                            [~,ti] = intersect(t,[sc,tc]);
                        end
                        
                        potential_trials = find(allocated_trials(ti)==0);
                        allocated_trials(ti(potential_trials(1:nr_trials))) = 1;
                        
                        st_idx(ti(potential_trials(1:nr_trials))) = true;
                    end
                end
            end
            
            curr_sc = unique(condition_master.stim_class(st_idx));
            curr_tc = unique(condition_master.task_class(st_idx));
            curr_tc = curr_tc(curr_tc~=1); % exclude fixation
            
            for jj = 1:length(curr_sc)
                for mm = 1:length(curr_tc)
                    
                    % check if we have such a crossing
                    if sum(ismember(condition_master.task_class(st_idx), curr_tc(mm)) & ismember(condition_master.stim_class(st_idx), curr_sc(jj)))>0
                        
                        % if so, then set the number of expected button presses
                        % **** !!! careful this is hardcoded stuff !!! *****
                        if ismember(curr_tc(mm),[2,4,5,6,7]) % cd, pc, wm, ltm, img have 2 response options
                            nr_responses = 2;
                        elseif ismember(curr_tc(mm),[3,8,10]) % scc, what, how have 4 response options
                            nr_responses = 4;
                        elseif ismember(curr_tc(mm),[9]) % where has 3 response options
                            nr_responses = 3;
                        end
                        
                        resp = condition_master.correct_response( st_idx & ...
                            condition_master.stim_class==curr_sc(jj) & ...
                            condition_master.task_class==curr_tc(mm));
                        n = histcounts(resp,1:(nr_responses+1));
                        
                        % For WHAT/HOW tasks we go by category info, because we combine
                        % object and foods into one button press..
                        if ismember(curr_sc(jj),[4,5]) && ismember(curr_tc(mm),8) % **** !!! careful this is hardcoded stuff !!! *****
                            supercat = condition_master.super_cat(st_idx & condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                            if curr_sc(jj) == 4 % objects
                                cued0  = condition_master.is_cued(st_idx & condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
                                n0     = histcounts([supercat(cued0==1,1);supercat(cued0==2,2)],1:6);
                            elseif curr_sc(jj) == 5 % scenes, no left/right, only center
                                n0     = histcounts(supercat,1:6);
                            end
                            n1 = [n0([1,2]), n0(3)+n0(4), n0(5)];  % we combine (3) objects and (4) food into a single button press (the ring finger)
                            
                            assert(all(n==n1))
                            
                        elseif ismember(curr_sc(jj),4) && ismember(curr_tc(mm),4) % PC-OBJ
                            resp2 = resp(condition_master.is_objectcatch( st_idx & ...
                                condition_master.stim_class==curr_sc(jj) & ...
                                condition_master.task_class==curr_tc(mm))==0); % only check non-obj catch  trials
                            n2 = histcounts(resp2,1:(nr_responses+1));
                            
                            if all(diff(n2)==0)
                                % ideally all responses are balanced
                            elseif all(diff(n)~=0) % but if not, we assume we have an uneven nr of trials due to objectcatch trials.
                                if mod(length(resp2),2)==1 % if we have an uneven nr of trials after we removed objectcatch trials
                                    if ~(diff(n2)<=1) % we allow for a difference of one, but not more
                                        error('[%s]: Button response counts diverge more than 1 and are considered unbalanced across options! Please rerun vcd_createConditions.m',mfilename);
                                    end
                                else % assume we have an even nr of trials after we removed objectcatch trials
                                    if all(diff(n2)~=n0) % we don't allow for a difference of one
                                        error('[%s]: Uneven nr of button responses! Please rerun vcd_createConditions.m',mfilename);
                                    end
                                end
                            end
                            
                        elseif ismember(curr_sc(jj),5) && ismember(curr_tc(mm),9) % NS-WHERE
                            subcat = condition_master.sub_cat(st_idx & condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                            n0     = histcounts(subcat,1:(nr_responses+1));
                            assert(all(n==n0))
                            
                        elseif ismember(curr_sc(jj),[4,5]) && ismember(curr_tc(mm),10) % OBJ-HOW & NS-HOW
                            affordcat = condition_master.affordance_cat(st_idx & condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                            if curr_sc(jj) == 4 % OBJ-HOW
                                cued0 = condition_master.is_cued(st_idx & condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
                                n1    = histcounts([affordcat(cued0==1,1);affordcat(cued0==2,2)],1:(nr_responses+1));
                                
                            elseif curr_sc(jj) == 5 % NS-HOW
                                n1     = histcounts(affordcat,1:(nr_responses+1));
                            end
                            
                            assert(all(n==n1))
                            
                        elseif ismember(curr_sc(jj),99) || ismember(curr_tc(mm),3) % scc-all
                            % Check if we sample all core images across trials
                            stmclass0           = condition_master.stim_class_name(st_idx & condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                            cued0               = condition_master.is_cued(st_idx & condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
                            stmclass0_cued      = [stmclass0(cued0==1,1);stmclass0(cued0==2,2)];
                            [~,stmclass_cued_i] = ismember(stmclass0_cued,params.exp.stimclassnames([1,3,2,4]));
                            
                            n1 = histcounts(stmclass_cued_i,1:(nr_responses+1));
                            
                            assert(all(n==n1))
                            
                        elseif all(n==0) && curr_tc(mm)~=2 && ~ismember(curr_tc(mm),[8,9,10])
                            error('[%s]: No correct responses found?!',mfilename);
                        else
                            assert(all(diff(n)==0)); % assert balanced button presses
                        end
                    end
                end
            end
        end
    end
end

return