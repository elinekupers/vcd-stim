function condition_master = vcd_determineContrastDecrementChangeTrials(params, condition_master)
% VCD bookkeeping function to determine which trials in a CD-block will
% have a contrast decrement and what time this will occur during the
% stimulus on period.
%
%   [condition_master] = vcd_determineContrastDecrementChangeTrials(condition_master)
%
% Steps:
% - Find CD blocks across all runs and sessions.
% - Calculate the number of trials across all sessions that will have a CD
%   change. Probability of CD change is determined by params.exp.trial.cd.prob_change
% - Make a [yes/no] CD change vector.
% - Apply yes/no change vector to cued side. For the uncued side, we do not
%   change contrast.
% - Ensure the nr of changes is equally distributed across cued sides
% - Get onset time of cd modulation within a trial
%
% The contrast temporal modulation function has 2 frames of high contrast,
% before it dips to a lower value. When we add the onset time, it will be
% in relative time frames from the onset of the stimulus, and refer to the
% first time frame that the contrast has changed, i.e., the third time
% frame of the defined contrast temporal modulation function.

% do it separately for each session
for st = unique(condition_master.session_type)'
    
    % if we have more than 1 session type, then update global counters separately.
    session_idx = find(condition_master.session_type==st);
    
    condition_master0 = condition_master(session_idx,:);
    
    % Get CD trials
    cd_trials_idx  = condition_master0.task_class == 2;
    cd_trials      = condition_master0(cd_trials_idx,:);
    
    % What side is cued for each CD trial?
    cd_blocks_cued_side = cd_trials.is_cued;
    
    % what is the stimclass?
    tmp_cued_side = cd_blocks_cued_side;
    tmp_cued_side(tmp_cued_side==3) = 1;
    cd_blocks_stimclass = cd_trials.stim_class;
    
    % There should be no nans for "cd_blocks_cued_side" (we let this constraint
    % loose for demo sessions)
    if ~params.is_demo
        assert(isequal(size(cd_trials,1),length(~isnan(cd_blocks_cued_side))))
    end
    
    % How many cued trials across all CD blocks do we have?
    nr_cued_cd_trials = length(cd_blocks_cued_side);
    
    % How often do we expect a trial with a CD change?
    expected_nr_cd_trials = round((nr_cued_cd_trials*params.exp.trial.cd.prob_change)/2)*2;
    expected_nr_cd_trials_per_cue = round(histcounts(cd_blocks_cued_side)*params.exp.trial.cd.prob_change);
    % If the rounding of separate left/right/neutral cued trials results
    % in 1 less CD+ trial than the overal expected nr of CD+ trials across
    % all cueing conditions, then we remove one trial from the overall
    % expected nr of CD+ trials
    if sum(expected_nr_cd_trials_per_cue) < expected_nr_cd_trials
       if isequal(sum(expected_nr_cd_trials_per_cue)+1, expected_nr_cd_trials)
           expected_nr_cd_trials = expected_nr_cd_trials-1;
       end
    elseif sum(expected_nr_cd_trials_per_cue) > expected_nr_cd_trials
        if isequal(sum(expected_nr_cd_trials_per_cue)-1, expected_nr_cd_trials)
            expected_nr_cd_trials = expected_nr_cd_trials+1;
        end
    end
    assert(isequal(sum(expected_nr_cd_trials_per_cue),expected_nr_cd_trials));
    
    %%% CHECK STIM CLASS
    min_nr_cd_trials_per_stimclass = zeros(1,5);
    min_nr_cd_trials_per_stimclass(unique(cd_trials.stim_class)') = floor(expected_nr_cd_trials/length(unique(cd_trials.stim_class)));
    
    % Check how many conditions we have that we want to check for balancing:
    % 5 stim classes, 3 cue types, 3 contrast levels, 3 coherence levels, 5 super ordinate categories
    nConditions = length(params.exp.stimclassnames)+length(unique(cd_trials.is_cued))+ ...
        length(params.stim.gabor.contrast)+length(params.stim.rdk.dots_coherence)+length(params.stim.obj.super_cat);
    
    cue_idx = 1:length(unique(cd_trials.is_cued));
    sc_idx  =  (max(cue_idx)+1):(max(cue_idx)+length(params.exp.stimclassnames));
    con_idx = (max(sc_idx)+1):(max(sc_idx)+length(params.stim.gabor.contrast));
    coh_idx = (max(con_idx)+1):(max(con_idx)+length(params.stim.rdk.dots_coherence));
    sup_idx = (max(coh_idx)+1):(max(coh_idx)+length(params.stim.obj.super_cat));
    
    cond_mat = zeros(size(cd_trials,1),nConditions);
    for ii = cue_idx
        cond_mat(:,ii) = double(cd_trials.is_cued==ii);
    end
    for ii = sc_idx
        cond_mat(:,ii) = double(cd_trials.stim_class==find(ii==sc_idx));
    end
    cue_type = mod(cd_trials.is_cued-1,2)+1;
    
    base_column = length(params.exp.stimclassnames)+length(unique(cd_trials.is_cued));
    for cc = 1:size(cue_type,1)
        
        switch cd_trials.stim_class(cc)
            case 1
                cond_mat(cc, + base_column + ...
                    find(cd_trials.contrast(cc,cue_type(cc))==params.stim.gabor.contrast)) = 1;
            case 2
                cond_mat(cc, + base_column + length(params.stim.gabor.contrast) + ...
                    find(cd_trials.rdk_coherence(cc,cue_type(cc))==params.stim.rdk.dots_coherence)) = 1;
            case {4,5}
                cond_mat(cc, + base_column + length(params.stim.gabor.contrast) + length(params.stim.rdk.dots_coherence) + ...
                    cd_trials.super_cat(cc,cue_type(cc))) = 1;
        end
    end
    
    % Do some voodoo to balance the nr of cd changes across
    % spatial cueing conditions, stim feature levels..
    % basically we resample trials until we find the sample we like
    
    fprintf('[%s]: Sample cd trials',mfilename)
    max_attempts = 200000;
    nr_attempts = 0;
    
    selected_cd_trials = [];
    
    %%%%% select cd trials %%%%%
    for ct = unique(cd_trials.is_cued)'
        
        % Get trials with the given spatial cue
        curr_cued_trials = cond_mat(:,ct)==1;
        
        % Reset flags
        reshuffle_trials = zeros(1,length(min_nr_cd_trials_per_stimclass));
        
        if ismember(ct,[1,2])
            sampled_stimclass = shuffle_concat([1:sum(unique(cd_trials.stim_class)<5)], ceil(expected_nr_cd_trials_per_cue(ct)/sum(unique(cd_trials.stim_class)<5)));
            sampled_stimclass = sampled_stimclass(1:expected_nr_cd_trials_per_cue(ct));
            
            if ct == 2
                curr_stimclass_distribution = histcounts(cd_trials.stim_class(selected_cd_trials),[1:5]);
                sampled_distribution        = histcounts(sampled_stimclass,[1:5]);
                sample_not_ok = any((sampled_distribution + curr_stimclass_distribution) < min_nr_cd_trials_per_stimclass([1:sum(unique(cd_trials.stim_class)<5)]));
                while sample_not_ok
                    sampled_stimclass = shuffle_concat([1:sum(unique(cd_trials.stim_class)<5)], ceil(expected_nr_cd_trials_per_cue(ct)/sum(unique(cd_trials.stim_class)<5)));
                    sampled_stimclass = sampled_stimclass(1:expected_nr_cd_trials_per_cue(ct));
                    sampled_distribution        = histcounts(sampled_stimclass,[1:5]);
                    sample_not_ok = any((sampled_distribution + curr_stimclass_distribution) < min_nr_cd_trials_per_stimclass([1:sum(unique(cd_trials.stim_class)<5)]));
                end
            end
            
            for sc = 1:sum(unique(cd_trials.stim_class)<5)
                
                while 1
                    nr_attempts = nr_attempts+1;
                    
                    % Sample cued trials for given stimclass
                    sampled_trials = datasample(find(curr_cued_trials & cd_trials.stim_class==sc), sum(sampled_stimclass==sc), 'Replace',false); % sample without replacement
                    
                    stimclass_trials = cond_mat(sampled_trials,sc_idx(sc));
                    
                    % Check constraints for each stimulus class
                    if sc == 1 % Gabors
                        % check contrast levels if we have gabor trials
                        nr_trials = sum(cond_mat(sampled_trials, con_idx),1);
                        % if we sampled fewer trials than nr of contrast levels
                        if length(sampled_trials) < length(params.stim.gabor.contrast)
                            % then make sure we sample as many different
                            % contrast levels as we have trials
                            constraint = zeros(1,length(params.stim.gabor.contrast));
                            constraint(1:length(sampled_trials)) = ones(1,length(sampled_trials));
                            if  any(sort(nr_trials) < sort(constraint))
                                reshuffle_trials(sc) = true;
                            else
                                reshuffle_trials(sc) = false;
                            end
                        else % otherwise we will check against the expected number of contrast levels
                            constraint = repmat(floor(min_nr_cd_trials_per_stimclass(sc)/length(params.stim.gabor.contrast)),1,length(params.stim.gabor.contrast));
                            if  any(nr_trials < constraint)
                                reshuffle_trials(sc) = true;
                            else
                                reshuffle_trials(sc) = false;
                            end
                        end
                        
                    elseif sc == 2 % RDKs
                        % check coherence levels if we have rdk trials
                        nr_trials = sum(cond_mat(sampled_trials(stimclass_trials==1), coh_idx),1);
                        % if we sampled fewer trials than nr of coherence levels
                        if length(sampled_trials) < length(params.stim.rdk.dots_coherence)
                            % then make sure we sample as many different
                            % coherence levels as we have trials
                            constraint = zeros(1,length(params.stim.rdk.dots_coherence));
                            constraint(1:length(sampled_trials)) = ones(1,length(sampled_trials));
                            if  any(sort(nr_trials) < sort(constraint))
                                reshuffle_trials(sc) = true;
                            else
                                reshuffle_trials(sc) = false;
                            end
                        else % otherwise we will check against the expected number of coherence levels
                            constraint = repmat(floor(min_nr_cd_trials_per_stimclass(sc)/length(params.stim.rdk.dots_coherence)),1,length(params.stim.rdk.dots_coherence));
                            if any(nr_trials < constraint)
                                reshuffle_trials(sc) = true;
                            else
                                reshuffle_trials(sc) = false;
                            end
                        end
                        
                    elseif sc == 3 % single dot
                        % do nothing
                    elseif sc == 4 % objects
                        % check super ordinate category levels if we have obj trials
                        nr_trials = sum(cond_mat(sampled_trials(stimclass_trials==1), sup_idx),1);
                        % if we sampled fewer trials than nr of superordinate categories
                        if length(sampled_trials) < length(params.stim.obj.super_cat)
                            % then make sure we sample as many different
                            % supercategories as we have trials
                            constraint = zeros(1,length(params.stim.obj.super_cat));
                            constraint(1:length(sampled_trials)) = ones(1,length(sampled_trials));
                            if  any(sort(nr_trials) < sort(constraint))
                                reshuffle_trials(sc) = true;
                            else
                                reshuffle_trials(sc) = false;
                            end
                        else % otherwise we will check against the expected number of superordinate categories
                            constraint = repmat(floor(min_nr_cd_trials_per_stimclass(sc)/length(params.stim.obj.super_cat)),1,length(params.stim.obj.super_cat));
                            if  any(nr_trials < constraint)
                                reshuffle_trials(sc) = true;
                            else
                                reshuffle_trials(sc) = false;
                            end
                        end
                    end
                    if sum(reshuffle_trials(:))==0
                        selected_cd_trials = cat(1,selected_cd_trials,sampled_trials);
                        break;
                    else
                        if nr_attempts > max_attempts
                            error('\n[%s]: Can''t reach a solution for CD trial selection after %d attempts!',mfilename, max_attempts)
                        end
                    end
                end
            end
            
        elseif ct == 3 % scenes % NS
            sc = 5;
            
            while 1
                nr_attempts = nr_attempts+1;
                
                % Sample neutral cued trials
                sampled_stimclass = ones(1,expected_nr_cd_trials_per_cue(ct));
                
                % Sample cued trials for given stimclass
                sampled_trials = datasample(find(curr_cued_trials & cd_trials.stim_class==sc), sum(sampled_stimclass), 'Replace',false); % sample without replacement
                
                stimclass_trials = cond_mat(sampled_trials,sc_idx(sc));
                
                % check super ordinate category levels if we have ns trials
                nr_trials = sum(cond_mat(sampled_trials(stimclass_trials==1), sup_idx),1);
                % if we sampled fewer trials than nr of superordinate categories
                if length(sampled_trials) < length(params.stim.ns.super_cat)
                    % then make sure we sample as many different
                    % supercategories as we have trials
                    constraint = zeros(1,length(params.stim.ns.super_cat));
                    constraint(1:length(sampled_trials)) = ones(1,length(sampled_trials));
                    if  any(sort(nr_trials) < sort(constraint))
                        reshuffle_trials(sc) = true;
                    else
                        reshuffle_trials(sc) = false;
                    end
                else % otherwise we will check against the expected number of superordinate categories
                    constraint = repmat(floor(min_nr_cd_trials_per_stimclass(sc)/length(params.stim.ns.super_cat)),1,length(params.stim.ns.super_cat));
                    if any(nr_trials < constraint)
                        reshuffle_trials(sc) = true;
                    else
                        reshuffle_trials(sc) = false;
                    end
                end
                
                if sum(reshuffle_trials(:))==0
                    selected_cd_trials = cat(1,selected_cd_trials,sampled_trials);
                    break;
                else
                    if nr_attempts > max_attempts
                        error('\n[%s]: Can''t reach a solution for CD trial selection after %d attempts!',mfilename, max_attempts)
                    end
                end
            end
        end
    end
    
    % Make sure we selected the nr of trials we expect for the given probability.
    assert(isequal(expected_nr_cd_trials,length(selected_cd_trials)));
    
    cd_trials_mat = cond_mat(sort(selected_cd_trials),:);
    cond_cnt      = sum(cd_trials_mat,1);
    
    % check cueing left/right
    assert(isequal(cond_cnt(1),cond_cnt(2))) % nr of trials should be the same
    
    % check stim classes
    assert(all(cond_cnt(sc_idx) >= min_nr_cd_trials_per_stimclass))
    
    % check contrast levels if we have gabor trials
    if (cond_cnt(sc_idx(1))/2)>2  && any(cond_cnt(con_idx) < 1) % divide nr of trials by two to account for independent sampling of left/right cueing condition
        error('[%s]: welp, something went wrong in the sampling', mfilename)
    end
    
    % check coherence levels if we have rdk trials
    if (cond_cnt(sc_idx(2))/2)>2 && any(cond_cnt(coh_idx) < 1) % divide nr of trials by two to account for independent sampling of left/right cueing condition
        error('[%s]: welp, something went wrong in the sampling', mfilename)
    end
    
    % check super ordinate category levels if we have obj trials
    if (cond_cnt(sc_idx(4))/2)>2 && any(cond_cnt(sup_idx) < 1) % divide nr of trials by two to account for independent sampling of left/right cueing condition
        error('[%s]: welp, something went wrong in the sampling', mfilename)
    end
    
    % check super ordinate category levels if we have ns trials
    if (cond_cnt(sc_idx(5))/2)>2 && any(cond_cnt(sup_idx) < 1) % divide nr of trials by two to account for independent sampling of left/right cueing condition
        error('[%s]: welp, something went wrong in the sampling', mfilename)
    end
    
    fprintf('Done!\n')
    
    
    % Get onset time frame of the cd modulation within a trial
    % (the modulation function has 2 frames of high contrast, prior to
    % dip). The onset time refers to the time the actual contrast dip
    % happens (so after the second time frames finishes and the third
    % starts).
    correct_response = 2.*ones(nr_cued_cd_trials,1); % all no's to begin with
    correct_response(selected_cd_trials)=1; % 1=yes change, 2=no change
    
    cd_start = NaN(nr_cued_cd_trials,1);
    for rpt = 1:length(selected_cd_trials)
        cd_start(selected_cd_trials(rpt)) = feval(params.stim.cd.cdsoafun);
    end
    
    % insert cd onset and correct response into the condition_master table
    condition_master0.cd_start(cd_trials_idx)         = cd_start;
    condition_master0.correct_response(cd_trials_idx) = correct_response;
    
    % update condition name and numbr for cued sides with cd contrast
    when_idx = find(cd_trials_idx);
    when_cd_trials_sub = when_idx(selected_cd_trials);
    for ii = 1:length(selected_cd_trials)
        condition_master0.condition_name{when_cd_trials_sub(ii),mod(condition_master0.is_cued(when_cd_trials_sub(ii))-1,2)+1} = ...
            catcell(2,[condition_master0.condition_name(when_cd_trials_sub(ii),mod(condition_master0.is_cued(when_cd_trials_sub(ii))-1,2)+1),{'+'}]);
        condition_master0.condition_nr(when_cd_trials_sub(ii),mod(condition_master0.is_cued(when_cd_trials_sub(ii))-1,2)+1) = ...
            vcd_conditionName2Number(condition_master0.condition_name{when_cd_trials_sub(ii),mod(condition_master0.is_cued(when_cd_trials_sub(ii))-1,2)+1});
    end
    
    % check expected probability to nr of cd trials with a change
    assert(isequal(expected_nr_cd_trials,sum(~isnan(condition_master0.cd_start))))
    assert(isequal(expected_nr_cd_trials,sum(condition_master0.correct_response(~isnan(condition_master0.cd_start))==1)))
    
    % combine sessions
    condition_master(session_idx,:) = condition_master0;
end

% check nr of changes across all cd trials
assert(isequal(sum(condition_master.correct_response(condition_master.task_class==2 & condition_master.is_cued>0)==1),sum(~isnan(condition_master.cd_start))));
assert(isequal(sum(condition_master.correct_response(condition_master.task_class==2 & condition_master.is_cued>0)==2),sum(isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued>0)))));

% check nr of changes across left/right cueing conditions
if any(ismember(condition_master.is_cued(condition_master.task_class==2 & condition_master.is_cued>0),[1,2]))
    assert(isequal(sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==1))),...
        sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==2)))));
end
if any(condition_master.is_cued(condition_master.task_class==2 & condition_master.is_cued>0)==3)
    assert(isequal(sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==3))),...
         sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==3)))))
end




