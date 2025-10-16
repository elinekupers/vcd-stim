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
    
    % Get catch trials
    catch_idx = condition_master0.is_catch==1;
    
    % Get CD trials
    cd_trials_idx  = condition_master0.task_class == 2 & ~catch_idx;
    cd_trials      = condition_master0(cd_trials_idx,:);
    
    % What side is cued for each CD trial?
    cd_blocks_cued_side = cd_trials.is_cued;
    
    % There should be no nans for "cd_blocks_cued_side" (we let this constraint
    % loose for demo sessions)
    if ~params.is_demo
        assert(isequal(size(cd_trials,1),length(~isnan(cd_blocks_cued_side))))
    end
    
    % How many cued trials across all CD blocks do we have?
    nr_cued_cd_trials = length(cd_blocks_cued_side);
    
    % How often do we expect a trial with a CD change?
    expected_nr_cd_trials = round(nr_cued_cd_trials*params.exp.trial.cd.prob_change);
    expected_nr_cd_trials_per_cue = round(histcounts(cd_blocks_cued_side)*params.exp.trial.cd.prob_change);
    
    % If the rounding of separate left/right/neutral cued trials results
    % in 1 less/more CD+ trial than the overal expected nr of CD+ trials
    % across all cueing conditions, then we add/remove one trial
    % to/from expected nr of CD+ trials across spatial cues. Here we pick
    % the cue type with the most/least nr of trials.
    if sum(expected_nr_cd_trials_per_cue) < expected_nr_cd_trials
        if isequal(sum(expected_nr_cd_trials_per_cue)+1, expected_nr_cd_trials)
            if exist('add_me','var') && st == 2 % don't add a trial to the same cue type twice
                [~,mi] = min(expected_nr_cd_trials_per_cue);
                tmp    = intersect(mi, setdiff(1:length(expected_nr_cd_trials_per_cue),add_me));
                add_me = datasample(tmp,1);
            else
                [~,mi] = min(expected_nr_cd_trials_per_cue);
                add_me = randi(mi,1);
            end
            expected_nr_cd_trials_per_cue(add_me) = expected_nr_cd_trials_per_cue(add_me)+1;
        end
    elseif sum(expected_nr_cd_trials_per_cue) > expected_nr_cd_trials
        if isequal(sum(expected_nr_cd_trials_per_cue)-1, expected_nr_cd_trials)
            if exist('remove_me','var') && st == 2 % don't remove a trial from the same cue type twice
                [~,ma]    = max(expected_nr_cd_trials_per_cue);
                tmp       = intersect(ma, setdiff(1:length(expected_nr_cd_trials_per_cue),remove_me));
                remove_me = datasample(tmp,1);
            else
                [~,ma]    = max(expected_nr_cd_trials_per_cue);
                remove_me = randi(ma,1);
            end
            expected_nr_cd_trials_per_cue(remove_me) = expected_nr_cd_trials_per_cue(remove_me)-1;
        end
    end
    assert(isequal(sum(expected_nr_cd_trials_per_cue),expected_nr_cd_trials));
    
    % How many CD+ trials per stimulus class?
    min_nr_cd_trials_per_stimclass = zeros(1,5);
    unique_stim_class = unique(cd_trials.stim_class)';
    min_nr_cd_trials_per_stimclass(5) = floor(sum(expected_nr_cd_trials_per_cue(3))/sum(ismember(unique_stim_class,5)));
    cd_l_r = sum(expected_nr_cd_trials_per_cue(1:2));
    if mod(cd_l_r/sum(ismember(unique_stim_class,[1:4])),1)~=0
        % if we deal with uneven nr of trials, we distribute the left over trials across
        % the stimulus classes
        min_nr_cd_trials_per_stimclass(unique_stim_class(unique_stim_class<5)) = floor(cd_l_r/sum(ismember(unique_stim_class,[1:4])));
        remainder = rem(cd_l_r,sum(ismember(unique_stim_class,[1:4])));
        bonus_idx = randsample(unique_stim_class(ismember(unique_stim_class,[1:4])),remainder);
        min_nr_cd_trials_per_stimclass(bonus_idx)=min_nr_cd_trials_per_stimclass(bonus_idx)+1;
    else
        min_nr_cd_trials_per_stimclass(ismember(unique_stim_class,[1:4])) = cd_l_r/sum(ismember(unique_stim_class,[1:4]));
    end
    assert(isequal(sum(min_nr_cd_trials_per_stimclass),sum(expected_nr_cd_trials_per_cue)))
    assert(isequal(sum(expected_nr_cd_trials_per_cue),expected_nr_cd_trials));
    
    % Check how many conditions we have that we want to check for balancing:
    % 5 stim classes, 3 cue types, 3 contrast levels, 3 coherence levels, 5 super ordinate categories
    nConditions = length(params.exp.stimclassnames)+length(unique(cd_trials.is_cued))+ ...
        length(params.stim.gabor.contrast)+length(params.stim.rdk.dots_coherence)+length(params.stim.obj.super_cat);
    
    % Now we will create a NxM design matrix, where N are the rows (nr of
    % CD+ trials), and M are the columns, which represent spatial cueing
    % direction (1=left, 2=right, 3=neutral), stimulus class (4=GBR, 5=RDK,
    % 6=DOT, 7=OBJ, 8=NS), Gabor contrast levels (9=low, 10=medium,
    % 11=high), RDK motion direction (12=low, 13=medium, 14=high), OBJ and
    % NS superordinate semantic category (15=human, 16=animal, 17=object,
    % 18=food, 19=place).
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
                cond_mat(cc, base_column + ...
                    find(cd_trials.contrast(cc,cue_type(cc))==params.stim.gabor.contrast)) = 1;
            case 2
                cond_mat(cc, base_column + length(params.stim.gabor.contrast) + ...
                    find(cd_trials.rdk_coherence(cc,cue_type(cc))==params.stim.rdk.dots_coherence)) = 1;
            case {4,5}
                cond_mat(cc, base_column + length(params.stim.gabor.contrast) + length(params.stim.rdk.dots_coherence) + ...
                    cd_trials.super_cat(cc,cue_type(cc))) = 1;
        end
    end
    
    % Do some voodoo to balance the nr of cd changes across
    % spatial cueing conditions, stim feature levels..
    % basically we resample trials until we find the sample we like
    
    fprintf('[%s]: Sample cd trials for session type %d.',mfilename,st)
    max_attempts = 50000;
    max_ct2_attempts = 100000; % left cue (cue type 2)
    
    % We will shuffle trial order based on stimulus class and randomly select CD+ trials.
    sample_not_ok = true;
    while sample_not_ok
        
        % Reset flags, counter, and selected trials
        nr_attempts = 0; 
        ct2_shuffle_attempts = 0;

        if ismember(nr_attempts,round(linspace(1,max_attempts,20))); fprintf('.'); end
        selected_cd_trials = [];
        reshuffle_trials = zeros(1,length(min_nr_cd_trials_per_stimclass));
        
        %%%%% select cd trials %%%%%
        for ct = unique(cd_trials.is_cued)'
            
            switch ct
                
                case {1,2} % GBR-CD, RDK-CD, DOT-CD, OBJ-CD
                    
                    % Get trials with the given spatial cue
                    curr_cued_trials = cond_mat(:,ct)==1;
                    
                    if ct == 1
                        % Reset flags for classic stim, which have left/right spatial cues
                        reshuffle_trials(1:4) = 0;
                        sample_not_ok = false;
                    end
                    
                    % Sample trials, get their stimulus class
                    sc_to_sample = unique(cd_trials.stim_class)';
                    sc_to_sample = sc_to_sample(sc_to_sample<5); % only classic stim here
                    
                    sampled_stimclass = shuffle_concat(sc_to_sample, ceil(expected_nr_cd_trials_per_cue(ct)/length(sc_to_sample))); % shuffle order of sampled stimulus classes
                    sampled_stimclass = sampled_stimclass(1:expected_nr_cd_trials_per_cue(ct)); % truncate in case we oversampled
                    
                    % After selecting left cued trials, we will select right
                    % cued trials and make sure that the combined set of trials
                    % have a distribution that is equally distributed across
                    % the 4 classic stimulus classes
                    if ct == 2
                        % try finding a shuffle that works with the sample for
                        % left cued stimuli..
                        ct2_shuffle_attempts = 0;
                        if ~isempty(selected_cd_trials)
                            while 1
                                ct2_shuffle_attempts = ct2_shuffle_attempts+1;
                                
                                curr_stimclass_distribution = histcounts(cd_trials.stim_class(selected_cd_trials),[1:5]);       % get CD trials per stimulus class
                                left_to_sample =  min_nr_cd_trials_per_stimclass(1:length(curr_stimclass_distribution)) - curr_stimclass_distribution;
                                if expected_nr_cd_trials_per_cue(ct) <= length(sc_to_sample) % if we expect as many trials as we have stimulus classes
                                    sampled_stimclass = repmat(find(left_to_sample),1,ceil(sum(left_to_sample)/length(find(left_to_sample))));
                                else
                                   sampled_stimclass     = shuffle_concat(repelem(sc_to_sample,left_to_sample(sc_to_sample)),1); % shuffle the stimulus class samples
                                   
                                end
                                sampled_stimclass     = sampled_stimclass(1:expected_nr_cd_trials_per_cue(ct));
                                sampled_distribution  = histcounts(sampled_stimclass,[1:5]);
                                tmp = (sampled_distribution + curr_stimclass_distribution);
                                sample_not_ok         = any( tmp(sc_to_sample) ~= min_nr_cd_trials_per_stimclass(sc_to_sample) );
                                
                                if ~sample_not_ok
                                    break;
                                elseif ct2_shuffle_attempts > max_ct2_attempts
                                    error('[%s]: Reached max nr of attempts!',mfilename)
                                end
                            end
                        end
                    end
                    
                    if ~sample_not_ok % if sampled trials are ok, we continue with our checks
                        
                        for sc = sc_to_sample % only classic stim here
                            
                            while 1
                                nr_attempts = nr_attempts+1;
                                % if the session includes this stimulusclass
                                if any(cd_trials.stim_class==sc)
                                    
                                    % Sample cued trials for given stimclass
                                    sampled_trials = datasample(find(curr_cued_trials & cd_trials.stim_class==sc), sum(sampled_stimclass==sc), 'Replace',false); % sample without replacement
                                    
                                    stimclass_trials = cond_mat(sampled_trials,3+sc); % columns: spatial cueing direction (1=left, 2=right, 3=neutral), stimulus class (4=GBR, 5=RDK, 6=DOT, 7=OBJ, 8=NS)
                                    assert(all(stimclass_trials==1))
                                    assert(all(cond_mat(sampled_trials,ct)==1))
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
                                            constraint = repmat(floor(length(sampled_trials)/length(params.stim.gabor.contrast)),1,length(params.stim.gabor.contrast));
                                            %                                         constraint = repmat(floor(min_nr_cd_trials_per_stimclass(sc)/length(params.stim.gabor.contrast)),1,length(params.stim.gabor.contrast));
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
                                            constraint = repmat(floor(length(sampled_trials)/length(params.stim.rdk.dots_coherence)),1,length(params.stim.rdk.dots_coherence));
                                            %                                         constraint = repmat(floor(min_nr_cd_trials_per_stimclass(sc)/length(params.stim.rdk.dots_coherence)),1,length(params.stim.rdk.dots_coherence));
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
                                            constraint = repmat(floor(length(sampled_trials)/length(params.stim.obj.super_cat)),1,length(params.stim.obj.super_cat));
                                            %                                         constraint = repmat(floor(min_nr_cd_trials_per_stimclass(sc)/length(params.stim.obj.super_cat)),1,length(params.stim.obj.super_cat));
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
                                else
                                    break;
                                end
                            end
                        end
                    end
                    
                case 3 % NS-CD

                    if ~sample_not_ok
                        
                        % Get trials with the given spatial cue
                        curr_cued_trials = cond_mat(:,ct)==1;
                        
                        % Get stim class
                        sc = find(ismember(params.exp.stimclassnames,'ns'));
                        
                        while 1
                            nr_attempts = nr_attempts+1;
                            
                            % Sample neutral cued trials
                            sampled_stimclass = ones(1,expected_nr_cd_trials_per_cue(ct));
                            
                            % Sample cued trials for given stimclass
                            sampled_trials = datasample(find(curr_cued_trials & cd_trials.stim_class==sc), sum(sampled_stimclass), 'Replace',false); % sample without replacement
                            
                            stimclass_trials = cond_mat(sampled_trials,8); % col 8 is NS
                            
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
                                sample_not_ok = false;
                                break;
                            else
                                if nr_attempts > max_attempts
                                    error('\n[%s]: Can''t reach a solution for CD trial selection after %d attempts!',mfilename, max_attempts)
                                end
                            end
                        end
                    end
            end
        end
    end
    
    
    % Make sure we selected the nr of trials we expect for the given probability.
    assert(isequal(expected_nr_cd_trials,length(selected_cd_trials)));
    
    % Collapse the matrix distribution across trials
    cd_trials_mat = cond_mat(sort(selected_cd_trials),:);
    cond_cnt      = sum(cd_trials_mat,1);
    
    % check cueing left/right
    if isequal(cond_cnt(1),cond_cnt(2))
        assert(isequal(cond_cnt(1),cond_cnt(2))) % nr of trials should be the same
    elseif any(mod(cond_cnt([1,2]),2)==1)
        assert(isequal(abs(diff([cond_cnt(1),cond_cnt(2)])),1)) % nr of trials for left vs right cue should differ by 1 due to uneven nr of trials
    end
    
    % check stim classes
    assert(all(cond_cnt(3+sc_to_sample) >= min_nr_cd_trials_per_stimclass(sc_to_sample)))
    
    % check contrast levels if we have gabor trials
    if (cond_cnt(4)/2)>2  && any(cond_cnt(con_idx) < 1) % divide nr of trials by two to account for independent sampling of left/right cueing condition
        error('[%s]: welp, something went wrong in the sampling', mfilename)
    end
    
    % check coherence levels if we have rdk trials
    if (cond_cnt(5)/2)>2 && any(cond_cnt(coh_idx) < 1) % divide nr of trials by two to account for independent sampling of left/right cueing condition
        error('[%s]: welp, something went wrong in the sampling', mfilename)
    end
    
    % check super ordinate category levels if we have obj trials
    if (cond_cnt(7)/2)>2 && any(cond_cnt(sup_idx) < 1) % divide nr of trials by two to account for independent sampling of left/right cueing condition
        error('[%s]: welp, something went wrong in the sampling', mfilename)
    end
    
    % check super ordinate category levels if we have ns trials
    if (cond_cnt(8)/2)>2 && any(cond_cnt(sup_idx) < 1) % divide nr of trials by two to account for independent sampling of left/right cueing condition
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
    condition_master0.cd_start(cd_trials_idx & ~catch_idx)         = cd_start;
    condition_master0.correct_response(cd_trials_idx & ~catch_idx) = correct_response;
    
    % update condition name and number for cued sides with cd contrast
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
    
    % check left/right cueing    
    diff_l_r = abs(diff([sum(condition_master0.is_cued(~isnan(condition_master0.cd_start))==1),sum(condition_master0.is_cued(~isnan(condition_master0.cd_start))==2)]));
    assert( diff_l_r <= 1)
    
    % combine sessions
    condition_master(session_idx,:) = condition_master0;
end

% check nr of changes across all cd trials
assert(isequal(sum(condition_master.correct_response(condition_master.task_class==2 & condition_master.is_cued>0 & condition_master.is_catch==0)==1),sum(~isnan(condition_master.cd_start(condition_master.is_catch==0)))));
assert(isequal(sum(condition_master.correct_response(condition_master.task_class==2 & condition_master.is_cued>0 & condition_master.is_catch==0)==2),sum(isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued>0 & condition_master.is_catch==0)))));

% check nr of changes across left/right cueing conditions
if any(ismember(condition_master.is_cued(condition_master.task_class==2 & condition_master.is_cued>0),[1,2]))
    leftcued_cd_trials  = sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==1 & condition_master.is_catch==0)));
    rightcued_cd_trials = sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==2 & condition_master.is_catch==0)));
    if ~isequal(leftcued_cd_trials,rightcued_cd_trials)
        assert(abs(diff([leftcued_cd_trials,rightcued_cd_trials]))<=2);
    end
end

if any(condition_master.is_cued(condition_master.task_class==2 & condition_master.is_cued>0)==3)
    assert(isequal(sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==3))),...
        sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==3)))))
end

% check if the number CD+ trials matches with our probability
assert( ismember(sum(~isnan(condition_master.cd_start(condition_master.is_catch==0))),round(sum(condition_master.task_class(condition_master.is_catch==0)==2).*params.exp.trial.cd.prob_change)+[-1:1]) ); % we allow rounding error of 1 trial.

return


