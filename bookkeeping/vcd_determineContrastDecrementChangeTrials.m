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

% Get CD trials
cd_trials_idx  = condition_master.task_class == 2;
cd_trials      = condition_master(cd_trials_idx,:);

% What side is cued for each CD trial?
cd_blocks_cued_side = cd_trials.is_cued;

% what is the stimclass?
tmp_cued_side = cd_blocks_cued_side;
tmp_cued_side(tmp_cued_side==3) = 1;
cd_blocks_stimclass = cd_trials.stim_class;

% There should be no nans for "cd_blocks_cued_side"
assert(isequal(size(cd_trials,1),length(~isnan(cd_blocks_cued_side))))

% How many cued trials across all CD blocks do we have?
nr_cued_cd_trials = length(cd_blocks_cued_side);

% How often do we expect a trial with a CD change?
expected_nr_cd_trials = nr_cued_cd_trials*params.exp.trial.cd.prob_change;
expected_nr_cd_trials_per_cue = round(histcounts(cd_blocks_cued_side)*params.exp.trial.cd.prob_change);
assert(isequal(sum(expected_nr_cd_trials_per_cue),expected_nr_cd_trials));

%%% CHECK STIM CLASS
min_nr_cd_trials_per_stimclass = zeros(1,5);
min_nr_cd_trials_per_stimclass(unique(cd_trials.stim_class)') = floor(expected_nr_cd_trials/length(unique(cd_trials.stim_class)));

% Do some voodoo forcefully trying to balance the nr of cd changes across
% spatial cueing conditions, stim feature levels..
fprintf('[%s]: Sample cd trials',mfilename)
while 1
    fprintf('.')
    reshuffle_trials = [];
    
    % select cd trials
    when_cd_trials = sort(randsample(nr_cued_cd_trials,expected_nr_cd_trials)); % default is without replacement
    
    % Make sure we selected the nr of trials we expect for the given probability.
    assert(isequal(expected_nr_cd_trials,length(when_cd_trials)));
    
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
    for ii = 1:size(cd_trials,1)
        cond_mat(ii, cd_trials.is_cued(ii)) = 1;
        cond_mat(ii, length(unique(cd_trials.is_cued))+cd_trials.stim_class(ii)) = 1;
        
        cc = mod(cd_trials.is_cued(ii)-1,2)+1;
        
        switch cd_trials.stim_class(ii)
            case 1
                cond_mat(ii, length(params.exp.stimclassnames)+length(unique(cd_trials.is_cued))+ ...
                    find(cd_trials.contrast(ii,cc)==params.stim.gabor.contrast)) = 1;
            case 2
                cond_mat(ii, length(params.exp.stimclassnames)+length(unique(cd_trials.is_cued))+ length(params.stim.gabor.contrast) +...
                    find(cd_trials.rdk_coherence(ii,cc)==params.stim.rdk.dots_coherence)) = 1;
            case {4,5}
                cond_mat(ii, length(params.exp.stimclassnames)+length(unique(cd_trials.is_cued))+ length(params.stim.gabor.contrast) +...
                    length(params.stim.rdk.dots_coherence) + ...
                    cd_trials.super_cat(ii,cc)) = 1;
        end
        
    end
    

    
    selected_cd_trials = cond_mat(when_cd_trials,:);
    
    cond_cnt = sum(selected_cd_trials,1);
    
    % check cueing left/right
    if diff(cond_cnt([1,2]))~=0
        if cond_cnt(1)>cond_cnt(2) % if we have more left than right
            reshuffle_trials(1) = true;
        elseif cond_cnt(2)>cond_cnt(1) % if we have more right than left
            reshuffle_trials(1) = true;
        end
    end
    
    % check stim classes
    if any(cond_cnt(sc_idx) < min_nr_cd_trials_per_stimclass)
        reshuffle_trials(2) = true;
    end
    
    
    % check contrast levels if we have gabor trials
    if cond_cnt(sc_idx(1))>1 && any(cond_cnt(con_idx) >= 2)
        reshuffle_trials(3) = true;       
    end
    
    % check coherence levels if we have rdk trials
    if cond_cnt(sc_idx(2))>1 && any(cond_cnt(coh_idx) >= 2) 
        reshuffle_trials(4) = true;
    end
    
    
    % check super ordinate category levels if we have obj trials
    if cond_cnt(sc_idx(4))>1 && any(cond_cnt(sup_idx) >= 2)
        reshuffle_trials(5) = true;      
    end
    
    % check super ordinate category levels if we have ns trials
    if cond_cnt(sc_idx(5))>1 && any(cond_cnt(sup_idx) >= 2)
        reshuffle_trials(6) = true;      
    end
    
    if sum(reshuffle_trials)==0
        break;
    end
    
end
fprintf('Done!\n')


% Get onset time frame of the cd modulation within a trial
% (the modulation function has 2 frames of high contrast, prior to
% dip). The onset time refers to the time the actual contrast dip
% happens (so after the second time frames finishes and the third
% starts).
correct_response = 2.*ones(nr_cued_cd_trials,1); % all no's to begin with
correct_response(when_cd_trials)=1; % 1=yes change, 2=no change

cd_start = NaN(nr_cued_cd_trials,1);
for rpt = 1:length(when_cd_trials)
    cd_start(when_cd_trials(rpt)) = feval(params.stim.cd.cdsoafun);
end

% insert cd onset and correct response into the condition_master table
condition_master.cd_start(cd_trials_idx)         = cd_start;
condition_master.correct_response(cd_trials_idx) = correct_response;

% check nr of changes across all cd trials
assert(isequal(sum(correct_response==1),sum(~isnan(condition_master.cd_start))));
assert(isequal(sum(correct_response==2),sum(isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued>0)))));

% check nr of changes across left/right cueing conditions
if any(ismember(cd_blocks_cued_side,[1,2]))
    assert(isequal(sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==1))),...
        sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==2)))));
end
if any(cd_blocks_cued_side==3)
    assert(isequal(sum(~isnan(condition_master.cd_start(condition_master.task_class == 2 & condition_master.is_cued==3))),cond_cnt(3)));
end

% check expected probability to nr of cd trials with a change
assert(isequal(expected_nr_cd_trials,sum(~isnan(condition_master.cd_start))))
assert(isequal(expected_nr_cd_trials,sum(condition_master.correct_response(~isnan(condition_master.cd_start))==1)))


