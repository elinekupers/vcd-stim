function condition_master = vcd_getTrialRepeatNr(condition_master)

% Get catch trials
catch_trials = condition_master.is_catch==1; 

% Get unique conditions and repeat nr for each trial
unique_cond_nrs = NaN(size(condition_master,1),2);

% keep left/right cols separate
unique_cond_nrs_no_catch = unique(condition_master.condition_nr(~catch_trials,:),'rows','stable');

% Check if stim left and right have overlapping conditions
col1    = unique_cond_nrs_no_catch(~isnan(unique_cond_nrs_no_catch(:,1))>0,1);
col2    = unique_cond_nrs_no_catch(~isnan(unique_cond_nrs_no_catch(:,2))>0,2);
[bi,ci] = intersect(col1,col2);

% assert that condition numbers for non-catch trials do not overlap between
% left/right stimulus locations (col1 and col2)
assert(isempty(bi)); assert(isempty(ci));

% include catch trials (which are the same for left and right)
unique_cond_nrs_with_catch = unique(condition_master.condition_nr,'rows','stable');

repeat_nr = NaN(size(condition_master,1),2);
for uc = 1:size(unique_cond_nrs_with_catch,1)
    if any(isnan(unique_cond_nrs_with_catch(uc,:)))
        assert(find(isnan(unique_cond_nrs_with_catch(uc,:)))==2)
        ai = find(ismember(condition_master.condition_nr(:,1),unique_cond_nrs_with_catch(uc,1)));
        assert(isequalwithequalnans(condition_master.condition_nr(ai,2), nan(length(ai),1)))
    else
        ai = find(sum(ismember(condition_master.condition_nr,unique_cond_nrs_with_catch(uc,:)),2)==2);
    end
    repeat_nr(ai) = [1:length(ai)]';
end

% update repeat_nr col
condition_master.repeat_nr = repeat_nr;

return