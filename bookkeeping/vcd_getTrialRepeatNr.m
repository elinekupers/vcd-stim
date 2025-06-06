function condition_master = vcd_getTrialRepeatNr(condition_master)


% Get unique conditions and repeat nr for each trial
unique_cond_nrs = NaN(size(condition_master,1),2);

for cols = [1,2] % keep left/right cols separate for now
    [uniquecond,uniquerow] = unique(condition_master.condition_nr(:,cols),'stable');
    unique_cond_nrs(uniquerow,cols) = uniquecond;
end

% Check if stim left and right have overlapping conditions
col1    = unique_cond_nrs(~isnan(unique_cond_nrs(:,1))>0,1);
col2    = unique_cond_nrs(~isnan(unique_cond_nrs(:,2))>0,2);
[bi,ci] = intersect(col1,col2);
assert(isempty(bi)); assert(isempty(ci));
unique_cond_nrs2 = unique([col1;col2]);

repeat_nr = NaN(size(condition_master,1),2);
for uc = 1:size(unique_cond_nrs,1)
    for col = [1,2]
        ai = find(unique_cond_nrs(uc,col)==condition_master.condition_nr(:,col));
        repeat_nr(ai, col) = 1:length(ai);
    end
end

% update repeat_nr col
condition_master.repeat_nr = repeat_nr;

return