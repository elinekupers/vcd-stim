function [trial_order, swap_ok] = vcd_doubleDotAngleCheck(min_angle, dot_idx, orient_dirs, spatial_cues)
% VCD function to check if trials within a stimulus block contain single dot
% stimuli on the left and right side that are too close to one another.  
% If they are closer than provided threshold by <min_angle> than we try to
% see if there is another trial with two single dot stimuli in the same 
% stimulus block, with the same spatial cue direction, and if swapping one 
% of the stimuli (left or right) will prevent the two single dots from being 
% too close. 
%
% If swapping dot stimuli across trials won't solve the issue, or if there
% are no "double" dot trials with the same spatial cue, or if there are no
% trials with dot stimuli at all, this function will return an empty
% vector.
%
% INPUTS:
%   dot_idx      :  (logical) which trials contain single dot stimuli 
%                       (dims: n trials x 2 spatial locations)
%                        (col1 = left, col2 = right). Can be empty.
%   orient_dirs  :  (double) tilt/motion direction/facing direction of stimuli 
%                       (expected dims: n trials x 2 spatial locations (col1 = left, col2 = right)
%   spatial_cues :  (double) which side is cued by the central covert spatial 
%                       attention cue (1 = left, 2 = right, 3 = both sides). 
%                       (expected dims: n trials x 1)
%   min_angle    :  (double) minimum distance required between left and right dot 
%                       (in angular degrees, where 0 = 12 o'clock, 90 = 3 o'clock) 
%                       this info is stored in params.stim.dot.min_ang_distance_dot.
%
% OUTPUTS:
%   trial_order  :  (double) reorded stimulus indices that would fix the 
%                   issue of a trial with two single dot stimuli that are
%                   too close to one another. If new_order is empty, there %
%                   is no solution. (n trials x 2 spatial locations) (col1 = left, col2 = right)
%   new_orient   :  (double) reorded stimulus orientation in degrees. 
%                   (n trials x 2 spatial locations) (col1 = left, col2 = right)
%   swap_ok      :  (double) if a swap is needed and possible (swap_ok=1),
%                   if a swap is needed, but wasn't possible (swap_ok=-1),
%                   if no swap was needed (swap_ok= 0)
%
% Examples:
% [trial_order, swap_ok] = vcd_doubleDotAngleCheck(45, [2,4], [90 90 30 60; 180 60 180 40]', [1 2 1 2])
% [trial_order, swap_ok] = vcd_doubleDotAngleCheck(45, [2,4], [90 90 30 60; 180 60 180 40]', [1 2 1 2])
%
% Written by Eline Kupers @ UMN (Jan 2026)

%% Check inputs and set trial order

% create index of stimulus order
trial_order = repmat([1:size(orient_dirs,1)],[2 1])'; % nr indicates trial index, col1 = stim left, col2 = stim right

% check if trials with single dot stimuli share the same spatial cue: 
n2 = histcounts(spatial_cues(dot_idx),[1:3]); % [1:3] --> histcount bins = [0,1,2]
    
% Caclulate the difference in angle between to dots
diff_orient = abs(circulardiff(orient_dirs(dot_idx,1),orient_dirs(dot_idx,2),360));

% Do we need to swap?
if any(diff_orient <= min_angle) % find the orientations that are too small
    swap_bool = true; % we should swap
else
    swap_bool = false; % we should NOT swap, not sure why we are even here..
end

swap_ok = [];

% If we need to swap could possibly do so
if swap_bool && any(n2 > 1)

    for jj = 1:length(n2) % go over left/right spatial cues (1=left, 2=right)
        
        % find single dot trials with the same spatial cues
        same_cue_idx = dot_idx(ismember(spatial_cues(dot_idx),jj));
        
        if ~isempty(same_cue_idx)
            
            if length(same_cue_idx)==2 % if we have two "double dot" scc/ltm trials, then:
                
                % check if we can do something
                if any(find(diff_orient(ismember(dot_idx,same_cue_idx)) <= min_angle))
                                
                    new_order = same_cue_idx([2,1]); % swap the two indices

                    % make a copy of the current dot orientations
                    tmp0 = orient_dirs;

                    % shuffle the right side, keep the left side the same
                    tmp0(same_cue_idx,2) = tmp0(new_order,2);

                    % recalculate difference angles
                    tmp_diff_orient = abs(circulardiff(tmp0(:,1),tmp0(:,2),360));

                    % check if this solved anything
                    if all(tmp_diff_orient(dot_idx(ismember(dot_idx,same_cue_idx)))>min_angle)
                        % if it did, update indices in current list and master list
                        trial_order(same_cue_idx,2) = new_order;
                        swap_ok(jj) = 1;
                    else
                        swap_ok(jj) = -1;
                    end
                else
                    swap_ok(jj) = 0; % no need to swap
                end

            elseif length(same_cue_idx)>2 % if we have more than two "double dot" scc or ltm trials, then:
            
                % find the orientations that are too small
                swap_me   = find(diff_orient(ismember(dot_idx,same_cue_idx))<=min_angle);

                % set apart the other orientations
                other_ori = setdiff([1:length(same_cue_idx)], swap_me); 

                if isempty(other_ori)
                    % shuffle order
                    new_order = same_cue_idx(shuffle_concat(swap_me,1)); 
                else
                    % place orientation that needs to be swapped first and add the rest
                    new_order = same_cue_idx([sort(cat(1,swap_me,other_ori(1)),'descend'); other_ori(2:end)']); 
                end
            
                % make a copy of the current dot orientations
                tmp0 = orient_dirs;

                % shuffle the right side, keep the left side the same
                tmp0(same_cue_idx,2) = tmp0(new_order,2);

                % recalculate difference angles
                tmp_diff_orient = abs(circulardiff(tmp0(:,1),tmp0(:,2),360));

                % check if this solved anything
                if all(tmp_diff_orient(dot_idx)>min_angle)
                    % if it did, update indices in current list and master list
                    trial_order(same_cue_idx,2) = new_order;
                    swap_ok(jj) = 1;
                else
                    swap_ok(jj) = -1;
                end
            
            
            elseif length(same_cue_idx)==1 % if we have only one double dot trial, then we cannot swap because we want to keep the spatial cue consistent
                % do nothing but setting swap_ok to -1
                swap_ok(jj)   = -1;
            end
        else
            % do nothing but setting swap_ok to 0
            swap_ok(jj) = 0;
        end
    end
    
    % Consolidate output
    if isequal(swap_ok,[1, 1]) || isequal(sort(swap_ok),[0, 1]) || isequal(sort(swap_ok),[-1, 1]) % if a swap for at least one spatial cue worked, we call it a success
        swap_ok = 1;
    elseif isequal(swap_ok,[0, 0])  % no swap was needed for both spatial cues
        swap_ok = 0;
    elseif isequal(swap_ok,[-1, -1]) % if swap was needed for both spatial cues, but neither worked out.
        swap_ok = -1;
    elseif isequal(sort(swap_ok),[-1, 0]) % if swap was needed for both spatial cues, but neither worked out.
        swap_ok = -1;
    else
        error('wtf')
    end
    
elseif swap_bool && all(n2 <= 1)
    swap_ok = -1;
    
elseif ~swap_bool 
    % do nothing
    swap_ok = 0;
else
    error('wtf')
end



return


