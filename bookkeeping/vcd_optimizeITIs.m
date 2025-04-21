function optimized_ITIs = vcd_optimizeITIs(block_dur, trial_dur, itis, nr_trials)
% VCD function to find the best combination of ITIs to allocate between
% trials of a given block. 
% For example, if we have 8 trials (each 4.2 s) and we want 7 ITIs sampled 
% between 0.2:0.1:1.6 seconds, that will give us a max block duration of 42 s.
% This function will find those 7 ITIs that exactly sum to 8.4 s.
%
%    optimized_ITIs = vcd_optimizeITIs(max_dur, itis, nr_trials)
%
% INPUTS:
%  block_dur                 : (double) duration of block (in frames)
%  trial_dur                 : (double) duration of single trial (in frames)
%  itis                      : (double, length 1xN) ITIs to optimize for
%                               (in nr of frames)
%  nr_trials                 : (double) nr of trials within a block.
% 
% Example:
% optimized_ITIs = vcd_optimizeITIs(42,4.2,[6:2:48], 8)

% Calculate how much time we need to account for
dur_to_optimize = block_dur - sum(trial_dur*nr_trials);

% get the possible combinations of requested itis, given all itis
C = nchoosek(itis,nr_trials-1);

% get the sum of each combination of itis
total_dur = sum(C,2);

% Find those equal to the time we want
optimized_ITIs = C((total_dur == dur_to_optimize),:);

if isempty(optimized_ITIs)
    error('[%s]: Cannot find a combination of ITIs that works!!')
end
    
% randomly select one combinaton, if we have multiple
if size(optimized_ITIs,1)>1
    optimized_ITIs = optimized_ITIs(randi(size(optimized_ITIs,1),1),:);
end

% shuffle the order of the possible iti combinations
optimized_ITIs = shuffle_concat(optimized_ITIs,1);

return