function [optimized_IBIs, preblank_to_add, postblank_to_add] = vcd_optimizeIBIs(run_dur, block_dur, ibis, nr_blocks, prepost_blank)
% VCD function to find the best combination of IBIs to allocate between
% blocks of a given run. 
%
% We want to sample IBIs such that the run duration will be exactly 363.2
% s, regardless of run type and nr of blocks in a run. This function will
% find those IBIs to exactly sum to the desired block length.
%
%    optimized_IBIs = vcd_optimizeIBIs(block_dur, trial_dur, itis, nr_trials)
%
% INPUTS:
%  run_dur                   : (double) duration of single run (in frames)
%  block_dur                 : (double) duration of blocks (in frames), can
%                               be integer or vector of integers. 
%  itis                      : (double, length 1xN) IBIs to optimize for
%                               (in nr of frames)
%  nr_blocks                 : (double) nr of trials within a block.
% 
% Example:
% [optimized_IBIs, preblank_to_add, postblank_to_add] = vcd_optimizeIBIs(300, 42, [5:1:9], 7, 0)


if numel(block_dur) == nr_blocks
    dur_to_optimize = run_dur - prepost_blank - sum(block_dur);
else
    dur_to_optimize = run_dur - prepost_blank - sum(block_dur*nr_blocks);
end

max_possible_IBI_dur = max(ibis)*(nr_blocks-1);
if dur_to_optimize > max_possible_IBI_dur
    warning('[%s]: IBIs cannot account for the total run duration. Will use max IBI and increase pre- and post-blank duration', mfilename)
    diff_dur = dur_to_optimize-max_possible_IBI_dur;
    preblank_to_add  = floor(diff_dur/2);
    postblank_to_add = ceil(diff_dur/2);

    optimized_IBIs = max(ibis)*ones(1,nr_blocks-1);
else
    % get the possible combinations of requested itis, given all itis
    C = nchoosek(repmat(ibis,1,3),nr_blocks-1);
    
    % get the sum of each combination of itis
    total_dur = sum(C,2);
    
    % Find those equal to the time we want
    optimized_IBIs = C((total_dur == dur_to_optimize),:);
    
    if isempty(optimized_IBIs)
        error('[%s]: Cannot find a combination of IBIs that works!!')
    end
    
    % randomly select one combinaton, if we have multiple
    if size(optimized_IBIs,1)>1
        optimized_IBIs = optimized_IBIs(randi(size(optimized_IBIs,1),1),:);
    end
    
    % shuffle the order of the possible iti combinations
    optimized_IBIs = shuffle_concat(optimized_IBIs,1);
    
    preblank_to_add = 0;
    postblank_to_add = 0;
end

return