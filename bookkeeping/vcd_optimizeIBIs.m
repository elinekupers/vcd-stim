function [optimized_IBIs, postblank_to_add] = vcd_optimizeIBIs(run_dur, block_dur, ibis, nr_blocks, prepost_blank)
% VCD function to find the best combination of IBIs to allocate between
% blocks of a given run. 
%
% We want to sample IBIs such that the run duration will be exactly 363.2
% s, regardless of run type and nr of blocks in a run. This function will
% find those IBIs from the provided list of possible IBIs to exactly sum to
% the desired block length.
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
    postblank_to_add = dur_to_optimize-max_possible_IBI_dur;
    dur_to_optimize  = max_possible_IBI_dur;
    fprintf('[%s]: **** WARNING START ****\n', mfilename)
    fprintf('[%s]: IBIs cannot account for the total run duration. Will increase post-blank duration by %d time frames.\n', mfilename,postblank_to_add)
    fprintf('[%s]: **** WARNING END ****\n', mfilename)
    if postblank_to_add > 3600
        error('[%s]: **** We are adding more than a minute of blank?! That doesn''t seem right.. ****\n', mfilename)
    end
else
    postblank_to_add = 0;
end

% EK HACK
if mod(dur_to_optimize,30)>0 % want IBIs to be in increments of 0.5 sec
    frames_to_add = 30 - mod(dur_to_optimize,30);
    dur_to_optimize = dur_to_optimize + frames_to_add;
end
    
% Find a distributed combination of possible IBIs
% inputrs: totalamt,validoptions,numbins,mode,numlookback
distribution_mode = 0; % uniformly sampled bins
n_iterations      = 100; % when do we decide to give up before we find a solution
f = distributewithconstraints(dur_to_optimize, ibis, nr_blocks-1, distribution_mode, n_iterations);

if isempty(f)
    error('[%s]: Cannot find a combination of IBIs that works!!')
else
    % shuffle the order of the possible iti combinations
    optimized_IBIs = shuffle_concat(f,1);
end



return