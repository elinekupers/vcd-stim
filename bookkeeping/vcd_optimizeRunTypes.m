function optimized_runs = vcd_optimizeRunTypes(run_types, nr_blocks)
% VCD function to find the smallest number of runs to allocate stimulus blocks, 
% given the run-types available and run-type block contraints. For example,
% if we have 100 blocks, 25 double-stim and 75 single-stim, how many runs
% do we need to allocate those blocks, if we can only use two types of runs:
% one that has space for 4 single and 2 double blocks, and one that has
% space for 7 single and 0 double blocks.
%
%    optimized_runs = vcd_optimizeRunTypes(params, run_types_to_optimize, nr_blocks)
%
% INPUTS:
%  run_types  : (int) 2xN nr of [single;double] blocks allowed per run.
%  nr_blocks  : (double) nr of [single-stim, double-stim] blocks to
%                optimize for. Can be a non-integer.
% 
% Example:
% vcd_optimizeRunTypes([6, 4, 2, 1; 1, 2, 4, 5], [44,33])
%
% Written by E Kupers @ UMN [2025/04]


% Objective function: minimize total number of runs
f_obj = [ones(1,size(run_types,2))];

% Lower bound is 0 runs
lb  = zeros(size(run_types,2),1);

% minimize nr of runs
x = linprog(f_obj,[],[],run_types,nr_blocks,lb);
x = ceil(x); % can't have half runs

assert(isequal(sum((x'.*run_types),2)>=nr_blocks',[1;1]))

% Tell the user
optimized_runs = x;



