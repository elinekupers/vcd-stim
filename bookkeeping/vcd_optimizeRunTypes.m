function optimized_runs = vcd_optimizeRunTypes(params, run_types_to_optimize, nr_blocks)
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
%  params                    : struct with parameters 
%       exp.run.run_type1       :  single-stim, double-stim blocks within
%                                   runtype 1 (for example: [7; 0])
%       exp.run.run_type2       :  single-stim, double-stim blocks within
%                                   runtype 2  (for example: [4; 2])
%       exp.run.run_type3       :  single-stim, double-stim blocks within
%                                   runtype 3  (for example: [1; 4])
%       exp.run.run_type4       :  single-stim, double-stim blocks within
%                                   runtype 4  (for example: [0; 5])
%  run_types_to_optimize     : (int) index for run type to use, e.g., [1:3]
%                               indicates only use runtypes 1, 2, and 3
%  nr_blocks                 : (double) nr of [single-stim, double-stim]
%                               blocks to optimize for. Can be a non-integer.
% 
% Example:
% vcd_optimizeRunTypes(params, [1:3], [233,259])

all_run_types = cat(2, params.exp.run.run_type1, ... [7; 0];
    params.exp.run.run_type2, ...% = [4; 2];
    params.exp.run.run_type3, ...% = [1; 4];
    params.exp.run.run_type4);   % = [0; 5]; 
    
run_types = all_run_types(:,run_types_to_optimize);

% Objective function: minimize total number of runs
f_obj = [ones(1,size(run_types,2))];

lb  = zeros(size(run_types,2),1);

% minimize nr of runs
x = linprog(f_obj,[],[],run_types,nr_blocks,lb);
x = ceil(x);

assert(isequal(sum((x'.*run_types),2)>=nr_blocks',[1;1]))

optimized_runs = x;



