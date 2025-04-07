function optimized_runs = vcd_optimizeRunTypes(params, run_types_to_optimize, nr_blocks)

% nr_blocks     : nr of [single-stim, double-stim] blocks to optimize for
% 
% Example:
%
% nr_blocks = [233,     259];
% nr_blocks = [931.5,  1035];
% nr_blocks = [16.2,   27.5];
% vcd_optimizeRunTypes(nr_blocks)

all_run_types = cat(2, params.exp.run.run_typeA, ... [7; 0]; % single-stim, double-stim blocks within a run --> counts to ~ 330 s
    params.exp.run.run_typeB, ...% = [4, 2];
    params.exp.run.run_typeC, ...% = [1, 4];
    params.exp.run.run_typeD); % = [0, 5]; 
    
run_types = all_run_types(:,run_types_to_optimize);

% Objective function: minimize total number of runs
f_obj = [ones(1,size(run_types,2))];

lb  = zeros(size(run_types,2),1);

% minimize nr of runs
x = linprog(f_obj,[],[],run_types,nr_blocks,lb);
x = ceil(x);

assert(isequal(sum((x'.*run_types),2)>=nr_blocks',[1;1]))

optimized_runs = x;



