function p = vcd_getTaskParams()

% Preallocate space
p = struct('session',[],'run',[],'miniblock',[],'trial', []);

%% Define big stim-task crossing table
p.stimClassLabels = {'gabor','rdk','dot','cobj','ns'};
p.taskClassLabels = {'FIX','CD','SCC','PC','WM','LTM','IMG','WHAT','WHERE','ACT'};

p.crossings = false(length(p.stimClassLabels),length(p.taskClassLabels));

% Set Classic block
p.crossings(1:5,1:7) = true;
p.crossings(5,3) = false;

% Set Naturalistic tail
p.crossings(4:5,8:10) = true;
p.crossings(4,9) = false;

% Set session params
p.session.task_start = [1,1,1,1,1,4,4,1,1,1];
p.sesion.runs_per_session = 12;

% Set general miniblock params
p.miniblock.n_trials_single_epoch = 8;
p.miniblock.n_trials_double_epoch = 4;

p.run.n_single_epoch_miniblocks = 3;
p.run.n_double_epoch_miniblocks = 3;

p.sesion.miniblocks_per_session = p.run.n_single_epoch_miniblocks + p.run.n_double_epoch_miniblocks;


% Set trial design params
p.trial.single_epoch_tasks = logical([1 1 1 1 0 0 0 1 1 1]);
p.trial.double_epoch_tasks = ~p.trial.single_epoch_tasks;
p.trial.stim_LR_loc_cue    = logical([1 1 1 1 0]);
p.trial.nr_unique_repeats  = 10;





