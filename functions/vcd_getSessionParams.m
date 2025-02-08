function p = vcd_getSessionParams()

% Preallocate space
p = struct('session',[],'run',[],'miniblock',[],'trial', []);

%% Define big stim-task crossing table
p.stimClassLabels = {'gabor','rdk','dot','cobj','ns'};
p.taskClassLabels = {'fix','cd','scc','pc','wm','ltm','img','what','where','how'};

p.crossings = false(length(p.stimClassLabels),length(p.taskClassLabels));

p.stimTaskLabels = cell(length(p.taskClassLabels),length(p.stimClassLabels));
for row = 1:size(p.crossings,1)
    for col = 1:size(p.crossings,2)
        p.stimTaskLabels{col,row} = sprintf('%s-%s',lower(p.taskClassLabels{col}),lower(p.stimClassLabels{row}));   
    end
end

% Set Classic block
p.crossings(1:5,1:7) = true;
p.crossings(5,3) = false;

% Set Naturalistic tail
p.crossings(4:5,8:10) = true;
p.crossings(4,9) = false;

p.stimTaskLabels = p.stimTaskLabels(p.crossings');

% General exp params
p.n_unique_trial_repeats = 4;

% Set session params
p.session.task_start = [1,1,1,1,1,6,6,1,1,1];
p.session.runs_per_session = 10;

% Set miniblock params
p.miniblock.n_trials_single_epoch = 8;
p.miniblock.n_trials_double_epoch = 4;

% Set trial params
p.trial.single_epoch_tasks = logical([1 1 1 1 0 0 0 1 1 1]);
p.trial.double_epoch_tasks = ~p.trial.single_epoch_tasks;
p.trial.stim_LR_loc_cue    = logical([1 1 1 1 0]);
% p.trial.nr_unique_repeats  = 10; % <-- what's this for?


% Infer run and session params
p.run.n_single_epoch_miniblocks = 3;
p.run.n_double_epoch_miniblocks = 3;
p.run.miniblocks_per_run = p.run.n_single_epoch_miniblocks + p.run.n_double_epoch_miniblocks;








