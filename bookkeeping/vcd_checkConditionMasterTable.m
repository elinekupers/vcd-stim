function [] = vcd_checkConditionMasterTable(condition_master)
% VCD bookkeeping function to check basic numbers about the VCD core
% experiment. 
%   
%    [] = vcd_checkConditionMasterTable(condition_master)
%
% We want to know checks like:
% -What is the total number of trials per run?
% -What is the total number of trials and blocks per stimulus class
% -What is the total number of trials and blocks per task class?
% -What is the total number of unique crossing across a session?
% -What is the total number of unique stimuli shown across a session?
% -How many distinct conditions do we have? (e.g., GBR-0001-L-CUED-PC, etc.),
% -How many repetitions do we have per distinct condition?
% -Does each stimulus class have at least one fixation block?
% -How many left/right/neutral spatial cues do we have within a block?
% -Are the spatial cue L/R balanced across blocks/runs/entire session for each crossing involving classic stimuli?
% -How many contrast levels, rdk coherence levels, dot angles, object & ns super/basic/sub categories do we sample across a session/run/block?
% -Do we sample all three coherence levels in an RDK block?
% -Do we sample all three contrast levels within a GBR block?
% -Do we sample all at least 4 superordinate categories within an OBJ or NS block?
% -To what extent are button presses balanced for each crossing?
% -Does this hold across blocks/runs/the entire session?
%
% INPUT:
% * condition_master    : (table) Table defining all trials across all 
%                         conditions for either the behavioral or MRI VCD 
%                         version of the experiment. Rows represent the
%                         individual trials, columns represent different
%                         aspects of the experiment/condition/trial/stimulus.
% OUTPUT:
% * None
%
% Written by E Kupers @ UMN 2025/06

% Bird's eye view stats
total_nr_of_sessions = unique(condition_master.session_nr);
total_nr_of_runs     = max(condition_master.global_run_nr);
total_nr_of_blocks   = max(condition_master.global_block_nr);
total_nr_of_trials   = max(condition_master.global_trial_nr);
total_nr_of_taskclasses = max(condition_master.global_trial_nr);

% assert(total_nr_of_runs


% More specific stats
total_nr_of_unique_conditions = length(unique(condition_master.condition_nr(~isnan(condition_master.condition_nr))));


% Stats per run:
[~,run_start] = unique(condition_master.global_run_nr);
run_start = cat(1,run_start,size(condition_master,1));

trials_per_run = zeros(1,total_nr_of_sessions);
for rr = 1:length(run_start)-1
    trials_per_run(rr) = length(condition_master.trial_nr(run_start(rr):(run_start(rr+1)-1)));
end

trials_per_taskclass = zeros(1,total_nr_of_sessions);
for rr = 1:length(run_start)-1
    trials_per_run(rr) = length(condition_master.trial_nr(run_start(rr):(run_start(rr+1)-1)));
end
   

% Tell the user!
fprintf('[%s]: Stats of condition_master:\n',mfilename) 
fprintf('Total number of sessions: \t\t%03d\n', total_nr_of_sessions);
fprintf('Total number of runs: \t\t\t%03d\n', total_nr_of_runs);
fprintf('Total number of block: \t\t\t%03d\n', total_nr_of_blocks);
fprintf('Total number of trials: \t\t%03d\n', total_nr_of_trials);
fprintf('Total number of unique conditions: \t%03d\n', total_nr_of_unique_conditions);
fprintf('Total number of trials for runs 1-%d: \t%s \n', total_nr_of_runs,num2str(trials_per_run))

fprintf('Total number of trials and blocks per task class?', 

