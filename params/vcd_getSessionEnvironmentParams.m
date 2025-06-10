function [all_sessions,session_types,n_runs_per_session,min_run_dur, ...
          total_run_dur,actual_task_run_dur, IBI, nr_session_types, ...
          preblank_run_dur, postblank_run_dur, unique_trial_repeats, ...
          nr_blocks_per_run,catch_trial_flag] = ...
            vcd_getSessionEnvironmentParams(params, session_env)

%% WRITE ME DESCRIPTION OF OUTPUTS

% all_sessions          : 
% session_types         : 
% n_runs_per_session    : 
% min_run_dur           : 
% total_run_dur         : 
% actual_task_run_dur   :  nr of time frames a subject actually spends
%                           doing the experiment: only stimulus blocks and
%                           IBI duration.                    
% IBI                   : 
% nr_session_types      : 
% preblank_run_dur      : 
% postblank_run_dur     : 
% unique_trial_repeats  : 
% nr_blocks_per_run     : 
% catch_trial_flag      : 


%%
% Get params for session environment (MRI or Behavior)
if strcmp(session_env,'MRI')
    
    % Concatenate wide and deep sessions
    all_sessions         = cat(3, params.exp.session.mri.wide.ses_blocks, params.exp.session.mri.deep.ses_blocks);
    session_types        = cat(1, params.exp.session.mri.wide.session_types, params.exp.session.mri.deep.session_types);
    nr_session_types     = 2;
    
    n_runs_per_session   = cat(1, params.exp.session.mri.wide.n_runs_per_session, params.exp.session.mri.deep.n_runs_per_session);
    nr_blocks_per_run    = cat(1, params.exp.session.mri.wide.nr_blocks_per_run, params.exp.session.mri.deep.nr_blocks_per_run);
    unique_trial_repeats = params.exp.n_unique_trial_repeats_mri;
    catch_trial_flag     = true;
    
    min_run_dur          = params.exp.run.min_run_dur_MRI;
    total_run_dur        = params.exp.run.total_run_dur_MRI;
    actual_task_run_dur  = params.exp.run.actual_task_dur_MRI;
    IBI                  = params.exp.block.IBI_MRI;
    preblank_run_dur     = params.exp.run.pre_blank_dur_MRI;
    postblank_run_dur    = params.exp.run.post_blank_dur_MRI;

elseif strcmp(session_env,'BEHAVIOR')
    all_sessions         = params.exp.session.behavior.ses_blocks;
    session_types        = params.exp.session.behavior.session_types;
    nr_session_types     = 1;
    
    n_runs_per_session   = params.exp.session.behavior.n_runs_per_session;
    nr_blocks_per_run    = params.exp.session.behavior.nr_blocks_per_run;
    unique_trial_repeats = params.exp.n_unique_trial_repeats_behavior;
    catch_trial_flag     = false;
    
    min_run_dur          = params.exp.run.min_run_dur_BEHAVIOR;
    total_run_dur        = params.exp.run.total_run_dur_BEHAVIOR;
    actual_task_run_dur  = params.exp.run.actual_task_dur_BEHAVIOR;
    IBI                  = params.exp.block.IBI_BEHAVIOR;
    preblank_run_dur     = params.exp.run.pre_blank_dur_BEHAVIOR;
    postblank_run_dur    = params.exp.run.post_blank_dur_BEHAVIOR;

else 
    error('[%s]: Session type doesn''t exist! Choose: "MRI" or "BEHAVIOR".',mfilename)
end
    
