function [all_sessions,session_types,n_runs_per_session,min_run_dur, ...
          total_run_dur,actual_task_run_dur, IBI, nr_session_types, ...
          preblank_run_dur, postblank_run_dur] = ...
            vcd_getSessionEnvironmentParams(params, session_env)

% actual_task_run_dur   :  nr of presentation frames we actually spend doing the experiment

% Get params for session environment (MRI or Behavior)
if strcmp(session_env,'MRI')
    
    % Concatenate wide and deep sessions
    all_sessions         = cat(3, params.exp.session.wide.ses_blocks, params.exp.session.deep.ses_blocks);
    session_types        = cat(1, params.exp.session.mri.wide.session_types, params.exp.session.mri.deep.session_types);
    n_runs_per_session   = cat(1, params.exp.session.mri.wide.n_runs_per_session, params.exp.session.mri.deep.n_runs_per_session);
    
    min_run_dur          = params.exp.run.min_run_dur_MRI;
    total_run_dur        = params.exp.run.total_run_dur_MRI;
    actual_task_run_dur  = params.exp.run.actual_task_dur_MRI;
    IBI                  = params.exp.block.IBI_MRI;
    nr_session_types     = 2;
    preblank_run_dur     = params.exp.run.pre_blank_dur_MRI;
    postblank_run_dur    = params.exp.run.post_blank_dur_MRI;
    
    
elseif strcmp(session_env,'BEHAVIOR')
    all_sessions        = params.exp.session.behavior.ses_blocks;
    session_types       = params.exp.session.behavior.session_types;
    n_runs_per_session  = params.exp.session.behavior.n_runs_per_session;
    min_run_dur         = params.exp.run.min_run_dur_BEHAVIOR;
    total_run_dur       = params.exp.run.total_run_dur_BEHAVIOR;
    actual_task_run_dur = params.exp.run.actual_task_dur_BEHAVIOR;
    IBI                 = params.exp.block.IBI_BEHAVIOR;
    nr_session_types    = 1;
    preblank_run_dur    = params.exp.run.pre_blank_dur_BEHAVIOR;
    postblank_run_dur   = params.exp.run.post_blank_dur_BEHAVIOR;
end



