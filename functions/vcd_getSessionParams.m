function exp_session = vcd_getSessionParams(p,load_params,store_params)

if ~exist('store_params','var') || isempty(store_params)
    store_params = true;
end

if ~exist('load_params','var') || isempty(load_params)
    load_params = true;
end

if load_params
    if ~isfield(p,'disp')
        warning('[%s]: Can''t find display name! Will use 7TAS_BOLDSCREEN32 instead', mfilename)
        dispname = '7TAS_BOLDSCREEN32';
    else
        dispname = p.disp.name;
    end
    
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('exp_session_%s*.mat',dispname)));
    if ~isempty(d)
        if length(d) > 1
            warning('[%s]: Multiple .mat files! Will pick the most recent one', mfilename);
        end
        load(fullfile(d(end).folder,d(end).name),'exp_session');
    else
        error('[%s]: Can''t find experiment session params file!', mfilename)
    end
else
    

    
    % Preallocate space
    exp_session = struct('session',[],'run',[],'miniblock',[],'trial', []);
    
    %% Define big stim-task crossing table
    exp_session.stimClassLabels = {'gabor','rdk','dot','cobj','ns'};
    exp_session.taskClassLabels = {'fix','cd','scc','pc','wm','ltm','img','what','where','how'};
    
    exp_session.crossings = false(length(exp_session.stimClassLabels),length(exp_session.taskClassLabels));
    
    exp_session.stimTaskLabels = cell(length(exp_session.taskClassLabels),length(exp_session.stimClassLabels));
    for row = 1:size(exp_session.crossings,1)
        for col = 1:size(exp_session.crossings,2)
            exp_session.stimTaskLabels{col,row} = sprintf('%s-%s',lower(exp_session.taskClassLabels{col}),lower(exp_session.stimClassLabels{row}));
        end
    end
    
    % Set Classic block
    exp_session.crossings(1:5,1:7) = true;
    exp_session.crossings(5,3) = false;
    
    % Set Naturalistic tail
    exp_session.crossings(4:5,8:10) = true;
    exp_session.crossings(4,9) = false;
    
    exp_session.stimTaskLabels = exp_session.stimTaskLabels(exp_session.crossings');
    
    % General exp params
    exp_session.n_unique_trial_repeats = 4;
    exp_session.n_sessions             = 1;%30;
    exp_session.n_runs_per_session     = 10;
    exp_session.TR                     = 1.6; % seconds
    exp_session.total_subjects         = 3;
    
    % %%%% SESSION %%%%
    exp_session.session.wideSessions     = [1,2]; % Sessions dedicated to wide subject sampling
    exp_session.session.baselineSessions = [1:6]; % Sessions dedicated to establish baseline
    % response for LTM/IMAG (=2 repeats of unique image conditions)
    exp_session.session.task_start       = [1,1,1,1,1,7,7,1,1,1]; % When do we start sampling the tasks (LTM/IMG have later starts)
    
    % %%%% MINIBLOCK %%%%
    % general
    exp_session.miniblock.n_trials_single_epoch = 8;
    exp_session.miniblock.n_trials_double_epoch = 4;
    exp_session.miniblock.response_ID    = 93; % 
    exp_session.miniblock.trial_start_ID = 94; % 
    exp_session.miniblock.spatial_cue_ID = 95; % 
    exp_session.miniblock.delay_ID       = 96; % 
    exp_session.miniblock.task_cue_ID    = 97; % 
    exp_session.miniblock.ITI_ID         = 98; %
    exp_session.miniblock.IBI_ID         = 99; %
    
    % check if these ITI/IBI IDs do not already exist
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.miniblock.task_cue_ID)));
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.miniblock.ITI_ID)));
    assert(isempty(intersect([1:length(exp_session.stimTaskLabels)],exp_session.miniblock.IBI_ID)));
    
    % timing
    exp_session.miniblock.task_cue_dur        = p.stim.fps*60; % 2.0 seconds
    exp_session.miniblock.IBI                 = p.stim.fps*linspace(150,270,5); % [5:1:9] seconds Inter-block interval -- uniformly sample between [min,max]
    
    
    % %%%% RUN %%%%
    % general
    exp_session.run.n_single_epoch_miniblocks = 3;
    exp_session.run.n_double_epoch_miniblocks = 3;
    exp_session.run.miniblocks_per_run = exp_session.run.n_single_epoch_miniblocks + exp_session.run.n_double_epoch_miniblocks;
    
    % timing
    exp_session.run.pre_blank_dur     = p.stim.fps*330; % 11 s
    exp_session.run.post_blank_dur    = p.stim.fps*330; % 11s
    exp_session.run.total_run_dur     = p.stim.fps*9936; % 9936 fps or 207 TRs or 331.2 s
    exp_session.run.actual_task_dur = exp_session.run.total_run_dur - exp_session.run.pre_blank_dur - exp_session.run.post_blank_dur; %s

    
    % %%%% TRIAL %%%%
    % general
    exp_session.trial.single_epoch_tasks = logical([1 1 1 1 0 0 0 1 1 1]);
    exp_session.trial.double_epoch_tasks = ~exp_session.trial.single_epoch_tasks;
    exp_session.trial.stim_LR_loc_cue    = logical([1 1 1 1 0]);
    
    % timing
    exp_session.trial.start_cue_dur       = p.stim.fps*12; % 0.4 seconds (thickening of dot rim)
    exp_session.trial.spatial_cue_dur     = p.stim.fps*24; % 0.8 seconds
    exp_session.trial.stim_array_dur      = p.stim.fps*60; % 2.0 seconds
    exp_session.trial.response_win_dur    = p.stim.fps*30; % 1.0 seconds
    % p.trial.end_cue_dur         = 0.4; % seconds  (thinning of dot rim)
    exp_session.trial.ITI                 = p.stim.fps.*[6:6:48]; % 0.2:0.2:1.6 seconds (thinning of dot rim)
    exp_session.trial.delay_dur           = p.stim.fps*240; % 8.0 seconds
    
    exp_session.trial.single_epoch_dur   = ...
        sum([exp_session.trial.start_cue_dur,... % seconds
        exp_session.trial.spatial_cue_dur, ...
        exp_session.trial.stim_array_dur, ...
        exp_session.trial.response_win_dur]);
    
    exp_session.trial.double_epoch_dur   = ...
        sum([exp_session.trial.start_cue_dur,... % seconds
        exp_session.trial.spatial_cue_dur, ...
        exp_session.trial.stim_array_dur, ...
        exp_session.trial.delay_dur, ...
        exp_session.trial.stim_array_dur, ...
        exp_session.trial.response_win_dur]);
    
    assert( nearZero(mod(exp_session.trial.single_epoch_dur / p.stim.fps,1),tolerance))
    
    %%
    
    exp_session.ses_blocks = zeros(size(exp_session.crossings,1),size(exp_session.crossings,2),exp_session.n_sessions);
    
    for ses = 1:exp_session.n_sessions
        
        for ii = 1:size(exp_session.stimTaskLabels)
            
            curr_cross = exp_session.stimTaskLabels{ii};
            
            sc = ~cellfun(@isempty, regexp(curr_cross,exp_session.stimClassLabels, 'match','once'));
            tc = ~cellfun(@isempty, regexp(curr_cross,exp_session.taskClassLabels, 'match','once'));
            
            switch curr_cross
                %%% FIXATION %%%
                case {'fix-gabor','fix-rdk','fix-dot','fix-cobj','fix-ns'}
                    
                    exp_session.ses_blocks(sc,tc,ses) = 1;
                    
                    %%%%%%%%%%%%%%%%%%%%% %%% GABORS %%%
                case 'cd-gabor'
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'scc-gabor'
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'pc-gabor'
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'wm-gabor'
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 3;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'ltm-gabor'
                    if ses >= exp_session.session.task_start(tc)
                        if mod(ses,2)==1 % uneven sessions
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        elseif mod(ses,2)==0 % even sessions
                            exp_session.ses_blocks(sc,tc,ses) = 3;
                        end
                    else % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0; % skip
                    end
                    
                case 'img-gabor'
                    if ses >= exp_session.session.task_start(tc)
                        if mod(ses,2)==1 % uneven sessions
                            exp_session.ses_blocks(sc,tc,ses) = 3;
                        elseif mod(ses,2)==0 % even sessions
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        end
                    else % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0; % skip
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%% %%% RDK %%%
                case 'cd-rdk'
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'scc-rdk'
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'pc-rdk'
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'wm-rdk'
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 3;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                    
                case 'ltm-rdk'
                    if ses >= exp_session.session.task_start(tc)
                        if mod(ses,2)==1 % uneven sessions
                            exp_session.ses_blocks(sc,tc,ses) = 3;
                        elseif mod(ses,2)==0 % even sessions
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        end
                    else  % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0; % skip
                    end
                    
                case 'img-rdk'
                    
                    if ses >= exp_session.session.task_start(tc)
                        if mod(ses,2)==1 % uneven sessions
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        elseif mod(ses,2)==0 % even sessions
                            exp_session.ses_blocks(sc,tc,ses) = 3;
                        end
                    else  % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0; % skip
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%% %%% SIMPLE DOT %%%
                case 'cd-dot'
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        else
                            if mod(ses,2)==1 % uneven sessions
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            elseif mod(ses,2)==0 % even sessions
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            end
                        end
                        
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'scc-dot'
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        else
                            if mod(ses,2)==1 % uneven sessions
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            elseif mod(ses,2)==0 % even sessions
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                        end
                        
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'pc-dot'
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        else
                            if mod(ses,2)==1 % uneven sessions
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            elseif mod(ses,2)==0 % even sessions
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            end
                        end
                        
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'wm-dot'
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'ltm-dot'
                    if ses >= exp_session.session.task_start(tc)
                        doubleSessions = [10:4:exp_session.n_sessions];
                        
                        if any(ses == doubleSessions) % every fourth session has 2 blocks
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        elseif ~any(ses == doubleSessions) % otherwise it is 1
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                    else % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'img-dot'
                    if ses >= exp_session.session.task_start(tc)
                        doubleSessions = [7:4:exp_session.n_sessions];
                        
                        if any(ses == doubleSessions) % every fourth session has 2 blocks
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        elseif ~any(ses == doubleSessions) % otherwise it is 1
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                    else % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                    
                    %%%%%%%%%%%%%%%%%%%%% %%% COMPLEX OBJECTS %%%
                case 'cd-cobj'
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        else
                            if mod(ses,2)==1 % uneven sessions
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            elseif mod(ses,2)==0 % even sessions
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            end
                        end
                        
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'scc-cobj'
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        else
                            if mod(ses,2)==1 % uneven sessions
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            elseif mod(ses,2)==0 % even sessions
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                        end
                        
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'pc-cobj'
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        else
                            if mod(ses,2)==1 % uneven sessions
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            elseif mod(ses,2)==0 % even sessions
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            end
                        end
                        
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                    
                case 'wm-cobj'
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'ltm-cobj'
                    if ses >= exp_session.session.task_start(tc)
                        exp_session.ses_blocks(sc,tc,ses) = 1;
                    else % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'img-cobj'
                    if ses >= exp_session.session.task_start(tc)
                        exp_session.ses_blocks(sc,tc,ses) = 1;
                    else % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'what-cobj'
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        else
                            % First 13 sesions we have 1 miniblock per
                            % session
                            if ses < 13
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            else
                                % after that we alternate between 1/0 for even/uneven
                                if mod(ses,2) == 1 % uneven sessions
                                    exp_session.ses_blocks(sc,tc,ses) = 0;
                                elseif mod(ses,2) == 0 % even sessions
                                    exp_session.ses_blocks(sc,tc,ses) = 1;
                                end
                                
                            end
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'how-cobj'
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        else
                            % First 13 sesions we have 1 miniblock per
                            % session
                            if ses < 13
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            else
                                % after that we alternate between 0/1 for even/uneven
                                if mod(ses,2) == 1 % uneven sessions
                                    exp_session.ses_blocks(sc,tc,ses) = 1;
                                elseif mod(ses,2) == 0 % even sessions
                                    exp_session.ses_blocks(sc,tc,ses) = 0;
                                end
                                
                            end
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%% %%% NATURAL SCENES %%%
                case 'cd-ns'
                    skipSessions = [9:2:exp_session.n_sessions];
                    
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            if any(ses==exp_session.session.baselineSessions)
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                            
                            if any(ses == skipSessions) % every third session skips a block
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            elseif ~any(ses == skipSessions) % otherwise it is 1 block
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'pc-ns'
                    skipSessions = [7:2:exp_session.n_sessions];
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            if any(ses==exp_session.session.baselineSessions)
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                            
                            if any(ses == skipSessions) % every third session skips a block
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            elseif ~any(ses == skipSessions) % otherwise it is 1 block
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'wm-ns'
                    doubleSessions = [9:3:exp_session.n_sessions];
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 4;
                        else
                            if any(ses == exp_session.session.baselineSessions)
                                exp_session.ses_blocks(sc,tc,ses) = 2;
                            end
                            
                            if any(ses == doubleSessions) % every fourth session adds a block
                                exp_session.ses_blocks(sc,tc,ses) = 2;
                            elseif ~any(ses == doubleSessions) % otherwise it is 1 block
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                        end
                    else % shouldn't happen, but for completeness we add this
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'ltm-ns'
                    doubleSessions = [7:10, 12:2:exp_session.n_sessions];
                    
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses == doubleSessions) % every fourth session adds a block
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        elseif ~any(ses == doubleSessions) % otherwise it is 1 block
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                    else % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'img-ns'
                    doubleSessions = [7:11, 13:2:exp_session.n_sessions];
                    if ses >= exp_session.session.task_start(tc)
                        
                        if any(ses == doubleSessions) % every fourth session adds a block
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        elseif ~any(ses == doubleSessions) % otherwise it is 1 block
                            exp_session.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                    else % skip until session 7
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'what-ns'
                    skipSessions = [7:2:exp_session.n_sessions,exp_session.n_sessions];
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            if any(ses == exp_session.session.baselineSessions)
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                            
                            if any(ses == skipSessions) % every third session skips a block
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            elseif ~any(ses == skipSessions) % otherwise it is 1 block
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'where-ns'
                    skipSessions = [7,8:2:exp_session.n_sessions];
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            if any(ses == exp_session.session.baselineSessions)
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                            
                            if any(ses == skipSessions) % every third session skips a block
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            elseif ~any(ses == skipSessions) % otherwise it is 1 block
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
                    
                case 'how-ns'
                    skipSessions = [7:2:exp_session.n_sessions,exp_session.n_sessions];
                    if ses >= exp_session.session.task_start(tc)
                        if any(ses==exp_session.session.wideSessions)
                            exp_session.ses_blocks(sc,tc,ses) = 2;
                        else
                            if any(ses == exp_session.session.baselineSessions)
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                            
                            if any(ses == skipSessions) % every third session skips a block
                                exp_session.ses_blocks(sc,tc,ses) = 0;
                            elseif ~any(ses == skipSessions) % otherwise it is 1 block
                                exp_session.ses_blocks(sc,tc,ses) = 1;
                            end
                        end
                    else % shouldn't occur, but for completeness..
                        exp_session.ses_blocks(sc,tc,ses) = 0;
                    end
            end
            
            
        end
        
    end

    if store_params
        fprintf('[%s]:Storing session data..\n',mfilename)
        saveDir = fileparts(fullfile(p.stim.fix.infofile));
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir, sprintf('exp_session_%s.mat',datestr(now,30))),'exp_session','-v7.3');
    end
end


return



