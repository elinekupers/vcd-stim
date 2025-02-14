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
p.n_sessions             = 1;%30;
p.n_runs_per_session     = 10;
p.TR                     = 1.6; % seconds
p.total_subjects         = 3;

% %%%% SESSION %%%%
p.session.wideSessions     = [1,2]; % Sessions dedicated to wide subject sampling
p.session.baselineSessions = [1:6]; % Sessions dedicated to establish baseline
% response for LTM/IMAG (=2 repeats of unique image conditions)
p.session.task_start       = [1,1,1,1,1,7,7,1,1,1]; % When do we start sampling the tasks (LTM/IMG have later starts)

% %%%% MINIBLOCK %%%%
% general
p.miniblock.n_trials_single_epoch = 8;
p.miniblock.n_trials_double_epoch = 4;
p.miniblock.task_cue_ID = 97; % pick a random high nr
p.miniblock.ITI_ID      = 98; % 
p.miniblock.IBI_ID      = 99; % 

% check if these ITI/IBI IDs do not already exist
assert(isempty(intersect([1:length(p.stimTaskLabels)],p.miniblock.task_cue_ID)));
assert(isempty(intersect([1:length(p.stimTaskLabels)],p.miniblock.ITI_ID)));
assert(isempty(intersect([1:length(p.stimTaskLabels)],p.miniblock.IBI_ID)));



% timing
p.miniblock.task_cue_dur          = 2.0; % seconds
p.miniblock.IBI                 = [5,9]; % seconds Inter-block interval -- uniformly sample between [min,max]


% %%%% RUN %%%%
% general
p.run.n_single_epoch_miniblocks = 3;
p.run.n_double_epoch_miniblocks = 3;
p.run.miniblocks_per_run = p.run.n_single_epoch_miniblocks + p.run.n_double_epoch_miniblocks;

% timing
p.run.pre_blank_dur     = 11; % s
p.run.post_blank_dur    = 11; % s
p.run.total_run_dur     = 60*5.5; % s
p.run.actual_task_dur = p.run.total_run_dur - p.run.pre_blank_dur - p.run.post_blank_dur; %s



% %%%% TRIAL %%%%
% general
p.trial.single_epoch_tasks = logical([1 1 1 1 0 0 0 1 1 1]);
p.trial.double_epoch_tasks = ~p.trial.single_epoch_tasks;
p.trial.stim_LR_loc_cue    = logical([1 1 1 1 0]);

% timing
p.trial.start_cue_dur       = 0.4; % seconds (thickening of dot rim)
p.trial.spatial_cue_dur     = 0.8; % seconds
p.trial.stim_array_dur      = 2.0; % seconds
p.trial.response_win_dur    = 1.0; % seconds
% p.trial.end_cue_dur         = 0.4; % seconds  (thinning of dot rim)
p.trial.ITI                 = [0.2:0.2:1.6]; % seconds (thinning of dot rim)
p.trial.delay_dur           = 8.0; % seconds

p.trial.single_epoch_dur   = ...
    sum([p.trial.start_cue_dur,... % seconds
        p.trial.spatial_cue_dur, ...
        p.trial.stim_array_dur, ...
        p.trial.response_win_dur]);

p.trial.double_epoch_dur   = ...
    sum([p.trial.start_cue_dur,... % seconds
        p.trial.spatial_cue_dur, ...
        p.trial.stim_array_dur, ...
        p.trial.delay_dur, ...
        p.trial.stim_array_dur, ...
        p.trial.response_win_dur]);



%%

p.ses_blocks = zeros(size(p.crossings,1),size(p.crossings,2),p.n_sessions);

for ses = 1:p.n_sessions
    
    for ii = 1:size(p.stimTaskLabels)
        
        curr_cross = p.stimTaskLabels{ii};
        
        sc = ~cellfun(@isempty, regexp(curr_cross,p.stimClassLabels, 'match','once'));
        tc = ~cellfun(@isempty, regexp(curr_cross,p.taskClassLabels, 'match','once'));
        
        switch curr_cross
            %%% FIXATION %%%
            case {'fix-gabor','fix-rdk','fix-dot','fix-cobj','fix-ns'}
                
                p.ses_blocks(sc,tc,ses) = 1;
                
                %%%%%%%%%%%%%%%%%%%%% %%% GABORS %%%
            case 'cd-gabor'
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'scc-gabor'
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'pc-gabor'
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'wm-gabor'
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 3;
                    else
                        p.ses_blocks(sc,tc,ses) = 2;
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'ltm-gabor'
                if ses >= p.session.task_start(tc)
                    if mod(ses,2)==1 % uneven sessions
                        p.ses_blocks(sc,tc,ses) = 2;
                    elseif mod(ses,2)==0 % even sessions
                        p.ses_blocks(sc,tc,ses) = 3;
                    end
                else % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0; % skip
                end
                
            case 'img-gabor'
                if ses >= p.session.task_start(tc)
                    if mod(ses,2)==1 % uneven sessions
                        p.ses_blocks(sc,tc,ses) = 3;
                    elseif mod(ses,2)==0 % even sessions
                        p.ses_blocks(sc,tc,ses) = 2;
                    end
                else % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0; % skip
                end
                
                %%%%%%%%%%%%%%%%%%%%% %%% RDK %%%
            case 'cd-rdk'
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'scc-rdk'
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'pc-rdk'
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'wm-rdk'
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 3;
                    else
                        p.ses_blocks(sc,tc,ses) = 2;
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
                
            case 'ltm-rdk'
                if ses >= p.session.task_start(tc)
                    if mod(ses,2)==1 % uneven sessions
                        p.ses_blocks(sc,tc,ses) = 3;
                    elseif mod(ses,2)==0 % even sessions
                        p.ses_blocks(sc,tc,ses) = 2;
                    end
                else  % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0; % skip
                end
                
            case 'img-rdk'
                
                if ses >= p.session.task_start(tc)
                    if mod(ses,2)==1 % uneven sessions
                        p.ses_blocks(sc,tc,ses) = 2;
                    elseif mod(ses,2)==0 % even sessions
                        p.ses_blocks(sc,tc,ses) = 3;
                    end
                else  % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0; % skip
                end
                
                %%%%%%%%%%%%%%%%%%%%% %%% SIMPLE DOT %%%
            case 'cd-dot'
                if ses >= p.session.task_start(tc)
                    
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 1;
                    else
                        if mod(ses,2)==1 % uneven sessions
                            p.ses_blocks(sc,tc,ses) = 1;
                        elseif mod(ses,2)==0 % even sessions
                            p.ses_blocks(sc,tc,ses) = 0;
                        end
                    end
                    
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'scc-dot'
                if ses >= p.session.task_start(tc)
                    
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 1;
                    else
                        if mod(ses,2)==1 % uneven sessions
                            p.ses_blocks(sc,tc,ses) = 0;
                        elseif mod(ses,2)==0 % even sessions
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                    end
                    
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'pc-dot'
                if ses >= p.session.task_start(tc)
                    
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 1;
                    else
                        if mod(ses,2)==1 % uneven sessions
                            p.ses_blocks(sc,tc,ses) = 1;
                        elseif mod(ses,2)==0 % even sessions
                            p.ses_blocks(sc,tc,ses) = 0;
                        end
                    end
                    
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'wm-dot'
                if ses >= p.session.task_start(tc)
                    
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                    
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'ltm-dot'
                if ses >= p.session.task_start(tc)
                    doubleSessions = [10:4:p.n_sessions];
                    
                    if any(ses == doubleSessions) % every fourth session has 2 blocks
                        p.ses_blocks(sc,tc,ses) = 2;
                    elseif ~any(ses == doubleSessions) % otherwise it is 1
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                    
                else % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'img-dot'
                if ses >= p.session.task_start(tc)
                    doubleSessions = [7:4:p.n_sessions];
                    
                    if any(ses == doubleSessions) % every fourth session has 2 blocks
                        p.ses_blocks(sc,tc,ses) = 2;
                    elseif ~any(ses == doubleSessions) % otherwise it is 1
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                    
                else % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
                
                %%%%%%%%%%%%%%%%%%%%% %%% COMPLEX OBJECTS %%%
            case 'cd-cobj'
                if ses >= p.session.task_start(tc)
                    
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 1;
                    else
                        if mod(ses,2)==1 % uneven sessions
                            p.ses_blocks(sc,tc,ses) = 1;
                        elseif mod(ses,2)==0 % even sessions
                            p.ses_blocks(sc,tc,ses) = 0;
                        end
                    end
                    
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'scc-cobj'
                if ses >= p.session.task_start(tc)
                    
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 1;
                    else
                        if mod(ses,2)==1 % uneven sessions
                            p.ses_blocks(sc,tc,ses) = 0;
                        elseif mod(ses,2)==0 % even sessions
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                    end
                    
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'pc-cobj'
                if ses >= p.session.task_start(tc)
                    
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 1;
                    else
                        if mod(ses,2)==1 % uneven sessions
                            p.ses_blocks(sc,tc,ses) = 1;
                        elseif mod(ses,2)==0 % even sessions
                            p.ses_blocks(sc,tc,ses) = 0;
                        end
                    end
                    
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
                
            case 'wm-cobj'
                if ses >= p.session.task_start(tc)
                    
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                    
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'ltm-cobj'
                if ses >= p.session.task_start(tc)
                    p.ses_blocks(sc,tc,ses) = 1;
                else % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'img-cobj'
                if ses >= p.session.task_start(tc)
                    p.ses_blocks(sc,tc,ses) = 1;
                else % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'what-cobj'
                if ses >= p.session.task_start(tc)
                    
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 1;
                    else
                        % First 13 sesions we have 1 miniblock per
                        % session
                        if ses < 13
                            p.ses_blocks(sc,tc,ses) = 1;
                        else
                            % after that we alternate between 1/0 for even/uneven
                            if mod(ses,2) == 1 % uneven sessions
                                p.ses_blocks(sc,tc,ses) = 0;
                            elseif mod(ses,2) == 0 % even sessions
                                p.ses_blocks(sc,tc,ses) = 1;
                            end
                            
                        end
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'how-cobj'
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 1;
                    else
                        % First 13 sesions we have 1 miniblock per
                        % session
                        if ses < 13
                            p.ses_blocks(sc,tc,ses) = 1;
                        else
                            % after that we alternate between 0/1 for even/uneven
                            if mod(ses,2) == 1 % uneven sessions
                                p.ses_blocks(sc,tc,ses) = 1;
                            elseif mod(ses,2) == 0 % even sessions
                                p.ses_blocks(sc,tc,ses) = 0;
                            end
                            
                        end
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
                %%%%%%%%%%%%%%%%%%%%% %%% NATURAL SCENES %%%
            case 'cd-ns'
                skipSessions = [9:2:p.n_sessions];
                
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        if any(ses==p.session.baselineSessions)
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                        if any(ses == skipSessions) % every third session skips a block
                            p.ses_blocks(sc,tc,ses) = 0;
                        elseif ~any(ses == skipSessions) % otherwise it is 1 block
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'pc-ns'
                skipSessions = [7:2:p.n_sessions];
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        if any(ses==p.session.baselineSessions)
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                        if any(ses == skipSessions) % every third session skips a block
                            p.ses_blocks(sc,tc,ses) = 0;
                        elseif ~any(ses == skipSessions) % otherwise it is 1 block
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'wm-ns'
                doubleSessions = [9:3:p.n_sessions];
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 4;
                    else
                        if any(ses == p.session.baselineSessions)
                            p.ses_blocks(sc,tc,ses) = 2;
                        end
                        
                        if any(ses == doubleSessions) % every fourth session adds a block
                            p.ses_blocks(sc,tc,ses) = 2;
                        elseif ~any(ses == doubleSessions) % otherwise it is 1 block
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                    end
                else % shouldn't happen, but for completeness we add this
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'ltm-ns'
                doubleSessions = [7:10, 12:2:p.n_sessions];
                
                if ses >= p.session.task_start(tc)
                    
                    if any(ses == doubleSessions) % every fourth session adds a block
                        p.ses_blocks(sc,tc,ses) = 2;
                    elseif ~any(ses == doubleSessions) % otherwise it is 1 block
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                    
                else % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'img-ns'
                doubleSessions = [7:11, 13:2:p.n_sessions];
                if ses >= p.session.task_start(tc)
                    
                    if any(ses == doubleSessions) % every fourth session adds a block
                        p.ses_blocks(sc,tc,ses) = 2;
                    elseif ~any(ses == doubleSessions) % otherwise it is 1 block
                        p.ses_blocks(sc,tc,ses) = 1;
                    end
                    
                else % skip until session 7
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'what-ns'
                skipSessions = [7:2:p.n_sessions,p.n_sessions];
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        if any(ses == p.session.baselineSessions)
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                        if any(ses == skipSessions) % every third session skips a block
                            p.ses_blocks(sc,tc,ses) = 0;
                        elseif ~any(ses == skipSessions) % otherwise it is 1 block
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'where-ns'
                skipSessions = [7,8:2:p.n_sessions];
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        if any(ses == p.session.baselineSessions)
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                        if any(ses == skipSessions) % every third session skips a block
                            p.ses_blocks(sc,tc,ses) = 0;
                        elseif ~any(ses == skipSessions) % otherwise it is 1 block
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
                
            case 'how-ns'
                skipSessions = [7:2:p.n_sessions,p.n_sessions];
                if ses >= p.session.task_start(tc)
                    if any(ses==p.session.wideSessions)
                        p.ses_blocks(sc,tc,ses) = 2;
                    else
                        if any(ses == p.session.baselineSessions)
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                        
                        if any(ses == skipSessions) % every third session skips a block
                            p.ses_blocks(sc,tc,ses) = 0;
                        elseif ~any(ses == skipSessions) % otherwise it is 1 block
                            p.ses_blocks(sc,tc,ses) = 1;
                        end
                    end
                else % shouldn't occur, but for completeness..
                    p.ses_blocks(sc,tc,ses) = 0;
                end
        end
        
        
    end
    
end



return



