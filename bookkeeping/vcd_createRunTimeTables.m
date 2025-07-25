function [time_table_master, all_run_frames] = vcd_createRunTimeTables(params, varargin)
% VCD function to add trial events to condition master, turning it into a
% time table. This function also cleans some table columns that are now obsolete.
%
%   [time_table_master, run_frames] = vcd_createRunTimeTables(params,condition_master, env_type, varargin)
%
% INPUTS:
%  params                : 	(struct) parameter struct needed to get subject
%                           nrs and run type params.
%  [condition_master]    :  (table) condition master table with subject
%                           session,runs,blocks,trials, but no individual
%                           trial events.
%  [env_type]         :  (char) Are we dealing with MRI or BEHAVIOR sessions
%                           Default: ''
%  [saveDir]             :  (char; optional) Where do we store
%                           time_table_master? If saveDir is not defined or 
%                           left empty, we will store it under
%                           fullfile(vcd_rootPath,'data',[subj_id])
%  [subj_id]              :  (char) fixed string added to filename when
%                           storing the time_table_master file. If left 
%                           empty, subj_id will be set to 'subj000'.
%
% OUTPUTS:
%   time_table_master   :  (struct) updated condition master table with single trial events.
%     {'session_nr'          } : (double) session number, ranging from
%                                 1-27 for MRI and there is only 1
%                                 behavioral session. Every row of run_nr
%                                 is defined and there are no NaNs in this
%                                 column.
%     {'session_type'        } : (double) session type, can be either 1 for
%                                 version "A" or 2 for version "B". Session
%                                 type is only relevant for the first
%                                 (wide) and the last (deep) MRI session,
%                                 where the stimulus blocks are different
%                                 for the two session versions. All other
%                                 sessions only have one session_type and
%                                 are labeled as 1. Every row of run_nr is
%                                 defined and there are no NaNs in this
%                                 column.
%     {'run_nr'              } : (double) local run number within a session
%                                 (resets every session), ranging from 1-12
%                                 for behavioral experiment, and 1-10 for
%                                 MRI experiment. Every row of run_nr is
%                                 defined and there are no NaNs in this
%                                 column.
%     {'block_nr'            } : (double) local block number within a run
%                                 (resets every run), ranging from 1-7 or
%                                 999 or NaN, where: - 1-7 is used to
%                                 indicate a regular stimulus block. - 999
%                                 is used to indicate an eyetracking block.
%                                 NaNs are used for events outside
%                                 eyetracking of stimulus blocks: pre- or
%                                 post-run blank rest periods and
%                                 inter-block intervals (IBIs).
%     {'trial_nr'            } : (double) local trial number within a block
%                                 (resets every block), ranging from 1-4 or
%                                 1-8 depending if it is a single- or
%                                 double-stimulus presentation block. This
%                                 column contains NaN for events outside
%                                 VCD trials, including for the eyetracking
%                                 block, pre- or post-run blank rest
%                                 periods, inter-block intervals (IBIs),
%                                 and inter-trial intervals (ITIs).
%     {'global_run_nr'       } : (double) global run number across all sessions
%                                 for a subject. Numbers continue counting across sessions.
%                                 Note that if a session has two session types, we the runs separately
%                                 for each session type. For example, global run nrs for session 1A 
%                                 count 1-10 and global run nrs for session 1B count 1-10.
%     {'global_block_nr'     } : (double) global block number across all sessions for a subject
%     {'global_trial_nr'     } : (double) global trial number across all sessions for a subject.
%     {'condition_name'      } : (char) condition name see vcd_getConditionNames.
%     {'condition_nr'        } : (double) condition number (ranging from 1-1331) associated with condition name
%     {'stim_class'          } : (double) stimulus class number [1-5]: GBR, RDK, DOT, OBJ, NS
%     {'task_class'          } : (double) task class number [1-10]: FIX, CD, SCC, PC, WM, LTM, IMG, WHAT, WHERE, HOW
%     {'stim_nr_left'        } : (double) unique stimulus number for left spatial stimulus location or center location in case of NS
%     {'stim_nr_right'       } : (double) unique stimulus number for left spatial stimulus location, left empty in case of NS
%     {'is_cued'             } : (double) if stimulus is cued left (1), right (2), or neutral (3)
%     {'is_catch'            } : (double) if stimulus is a catch trial (1) or not (NaN)
%     {'is_objectcatch'      } : (double) if a PC-OBJ trial has a catch rotation (1) or not (NaN)
%     {'trial_type'          } : (double) if it is a single-stimulus or double-stimulus presentation trial
%     {'crossing_nr'         } : (double) stimulus-task crossing (1:32), see params.exp.stimtaskcrossings
%     {'event_start'         } : (double) event onset in 33 ms frames
%     {'event_dur'           } : (double) event duration in 33 ms frames
%     {'event_end'           } : (double) event end in 33 ms frames
%     {'event_id'            } : (double) event ID (e.g., 91 for stimulus)
%     {'event_name'          } : (cell w str) same as event ID but human readable
%     {'stim_class_name'     } : (cell w str) same as stim_class but human readable
%     {'task_class_name'     } : (cell w str) same as task_class but human readable
%     {'orient_dir'          } : (double) Gabor tilt orientation, RDK motion direction, dot spatial location angle, object facing direction in degrees
%     {'contrast'            } : (double) stimulus michelson contrast (fraction), 1 = 100%
%     {'gbr_phase'           } : (double) Gabor stimulus phase (in degrees) or NaN for other stimulus classes
%     {'rdk_coherence'       } : (double) RDK stimulus coherence (in fraction of dots) or NaN for other stimulus classes
%     {'super_cat'           } : (double) Object/Scene superordinate semantic category or NaN for other stimulus classes
%     {'basic_cat'           } : (double) Object/Scene basic semantic category or NaN for other stimulus classes
%     {'sub_cat'             } : (double) Object/Scene subordinate semantic category or NaN for other stimulus classes
%     {'super_cat_name'      } : (cell w str) same as super_cat but human readable
%     {'basic_cat_name'      } : (cell w str) same as basic_cat but human readable
%     {'sub_cat_name'        } : (cell w str) same as sub_cat but human readable
%     {'is_in_img'           }
%     {'stim2_delta'         } : (double) WM or IMG quiz image change. For
%                                 example, for WM it is the difference from reference for gabor degrees tilt, rdk motion direction,
%                                 dot angle, object facing direction, or
%                                 type of change in scene.
%     {'stim2_im_nr'         } : (cell w str or double) outcome of 'stim2_delta'. For
%                                 example, for -8 + 18 = 10 degrees gabor
%                                 tilt for quiz image after delay.
%     {'ltm_stim_pair'       } : (double) Long term memory paired stimulus
%     {'is_lure'             } : (double) if LTM quiz image is lure (1) or not (0)
%     {'repeat_nr'           } : (double) how often this particular stimulus-task crossing has been repeated across sessions for a subject

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'           , @isstruct);
p0.addParameter('condition_master', [], @istable);
p0.addParameter('load_params'     , false, @islogical);
p0.addParameter('store_params'    , true, @islogical);
p0.addParameter('store_imgs'      , false, @islogical);
p0.addParameter('verbose'         , false, @islogical);
p0.addParameter('env_type'        , '', @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));
p0.addParameter('saveDir'         , fullfile(vcd_rootPath,'data'), @ischar);
p0.addParameter('subj_id'         , 'vcd_subj000', @ischar);

% Parse inputs
p0.parse(params,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

% Infer environment if we haven't set this parameter  
if ~exist('env_type','var') || isempty(env_type)
    if strcmp(params.disp.name,'7TAS_BOLDSCREEN32')
        env_type = 'MRI';
    elseif ismember(params.disp.name,{'CCNYU_VIEWPIXX3D','PPROOM_EIZOFLEXSCAN'}) 
        env_type = 'BEHAVIOR';
    elseif ismember(params.disp.name,{'KKOFFICE_AOCQ3277','EKHOME_ASUSVE247'}) 
        env_type = 'BEHAVIOR';
    end
end

%% Load params if requested and we can find the file
if load_params
    if ~isfield(params,'is_demo'), params.is_demo = false; end
    d = dir(fullfile(vcd_rootPath,'data', env_type, subj_id, sprintf('%s_time_table_master_%s%s*.mat',subj_id,choose(params.is_demo,'demo_',''),params.disp.name)));
    if isempty(d)
        error('[%s]: Can''t find time table file with subj_id: %s!\n', mfilename, subj_id)
    elseif ~isempty(d(end).name)
        if length(d) > 1
            warning('[%s]: Multiple trial .mat files! Will pick the most recent one.\n', mfilename);
        end
        load(fullfile(d(end).folder,d(end).name),'time_table_master');
    end
    
else
    % Create subject sessions
    if ~isfield(params,'is_demo'), params.is_demo = false; end
    
    % check if trial struct is already defined and load it if needed
    if ~exist('condition_master','var') || isempty(condition_master)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','data',env_type, subj_id, sprintf('%s_condition_master_%s*%s*.mat',subj_id, choose(params.is_demo,'demo_',''),params.disp.name)));
        
        if ~isempty(d(end).name)
            if length(d) > 1
                warning('[%s]: Multiple condition_master .mat files! Will pick the most recent one', mfilename);
            end
            a = load(fullfile(d(end).folder,d(end).name),'condition_master_shuffled');
            condition_master = a.condition_master_shuffled; clear a;
        else
            error('[%s]: Can''t find condition_master .mat files! Please check! You may have to run vcd_createConditions.m and/or vcd_createSessions.m', mfilename);
        end
    end
    
    % check if experimental parameters are already defined and load it if needed
    if ~isfield(params,'exp') || isempty(params.exp)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info','exp*.mat'));
        
        if length(d) == 1
            load(fullfile(d(end).folder,d(end).name),'exp');
            params.exp = exp;
        elseif length(d) > 1
            warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            load(fullfile(d(end).folder,d(end).name),'exp');
            params.exp = exp; clear exp;
        else
            warning('[%s]: Can''t find exp session .mat files! Will run vcd_getSessionParams.m', mfilename);
            params.exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
        end
    end
    
    % check if stimulus parameters are already defined and load it if needed
    if ~isfield(params,'stim') || isempty(params.stim)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info','stim*.mat'));
        
        if ~isempty(d)
            if length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            end
            load(fullfile(d(end).folder,d(end).name),'stim');
            params.stim = stim;
        else
            warning('[%s]: Can''t find stim session .mat files! Will run vcd_getStimParams.m', mfilename);
            params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);
        end
    end
end

%% Copy condition master input
condition_master0 = condition_master;

% check inputs
if (nargin ==2 && ~exist(env_type,'var')) || isempty(env_type)
    env_type = 'MRI';
else
    [all_sessions,all_sessiontypes,runs_per_session, ~, ...
          session_totalrundur,actual_task_run_dur, IBI_to_use, ~, ...
          session_preblankdur, session_postblankdur,~,~,~,block_durs] = ...
            vcd_getSessionEnvironmentParams(params, env_type);
end
if ~isfield(params,'is_demo'), params.is_demo = false; end

% Define event ID within trial for single and double epoch trial types
trial_ID_single_epoch = [params.exp.block.spatial_cue_ID, ...
    params.exp.block.pre_stim_blank_ID, ...
    params.exp.block.stim_epoch1_ID, ...
    params.exp.block.response_ID];
trial_ID_double_epoch = [params.exp.block.spatial_cue_ID, ...
    params.exp.block.pre_stim_blank_ID,...
    params.exp.block.stim_epoch1_ID, ...
    params.exp.block.delay_ID, ...
    params.exp.block.stim_epoch2_ID, ...
    params.exp.block.response_ID];

% pick an arbitrary large number of rows for new table
% (otherwise table uses 0 for missing values and we don't want that)
tbl_nrows = 1000; 

% Loop over sessions
for ses = 1:size(all_sessions,3)
    
    % Loop over session types (version A or B)
    for st = 1:size(all_sessions,4)
        
        % Check if the selected session type exists
        if ~isnan(all_sessiontypes(ses,st))
            
            fprintf('\nSESSION %03d %s..\n',ses,choose(st==1,'A','B'))
            
            clear run_time_table

            % check if the number of runs in table matches with how many runs we expect..
            curr_run_nrs = unique(condition_master.run_nr(~isnan(condition_master.run_nr) & condition_master.session_nr==ses & condition_master.session_type==st));
            assert(runs_per_session(ses,st)==length(curr_run_nrs)); 
            
            % Loop over runs
            for rr = curr_run_nrs'
                
                % copy condition order table headers and scrub content
                sz = [tbl_nrows size(condition_master(1,:),2)];
                time_table = vcd_preallocateNaNTable(sz(1), sz(2), condition_master(1,:), []);
                
                % add columns for event timing and id
                time_table.event_start = NaN(tbl_nrows,1);
                time_table.event_dur   = NaN(tbl_nrows,1);
                time_table.event_end   = NaN(tbl_nrows,1);
                time_table.event_id    = NaN(tbl_nrows,1);
                time_table.event_name  = repmat({''},tbl_nrows,1);
                time_table.nr_fix_changes = NaN(tbl_nrows,1);
                
                % Find corresponding trials (rows in table) and block nrs for
                % this particular subject, run, and session. The t_trial table
                % will be inserted into a bigger time table that includes other
                % event types, like iti, cues, etc.
                idx0      = (condition_master.session_nr==ses & ...
                             condition_master.session_type==st & ...
                             condition_master.run_nr==rr);
                assert(sum(idx0)>0)
                t_trial   = condition_master(idx0,:); % trials in this run
                block_nrs = t_trial.block_nr; % get all the block numbers
                
                % Reorder blocks and trials according to subject's specific 
                % block and trial order
                [block_nrs,block_nrs_idx] = sort(block_nrs,'ascend');
                [~,nr_unique_blocks]      = unique(block_nrs,'stable');
                t_trial                   = t_trial(block_nrs_idx,:);

                % Get block durations
                block_dur = block_durs(t_trial.trial_type(nr_unique_blocks));
                
                % reset frame counter and finished run boolean flag
                total_run_frames = 0;
                run_finished     = 0;
                
                if ~isempty(block_nrs)
                    % remove frames such that duration in seconds is an integer (we will add those later)
                    rounded_session_totalrundur = (floor(session_totalrundur/params.stim.presentationrate_hz)*params.stim.presentationrate_hz);
                    total_ses_dur = rounded_session_totalrundur - params.exp.block.total_eyetracking_block_dur - session_postblankdur - session_preblankdur;
                    assert(isequal(actual_task_run_dur,total_ses_dur));
                    
                    % predefine IBIs, make sure we don't go over the total
                    % run duration we want
                    while 1
                        % vcd_optimizeIBIs inputs are run_dur, block_dur, ibis, nr_blocks, prepost_blank
                        [ibis, postblank_to_add] = vcd_optimizeIBIs(total_ses_dur, block_dur, IBI_to_use, length(block_dur), 0, verbose); % prepost_blank is 0 because we already subtract this when defining total_ses_dur
                        if (total_ses_dur - sum(ibis) - sum(block_dur) - postblank_to_add) <= 0  
                            break;
                        end
                    end
                    
                    % for debug purposes
                    % fprintf('[%s]: Selected IBIs (in time frames): %s \n', mfilename, num2str(ibis))
                end
                
                % Add eyetracking block and pre-blank period
                if total_run_frames == 0

                    table_idx = 1;
                    % Add prefixation period (1s)
                    time_table.session_nr(table_idx)       = ses;
                    time_table.session_type(table_idx)     = st;
                    time_table.run_nr(table_idx)           = rr;
                    time_table.block_nr(table_idx)         = 999;
                    time_table.crossing_nr(table_idx)      = NaN;
                    time_table.stim_class(table_idx)       = NaN;
                    time_table.task_class(table_idx)       = NaN;
                    time_table.repeat_nr(table_idx)        = NaN;
                    time_table.trial_nr(table_idx)         = NaN;
                    time_table.global_run_nr(table_idx)    = t_trial.global_run_nr(1);
                    time_table.global_trial_nr(table_idx)  = NaN; 
                    time_table.global_block_nr(table_idx)  = NaN; 
                    time_table.correct_response(table_idx) = NaN;
                    time_table.is_objectcatch(table_idx)   = NaN;
                    time_table.event_start(table_idx)      = total_run_frames;
                    time_table.event_dur(table_idx)        = params.exp.block.eye_gaze_fix0;
                    time_table.event_end(table_idx)        = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)         = params.exp.block.eye_gaze_fix_ID;
                    time_table.event_name(table_idx)       = {'et_fix'};
                    total_run_frames = time_table.event_end(table_idx);
                    
                    % update time table idx
                    table_idx = table_idx+1;
                    
                    % Add saccade events (1.2s each)
                    sac_IDs = NaN(1,params.exp.block.nr_of_saccades);
                    insert_fix_target = randi([2,params.exp.block.nr_of_saccades-1],1);
                    sac_IDs(insert_fix_target) = params.exp.block.eye_gaze_sac_target_ID(1);
                    sac_IDs(isnan(sac_IDs)) = shuffle_concat(params.exp.block.eye_gaze_sac_target_ID(2:end),1);
                    
                    for sac = 1:params.exp.block.nr_of_saccades
                        
                        time_table.session_nr(table_idx)        = ses;
                        time_table.session_type(table_idx)      = st;
                        time_table.run_nr(table_idx)            = rr;
                        time_table.block_nr(table_idx)          = 999;
                        time_table.global_run_nr(table_idx)     = t_trial.global_run_nr(1);
                        time_table.global_trial_nr(table_idx)   = NaN; 
                        time_table.global_block_nr(table_idx)   = NaN;
                        time_table.correct_response(table_idx)  = NaN;
                        time_table.is_objectcatch(table_idx)    = NaN;
                        time_table.crossing_nr(table_idx)       = NaN;
                        time_table.trial_nr(table_idx)          = NaN;
                        time_table.stim_class(table_idx)        = NaN;
                        time_table.task_class(table_idx)        = NaN;
                        time_table.event_start(table_idx)       = total_run_frames;
                        time_table.event_dur(table_idx)         = params.exp.block.eye_gaze_sac_target;
                        time_table.event_end(table_idx)         = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                        time_table.event_id(table_idx)          = sac_IDs(sac);
                        time_table.event_name(table_idx)        = {'et_sac'};
                        total_run_frames = time_table.event_end(table_idx);
                        
                        % update time table idx
                        table_idx = table_idx+1;
                    end
                    
                    % Add post saccade fixation period (1s)
                    time_table.session_nr(table_idx)            = ses;
                    time_table.session_type(table_idx)          = st;
                    time_table.run_nr(table_idx)                = rr;
                    time_table.crossing_nr(table_idx)           = NaN;
                    time_table.global_run_nr(table_idx)         = t_trial.global_run_nr(1);
                    time_table.global_trial_nr(table_idx)       = NaN; 
                    time_table.global_block_nr(table_idx)       = NaN;
                    time_table.correct_response(table_idx)      = NaN;
                    time_table.is_objectcatch(table_idx)        = NaN;
                    time_table.block_nr(table_idx)              = 999;
                    time_table.trial_nr(table_idx)              = NaN;
                    time_table.stim_class(table_idx)            = NaN;
                    time_table.task_class(table_idx)            = NaN;
                    time_table.event_start(table_idx)           = total_run_frames;
                    time_table.event_dur(table_idx)             = params.exp.block.eye_gaze_fix1;
                    time_table.event_end(table_idx)             = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)              = params.exp.block.eye_gaze_fix_ID;
                    time_table.event_name(table_idx)            = {'et_fix'};
                    
                    total_run_frames = time_table.event_end(table_idx);

                    
                    % update time table idx
                    table_idx = table_idx+1;
                    
                    % Add pupil trial (3s black)
                    time_table.session_nr(table_idx)            = ses;
                    time_table.session_type(table_idx)          = st;
                    time_table.run_nr(table_idx)                = rr;
                    time_table.crossing_nr(table_idx)           = NaN;
                    time_table.global_run_nr(table_idx)         = t_trial.global_run_nr(1);
                    time_table.global_trial_nr(table_idx)       = NaN;
                    time_table.global_block_nr(table_idx)       = NaN;
                    time_table.correct_response(table_idx)      = NaN;
                    time_table.is_objectcatch(table_idx)        = NaN;
                    time_table.block_nr(table_idx)              = 999;
                    time_table.trial_nr(table_idx)              = NaN;
                    time_table.stim_class(table_idx)            = NaN;
                    time_table.task_class(table_idx)            = NaN;
                    time_table.event_start(table_idx)           = total_run_frames;
                    time_table.event_dur(table_idx)             = params.exp.block.eye_gaze_pupil_black;
                    time_table.event_end(table_idx)             = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)              = params.exp.block.eye_gaze_pupil_black_ID;
                    time_table.event_name(table_idx)            = {'et_pupil_black'};
                    
                    total_run_frames = time_table.event_end(table_idx);
                    
                    % update time table idx
                    table_idx = table_idx+1;
                    

                    
                    % Add pupil trial (1s white)
                    time_table.session_nr(table_idx)            = ses;
                    time_table.session_type(table_idx)          = st;
                    time_table.run_nr(table_idx)                = rr;
                    time_table.crossing_nr(table_idx)           = NaN;
                    time_table.global_run_nr(table_idx)         = t_trial.global_run_nr(1);
                    time_table.global_trial_nr(table_idx)       = NaN; 
                    time_table.global_block_nr(table_idx)       = NaN;
                    time_table.correct_response(table_idx)      = NaN;
                    time_table.is_objectcatch(table_idx)        = NaN;
                    time_table.block_nr(table_idx)              = 999;
                    time_table.trial_nr(table_idx)              = NaN;
                    time_table.stim_class(table_idx)            = NaN;
                    time_table.task_class(table_idx)            = NaN;
                    time_table.event_start(table_idx)           = total_run_frames;
                    time_table.event_dur(table_idx)             = params.exp.block.eye_gaze_pupil_white;
                    time_table.event_end(table_idx)             = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)              = params.exp.block.eye_gaze_pupil_white_ID;
                    time_table.event_name(table_idx)            = {'et_pupil_white'};
                    
                    total_run_frames = time_table.event_end(table_idx);
                    
                    % update time table idx
                    table_idx = table_idx+1;
                    

                    
                    % Add final rest period (3s gray)
                    time_table.session_nr(table_idx)            = ses;
                    time_table.session_type(table_idx)          = st;
                    time_table.run_nr(table_idx)                = rr;
                    time_table.crossing_nr(table_idx)           = NaN;
                    time_table.global_run_nr(table_idx)         = t_trial.global_run_nr(1);
                    time_table.global_trial_nr(table_idx)       = NaN;
                    time_table.global_block_nr(table_idx)       = NaN;
                    time_table.correct_response(table_idx)      = NaN;
                    time_table.is_objectcatch(table_idx)        = NaN;
                    time_table.block_nr(table_idx)              = 999;
                    time_table.trial_nr(table_idx)              = NaN;
                    time_table.stim_class(table_idx)            = NaN;
                    time_table.task_class(table_idx)            = NaN;
                    time_table.event_start(table_idx)           = total_run_frames;
                    time_table.event_dur(table_idx)             = params.exp.block.eye_gaze_fix2;
                    time_table.event_end(table_idx)             = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)              = params.exp.block.eye_gaze_fix_ID;
                    time_table.event_name(table_idx)            = {'et_fix'};
                    
                    total_run_frames = time_table.event_end(table_idx);
                    
                    % update time table idx
                    table_idx = table_idx+1;
                                        
                    % Add pre-blank period
                    time_table.event_start(table_idx)      = total_run_frames;
                    time_table.event_dur(table_idx)        = session_preblankdur;
                    time_table.event_end(table_idx)        = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)         = 0;
                    time_table.event_name(table_idx)       = {'pre-blank'};
                    time_table.session_nr(table_idx)       = ses;
                    time_table.session_type(table_idx)     = st;
                    time_table.run_nr(table_idx)           = rr;
                    time_table.trial_nr(table_idx)         = NaN;
                    time_table.block_nr(table_idx)         = NaN;
                    time_table.global_run_nr(table_idx)    = t_trial.global_run_nr(1);
                    time_table.global_trial_nr(table_idx)  = NaN;
                    time_table.global_block_nr(table_idx)  = NaN;
                    time_table.correct_response(table_idx) = NaN;
                    time_table.is_objectcatch(table_idx)   = NaN;
                    time_table.crossing_nr(table_idx)      = NaN;
                    time_table.stim_class(table_idx)       = NaN;
                    time_table.task_class(table_idx)       = NaN;
                    total_run_frames = time_table.event_end(table_idx);
                    
                    % update time table idx
                    table_idx = table_idx+1;
                end
                
                % reset block vector counter
                block_vec_idx = 1;
                
                while 1
                    
                    if isempty(t_trial)
                        run_finished = 1;
                    end
                    
                    if run_finished
                        break;
                    end
                    
                    curr_trial = t_trial.trial_nr(1);
                    
                    % get current block nr
                    curr_block = block_nrs(block_vec_idx);
                    curr_global_run_nr            = t_trial.global_run_nr(1);
                    curr_global_block_nr          = t_trial.global_block_nr(1);
                    curr_global_trial_nr          = t_trial.global_trial_nr(1);
                    curr_crossing_nr              = t_trial.crossing_nr(1);
                    curr_crossing_name            = t_trial.crossing_name(1);
                    curr_stim_class               = t_trial.stim_class(1);
                    curr_task_class               = t_trial.task_class(1);
                    curr_stim_class_name          = t_trial.stim_class_name(1,:);
                    curr_task_class_name          = t_trial.task_class_name(1);
                    curr_repeat_nr                = t_trial.repeat_nr(1);
                    curr_trial_type               = t_trial.trial_type(1);
                    curr_trial_nr                 = t_trial.trial_nr(1);
                    curr_spatial_cue              = t_trial.is_cued(1);
                    curr_correct_response         = t_trial.correct_response(1);
                    curr_cd_start                 = t_trial.cd_start(1,:);
                    curr_objectcatch              = t_trial.is_objectcatch(1);
                    
                    % first trial of the block has a task cue
                    if t_trial.trial_nr(1) == 1 || ...
                            (block_vec_idx == 1) || (curr_block ~= block_nrs(block_vec_idx-1)) % new block
                        
                        % Reset shuffled ITIs prior to block start to ensure we
                        % get the right block lenght.
                        if curr_trial_type == 1 % single stim presentation block
                            if params.is_demo
                                itis = shuffle_concat(params.exp.trial.demo.ITI_single_block,1);
                            else
                                itis = shuffle_concat(params.exp.trial.ITI_single_block,1);
                            end
                        elseif curr_trial_type == 2  % double stim presentation block
                            if params.is_demo
                                itis = shuffle_concat(params.exp.trial.demo.ITI_double_block,1);
                            else
                                itis   = shuffle_concat(params.exp.trial.ITI_double_block,1);
                            end
                        end
                        
                        time_table.event_start(table_idx)    = total_run_frames;
                        time_table.event_dur(table_idx)      = params.exp.trial.task_cue_dur;
                        time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                        time_table.event_id(table_idx)       = params.exp.block.task_cue_ID;
                        time_table.event_name(table_idx)     = {'task-cue'};
                        time_table.session_nr(table_idx)      = ses;
                        time_table.session_type(table_idx)    = st;
                        time_table.run_nr(table_idx)          = rr;
                        time_table.global_run_nr(table_idx)   = curr_global_run_nr;
                        time_table.global_block_nr(table_idx) = curr_global_block_nr;
                        time_table.global_trial_nr(table_idx) = curr_global_trial_nr;
                        time_table.block_nr(table_idx)        = curr_block;
                        time_table.trial_nr(table_idx)        = curr_trial_nr;
                        time_table.crossing_nr(table_idx)     = curr_crossing_nr;
                        time_table.crossing_name(table_idx)   = curr_crossing_name;
                        time_table.stim_class(table_idx)      = curr_stim_class;
                        time_table.task_class(table_idx)      = curr_task_class;
                        time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                        time_table.task_class_name(table_idx) = curr_task_class_name;
                        time_table.repeat_nr(table_idx)       = NaN;
                        time_table.trial_type(table_idx)      = curr_trial_type;
                        time_table.is_cued(table_idx)         = curr_spatial_cue;
                        time_table.correct_response(table_idx) = NaN;
                        time_table.is_objectcatch(table_idx)  = NaN;
                        
                        % update time/table idx
                        total_run_frames = time_table.event_end(table_idx);
                        table_idx = table_idx+1;
                        
                        
                        % add iti that lasts at least 1 second
                        time_table.event_start(table_idx)    = total_run_frames;
                        time_table.event_dur(table_idx)      = params.exp.trial.post_task_cue_ITI_dur;
                        time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                        time_table.event_id(table_idx)       = params.exp.block.post_task_cue_ITI_ID;
                        time_table.event_name(table_idx)     = {'post-task-ITI'};
                        time_table.session_nr(table_idx)      = ses;
                        time_table.session_type(table_idx)    = st;
                        time_table.run_nr(table_idx)          = rr;
                        time_table.global_run_nr(table_idx)   = curr_global_run_nr;
                        time_table.global_block_nr(table_idx) = curr_global_block_nr;
                        time_table.global_trial_nr(table_idx) = curr_global_trial_nr;
                        time_table.block_nr(table_idx)        = curr_block;
                        time_table.trial_nr(table_idx)        = curr_trial_nr;
                        time_table.crossing_nr(table_idx)     = curr_crossing_nr;
                        time_table.crossing_name(table_idx)   = curr_crossing_name;
                        time_table.stim_class(table_idx)      = curr_stim_class;
                        time_table.task_class(table_idx)      = curr_task_class;
                        time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                        time_table.task_class_name(table_idx) = curr_task_class_name;
                        time_table.repeat_nr(table_idx)       = NaN;
                        time_table.trial_type(table_idx)      = curr_trial_type;
                        time_table.is_cued(table_idx)         = curr_spatial_cue;
                        time_table.correct_response(table_idx) = NaN;
                        time_table.is_objectcatch(table_idx)  = NaN;
                        % update time/table idx
                        total_run_frames = time_table.event_end(table_idx);
                        table_idx = table_idx+1;
                    end
                    
                    % Get trial IDs of current trial
                    if t_trial.trial_type(1) ==1
                        trial_IDs = trial_ID_single_epoch;
                    else
                        trial_IDs = trial_ID_double_epoch;
                    end
                    
                    % Loop over trial events within a block
                    for id = 1:length(trial_IDs)
                        
                        % event start is the same for all IDs
                        time_table.event_start(table_idx) = total_run_frames;
                        
                        % add ses, run nr
                        time_table.session_nr(table_idx)      = ses;
                        time_table.session_type(table_idx)   = st;
                        time_table.run_nr(table_idx)          = rr;
                        time_table.global_run_nr(table_idx)   = curr_global_run_nr;
                        time_table.global_block_nr(table_idx) = curr_global_block_nr;
                        time_table.global_trial_nr(table_idx) = curr_global_trial_nr;
                        time_table.block_nr(table_idx)        = curr_block;
                        time_table.trial_nr(table_idx)        = curr_trial_nr;
                        time_table.crossing_nr(table_idx)     = curr_crossing_nr;
                        time_table.stim_class(table_idx)      = curr_stim_class;
                        time_table.task_class(table_idx)      = curr_task_class;
                        time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                        time_table.task_class_name(table_idx) = curr_task_class_name;
                        time_table.repeat_nr(table_idx)       = curr_repeat_nr;
                        time_table.trial_type(table_idx)      = curr_trial_type;
                        time_table.crossing_name(table_idx)      = curr_crossing_name;
                        time_table.correct_response(table_idx) = NaN;
                        time_table.is_objectcatch(table_idx)  = NaN;
                        
                        % Add individual trial events
                        switch trial_IDs(id)
                            
                            case params.exp.block.spatial_cue_ID  % 92 spatial cue
                                time_table.event_dur(table_idx)         = params.exp.trial.spatial_cue_dur;
                                time_table.event_end(table_idx)         = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                                time_table.event_id(table_idx)          = trial_IDs(id);
                                time_table.event_name(table_idx)        = {'spatial-cue'};
                                time_table.is_cued(table_idx)           = curr_spatial_cue;
                                time_table.repeat_nr(table_idx)         = NaN;
                                time_table.is_objectcatch(table_idx)    = NaN;
                                total_run_frames = time_table.event_end(table_idx);
                                table_idx   = table_idx+1;
                                
                            case params.exp.block.pre_stim_blank_ID % 93 blank time between spatial cue and stim onset
                                time_table.event_dur(table_idx)         = params.exp.trial.pre_stim_blank_dur;
                                time_table.event_end(table_idx)         = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                                time_table.event_id(table_idx)          = trial_IDs(id);
                                time_table.event_name(table_idx)        = {'pre-stim-blank'};
                                time_table.is_cued(table_idx)           = curr_spatial_cue;
                                time_table.repeat_nr(table_idx)         = NaN;
                                time_table.is_objectcatch(table_idx)    = NaN;
                                total_run_frames = time_table.event_end(table_idx);
                                table_idx   = table_idx+1;
                                
                            case params.exp.block.stim_epoch1_ID % 94 stim epoch 1 - left&right
                                tmp1 = t_trial(1,:);
                                tmp1.event_start    = total_run_frames(end);
                                tmp1.event_dur      = params.exp.trial.stim_array_dur;
                                tmp1.event_end      = tmp1.event_start + tmp1.event_dur;
                                tmp1.event_id       = trial_IDs(id);
                                tmp1.event_name     = {'stim1'};
                                if curr_task_class == 2 % cd 
                                    % convert relative time frame to stim
                                    % start, to relative time frame to run
                                    % start
                                    if ~isnan(tmp1.cd_start)
                                        tmp1.cd_start  = tmp1.cd_start + tmp1.event_start -1;
                                        assert( (tmp1.cd_start > tmp1.event_start) && (tmp1.cd_start < tmp1.event_end))
                                    else
                                        tmp1.cd_start  = NaN;
                                    end
                                end
                                tmp1.block_nr        = curr_block;
                                tmp1.session_nr      = ses;
                                tmp1.session_type    = st;
                                tmp1.run_nr          = rr;
                                tmp1.is_cued         = curr_spatial_cue;
                                tmp1.global_run_nr   = curr_global_run_nr;
                                tmp1.global_block_nr = curr_global_block_nr;
                                tmp1.global_trial_nr = curr_global_trial_nr;
                                tmp1.is_objectcatch  = curr_objectcatch;
                                tmp1.nr_fix_changes  = NaN;
                                time_table(table_idx,:) = tmp1;
                                total_run_frames = time_table.event_end(table_idx);
                                
                                table_idx   = table_idx+1;
                                
                                % if there is a second stim epoch, we wait with
                                % updating trial table. If single stim epoch,
                                % then we can delete trial
                                if t_trial.trial_type(1) == 1
                                    t_trial(1,:) = [];
                                end
                                
                            case params.exp.block.stim_epoch2_ID % 95: stim epoch 2 (after delay) - left&right
                                tmp2 = tmp1;
                                tmp2.event_start     = total_run_frames(end);
                                tmp2.event_dur       = params.exp.trial.stim_array_dur;
                                tmp2.event_end       = tmp2.event_start + tmp2.event_dur;
                                tmp2.event_id        = trial_IDs(id);
                                tmp2.event_name      = {'stim2'};
                                tmp2.block_nr        = curr_block;
                                tmp2.session_nr      = ses;
                                tmp2.session_type    = st;
                                tmp2.run_nr          = rr;
                                tmp2.global_run_nr   = curr_global_run_nr;
                                tmp2.global_block_nr = curr_global_block_nr;
                                tmp2.global_trial_nr = curr_global_trial_nr;
                                tmp2.is_cued         = curr_spatial_cue;
                                tmp2.correct_response = curr_correct_response;
                                tmp2.is_objectcatch  = curr_objectcatch;
                                tmp2.nr_fix_changes  = NaN;
                                time_table(table_idx,:) = tmp2;
                                
                                % reset correct_response to NaN for
                                % reference stim:
                                time_table.correct_response(table_idx-2,:) = NaN;
                                
                                total_run_frames = time_table.event_end(table_idx);
                                
                                table_idx = table_idx+1;
                                
                                % Now we update trial
                                t_trial(1,:) = [];
                                
                            case params.exp.block.response_ID % response
                                time_table.event_dur(table_idx)         = params.exp.trial.response_win_dur;
                                time_table.event_end(table_idx)         = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                                time_table.event_id(table_idx)          = params.exp.block.response_ID;
                                time_table.event_name(table_idx)        = {'response-cue'};
                                time_table.block_nr(table_idx)          = curr_block;
                                time_table.trial_nr(table_idx)          = curr_trial_nr;
                                time_table.crossing_nr(table_idx)       = curr_crossing_nr;
                                time_table.stim_class(table_idx)        = curr_stim_class;
                                time_table.task_class(table_idx)        = curr_task_class;
                                time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                                time_table.task_class_name(table_idx)   = curr_task_class_name;
                                time_table.repeat_nr(table_idx)         = NaN;
                                time_table.trial_type(table_idx)        = curr_trial_type;
                                time_table.crossing_name(table_idx)     = curr_crossing_name;
                                time_table.is_cued(table_idx)           = curr_spatial_cue;
                                time_table.correct_response(table_idx)  = NaN; % correct response is tied to last stimulus period
                                time_table.global_run_nr(table_idx)     = curr_global_run_nr;
                                time_table.global_block_nr(table_idx)   = curr_global_block_nr;
                                time_table.global_trial_nr(table_idx)   = curr_global_trial_nr;
                                time_table.is_objectcatch(table_idx)    = NaN;
                                
                                total_run_frames = time_table.event_end(table_idx);
                                
                                table_idx = table_idx+1;
                                
                            case params.exp.block.delay_ID % delay
                                time_table.event_dur(table_idx)         = params.exp.trial.delay_dur;
                                time_table.event_end(table_idx)         = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                                time_table.event_id(table_idx)          = trial_IDs(id);
                                time_table.event_name(table_idx)        = {'delay'};
                                time_table.block_nr(table_idx)          = curr_block;
                                time_table.trial_nr(table_idx)          = curr_trial_nr;
                                time_table.crossing_nr(table_idx)       = curr_crossing_nr;
                                time_table.stim_class(table_idx)        = curr_stim_class;
                                time_table.task_class(table_idx)        = curr_task_class;
                                time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                                time_table.task_class_name(table_idx)   = curr_task_class_name;
                                time_table.repeat_nr(table_idx)         = NaN;
                                time_table.trial_type(table_idx)        = curr_trial_type;
                                time_table.crossing_name(table_idx)     = curr_crossing_name;
                                time_table.is_cued(table_idx)           = curr_spatial_cue;
                                time_table.correct_response(table_idx)  = NaN; % correct response is tied to last stimulus period
                                time_table.global_run_nr(table_idx)     = curr_global_run_nr;
                                time_table.global_block_nr(table_idx)   = curr_global_block_nr;
                                time_table.global_trial_nr(table_idx)   = curr_global_trial_nr;
                                time_table.is_objectcatch(table_idx)    = NaN;
                                
                                total_run_frames = time_table.event_end(table_idx);
                                
                                table_idx = table_idx+1;
                        end
                        
                    end % trial IDs
                    
                    
                    %% ADD inter-trial, inter-block interval
                    if block_vec_idx == length(block_nrs) %  last trial of the last block
                        
                        run_finished = 1;
                        break;
                        
                    elseif block_vec_idx < length(block_nrs) % more trials to go
                        
                        next_block_id = block_nrs(block_vec_idx+1);
                        
                        if curr_block ~= next_block_id
                            % next trial is from different block (tt is already
                            % updated), so this is the last trial of the
                            % current block
                            
                            % add IBI
                            time_table.event_start(table_idx)  = total_run_frames;
                            time_table.event_start(table_idx)  = total_run_frames;
                            time_table.event_dur(table_idx)    = ibis(1);
                            time_table.event_name(table_idx)   = {'IBI'};
                            time_table.event_id(table_idx)     = params.exp.block.IBI_ID;
                            time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                            % add ses, block, run nr
                            time_table.session_nr(table_idx)       = ses;
                            time_table.session_type(table_idx)     = st;
                            time_table.run_nr(table_idx)           = rr;
                            time_table.block_nr(table_idx)         = NaN; % insert NaN for block nr
                            time_table.trial_nr(table_idx)         = NaN;
                            time_table.crossing_nr(table_idx)      = NaN; % insert NaN for crossing nr
                            time_table.stim_class(table_idx)       = NaN;
                            time_table.task_class(table_idx)       = NaN;
                            time_table.is_cued(table_idx)          = NaN;
                            time_table.correct_response(table_idx) = NaN;
                            time_table.repeat_nr(table_idx)        = NaN;
                            time_table.global_run_nr(table_idx)    = curr_global_run_nr;
                            time_table.global_block_nr(table_idx)  = NaN; % insert NaN for global block nr
                            time_table.global_trial_nr(table_idx)  = NaN;
                            time_table.is_objectcatch(table_idx)   = NaN;
                            
                            total_run_frames = time_table.event_end(table_idx);
                            
                            ibis(1) = [];
                            table_idx = table_idx+1;
                            
                        elseif curr_block == next_block_id
                            % next trial is from the same block (tt is already
                            % updated), so not last trial of the block, and we
                            % need to insert ITI
                            if isempty(itis)
                                if time_table.trial_type(table_idx-1) == 1
                                    if params.is_demo
                                        itis = shuffle_concat(params.exp.trial.demo.ITI_single_block,1);
                                    else
                                        itis = shuffle_concat(params.exp.trial.ITI_single_block,1);
                                    end
                                else
                                    if params.is_demo
                                        itis = shuffle_concat(params.exp.trial.demo.ITI_double_block,1);
                                    else
                                        itis = shuffle_concat(params.exp.trial.ITI_double_block,1);
                                    end
                                end
                                % OLD: itis = vcd_optimizeITIs(max_block_dur-params.exp.trial.task_cue_dur, max_trial_dur, itis_to_use, max_nr_trials);
                            end
                            
                            % add ITI
                            time_table.event_start(table_idx)  = total_run_frames;
                            time_table.event_dur(table_idx)    = itis(1); % grab the first one from the shuffled list
                            time_table.event_name(table_idx)   = {'ITI'};
                            time_table.event_id(table_idx)     = params.exp.block.ITI_ID;
                            time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                            time_table.session_nr(table_idx)        = ses;
                            time_table.session_type(table_idx)      = st;
                            time_table.run_nr(table_idx)            = rr;
                            time_table.block_nr(table_idx)          = curr_block;
                            time_table.trial_nr(table_idx)          = NaN;
                            time_table.crossing_nr(table_idx)       = NaN; % insert NaN for crossing nr
                            time_table.stim_class(table_idx)        = NaN;
                            time_table.task_class(table_idx)        = NaN;
                            time_table.is_cued(table_idx)           = NaN;
                            time_table.correct_response(table_idx)  = NaN;
                            time_table.repeat_nr(table_idx)         = NaN;
                            time_table.global_run_nr(table_idx)     = curr_global_run_nr;
                            time_table.global_block_nr(table_idx)   = curr_global_block_nr;
                            time_table.global_trial_nr(table_idx)   = NaN;
                            time_table.is_objectcatch(table_idx)    = NaN;
                            
                            total_run_frames = time_table.event_end(table_idx);
                            itis(1) = []; %  remove used ITI it
                            table_idx = table_idx+1;
                        end
                        
                        % update block idx
                        block_vec_idx = block_vec_idx+1;
                    end % ITI/IBI if statement
                    
                end % while loop
                
                % add post-run blank
                if run_finished || block_vec_idx > length(block_nrs) || (block_vec_idx+1)==length(block_nrs) % last trial of block and last block
                    
                    % remove last ITI
                    if strcmp(time_table.event_name(table_idx-1),'ITI')
                        time_table(table_idx-1,:) = [];
                        table_idx = table_idx-1;
                    end
                    
                    % session_totalrundur_run, total_ses_dur
                    % round out post blank period to finish on full TR
                    total_run_time2 = total_run_frames + session_postblankdur + postblank_to_add;
                    
                    if total_run_time2 == session_totalrundur
                        post_blank_dur = session_postblankdur + postblank_to_add;
                    elseif total_run_time2 < session_totalrundur
                        round_me_out = (session_totalrundur - total_run_time2);
                        post_blank_dur = session_postblankdur + postblank_to_add + round_me_out;
                    elseif total_run_time2 > session_totalrundur
                        if postblank_to_add  == 0
                            if params.is_demo % if it is a demo, then we don't care and shave off a little bit from the end
                                session_postblankdur = total_run_time2-session_totalrundur;
                                total_run_time2 = session_totalrundur;
                            end
                            error('[%s]: We have too many frames for this run!!')
                        elseif (total_run_frames + session_postblankdur) < session_totalrundur
                            left_over_blank = (session_totalrundur - (total_run_frames + session_postblankdur));
                            post_blank_dur  = session_postblankdur + left_over_blank;
                            assert(isequal((total_run_frames + post_blank_dur),session_totalrundur))
                        else
                            error('[%s]: Run duration doesn''t have the expected length',mfilename)
                        end
                    end
                    % add post_blank
                    time_table.event_start(table_idx)      = total_run_frames;
                    time_table.event_dur(table_idx)        = post_blank_dur;
                    time_table.event_name(table_idx)       = {'post-blank'};
                    time_table.event_id(table_idx)         = 0;
                    time_table.event_end(table_idx)        = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.session_nr(table_idx)       = ses;
                    time_table.session_type(table_idx)     = st;
                    time_table.run_nr(table_idx)           = rr;
                    time_table.block_nr(table_idx)         = NaN;
                    time_table.trial_nr(table_idx)         = NaN;
                    time_table.crossing_nr(table_idx)      = NaN;
                    time_table.stim_class(table_idx)       = NaN;
                    time_table.task_class(table_idx)       = NaN;
                    time_table.is_cued(table_idx)          = NaN;
                    time_table.repeat_nr(table_idx)        = NaN;
                    time_table.global_run_nr(table_idx)    = time_table.global_run_nr(table_idx-1);
                    time_table.global_block_nr(table_idx)  = NaN;
                    time_table.global_trial_nr(table_idx)  = NaN;
                    time_table.correct_response(table_idx) = NaN;
                    time_table.is_objectcatch(table_idx)   = NaN;
                    
                    total_run_frames = time_table.event_end(table_idx);
                    
                    assert(isequal(session_totalrundur, total_run_frames));
                    assert(nearZero(mod(total_run_frames,session_totalrundur)))
                    if strcmp(env_type,'MRI')
                        assert(nearZero(mod(total_run_frames,params.exp.TR)));
                    end
                    
                    run_finished = 0; % reset flag
                    
                    % trim time table
                    time_table2 = time_table;
                    if strcmp(time_table2.event_name(table_idx),{'post-blank'})
                        time_table2((table_idx+1):end,:) = [];
                    end
                    
                    %
                    % reorganize column order
                    col_order = {'session_nr','session_type','run_nr','block_nr','trial_nr',...
                        'global_run_nr','global_block_nr','global_trial_nr','condition_nr',...
                        'condition_name','stim_class','stim_class_name','task_class',...
                        'task_class_name','crossing_nr','crossing_name',...
                        'stim_nr_left','stim_nr_right',...
                        'is_cued','is_catch','is_objectcatch', 'correct_response',...
                        'event_start','event_dur', 'event_end','event_id','event_name', 'cd_start',...
                        'orient_dir','contrast','gbr_phase','rdk_coherence',...
                        'super_cat','super_cat_name','basic_cat','basic_cat_name',...
                        'sub_cat','sub_cat_name','affordance_cat','affordance_name',...
                        'stim2_im_nr', 'stim2_delta', 'stim2_orient_dir',...
                        'is_special_core','is_lure','repeat_nr','trial_type','nr_fix_changes'};
                    
                    % check if we missed any columns
                    assert(all(ismember(time_table2.Properties.VariableNames,col_order)));
                    
                    % apply the reordering
                    for col = 1:length(col_order)
                        new_col_order(col) = find(strcmp(time_table2.Properties.VariableNames,col_order(col)));
                    end
                    time_table2 = time_table2(:,new_col_order);
                    
                    % Clean up table
                    
                    % Set names to NaN if empty
                    time_table2.condition_name(cellfun(@isempty,time_table2.condition_name))   = {NaN};
                    time_table2.stim_class_name(cellfun(@isempty,time_table2.stim_class_name)) = {NaN};
                    time_table2.task_class_name(cellfun(@isempty,time_table2.task_class_name)) = {NaN};
                    time_table2.crossing_name(cellfun(@isempty,time_table2.crossing_name))     = {NaN};
                    time_table2.super_cat_name(cellfun(@isempty,time_table2.super_cat_name))   = {NaN};
                    time_table2.basic_cat_name(cellfun(@isempty,time_table2.basic_cat_name))   = {NaN};
                    time_table2.sub_cat_name(cellfun(@isempty,time_table2.sub_cat_name))       = {NaN};
                    time_table2.affordance_name(cellfun(@isempty,time_table2.affordance_name)) = {NaN};
                    
                    % set local and global trials/runs/sessions to NaN if zero
                    time_table2.run_nr(time_table2.run_nr==0) = NaN;
                    time_table2.block_nr(time_table2.block_nr==0) = NaN;
                    time_table2.trial_nr(time_table2.trial_nr==0) = NaN;
                    time_table2.global_run_nr(time_table2.global_run_nr==0) = NaN;
                    time_table2.global_block_nr(time_table2.global_block_nr==0) = NaN;
                    time_table2.global_trial_nr(time_table2.global_trial_nr==0) = NaN;
                    
                    % set stimclass/taskclass/crossingnr/is_cued to NaN if zero
                    time_table2.stim_class(time_table2.stim_class==0) = NaN;
                    time_table2.task_class(time_table2.task_class==0) = NaN;
                    time_table2.crossing_nr(time_table2.crossing_nr==0) = NaN;
                    time_table2.is_cued(time_table2.is_cued==0) = NaN;
                    
                    % concatenate time_table2 to run_time_table
                    if ~exist('run_time_table','var')
                        run_time_table = time_table2([],:); % create empty table with the same columns
                    end
                    run_time_table = cat(1, run_time_table, time_table2);
                    
                    clear time_table2
                end
            end
         
            % concatenate run_time_table to time_table_master
            if ~exist('time_table_master','var')
                time_table_master = run_time_table([],:);
            end
            time_table_master  = cat(1,time_table_master,run_time_table);
            
            if ~verLessThan('matlab', '9.6') % if we have an old MATLAB version, don't even bother making figures (this code is not compatible)
            
                if verbose

                    figure(99); clf; set(gcf, 'Position', [1 1 1200 500])
                    imagesc([time_table_master.run_nr,time_table_master.block_nr,time_table_master.trial_nr,time_table_master.crossing_nr]')
                    set(gca,'YTick',[1:4],'YTickLabel',{'run nr','block nr', 'trial nr', 'crossing nr'})
                    colormap(cmapturbo(32))
                    colorbar
                    set(gca,'CLim', [0 32]);
                    xlabel('events (in time)')
                    sgtitle(sprintf('Session %02d %s',ses, choose(st==1,'A','B')))

                    if store_imgs
                        saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master1_%s',env_type));
                        filename = sprintf('vcd_session%02d_%s_time_table_events.png', ses, choose(st==1,'A','B'));
                        print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                    end
                end
            end

        end % ~nan(session type)
    end % session type
end % session nr

%% AT LAST: We expand the "time_table_master" with the fixation sequence
% and onset of contrast dip, and correct button presses for FIX and CD
% task-crossings
[time_table_master, all_run_frames] = vcd_createRunFrames(params,time_table_master,env_type);

%% Store tables locally, if requested
if store_params
    if isempty(saveDir)
        saveDir = fullfile(vcd_rootPath,'data',env_type, subj_id);
    end
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    fname = sprintf('%stime_table_master_%s%s_%s.mat', [subj_id '_'], choose(params.is_demo,'demo_',''),params.disp.name, datestr(now,30));
    fprintf('[%s]: Storing time table master for subject in:\n',mfilename)
    fprintf('\t%s\n',fullfile(saveDir,fname))
    save(fullfile(saveDir, fname), 'time_table_master','all_run_frames')
end
    

return