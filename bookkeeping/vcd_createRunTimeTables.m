function time_table_master = vcd_createRunTimeTables(params,condition_master_in,session_type)
% VCD function to add trial events to condition master, turning it into a
% time table. This function also cleans some table columns that are now obsolete.
%
%   condition_master = vcd_allocateBlocksToRuns(params,condition_master_in)
%
% INPUTS:
%   params              : 	(struct) parameter struct needed to get subject
%                           nrs and run type params.
%  condition_master_in  :   (table) condition master table with subject
%                           session,runs,blocks,trials, but no individual
%                           trial events.
%  [session_type]       :   (optional) Are we dealing with MRI or BEHAVIOR sessions
%                           Default: 'MRI'
% OUTPUTS:
%   time_table_master   :  (struct) updated condition master table with
%                           single trial events.
%     {'session_nr'          }
%     {'run_nr'              }
%     {'block_nr'            }
%     {'trial_nr'}
%     {'condition_name'      }
%     {'stim_class'          } : (double) stimulus class number [1-5]: GBR, RDK, DOT, OBJ, NS
%     {'task_class'          } : (double) task class number [1-10]: FIX, CD, SCC, PC, WM, LTM, IMG, WHAT, WHERE, HOW
%     {'stim_nr_left'        } : (double) unique stimulus number for left spatial stimulus location or center location in case of NS
%     {'stim_nr_right'       } : (double) unique stimulus number for left spatial stimulus location, left empty in case of NS
%     {'is_cued'             } : (double) if stimulus is cued left (1), right (2), or neutral (3)
%     {'is_catch'            } : (bool) if stimulus is a catch trial (true) or not (false)
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
%     {'stim2_im_nr'        } : (cell w str or double) outcome of 'stim2_delta'. For
%                                 example, for -8 + 18 = 10 degrees gabor
%                                 tilt for quiz image after delay.
%     {'ltm_stim_pair'       } : (double) Long term memory paired stimulus
%     {'is_lure'             } : (bool) if LTM quiz image is lure (true) or not (false)
%     {'repeat_nr'           } : (double) how often this particular stimulus-task crossing has been repeated across sessions for a subject
%     {'global_block_nr'     } : (double) block number across all sessions for a subject


% Copy condition master input
condition_master = condition_master_in;

% check inputs
if (nargin ==2 && ~exist(session_type,'var')) || isempty(session_type)
    session_type = 'MRI';
elseif strcmp(session_type,'MRI')
    % Concatenate wide and deep sessions
    all_sessions     = cat(3, params.exp.session.wide.ses_blocks, params.exp.session.deep.ses_blocks);
    runs_per_session = cat(2, params.exp.session.mri.wide.n_runs_per_session,params.exp.session.mri.deep.n_runs_per_session);
    session_types    = cat(2, params.exp.session.mri.wide.session_types,params.exp.session.mri.deep.session_types);
    
    % Use longer plank periods
    session_preblankdur  = params.exp.run.pre_blank_dur_MRI;
    session_postblankdur = params.exp.run.post_blank_dur_MRI;
    session_totalrundur  = params.exp.run.total_run_dur_MRI;
    IBI_to_use = params.exp.block.IBI_MRI;
    
    
elseif strcmp(session_type,'BEHAVIOR')
    all_sessions         = params.exp.session.behavior.ses_blocks;
    runs_per_session     = params.exp.session.behavior.n_runs_per_session;
    session_types        = params.exp.session.behavior.session_types;
    
    % Use shorter plank periods
    session_preblankdur  = params.exp.run.pre_blank_dur_BEHAVIOR;
    session_postblankdur = params.exp.run.post_blank_dur_BEHAVIOR;
    session_totalrundur  = params.exp.run.total_run_dur_BEHAVIOR;
    IBI_to_use           = params.exp.block.IBI_BEHAVIOR;
end

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

% blocktype allocation
blocktype_allocation = cat(2,params.exp.run.run_type1, params.exp.run.run_type2,params.exp.run.run_type3,params.exp.run.run_type4);

% Preallocate space
time_table_master = [];
tbl_nrows    = 1000; % pick an arbitrary large number of rows (otherwise table uses 0 for missing values and we don't want that)
tr_in_frames = (params.exp.TR*params.stim.presentationrate_hz);

% preallocate session time table
session_time_table = [];

et_global_block_nr = 1;
et_global_trial_nr = 1;

% Loop over sessions
for ses = 1:size(all_sessions,3)
    
    for st = 1:size(all_sessions,4)
        if ~isnan(session_types(ses,st))
            fprintf('SESSION %s %03d %02d..\n',ses,session_types(st))
            
            % preallocate subject time table
            run_time_table = [];
            subj_run_nrs = unique(condition_master.run_nr(~isnan(condition_master.run_nr)));
            assert(isequal(runs_per_session(ses),length(subj_run_nrs)))
            
            % Loop over runs
            for rr = 1:runs_per_session(ses)
                
                % copy condition order table headers and scrub content
                sz = [tbl_nrows size(condition_master(1,:),2)];
                time_table = vcd_preallocateNaNTable(sz(1), sz(2), condition_master(1,:), []);
                
                % add column for event timing and id
                time_table.event_start = NaN(tbl_nrows,1);
                time_table.event_dur   = NaN(tbl_nrows,1);
                time_table.event_end   = NaN(tbl_nrows,1);
                time_table.event_id    = NaN(tbl_nrows,1);
                time_table.event_name  = repmat({''},tbl_nrows,1);
                time_table.global_block_nr = NaN(tbl_nrows,1);
                
                % Find corresponding trials (rows in table) and block nrs for
                % this particular subject, run, and session. The t_trial table
                % will be inserted into a bigger time table that includes other
                % event types, like iti, cues, etc.
                idx0      = (condition_master.session_nr==ses & condition_master.session_type==st & condition_master.run_nr==rr);
                assert(sum(idx0)>0)
                t_trial   = condition_master(idx0,:);
                block_nrs = t_trial.block_nr;
                
                % Reorder trials according to subject's specific order
                [block_nrs,block_nrs_idx] = sort(block_nrs,'ascend');
                t_trial = t_trial(block_nrs_idx,:);
                
                [~,nr_unique_blocks] = unique(block_nrs,'stable');
                tmp_durs = [params.exp.block.total_single_epoch_dur, params.exp.block.total_double_epoch_dur];
                block_dur = tmp_durs(t_trial.trial_type(nr_unique_blocks));
                
                [~,run_type_idx] = intersect(blocktype_allocation',[sum(t_trial.trial_type(nr_unique_blocks)==1);sum(t_trial.trial_type(nr_unique_blocks)==2)]','rows');
                if isempty(run_type_idx)
                    special_run_flag = true;
                else
                    special_run_flag = false;
                    expected_nr_of_blocktypes = blocktype_allocation(:,run_type_idx);
                end
                
                % reset counter
                total_run_frames = 0;
                run_finished     = 0;
                
                if ~isempty(block_nrs)
                    % remove frames such that duration in seconds is an integer (we will add those later)
                    rounded_session_totalrundur = (floor(session_totalrundur/params.stim.presentationrate_hz)*params.stim.presentationrate_hz);
                    total_ses_dur =  rounded_session_totalrundur - ...
                        params.exp.block.total_eyetracking_block_dur - session_postblankdur - session_preblankdur;
                    
                    % sneak in extra block for one behavioral run
                    if strcmp(session_type,'BEHAVIOR') && special_run_flag
                        total_ses_dur = total_ses_dur + params.exp.block.total_single_epoch_dur + min(IBI_to_use);
                        rounded_session_totalrundur = total_ses_dur + params.exp.block.total_eyetracking_block_dur + session_postblankdur + session_preblankdur;
                    end
                    
                    % predefine IBIs, make sure we don't go over the total
                    % run duration we want
                    while 1
                        % vcd_optimizeIBIs inputs are run_dur, block_dur, ibis, nr_blocks, prepost_blank
                        [ibis, postblank_to_add] = vcd_optimizeIBIs(total_ses_dur, block_dur, IBI_to_use, length(block_dur), 0); % prepost_blank is 0 because we already subtract this when defining total_ses_dur
                        if (rounded_session_totalrundur - sum(ibis) - sum(block_dur) - ...
                                postblank_to_add -  params.exp.block.total_eyetracking_block_dur - ...
                                session_postblankdur -  session_preblankdur) <= 0  
                            break;
                        end
                    end
                    fprintf('Selected IBIs (in time frames): \n')
                    disp(ibis)
                    fprintf('\n')
                end
                % Add eyetracking block and pre-blank period
                if total_run_frames == 0
                    
                    table_idx = 1;
                    % Add prefixation period (1s)
                    time_table.session_nr(table_idx)       = ses;
                    time_table.session_type(table_idx)     = st;
                    time_table.run_nr(table_idx)           = rr;
                    time_table.block_nr(table_idx)         = 999;
                    time_table.crossing_nr(table_idx)      = 999;
                    time_table.stim_class(table_idx)       = NaN;
                    time_table.task_class(table_idx)       = NaN;
                    time_table.repeat_nr(table_idx)        = t_trial.repeat_nr(1);
                    time_table.trial_nr(table_idx)         = 1;
                    time_table.global_trial_nr(table_idx)   = et_global_trial_nr;
                    time_table.global_block_nr(table_idx)   = et_global_block_nr;
                    time_table.global_session_nr(table_idx) = t_trial.global_session_nr(1);
                    time_table.correct_response(table_idx) = NaN;
                    time_table.event_start(table_idx)      = total_run_frames;
                    time_table.event_dur(table_idx)        = params.exp.block.eye_gaze_fix0;
                    time_table.event_end(table_idx)        = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)         = params.exp.block.eye_gaze_fix_ID;
                    time_table.event_name(table_idx)       = {'et_fix'};
                    
                    
                    total_run_frames = time_table.event_end(table_idx);
                    et_global_trial_nr = et_global_trial_nr+1;
                    
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
                        time_table.global_trial_nr(table_idx)   = et_global_trial_nr;
                        time_table.global_block_nr(table_idx)   = et_global_block_nr;
                        time_table.global_session_nr(table_idx) = t_trial.global_session_nr(1);
                        time_table.crossing_nr(table_idx)        = 999;
                        time_table.trial_nr(table_idx)          = 1+sac;
                        time_table.stim_class(table_idx)        = NaN;
                        time_table.task_class(table_idx)        = NaN;
                        time_table.event_start(table_idx)       = total_run_frames;
                        time_table.event_dur(table_idx)         = params.exp.block.eye_gaze_sac_target;
                        time_table.event_end(table_idx)         = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                        time_table.event_id(table_idx)          = sac_IDs(sac);
                        time_table.event_name(table_idx)        = {'et_sac'};
                        total_run_frames = time_table.event_end(table_idx);
                        
                        et_global_trial_nr = et_global_trial_nr+1;
                        
                        % update time table idx
                        table_idx = table_idx+1;
                    end
                    
                    % Add post saccade fixation period (1s)
                    time_table.run_nr(table_idx)                = rr;
                    time_table.session_nr(table_idx)            = ses;
                    time_table.session_type(table_idx)          = st;
                    time_table.crossing_nr(table_idx)           = 999;
                    time_table.global_trial_nr(table_idx)       = et_global_trial_nr;
                    time_table.global_block_nr(table_idx)   = et_global_block_nr;
                    time_table.global_session_nr(table_idx) = t_trial.global_session_nr(1);
                    time_table.block_nr(table_idx)              = 999;
                    time_table.trial_nr(table_idx)              = params.exp.block.nr_of_saccades+2;
                    time_table.stim_class(table_idx)            = NaN;
                    time_table.task_class(table_idx)            = NaN;
                    time_table.event_start(table_idx)           = total_run_frames;
                    time_table.event_dur(table_idx)             = params.exp.block.eye_gaze_fix1;
                    time_table.event_end(table_idx)             = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)              = params.exp.block.eye_gaze_fix_ID;
                    time_table.event_name(table_idx)            = {'et_fix'};
                    
                    total_run_frames = time_table.event_end(table_idx);
                    et_global_trial_nr = et_global_trial_nr+1;
                    
                    % update time table idx
                    table_idx = table_idx+1;
                    
                    % Add pupil trial (3s black)
                    time_table.run_nr(table_idx)                = rr;
                    time_table.session_nr(table_idx)            = ses;
                    time_table.session_type(table_idx)          = st;
                    time_table.crossing_nr(table_idx)           = 999;
                    time_table.global_trial_nr(table_idx)       = et_global_trial_nr;
                    time_table.global_block_nr(table_idx)       = et_global_block_nr;
                    time_table.global_session_nr(table_idx)     = t_trial.global_session_nr(1);
                    time_table.block_nr(table_idx)              = 999;
                    time_table.trial_nr(table_idx)              = params.exp.block.nr_of_saccades+3;
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
                    
                    et_global_trial_nr = et_global_trial_nr+1;
                    
                    % Add pupil trial (1s white)
                    time_table.run_nr(table_idx)                = rr;
                    time_table.session_nr(table_idx)            = ses;
                    time_table.session_type(table_idx)          = st;
                    time_table.crossing_nr(table_idx)           = 999;
                    time_table.global_trial_nr(table_idx)       = et_global_trial_nr;
                    time_table.global_block_nr(table_idx)       = et_global_block_nr;
                    time_table.global_session_nr(table_idx)     = t_trial.global_session_nr(1);
                    time_table.block_nr(table_idx)              = 999;
                    time_table.trial_nr(table_idx)              = params.exp.block.nr_of_saccades+4;
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
                    
                    et_global_trial_nr = et_global_trial_nr+1;
                    
                    % Add final rest period (3s gray)
                    time_table.run_nr(table_idx)                = rr;
                    time_table.session_nr(table_idx)            = ses;
                    time_table.session_type(table_idx)          = st;
                    time_table.crossing_nr(table_idx)           = 999;
                    time_table.global_trial_nr(table_idx)       = et_global_trial_nr;
                    time_table.global_block_nr(table_idx)       = et_global_block_nr;
                    time_table.global_session_nr(table_idx)     = t_trial.global_session_nr(1);
                    time_table.block_nr(table_idx)              = 999;
                    time_table.trial_nr(table_idx)              = params.exp.block.nr_of_saccades+5;
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
                    
                    % update global eyetracking counters
                    et_global_trial_nr = et_global_trial_nr+1;
                    et_global_block_nr = et_global_block_nr+1;
                    
                    
                    % Add pre-blank period
                    time_table.event_start(table_idx)    = total_run_frames;
                    time_table.event_dur(table_idx)      = session_preblankdur;
                    time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.event_id(table_idx)       = 0;
                    time_table.event_name(table_idx)     = {'pre-blank'};
                    
                    time_table.run_nr(table_idx)         = rr;
                    time_table.session_nr(table_idx)     = ses;
                    time_table.session_type(table_idx)   = st;
                    time_table.block_nr(table_idx)       = 0;
                    time_table.global_trial_nr(table_idx)   = t_trial.global_trial_nr(1);
                    time_table.global_block_nr(table_idx)   = t_trial.global_block_nr(1);
                    time_table.global_session_nr(table_idx) = t_trial.global_session_nr(1);
                    time_table.crossing_nr(table_idx)    = 0;
                    time_table.trial_nr(table_idx)       = 0;
                    time_table.crossing_nr(table_idx)    = 0;
                    time_table.stim_class(table_idx)     = 0;
                    time_table.task_class(table_idx)     = 0;
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
                    curr_global_session_nr        = t_trial.global_session_nr(1);
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
                    
                    % first trial of the block has a task cue
                    if t_trial.trial_nr(1) == 1 || ...
                            (block_vec_idx == 1) || (curr_block ~= block_nrs(block_vec_idx-1)) % new block
                        
                        % Reset shuffled ITIs prior to block start to ensure we
                        % get the right block lenght.
                        if curr_trial_type == 1 % single stim presentation block
                            itis  = shuffle_concat(params.exp.trial.ITI_single_block,1);
                        elseif curr_trial_type == 2  % double stim presentation block
                            itis   = shuffle_concat(params.exp.trial.ITI_double_block,1);
                        end
                        % itis = vcd_optimizeITIs(max_block_dur-params.exp.trial.task_cue_dur, max_trial_dur, itis_to_use, max_nr_trials);
                        
                        time_table.event_start(table_idx)    = total_run_frames;
                        time_table.event_dur(table_idx)      = params.exp.trial.task_cue_dur;
                        time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                        time_table.event_id(table_idx)       = params.exp.block.task_cue_ID;
                        time_table.event_name(table_idx)     = {'task-cue'};
                        time_table.session_nr(table_idx)      = ses;
                        time_table.session_type(table_idx)    = st;
                        time_table.run_nr(table_idx)          = rr;
                        time_table.global_session_nr(table_idx) = curr_global_session_nr;
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
                        time_table.repeat_nr(table_idx)       = curr_repeat_nr;
                        time_table.trial_type(table_idx)      = curr_trial_type;
                        time_table.correct_response(table_idx) = NaN;
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
                        time_table.global_session_nr(table_idx) = curr_global_session_nr;
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
                        time_table.repeat_nr(table_idx)       = curr_repeat_nr;
                        time_table.trial_type(table_idx)      = curr_trial_type;
                        time_table.correct_response(table_idx) = NaN;
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
                        time_table.global_session_nr(table_idx) = curr_global_session_nr;
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
                        
                        % Add individual trial events
                        switch trial_IDs(id)
                            
                            case params.exp.block.spatial_cue_ID  % 92 spatial cue
                                time_table.event_dur(table_idx)         = params.exp.trial.spatial_cue_dur;
                                time_table.event_end(table_idx)         = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                                time_table.event_id(table_idx)          = trial_IDs(id);
                                time_table.event_name(table_idx)        = {'spatial-cue'};
                                %                                 time_table.block_nr(table_idx)          = curr_block;
                                
                                %                                 time_table.trial_nr(table_idx)          = curr_trial_nr;
                                %                                 time_table.crossing_nr(table_idx)       = curr_crossing_nr;
                                %                                 time_table.stim_class(table_idx)        = curr_stim_class;
                                %                                 time_table.task_class(table_idx)        = curr_task_class;
                                %                                 time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                                %                                 time_table.task_class_name(table_idx)   = curr_task_class_name;
                                %                                 time_table.repeat_nr(table_idx)         = curr_repeat_nr;
                                %                                 time_table.trial_type(table_idx)      = curr_trial_type;
                                %                                 time_table.crossing_name(table_idx)      = curr_crossing_name;
                                time_table.is_cued(table_idx)          = curr_spatial_cue;
                                
                                total_run_frames = time_table.event_end(table_idx);
                                table_idx   = table_idx+1;
                                
                            case params.exp.block.pre_stim_blank_ID % 93 blank time between spatial cue and stim onset
                                time_table.event_dur(table_idx)         = params.exp.trial.pre_stim_blank_dur;
                                time_table.event_end(table_idx)         = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                                time_table.event_id(table_idx)          = trial_IDs(id);
                                time_table.event_name(table_idx)        = {'pre-stim-blank'};
                                time_table.is_cued(table_idx)           = curr_spatial_cue;
                                
                                total_run_frames = time_table.event_end(table_idx);
                                table_idx   = table_idx+1;
                                
                            case params.exp.block.stim_epoch1_ID % 94 stim epoch 1 - left&right
                                tmp1 = t_trial(1,:);
                                tmp1.event_start    = total_run_frames(end);
                                tmp1.event_dur      = params.exp.trial.stim_array_dur;
                                tmp1.event_end      = tmp1.event_start(1) + tmp1.event_dur(1);
                                tmp1.event_id       = trial_IDs(id);
                                tmp1.event_name     = {'stim1'};
                                if curr_task_class == 2 % cd
                                    for side = find(~isnan(curr_cd_start))
                                        if tmp1.cd_start(side)~=0
                                            tmp1.cd_start(side)  = tmp1.cd_start(side) + tmp1.event_start(1) -1;
                                        end
                                    end
                                end
                                tmp1.block_nr       = curr_block;
                                tmp1.session_nr     = ses;
                                tmp1.session_type   = st;
                                tmp1.run_nr         = rr;
                                tmp1.is_cued         = curr_spatial_cue;
                                tmp1.global_session_nr = curr_global_session_nr;
                                tmp1.global_block_nr = curr_global_block_nr;
                                tmp1.global_trial_nr = curr_global_trial_nr;
                                if tmp1.trial_type == 1 % we are dealing with a single stim presentation trial
                                    tmp1.correct_response = curr_correct_response;
                                else % we are dealing with a double stim presentation trial, and insert the correct response during stim2
                                    tmp1.correct_response = NaN;
                                end
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
                                tmp2.event_start    = total_run_frames(end);
                                tmp2.event_dur      = params.exp.trial.stim_array_dur;
                                tmp2.event_end      = tmp2.event_start(1) + tmp2.event_dur(1);
                                tmp2.event_id       = trial_IDs(id);
                                tmp2.event_name     = {'stim2'};
                                tmp2.stim_nr_left   = tmp1.stim_nr_left;
                                tmp2.stim_nr_right  = tmp1.stim_nr_right;
                                tmp2.block_nr       = curr_block;
                                tmp2.session_nr     = ses;
                                tmp2.session_type   = st;
                                tmp2.run_nr         = rr;
                                tmp2.global_session_nr = curr_global_session_nr;
                                tmp2.global_block_nr = curr_global_block_nr;
                                tmp2.global_trial_nr = curr_global_trial_nr;
                                tmp2.is_cued         = curr_spatial_cue;
                                tmp2.correct_response = curr_correct_response;
                                
                                time_table(table_idx,:) = tmp2;
                                
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
                                time_table.repeat_nr(table_idx)         = curr_repeat_nr;
                                time_table.trial_type(table_idx)        = curr_trial_type;
                                time_table.crossing_name(table_idx)     = curr_crossing_name;
                                time_table.is_cued(table_idx)           = curr_spatial_cue;
                                time_table.correct_response(table_idx)  = NaN;
                                time_table.global_session_nr(table_idx) = curr_global_session_nr;
                                time_table.global_block_nr(table_idx)   = curr_global_block_nr;
                                time_table.global_trial_nr(table_idx)   = curr_global_trial_nr;
                                
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
                                time_table.repeat_nr(table_idx)         = curr_repeat_nr;
                                time_table.trial_type(table_idx)        = curr_trial_type;
                                time_table.crossing_name(table_idx)     = curr_crossing_name;
                                time_table.is_cued(table_idx)           = curr_spatial_cue;
                                time_table.correct_response(table_idx)  = NaN;
                                time_table.global_session_nr(table_idx) = curr_global_session_nr;
                                time_table.global_block_nr(table_idx)   = curr_global_block_nr;
                                time_table.global_trial_nr(table_idx)   = curr_global_trial_nr;
                                
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
                            IBI = ibis(1);
                            %                             if ~nearZero(mod(time_table.event_start(table_idx)+IBI,1))
                            %                                 IBI = IBI + (1 - mod(time_table.event_start(table_idx)+IBI,1)); % round out to a full second
                            %                             end
                            time_table.event_start(table_idx)  = total_run_frames;
                            time_table.event_dur(table_idx)    = IBI;
                            time_table.event_name(table_idx)   = {'IBI'};
                            time_table.event_id(table_idx)     = params.exp.block.IBI_ID;
                            time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                            % add ses, block, run nr
                            time_table.session_nr(table_idx)       = ses;
                            time_table.session_type(table_idx)     = st;
                            time_table.run_nr(table_idx)           = rr;
                            time_table.block_nr(table_idx)         = 0;
                            time_table.trial_nr(table_idx)      = 0;
                            time_table.crossing_nr(table_idx)      = 0;
                            time_table.stim_class(table_idx)       = 0;
                            time_table.task_class(table_idx)       = 0;
                            time_table.is_cued(table_idx)          = 0;
                            time_table.correct_response(table_idx) = NaN;
                            time_table.global_session_nr(table_idx) = curr_global_session_nr;
                            time_table.global_block_nr(table_idx) = 0;
                            time_table.global_trial_nr(table_idx) = 0;
                            
                            total_run_frames = time_table.event_end(table_idx);
                            
                            ibis(1) = []; clear IBI;
                            table_idx = table_idx+1;
                            
                        elseif curr_block == next_block_id
                            % next trial is from the same block (tt is already
                            % updated), so not last trial of the block, and we
                            % need to insert ITI
                            if isempty(itis)
                                if time_table.trial_type(table_idx-1) == 1
                                    itis = shuffle_concat(params.exp.trial.ITI_single_block,1);
                                else
                                    itis = shuffle_concat(params.exp.trial.ITI_double_block,1);
                                end
                                %                                 itis = vcd_optimizeITIs(max_block_dur-params.exp.trial.task_cue_dur, max_trial_dur, itis_to_use, max_nr_trials);
                            end
                            
                            % add ITI
                            ITI =  itis(1); % grab the first one from the shuffled list
                            itis(1) = []; % and remove it
                            time_table.event_start(table_idx)  = total_run_frames;
                            time_table.event_dur(table_idx)    = ITI;
                            time_table.event_name(table_idx)   = {'ITI'};
                            time_table.event_id(table_idx)     = params.exp.block.ITI_ID;
                            time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                            time_table.session_nr(table_idx)        = ses;
                            time_table.session_type(table_idx)      = st;
                            time_table.run_nr(table_idx)            = rr;
                            time_table.block_nr(table_idx)          = curr_block;
                            time_table.trial_nr(table_idx)          = 0;
                            time_table.crossing_nr(table_idx)       = 0;
                            time_table.stim_class(table_idx)        = 0;
                            time_table.task_class(table_idx)        = 0;
                            time_table.is_cued(table_idx)           = 0;
                            time_table.correct_response(table_idx)  = NaN;
                            time_table.global_session_nr(table_idx) = curr_global_session_nr;
                            time_table.global_block_nr(table_idx)   = curr_global_block_nr;
                            time_table.global_trial_nr(table_idx)   = 0;
                            
                            total_run_frames = time_table.event_end(table_idx);
                            
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
                            if special_run_flag
                                postblank_to_add = session_postblankdur;
                            else
                                error('[%s]: We have too many frames for this run!!')
                            end
                        elseif (total_run_frames + session_postblankdur) < session_totalrundur
                            left_over_blank = (session_totalrundur - (total_run_frames + session_postblankdur));
                            post_blank_dur  = session_postblankdur + left_over_blank;
                            assert(isequal((total_run_frames + post_blank_dur),session_totalrundur))
                        end
                    end
                    % add post_blank
                    time_table.event_start(table_idx)  = total_run_frames;
                    time_table.event_dur(table_idx)    = post_blank_dur;
                    time_table.event_name(table_idx)   = {'post-blank'};
                    time_table.event_id(table_idx)     = 0;
                    time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx);
                    time_table.session_nr(table_idx)    = ses;
                    time_table.session_type(table_idx)     = st;
                    time_table.run_nr(table_idx)        = rr;
                    time_table.block_nr(table_idx)      = 0;
                    time_table.trial_nr(table_idx)      = 0;
                    time_table.crossing_nr(table_idx)   = 0;
                    time_table.stim_class(table_idx)    = 0;
                    time_table.task_class(table_idx)    = 0;
                    time_table.is_cued(table_idx)       = 0;
                    time_table.global_session_nr(table_idx) = curr_global_session_nr;
                    time_table.global_block_nr(table_idx) = 0;
                    time_table.global_trial_nr(table_idx) = 0;
                    time_table.correct_response(table_idx)  = NaN;
                    
                    total_run_frames = time_table.event_end(table_idx);
                    
                    if strcmp(session_type,'MRI')
                        assert(isequal(session_totalrundur, total_run_frames));
                        assert(nearZero(mod(total_run_frames,params.exp.TR)));
                        assert(nearZero(mod(total_run_frames,session_totalrundur)))
                    elseif strcmp(session_type,'BEHAVIOR')
                        %                         assert(isequal(session_totalrundur_run, total_run_frames));
                        %                         assert(total_run_frames <= session_totalrundur_run)
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
                        'global_session_nr','global_block_nr','global_trial_nr','condition_nr',...
                        'condition_name','stim_class','stim_class_name','task_class',...
                        'task_class_name','crossing_nr','crossing_name',...
                        'stim_nr_left','stim_nr_right',...
                        'is_cued','is_catch', 'correct_response',...
                        'event_start','event_dur', 'event_end','event_id','event_name', 'cd_start',...
                        'orient_dir','contrast','gbr_phase','rdk_coherence',...
                        'super_cat','super_cat_name','basic_cat','basic_cat_name',...
                        'sub_cat','sub_cat_name','affordance_cat','affordance_name',...
                        'stim2_im_nr', 'stim2_delta', 'stim2_orient_dir',...
                        'is_special_core','is_lure','repeat_nr','trial_type'};
                    
                    assert(all(ismember(time_table2.Properties.VariableNames,col_order)));
                    for col = 1:length(col_order)
                        new_col_order(col) = find(strcmp(time_table2.Properties.VariableNames,col_order(col)));
                    end
                    time_table2 = time_table2(:,new_col_order);
                    
                    time_table2.is_lure(isnan(time_table2.is_lure(:,1)),1)=0;
                    time_table2.is_lure(isnan(time_table2.is_lure(:,2)),2)=0;
                    time_table2.is_catch(isnan(time_table2.is_catch))=0;
                    time_table2.is_special_core(isnan(time_table2.is_special_core(:,1)),1)=0;
                    time_table2.is_special_core(isnan(time_table2.is_special_core(:,2)),2)=0;
                    time_table2.condition_name(cellfun(@isempty,time_table2.condition_name)) = {NaN};
                    time_table2.stim_class_name(cellfun(@isempty,time_table2.stim_class_name)) = {NaN};
                    time_table2.task_class_name(cellfun(@isempty,time_table2.task_class_name)) = {NaN};
                    time_table2.crossing_name(cellfun(@isempty,time_table2.crossing_name)) = {NaN};
                    time_table2.super_cat_name(cellfun(@isempty,time_table2.super_cat_name)) = {NaN};
                    time_table2.basic_cat_name(cellfun(@isempty,time_table2.basic_cat_name)) = {NaN};
                    time_table2.sub_cat_name(cellfun(@isempty,time_table2.sub_cat_name)) = {NaN};
                    time_table2.affordance_name(cellfun(@isempty,time_table2.affordance_name)) = {NaN};
                    
                    % add run to master
                    run_time_table = cat(1, run_time_table, time_table2);
                    
                    
                    clear time_table2
                end
            end
            
            session_time_table = cat(1, session_time_table, run_time_table);
            if params.verbose
%                 session_time_table.run_nr(isnan(session_time_table.run_nr))=0;
%                 session_time_table.block_nr(isnan(session_time_table.block_nr))=0;
%                 session_time_table.block_nr(session_time_table.block_nr == 999)=-1;
%                 session_time_table.trial_nr(isnan(session_time_table.trial_nr))=0;
%                 session_time_table.crossing_nr(isnan(session_time_table.crossing_nr))=0;
%                 session_time_table.crossing_nr(session_time_table.crossing_nr>990)=-1;
%                 session_time_table.block_nr( session_time_table.block_nr==99) = 0.5;
                
                figure(99); clf; set(gcf, 'Position', [1 1 1200 500])
                imagesc([session_time_table.run_nr,session_time_table.block_nr,session_time_table.trial_nr,session_time_table.crossing_nr]')
                set(gca,'YTick',[1:4],'YTickLabel',{'run nr','block nr', 'trial nr', 'crossing nr'})
                colormap(cmapturbo(max([session_time_table.block_nr',session_time_table.run_nr'])))
                colorbar
                set(gca,'CLim', [0 110]);
                xlabel('events (in time)')
                sgtitle(sprintf('Session %02d %s',ses, choose(unique(session_time_table.session_type),'A','B')))
                
                if params.store_imgs
                    saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master1_%s',session_type));
                    filename = sprintf('vcd_session%02d_%s_time_table_events.png', ses, choose(unique(session_time_table.session_type),'A','B'),rr);
                    print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                end
            end
        end
    end % session type
    time_table_master = cat(1,time_table_master,session_time_table);
end % session nr



if params.store_params
    fprintf('\n[%s]:Storing time table master..\n',mfilename)
    saveDir = fileparts(fullfile(params.stim.fix.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir,sprintf('time_table_master_%s_%s_%s.mat',params.disp.name, session_type, datestr(now,30))),'time_table_master')
end

return