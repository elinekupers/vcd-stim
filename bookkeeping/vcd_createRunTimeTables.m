function time_table_master = vcd_createRunTimeTables(params,session_type)
% VCD function to add trial events to condition master, turning it into a 
% time table. This function also cleans some table columns that are now obsolete.
% 
%   condition_master = vcd_allocateBlocksToRuns(params)
%
% INPUTS:
%   params              : 	(struct) parameter struct needed to get subject
%                           nrs and run type params. REQUIRES params.trials 
%                           to exist and contain condition_master v0.                   
%  [session_type]       :   (optional) Are we dealing with MRI or BEHAVIOR sessions
%                           Default: 'MRI'   
% OUTPUTS:
%   time_table_master   :  (struct) updated condition master table with
%                           single trial events.
%     {'subj_nr'             } : (double) subject number. We assign a unique integral number between ranging from 001-999 to subjects.                                 
%     {'session_nr'          } : (double) MRI or behavioral session nr. There are 
%     {'session_name'        }
%     {'run_nr'              }
%     {'block_nr'            }
%     {'block_local_trial_nr'}
%     {'cond_name'           }
%     {'stim_class'          }
%     {'task_class'          }
%     {'stim_nr_left'        }
%     {'stim_nr_right'       }
%     {'is_cued'             }
%     {'is_catch'            }
%     {'trial_type'          }
%     {'block_ID'            }
%     {'event_start'         }
%     {'event_dur'           }
%     {'event_end'           }
%     {'event_id'            }
%     {'event_name'          }
%     {'stim_class_name'     }
%     {'task_class_name'     }
%     {'orient_dir'          }
%     {'contrast'            }
%     {'gbr_phase'           }
%     {'rdk_coherence'       }
%     {'super_cat'           }
%     {'basic_cat'           }
%     {'sub_cat'             }
%     {'affordance_cat'      }
%     {'super_cat_name'      }
%     {'basic_cat_name'      }
%     {'sub_cat_name'        }
%     {'affordance_name'     }
%     {'is_in_img_ltm'       }
%     {'stim2_delta'         }
%     {'stim2'               }
%     {'img_txt_prompt'      }
%     {'ltm_stim_pair'       }
%     {'is_lure'             }
%     {'repeat_nr'           }
%     {'global_block_nr'     }

%     {'subj_nr'             } 
%     {'session_nr'          } 
%     {'session_name'        } 
%     {'run_nr'              } 
%     {'block_nr'            } 
%     {'block_local_trial_nr'} 
%     {'cond_name'           } 
%     {'stim_class'          } : (double) stimulus class number [1-5]: GBR, RDK, DOT, OBJ, NS
%     {'task_class'          } : (double) task class number [1-10]: FIX, CD, SCC, PC, WM, LTM, IMG, WHAT, WHERE, HOW
%     {'stim_nr_left'        } : (double) unique stimulus number for left spatial stimulus location or center location in case of NS
%     {'stim_nr_right'       } : (double) unique stimulus number for left spatial stimulus location, left empty in case of NS
%     {'is_cued'             } : (double) if stimulus is cued left (1), right (2), or neutral (3)
%     {'is_catch'            } : (bool) if stimulus is a catch trial (true) or not (false)
%     {'trial_type'          } : (double) if it is a single-stimulus or double-stimulus presentation trial
%     {'block_ID'            } : (double) stimulus-task crossing (1:32), see params.exp.stimtaskcrossings
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


% check inputs
if (nargin ==1 && ~exist(session_type,'var')) || isempty(session_type)
    session_type = 'MRI';
elseif strcmp(session_type,'MRI')
    % Concatenate wide and deep sessions
    all_sessions = cat(3, params.exp.session.wide.ses_blocks, params.exp.session.deep.ses_blocks);
    session_names = {'WIDE1A', 'WIDE1B'};
    for ii = params.exp.session.deep.session_nrs, session_names = cat(2,session_names{:},{sprintf('DEEP%03d',ii)}); end
    session_names{length(session_names)-1,length(session_names)} = {'WIDE1A', 'WIDE1B'};
    runs_per_session = cat(2,params.exp.session.wide.n_runs_per_session,params.exp.session.deep.n_runs_per_session);
    session_preblankdur = params.exp.run.pre_blank_dur_MRI;
    session_postblankdur = params.exp.run.post_blank_dur_MRI;
    session_totalrundur  = params.exp.run.total_run_dur_MRI;
    IBI_to_use = params.exp.block.IBI_MRI;
elseif strcmp(session_type,'BEHAVIOR')
    session_names = {'BEHAVIOR01'};
    all_sessions = params.exp.session.behavior.ses_blocks;
    runs_per_session = params.exp.session.behavior.n_runs_per_session;
    session_preblankdur = params.exp.run.pre_blank_dur_BEHAVIOR;
    session_postblankdur = params.exp.run.post_blank_dur_BEHAVIOR;
    session_totalrundur  = params.exp.run.total_run_dur_BEHAVIOR;
    IBI_to_use = repmat(params.exp.block.IBI_BEHAVIOR,1,5);
end

% Define event ID within trial for single and double epoch trial types
trial_ID_single_epoch = [params.exp.block.trial_start_ID, params.exp.block.spatial_cue_ID, params.exp.block.stim_epoch1_ID, params.exp.block.response_ID];
trial_ID_double_epoch = [params.exp.block.trial_start_ID, params.exp.block.spatial_cue_ID, params.exp.block.stim_epoch1_ID, params.exp.block.delay_ID, params.exp.block.stim_epoch2_ID, params.exp.block.response_ID];

% Preallocate space
time_table_master = [];
tbl_nrows = 1000; % pick an arbitrary large number of rows (otherwise table uses 0 for missing values and we don't want that)
tr_in_frames = (params.exp.TR*params.stim.presentationrate_hz);


% Loop over subjects
for sj = 1;%:params.exp.total_subjects
    fprintf('\nSUBJECT %03d..',sj)
    
    % preallocate session time table
    session_time_table = [];
    
    % Loop over sessions
    for ses = 1:size(all_sessions,3)
        fprintf('SESSION %s..\n',session_names{ses})
        
        % preallocate subject time table
        subj_time_table = [];
        
        % Loop over runs
        for rr = 1:runs_per_session(ses)
            
            % copy condition order table headers and scrub content
            sz = [tbl_nrows size(params.trials(1,:),2)];
            time_table = vcd_preallocateNaNTable(sz(1), sz(2), params.trials(1,:), []);

            % add column for event timing and id
            time_table.event_start = NaN(tbl_nrows,1);
            time_table.event_dur   = NaN(tbl_nrows,1);
            time_table.event_end   = NaN(tbl_nrows,1);
            time_table.event_id    = NaN(tbl_nrows,1);
            time_table.event_name  = repmat({''},tbl_nrows,1);
            time_table.subj_nr     = sj.*ones(tbl_nrows,1);
            time_table.subj_block_nr = NaN(tbl_nrows,1);
            time_table.global_block_nr = NaN(tbl_nrows,1);

            % Find corresponding trials (rows in table) and block nrs for
            % this particular subject, run, and session. The t_trial table
            % will be inserted into a bigger time table that includes other
            % event types, like iti, cues, etc.
            idx0      = (params.trials.session_nr==ses & params.trials.run_nr==rr);
            t_trial   = params.trials(idx0,:);
            block_nrs = t_trial.subj_block_nr(:,sj);
            
            % Reorder trials according to subject's specific order
            [block_nrs,block_nrs_idx] = sort(block_nrs,'ascend');
            t_trial = t_trial(block_nrs_idx,:);
            
            [~,tmp] = unique(block_nrs,'stable');
            tmp_durs = [params.exp.block.total_single_epoch_dur, params.exp.block.total_double_epoch_dur];
            block_dur = tmp_durs(t_trial.trial_type(tmp));
            
            % reset counter
            total_run_frames = 0;
            run_finished     = 0;
            
            % Add eyetracking block and pre-blank period 
            if total_run_frames == 0
                
                table_idx = 1;
                % Add prefixation period (1s)
                time_table.event_start(table_idx)    = total_run_frames;
                time_table.event_dur(table_idx)      = params.exp.block.eye_gaze_fix0;
                time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx)-1;
                time_table.event_id(table_idx)       = params.exp.block.eye_gaze_fix_ID;
                time_table.event_name(table_idx)     = {'et_fix'};
                time_table.run_nr(table_idx)         = rr;
                time_table.session_nr(table_idx)     = ses;
                time_table.session_name(table_idx)   = session_names(ses);
                time_table.block_nr(table_idx)       = 999;
                time_table.subj_block_nr(table_idx)  = 999;
                time_table.block_ID(table_idx)       = 999;
                time_table.block_local_trial_nr(table_idx) = 0;
                time_table.block_ID(table_idx)        = 0;
                time_table.stim_class(table_idx)      = NaN;
                time_table.task_class(table_idx)      = NaN;
                time_table.global_block_nr(table_idx) = ses;
                total_run_frames = time_table.event_end(table_idx) + 1;
                
                % update time table idx
                table_idx = table_idx+1;
                
                % Add saccade events (1.2s each)
                sac_IDs = NaN(1,params.exp.block.nr_of_saccades);
                insert_fix_target = randi([2,params.exp.block.nr_of_saccades-1],1);
                sac_IDs(insert_fix_target) = params.exp.block.eye_gaze_sac_target_ID(1);
                sac_IDs(isnan(sac_IDs)) = shuffle_concat(params.exp.block.eye_gaze_sac_target_ID(2:end),1);
                
                for sac = 1:params.exp.block.nr_of_saccades
                    time_table.event_start(table_idx)    = total_run_frames;
                    time_table.event_dur(table_idx)      = params.exp.block.eye_gaze_sac_target;
                    time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx)-1;
                    time_table.event_id(table_idx)       = sac_IDs(sac);
                    time_table.event_name(table_idx)     = {'et_sac'};
                    time_table.run_nr(table_idx)         = rr;
                    time_table.session_nr(table_idx)     = ses;
                    time_table.session_name(table_idx)   = session_names(ses);
                    time_table.block_nr(table_idx)       = 999;
                    time_table.subj_block_nr(table_idx)  = 999;
                    time_table.block_ID(table_idx)       = 999;
                    time_table.block_local_trial_nr(table_idx) = 0;
                    time_table.global_block_nr(table_idx) = ses;
                    time_table.block_ID(table_idx)        = 0;
                    time_table.stim_class(table_idx)      = NaN;
                    time_table.task_class(table_idx)      = NaN;
                    total_run_frames = time_table.event_end(table_idx) + 1;
                    
                    % update time table idx
                    table_idx = table_idx+1;
                end
                
                % Add post saccade fixation period (1s)
                time_table.event_start(table_idx)    = total_run_frames;
                time_table.event_dur(table_idx)      = params.exp.block.eye_gaze_fix1;
                time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx)-1;
                time_table.event_id(table_idx)       = params.exp.block.eye_gaze_fix_ID;
                time_table.event_name(table_idx)     = {'et_fix'};
                time_table.run_nr(table_idx)         = rr;
                time_table.session_nr(table_idx)     = ses;
                time_table.session_name(table_idx)   = session_names(ses);
                time_table.block_ID(table_idx)       = 999;
                time_table.subj_block_nr(table_idx)  = 999;
                time_table.block_nr(table_idx)       = 999;
                time_table.block_local_trial_nr(table_idx) = 0;
                time_table.block_ID(table_idx)        = 0;
                time_table.stim_class(table_idx)      = NaN;
                time_table.task_class(table_idx)      = NaN;
                time_table.global_block_nr(table_idx) = ses;
                total_run_frames = time_table.event_end(table_idx) + 1;
                
                % update time table idx
                table_idx = table_idx+1;

                % Add pupil trial (3s black, 1s white)
                for ppl = 1:length(params.exp.block.eye_gaze_pupil)
                    time_table.event_start(table_idx)    = total_run_frames;
                    time_table.event_dur(table_idx)      = params.exp.block.eye_gaze_pupil(ppl);
                    time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx)-1;
                    time_table.event_id(table_idx)       = params.exp.block.eye_gaze_pupil_ID;
                    time_table.event_name(table_idx)     = {['et_pupil' num2str(ppl)]};
                    time_table.run_nr(table_idx)         = rr;
                    time_table.session_nr(table_idx)     = ses;
                    time_table.session_name(table_idx)   = session_names(ses);
                    time_table.block_ID(table_idx)       = 999;
                    time_table.subj_block_nr(table_idx)  = 999;
                    time_table.block_nr(table_idx)       = 999;
                    time_table.block_local_trial_nr(table_idx) = 0;
                    time_table.global_block_nr(table_idx) = ses;
                    time_table.block_ID(table_idx)        = 0;
                    time_table.stim_class(table_idx)      = NaN;
                    time_table.task_class(table_idx)      = NaN;
                    total_run_frames = time_table.event_end(table_idx) + 1;
                     
                    % update time table idx
                    table_idx = table_idx+1;
                end
                
                % Add pre-blank period                 
                time_table.event_start(table_idx)    = total_run_frames;
                time_table.event_dur(table_idx)      = session_preblankdur;
                time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx)-1;
                time_table.event_id(table_idx)       = 0;
                time_table.event_name(table_idx)     = {'pre-blank'};
                
                time_table.run_nr(table_idx)         = rr;
                time_table.session_nr(table_idx)     = ses;
                time_table.session_name(table_idx)   = session_names(ses);
                time_table.block_nr(table_idx)       = 0;
                time_table.subj_block_nr(table_idx)  = 0;
                time_table.block_ID(table_idx)       = 0;
                time_table.block_local_trial_nr(table_idx) = 0;
                time_table.block_ID(table_idx)        = 0;
                time_table.stim_class(table_idx)      = 0;
                time_table.task_class(table_idx)      = 0;
                total_run_frames = time_table.event_end(table_idx) + 1;
                
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
                
                curr_trial = t_trial.block_local_trial_nr(1);
                
                % get current block nr
                curr_block = block_nrs(block_vec_idx);
                curr_global_session_block_nr  = t_trial.global_block_nr(1);
                curr_subj_within_run_block_nr = t_trial.subj_block_nr(1,sj);
                curr_block_ID                 = t_trial.block_ID(1);
                curr_trial_nr                 = t_trial.unique_trial_nr(1);
                curr_stim_class               = t_trial.stim_class(1);
                curr_task_class               = t_trial.task_class(1);
                curr_stim_class_name          = t_trial.stim_class_name(1,:);
                curr_task_class_name          = t_trial.task_class_name(1);
                curr_repeat_nr                = t_trial.repeat_nr(1); 
                curr_stim_class_unique_block_nr = t_trial.stim_class_unique_block_nr(1);
                curr_trial_type               = t_trial.trial_type(1);
                curr_block_name               = t_trial.block_name(1);
                curr_block_local_trial_nr     = t_trial.block_local_trial_nr(1);
                curr_spatial_cue              = t_trial.is_cued(1);

                
                
                % first trial of the block has a task cue
                if t_trial.block_local_trial_nr(1) == 1 || ...
                        (block_vec_idx == 1) || (curr_block ~= block_nrs(block_vec_idx-1)) % new block
                    time_table.event_start(table_idx)    = total_run_frames;
                    time_table.event_dur(table_idx)      = params.exp.block.task_cue_dur;
                    time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx) - 1;
                    time_table.event_id(table_idx)       = params.exp.block.task_cue_ID;
                    time_table.event_name(table_idx)     = {'task-cue'};
                    
                    % add ses, run nr, etc.
                    time_table.session_nr(table_idx)      = ses;
                    time_table.session_name(table_idx)    = session_names(ses);
                    time_table.run_nr(table_idx)          = rr;
                    time_table.block_nr(table_idx)        = curr_subj_within_run_block_nr;
                    time_table.global_block_nr(table_idx) = curr_global_session_block_nr;
                    time_table.block_nr(table_idx)        = curr_subj_within_run_block_nr;
                    time_table.subj_block_nr(table_idx)   = curr_subj_within_run_block_nr;
                    time_table.block_local_trial_nr(table_idx) = curr_block_local_trial_nr;
                    time_table.block_ID(table_idx)        = curr_block_ID;
                    time_table.unique_trial_nr(table_idx) = curr_trial_nr;
                    time_table.stim_class(table_idx)      = curr_stim_class;
                    time_table.task_class(table_idx)      = curr_task_class;
                    time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                    time_table.task_class_name(table_idx) = curr_task_class_name;
                    time_table.repeat_nr(table_idx)       = curr_repeat_nr;
                    time_table.stim_class_unique_block_nr(table_idx) = curr_stim_class_unique_block_nr;
                    time_table.trial_type(table_idx)      = curr_trial_type;
                    time_table.block_name(table_idx)      = curr_block_name;
                    
                    total_run_frames = time_table.event_end(table_idx) + 1;
                    
                    table_idx = table_idx+1;
                end
                
                % Get trial IDs of current trial
                if t_trial.trial_type(1) ==1
                    trial_IDs = trial_ID_single_epoch;
                else
                    trial_IDs = trial_ID_double_epoch;
                end
                
                % Loop over events within trial
                for id = 1:length(trial_IDs)
                    
                    % event start is the same for all IDs
                    time_table.event_start(table_idx) = total_run_frames;
                    
                    % add ses, run nr
                    time_table.session_nr(table_idx)      = ses; 
                    time_table.run_nr(table_idx)          = rr;
                    time_table.global_block_nr(table_idx) = curr_global_session_block_nr;
                    time_table.block_nr(table_idx)        = curr_subj_within_run_block_nr;
                    time_table.subj_block_nr(table_idx)   = curr_subj_within_run_block_nr;
                    time_table.session_name(table_idx)    = session_names(ses);
                    time_table.block_local_trial_nr(table_idx) = curr_block_local_trial_nr;
                    time_table.block_ID(table_idx)        = curr_block_ID;
                    time_table.unique_trial_nr(table_idx) = curr_trial_nr;
                    time_table.stim_class(table_idx)      = curr_stim_class;
                    time_table.task_class(table_idx)      = curr_task_class;
                    time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                    time_table.task_class_name(table_idx) = curr_task_class_name;
                    time_table.repeat_nr(table_idx)       = curr_repeat_nr;
                    time_table.stim_class_unique_block_nr(table_idx) = curr_stim_class_unique_block_nr;
                    time_table.trial_type(table_idx)      = curr_trial_type;
                    time_table.block_name(table_idx)      = curr_block_name;
                    
                    % Add individual trial events
                    switch trial_IDs(id)
                        
                        case params.exp.block.trial_start_ID % 94 trial start (dot thickening)
                            time_table.event_dur(table_idx)      = params.exp.trial.start_cue_dur;
                            time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx) - 1;
                            time_table.event_id(table_idx)       = trial_IDs(id);
                            time_table.event_name(table_idx)     = {'trial-start'};
                            time_table.block_local_trial_nr(table_idx) = curr_block_local_trial_nr;
                            time_table.block_ID(table_idx)        = curr_block_ID; 
                            time_table.unique_trial_nr(table_idx) = curr_trial_nr;
                            time_table.stim_class(table_idx)      = curr_stim_class;
                            time_table.task_class(table_idx)      = curr_task_class;
                            time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                            time_table.task_class_name(table_idx) = curr_task_class_name;
                            time_table.repeat_nr(table_idx)       = curr_repeat_nr;
                            time_table.stim_class_unique_block_nr(table_idx) = curr_stim_class_unique_block_nr;
                            time_table.trial_type(table_idx)      = curr_trial_type;
                            time_table.block_name(table_idx)      = curr_block_name;
                            time_table.is_cued(table_idx)          = curr_spatial_cue;
                            
                            total_run_frames = time_table.event_end(table_idx) + 1;
                            table_idx   = table_idx+1; 
                            
                        case params.exp.block.spatial_cue_ID  % 95 spatial cue
                            time_table.event_dur(table_idx)      = params.exp.trial.spatial_cue_dur;
                            time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx) - 1;
                            time_table.event_id(table_idx)       = trial_IDs(id);
                            time_table.event_name(table_idx)     = {'spatial-cue'};
                            time_table.block_nr(table_idx)       = curr_block;
                            time_table.subj_block_nr(table_idx)  = curr_block;
                            time_table.block_local_trial_nr(table_idx) = curr_block_local_trial_nr;
                            time_table.block_ID(table_idx)        = curr_block_ID; 
                            time_table.unique_trial_nr(table_idx) = curr_trial_nr;
                            time_table.stim_class(table_idx)      = curr_stim_class;
                            time_table.task_class(table_idx)      = curr_task_class;
                            time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                            time_table.task_class_name(table_idx) = curr_task_class_name;
                            time_table.repeat_nr(table_idx)       = curr_repeat_nr;
                            time_table.stim_class_unique_block_nr(table_idx) = curr_stim_class_unique_block_nr;
                            time_table.trial_type(table_idx)      = curr_trial_type;
                            time_table.block_name(table_idx)      = curr_block_name;
                            time_table.is_cued(table_idx)          = curr_spatial_cue;
                            total_run_frames = time_table.event_end(table_idx) + 1;
                            table_idx   = table_idx+1;
                            
                        case params.exp.block.stim_epoch1_ID % 91 stim epoch 1 - left&right                           
                            tmp1 = t_trial(1,:);
                            tmp1.event_start    = total_run_frames(end);
                            tmp1.event_dur      = params.exp.trial.stim_array_dur;
                            tmp1.event_end      = tmp1.event_start(1) + tmp1.event_dur(1) - 1;
                            tmp1.event_id       = trial_IDs(id); 
                            tmp1.event_name     = {'stim1'};
                            tmp1.subj_nr        = sj;
                            tmp1.subj_block_nr  = curr_block;
                            tmp1.block_nr       = curr_block;
                            tmp1.session_name   = session_names(ses);
                            tmp1.session_nr     = ses;
                            tmp1.run_nr         = rr;
                            tmp1.is_cued         = curr_spatial_cue;
                            tmp1.global_block_nr = curr_global_session_block_nr;

                            time_table(table_idx,:) = tmp1;
                            total_run_frames = time_table.event_end(table_idx) + 1;
                            
                            table_idx   = table_idx+1;
                            
                            % if there is a second stim epoch, we wait with
                            % updating trial table. If single stim epoch,
                            % then we can delete trial
                            if t_trial.trial_type(1) == 1 
                                t_trial(1,:) = [];
                            end
                            
                        case params.exp.block.stim_epoch2_ID % stim epoch 2 (after delay) - left&right
                                tmp2 = tmp1;
                                tmp2.event_start    = total_run_frames(end);
                                tmp2.event_dur      = params.exp.trial.stim_array_dur;
                                tmp2.event_end      = tmp2.event_start(1) + tmp2.event_dur(1) - 1;
                                tmp2.event_id       = trial_IDs(id);
                                tmp2.event_name     = {'stim2'};
                                tmp2.stim_nr_left   = NaN;
                                tmp2.stim_nr_right  = NaN;
                                tmp2.subj_nr        = sj;
                                tmp2.subj_block_nr  = curr_block;
                                tmp2.block_nr       = curr_block;
                                tmp2.session_name   = session_names(ses);
                                tmp2.session_nr     = ses;
                                tmp2.run_nr         = rr;
                                tmp2.global_block_nr = curr_global_session_block_nr;
                                tmp2.is_cued         = curr_spatial_cue;                    
                                time_table(table_idx,:) = tmp2;
                                
                                total_run_frames = time_table.event_end(table_idx) + 1;
                                
                                table_idx = table_idx+1;
                                
                                % Now we update trial
                                t_trial(1,:) = [];
                            
                        case params.exp.block.response_ID % response
                            time_table.event_dur(table_idx)      = params.exp.trial.response_win_dur;
                            time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx) - 1;
                            time_table.event_id(table_idx)       = params.exp.block.response_ID;
                            time_table.event_name(table_idx)     = {'response-cue'};
                            time_table.block_nr(table_idx)       = curr_block;
                            time_table.subj_block_nr(table_idx)  = curr_block;
                            time_table.block_local_trial_nr(table_idx) = curr_block_local_trial_nr;
                            time_table.block_ID(table_idx)        = curr_block_ID;
                            time_table.unique_trial_nr(table_idx) = curr_trial_nr;
                            time_table.stim_class(table_idx)      = curr_stim_class;
                            time_table.task_class(table_idx)      = curr_task_class;
                            time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                            time_table.task_class_name(table_idx) = curr_task_class_name;
                            time_table.repeat_nr(table_idx)       = curr_repeat_nr;
                            time_table.stim_class_unique_block_nr(table_idx) = curr_stim_class_unique_block_nr;
                            time_table.trial_type(table_idx)      = curr_trial_type;
                            time_table.block_name(table_idx)      = curr_block_name;
                            time_table.is_cued(table_idx)          = curr_spatial_cue;
                            
                            total_run_frames = time_table.event_end(table_idx) + 1;
                            
                            table_idx = table_idx+1;
                            
                        case params.exp.block.delay_ID % delay
                            time_table.event_dur(table_idx)      = params.exp.trial.delay_dur;
                            time_table.event_end(table_idx)      = time_table.event_start(table_idx) + time_table.event_dur(table_idx) - 1;
                            time_table.event_id(table_idx)       = trial_IDs(id);
                            time_table.event_name(table_idx)     = {'delay'};
                            time_table.block_nr(table_idx)       = curr_block;
                            time_table.subj_block_nr(table_idx)  = curr_block;
                            time_table.block_local_trial_nr(table_idx) = curr_block_local_trial_nr;
                            time_table.block_ID(table_idx)        = curr_block_ID;
                            time_table.unique_trial_nr(table_idx) = curr_trial_nr;
                            time_table.stim_class(table_idx)      = curr_stim_class;
                            time_table.task_class(table_idx)      = curr_task_class;
                            time_table.stim_class_name(table_idx,:) = curr_stim_class_name;
                            time_table.task_class_name(table_idx) = curr_task_class_name;
                            time_table.repeat_nr(table_idx)       = curr_repeat_nr;
                            time_table.stim_class_unique_block_nr(table_idx) = curr_stim_class_unique_block_nr;
                            time_table.trial_type(table_idx)      = curr_trial_type;
                            time_table.block_name(table_idx)      = curr_block_name;
                            time_table.is_cued(table_idx)          = curr_spatial_cue;
                            
                            total_run_frames = time_table.event_end(table_idx) + 1;
                            
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

                        if ~exist('ibis','var') || isempty(ibis)
                            if ~exist('postblank_to_add','var')
                                [ibis, postblank_to_add] = ...
                                    vcd_optimizeIBIs(session_totalrundur-params.exp.block.total_eyetracking_block_dur, ...
                                    block_dur, IBI_to_use, length(block_dur), sum([session_postblankdur,session_preblankdur]));
                            else
                                ibis = ...
                                    vcd_optimizeIBIs(session_totalrundur-params.exp.block.total_eyetracking_block_dur, ...
                                    block_dur, IBI_to_use, length(block_dur), sum([session_postblankdur,session_preblankdur]));
                            end
                        end
                        
                        % add IBI
                        time_table.event_start(table_idx)  = total_run_frames;
                        IBI = ibis(1);
                        if ~nearZero(mod(time_table.event_start(table_idx)+IBI,1))
                             IBI = IBI + (1 - mod(time_table.event_start(table_idx)+IBI,1)); % round out to a full second
                        end
                        time_table.event_start(table_idx)  = total_run_frames;
                        time_table.event_dur(table_idx)    = IBI;
                        time_table.event_name(table_idx)   = {'IBI'};
                        time_table.event_id(table_idx)     = params.exp.block.IBI_ID;
                        time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx) - 1;
                        % add ses, block, run nr
                        time_table.session_nr(table_idx) = ses;
                        time_table.session_name(table_idx)   = session_names(ses);
                        time_table.run_nr(table_idx) = rr;
                        time_table.block_nr(table_idx) = 0;
                        time_table.subj_block_nr(table_idx)   = 0;
                        time_table.block_ID(table_idx)        = 0;
                        time_table.stim_class(table_idx)      = 0;
                        time_table.task_class(table_idx)      = 0;
                        time_table.is_cued(table_idx)          = 0;

                        total_run_frames = time_table.event_end(table_idx) + 1;
                        
                        ibis(1) = [];
                        table_idx = table_idx+1;
                        
                    elseif curr_block == next_block_id  
                        % next trial is from the same block (tt is already
                        % updated), so not last trial of the block, and we
                        % need to insert ITI
                        
                        if curr_trial == 1 
                            % Reset shuffled ITIs prior to allocating to
                            % ensure the right block length
                            if time_table.trial_type(table_idx-2) == 1 % single stim epoch
                                max_nr_trials = params.exp.block.n_trials_single_epoch;
                                max_block_dur = params.exp.block.total_single_epoch_dur;
                                max_trial_dur = params.exp.trial.single_epoch_dur;
                            else % double stim epoch
                                max_nr_trials = params.exp.block.n_trials_double_epoch;
                                max_block_dur = params.exp.block.total_double_epoch_dur;
                                max_trial_dur = params.exp.trial.double_epoch_dur;
                            end
                            
                            itis = vcd_optimizeITIs(max_block_dur-params.exp.block.task_cue_dur, max_trial_dur, params.exp.trial.ITI, max_nr_trials);
                        end

                        % add ITI
                        ITI =  itis(1); % grab the first one from the shuffled list
                        itis(1) = []; % and remove it
                        time_table.event_start(table_idx)  = total_run_frames;
                        time_table.event_dur(table_idx)    = ITI;
                        time_table.event_name(table_idx)   = {'ITI'};
                        time_table.event_id(table_idx)     = params.exp.block.ITI_ID;
                        time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx) - 1;
                        % add ses, block, run nr
                        time_table.session_nr(table_idx)      = ses;
                        time_table.session_name(table_idx)    = session_names(ses);
                        time_table.run_nr(table_idx)          = rr;
                        time_table.subj_nr(table_idx)         = sj;
                        time_table.block_nr(table_idx)        = curr_block;
                        time_table.subj_block_nr(table_idx)   = curr_block;
                        time_table.block_ID(table_idx)        = 0;
                        time_table.stim_class(table_idx)      = 0;
                        time_table.task_class(table_idx)      = 0;
                        time_table.is_cued(table_idx)          = 0;
                        
                        total_run_frames = time_table.event_end(table_idx) + 1;
                        
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
                
                % round out post blank period to finish on full TR
                total_run_time2 = total_run_frames + session_postblankdur + postblank_to_add;
                
                if total_run_time2 == session_totalrundur
                    post_blank_dur = session_postblankdur + postblank_to_add + 1;
                elseif total_run_time2 < session_totalrundur
                    round_me_out = (ceil(total_run_time2/ tr_in_frames) - (total_run_time2/tr_in_frames)).*tr_in_frames;
                    post_blank_dur = session_postblankdur + postblank_to_add + round_me_out;
                elseif total_run_time2 > session_totalrundur
                    post_blank_dur = session_totalrundur - total_run_frames + 1;
                end
                % add post_blank
                time_table.event_start(table_idx)  = total_run_frames;
                time_table.event_dur(table_idx)    = post_blank_dur;
                time_table.event_name(table_idx)   = {'post-blank'};
                time_table.event_id(table_idx)     = 0;
                time_table.event_end(table_idx)    = time_table.event_start(table_idx) + time_table.event_dur(table_idx) - 1;
                % add ses, run nr
                time_table.session_nr(table_idx) = ses;
                time_table.session_name(table_idx)   = session_names(ses);
                time_table.run_nr(table_idx) = rr;
                time_table.subj_nr(table_idx) = sj;
                time_table.block_nr(table_idx) = 0;
                time_table.subj_block_nr(table_idx) = 0;
                time_table.block_ID(table_idx)        = 0;
                time_table.stim_class(table_idx)      = 0;
                time_table.task_class(table_idx)      = 0;
                time_table.is_cued(table_idx)          = 0;
                total_run_frames = time_table.event_end(table_idx);

%                 if strcmp(session_type,'MRI')
%                     assert(nearZero(mod(total_run_frames,params.exp.TR)));
%                     assert(nearZero(mod(total_run_frames,session_totalrundur)))
%                 elseif strcmp(session_type,'BEHAVIOR')
%                     assert(total_run_frames <= session_totalrundur)
%                 end
                
                run_finished = 0; % reset flag

                % trim time table
                time_table2 = time_table;
                if strcmp(time_table2.event_name(table_idx),{'post-blank'})
                    time_table2((table_idx+1):end,:) = [];
                end
            
                % remove obselete columns
                time_table2.unique_trial_nr = [];
                time_table2.stim_class_unique_block_nr = []; 
                time_table2.block_name = [];
                
                % reorganize column order
                col_order = {'subj_nr','session_nr','session_name','run_nr','block_nr','block_local_trial_nr',...
                    'cond_name','stim_class','task_class',...
                    'stim_nr_left','stim_nr_right','is_cued',...
                    'is_catch','trial_type',...
                    'block_ID','event_start','event_dur','event_end',...
                    'event_id','event_name', ...
                    'stim_class_name','task_class_name',...
                    'orient_dir','contrast','gbr_phase','rdk_coherence',...
                    'super_cat','basic_cat','sub_cat','affordance_cat','super_cat_name',...
                    'basic_cat_name','sub_cat_name','affordance_name','is_in_img_ltm','stim2_im_nr','stim2_delta','stim2_orient_dir',...
                    'is_lure','repeat_nr','global_block_nr'};
                assert(all(ismember(col_order,time_table2.Properties.VariableNames)));
                for col = 1:length(col_order)
                    new_col_order(col) = find(strcmp(time_table2.Properties.VariableNames,col_order(col)));
                end
                time_table2 = time_table2(:,new_col_order);
                
                time_table2.is_lure(isnan(time_table2.is_lure(:,1)),1)=0;
                time_table2.is_lure(isnan(time_table2.is_lure(:,2)),2)=0;
                time_table2.is_catch(isnan(time_table2.is_catch))=0;
                time_table2.is_in_img_ltm(isnan(time_table2.is_in_img_ltm(:,1)),1)=0;
                time_table2.is_in_img_ltm(isnan(time_table2.is_in_img_ltm(:,2)),2)=0;
                % add run to master
                subj_time_table = cat(1, subj_time_table, time_table2);
               
                if params.verbose
                    time_table2.run_nr(isnan(time_table2.run_nr))=0;
                    time_table2.block_nr(isnan(time_table2.block_nr))=0;
                    time_table2.block_local_trial_nr(isnan(time_table2.block_local_trial_nr))=0;
                    time_table2.block_ID(isnan(time_table2.block_ID))=0;
                    
                    time_table2.block_nr( time_table2.block_nr==99) = 0.5;
                    
                    figure(99); clf; set(gcf, 'Position', [1 1 1200 500])
                    subplot(211); imagesc([time_table2.run_nr,time_table2.block_nr,time_table2.block_local_trial_nr]')
                    set(gca,'YTick',[1:3],'YTickLabel',{'run nr','block nr', 'trial nr'})
                    subplot(212); imagesc(time_table2.block_ID')
                    set(gca,'YTick',[1],'YTickLabel',{'block ID'})
                    colormap(cmapturbo(40))
                    xlabel('events (in time)')
                    sgtitle(sprintf('subj %03d, ses %02d',sj,ses))
                    
                    if params.store_imgs
                        saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master1_%s',session_type));
                        filename = sprintf('vcd_session%02d_subj%03d_run%02d_event_table.png', ses,sj,rr);
                        print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                    end
                end
                clear time_table2
            end
        end
        
        session_time_table = cat(1, session_time_table, subj_time_table);

    end
    time_table_master = cat(1,time_table_master,session_time_table);
end

if params.store_params
    fprintf('\n[%s]:Storing subject session data..\n',mfilename)
    saveDir = fileparts(fullfile(params.stim.fix.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir,sprintf('time_table_master_%s_%s_%s.mat',params.disp.name, session_type, datestr(now,30))),'time_table_master')
end

return