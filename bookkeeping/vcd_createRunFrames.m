function [time_table_master,all_run_frames] = vcd_createRunFrames(params, time_table_master, env_type, varargin)
% VCD function to create an expanded version of the time_table_master at 
% the granularity of individual presentation frames. 
% 
%    [time_table_master,all_run_frames] = vcd_createRunFrames(params, time_table_master, env_type)
%
% This function takes the time_table_master and does two things:
% First, it expands several columns at the granularity of a single time 
% frame, where the frame rate defined by params.stim.presentationrate_hz, 
% see vcd_getStimulusParams.m). In the vcd-core experiment, this
% presentation rate is set to 60 Hz. This means that 1 one-second stimulus
% ON period will be coded as 60 rows in the all_run_frames table.
% 
% For PProom: time frames range from t=0 until t=22559.
% For demo, time frames range from t=0 until t=14999.
%
% Note that when we define time frames in the time_table_master, we use 
% 1-based indexing. This means that time_table_master.event_start = n
% corresponds to run_frames row n+1. One can refer to the time_table_master
% time in the run_frames.timing column.
% 
% For the event_end column in the time_table_master, we code the end of 
% last time frame (which is also the onset of the next frame, i.e., you are
% technically at a point in time in between monitor refreshes). 
% This means that one needs to subtract 1 time frame from time_table_event
% end, but given the 1-based indexing for the time_table_master,
% time_table_master.event_end = n, corresponds to run_frames row n as we 
% add one frame and subtract one frame, which cancels eachother out.
% 
% all_run_frames: 
%         {'session_nr'          }
%         {'session_type'        }
%         {'run_nr'              }
%         {'fix_abs_lum'         }
%         {'fix_correct_response'}
%         {'timing'              }
%         {'block_nr'            }
%         {'global_block_nr'     }
%         {'trial_nr'            }
%         {'global_trial_nr'     }
%         {'frame_event_nr'      }
%         {'is_cued'             }
%         {'is_catch'            }
%         {'is_objectcatch'      }
%         {'crossingIDs'         }
%         {'frame_im_nr'         }
%         {'contrast'            }
%         {'button_response_cd'  }


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'             , @isstruct);
p0.addRequired('time_table_master'  , @istable);
p0.addRequired('env_type'	        , @(x) ismember(x,{'MRI','BEHAVIOR'}));
p0.addParameter('store_params'      , true,  @islogical);
p0.addParameter('verbose'           , false, @islogical);
p0.addParameter('load_params'       , false, @islogical);
p0.addParameter('store_imgs'        , false, @islogical);


% Parse inputs
p0.parse(params, time_table_master, env_type, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0


% get session environment params
[~,session_types,~,~,~,~, ~, nr_session_types ] = vcd_getSessionEnvironmentParams(params, env_type);
        

% Preallocate space for generated subject run frames and updated time table master
session_nrs  = unique(time_table_master.session_nr);
run_nrs      = unique(time_table_master.run_nr);

% time frames at initial contrast level (next one will have dip)
pre_onset_time_frames = sum(params.stim.cd.t_cmodfun==1); 

all_run_frames     = [];
time_table_master2 = [];

% loop over sessions
for ses = 1:length(session_nrs)
    
    for st = 1:nr_session_types
        
        if ~isnan(session_types(ses,st))
            % loop over runs
            for rr = 1:length(run_nrs)
                
                % grab run trials from time_table_master
                this_run = time_table_master(time_table_master.session_nr==ses & ...
                    time_table_master.session_type==st & ...
                    time_table_master.run_nr==rr,:);
                
                % check if the last event is a post-blank
                assert(strcmp(this_run.event_name(end),'post-blank'))
                
                 % get run duration (in frames), minus one to get nr of
                 % frames
                run_dur = length(this_run.event_start(1):(this_run.event_end(end)-1));
                assert(isequal(run_dur,this_run.event_end(end)))
                
                % ensure a run is of reasonable length, and is not longer than 10 min
                assert((run_dur*params.stim.presentationrate_hz)/3600 < 600)
                
                % define time frames from t=0 until t=22559.
                % Note that the time_table_master will have the event_end of time frame
                % t=22559 as 22560, given that this is the end of last time frame (and
                % technically the start of the next time frame, if there was any..)
                t_frame = 0:(run_dur-1); 
                nr_frames_per_run = length(t_frame);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%% GENERATE FIXATION SEQUENCE %%%%%%%%%%%%%%%
                %
                % Note 1: every subject will get the same order.
                % Note 2: We freeze the fixation twinkle during rest periods.
                
                % Get rest / blank periods to know when to freeze the
                % fixation luminance sequence. 
                % 0: 'pre-blank', 90:'task-cue', 99:'IBI', 999: eyetracking block
                % Note that vcd_createFixationSequence will convert time
                % onset and offset to frame indices by blank_onset +1 and
                % blank_offset - 1;
                blank_onset = this_run.event_start(this_run.block_nr==999 | this_run.event_id==0 | this_run.event_id==90 | this_run.event_id==99); 
                blank_dur   = this_run.event_dur(this_run.block_nr==999 | this_run.event_id==0 | this_run.event_id==90 | this_run.event_id==99); % includes postblank period 
                
                % Get fixation sequence
                % fix_matrix is a matrix with dims: time frames x 4, where 
                % Column 1: time points (in frames)
                % Column 2: absolute luminance values (between 0 and 255)
                % Column 3: relative change in luminance values compared to the previous time point
                % Column 4: correct button press associated with the relative luminance change (1=brighter, 2=dimmer).
                fix_matrix = vcd_createFixationSequence(params,params.stim.fix.fixsoafun, run_dur, blank_onset, blank_dur); 
                
                % remove button responses from frozen fixation periods;
                if sum(fix_matrix(:,4)==0)>0
                    fix_matrix(fix_matrix(:,4)==0,4) = NaN;
                end
                
                % add fixation sequence to run_frames table
                run_frames = table();
                run_frames.session_nr           = repmat(ses,nr_frames_per_run,1);
                run_frames.session_type         = repmat(session_types(ses,st),nr_frames_per_run,1);
                run_frames.run_nr               = repmat(rr,nr_frames_per_run,1);
                run_frames.fix_abs_lum          = fix_matrix(:,2);
                run_frames.fix_correct_response = fix_matrix(:,4);
                run_frames.timing               = single(t_frame');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if ~verLessThan('matlab', '9.6')
                    % visualize fixation sequence
                    if verbose
                        makeprettyfigures;
                        
                        figure(101); clf; set(gcf,'Position',[1,1,1400,444])
                        sgtitle(sprintf('Fixation sequence: session %02d %s, run %02d', session_nrs(ses),choose(session_types(ses,st)==1,'A','B'),run_nrs(rr)))
                        subplot(211);
                        plot(fix_matrix(:,1)./params.stim.presentationrate_hz,run_frames.fix_abs_lum,'ko-');
                        xlabel('Time (s)'); ylabel('dot luminance');
                        ylim([0 255]); xlim([0, length(fix_matrix(:,1))/params.stim.presentationrate_hz])
                        title('fixation dot luminance sequence')
                        
                        subplot(212);
                        plot(fix_matrix(:,1)./params.stim.presentationrate_hz,  run_frames.fix_correct_response ,'ko-');
                        title('fixation dot luminance rel diff')
                        set(gca,'YTick', [0,1,2], 'YTickLabel', {'No Change','Brighter','Dimmer'})
                        xlabel('Time (s)');
                        xlim([0, length(fix_matrix(:,1))/params.stim.presentationrate_hz])
                        
                        if store_imgs
                            saveFigsFolder = fullfile(vcd_rootPath,'figs');
                            filename = sprintf('vcd_session%02d_%s_run%02d_fix_sequence.png', session_nrs(ses),choose(session_types(ses,st)==1,'A','B'),run_nrs(rr));
                            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                        end
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%% Expand run frames table %%%%%%%%%%%%%%%%%%
                % this makes it easier to loop through the frames when we
                % start flipping the screen
                all_events = []; all_cued = []; all_catch = []; all_crossings = []; all_objectcatch = [];
                all_block_nrs = []; all_trial_nrs = []; all_global_block_nrs = []; all_global_trial_nrs = [];
                for jj = 1:length(this_run.event_id)
                    if ~isnan(this_run.event_dur(jj)) && this_run.event_dur(jj)~=0
                        dur = this_run.event_dur(jj);
                        all_events      = cat(1, all_events,      repmat(this_run.event_id(jj),        dur,1));
                        all_cued        = cat(1, all_cued,        repmat(this_run.is_cued(jj),         dur,1));
                        all_catch       = cat(1, all_catch,       repmat(this_run.is_catch(jj),        dur,1));
                        all_objectcatch = cat(1, all_objectcatch, repmat(this_run.is_objectcatch(jj),  dur,1));
                        all_crossings   = cat(1, all_crossings,   repmat(this_run.crossing_nr(jj),     dur,1));
                        all_block_nrs   = cat(1, all_block_nrs,   repmat(this_run.block_nr(jj),        dur,1));
                        all_trial_nrs   = cat(1, all_trial_nrs,   repmat(this_run.trial_nr(jj),        dur,1));
                        all_global_trial_nrs = cat(1, all_global_trial_nrs, repmat(this_run.global_trial_nr(jj), dur,1));
                        all_global_block_nrs = cat(1, all_global_block_nrs, repmat(this_run.global_block_nr(jj), dur,1));
                    end
                end
                % If numeric vector: Change zero's to NaN
                all_cued(all_cued==0) = NaN;
                
                % If logical vector: Change  NaN's to zeros
                all_catch(isnan(all_catch))             = 0;
                all_objectcatch(isnan(all_objectcatch)) = 0;
                
                % Add events to run_frames struct
                run_frames.block_nr         = single(all_block_nrs);   clear all_block_nrs
                run_frames.global_block_nr  = single(all_global_block_nrs);   clear all_global_block_nrs
                run_frames.trial_nr         = single(all_trial_nrs);    clear all_trial_nrs
                run_frames.global_trial_nr  = single(all_global_trial_nrs);   clear all_global_trial_nrs
                run_frames.frame_event_nr   = single(all_events);       clear all_events
                run_frames.is_cued          = single(all_cued);         clear all_cued;
                run_frames.is_catch         = logical(all_catch);       clear all_catch;
                run_frames.is_objectcatch   = logical(all_objectcatch); clear all_objectcatch;
                run_frames.crossingIDs      = single(all_crossings);    clear all_crossings
                
                
                % use crossingIDs to remove any fixation related button responses outside fixation task blocks
                fix_crossing_nr = unique(this_run.crossing_nr(find(ismember(this_run.crossing_nr,find(~cellfun(@isempty, regexp(params.exp.crossingnames,'fix')))))));
                fix_block_on = [];
                for fc = 1:length(fix_crossing_nr)
                    fix_block_start = min(this_run.event_start(this_run.crossing_nr == fix_crossing_nr(fc)))+1; % note +1 for frame indexing
                    fix_block_end   = max(this_run.event_end(this_run.crossing_nr == fix_crossing_nr(fc))); % note -1 +1 cancel eachother to align with end of frame
                    fix_block_on    = cat(1,fix_block_on,fix_block_start:fix_block_end);
                end
                fix_block_off   = setdiff([1:run_dur],fix_block_on);
                assert(isempty(intersect(fix_block_on,fix_block_off)))
                run_frames.fix_correct_response(fix_block_off) = NaN;

                % find stimulus and event frames
                stim_idx    = ismember(this_run.event_id, [params.exp.block.stim_epoch1_ID, params.exp.block.stim_epoch2_ID, ...
                                params.exp.block.eye_gaze_fix_ID, params.exp.block.eye_gaze_sac_target_ID, ...
                                params.exp.block.eye_gaze_pupil_black_ID, params.exp.block.eye_gaze_pupil_white_ID]);
                stim_events = this_run.event_id(stim_idx);
                stim_row    = find(stim_idx);
                
                % Preallocate space for nr of fixation changes, frame nrs,
                % contrast levels, button_response related to contrast dip/
                % we use single class to save memory
                run_frames.frame_im_nr        = NaN(run_dur,2);
                run_frames.contrast           = ones(run_dur,2,'single'); % relative contrast (1 = mean contrast of image, <1 mean dip in contrast)
                run_frames.button_response_cd = NaN(run_dur,1);
                
                for ii = 1:length(stim_row)
                    
                    % Get frame nrs
                    nframes       = round(this_run.event_dur(stim_row(ii)));
                    timetable_start = this_run.event_start(stim_row(ii));
                    frametable_start = find(timetable_start==run_frames.timing);
                    
                    curr_frames   = frametable_start:(frametable_start+nframes-1);
                    
                    % Get unique stimulus image nrs fpr every frame
                    nsides = find(~cell2mat(cellfun(@(x) any(isnan(x)), this_run.stim_class_name(stim_row(ii),:), 'UniformOutput',false)));
                    
                    if ~isempty(nsides)
                        if isnan(this_run.stim_nr_right(stim_row(ii)))
                            nsides = 1;
                        end
                        
                        for side = nsides
                            if strcmp(this_run.event_name(stim_row(ii),1),'stim1')
                                if side == 1
                                    unique_im = this_run.stim_nr_left(stim_row(ii));
                                elseif side == 2
                                    unique_im = this_run.stim_nr_right(stim_row(ii));
                                end
                            elseif strcmp(this_run.event_name(stim_row(ii),1),'stim2')
                                unique_im = this_run.stim2_im_nr(stim_row(ii),side);
                            else
                                error('[%s]: Expecting a unique image but found none..',mfilename);
                            end
                            
                            % fill in unique_im nr
                            run_frames.frame_im_nr(curr_frames,side) = repmat(unique_im, length(curr_frames),1);
                        end
                        
                        
                        % ADD CD in frames
                        % only for STIMULI:  94 = Stim interval 1, 95 = Stim interval 2,
                        if (stim_events(ii) == params.exp.block.stim_epoch1_ID || stim_events(ii) == params.exp.block.stim_epoch2_ID)
                            
                            % NOT A CATCH TRIAL
                            if ~this_run.is_catch(stim_row(ii))
                                
                                % IF CONTRAST DECREMENT TASK BLOCK
                                if strcmp(this_run.task_class_name(stim_row(ii)),'cd')
                                    % 20% change we will actually apply the contrast
                                    % decrement change to the cued stimulus
                                    % (uncued stimulus will never change
                                    % contrast).
                                    c_onset   = this_run.cd_start(stim_row(ii)) + 1; % add one for frame indexing
                                    c_offset  = this_run.event_end(stim_row(ii)); % no additions: we add one for frame indexing, and remove one to align with stim offset
                                    if ~isnan(c_onset) && c_onset~=0
                                        cd_cued_side    = mod(this_run.is_cued(stim_row(ii))-1,2)+1;
                                        c_onset_support = c_onset - pre_onset_time_frames; % we shift the support function such that the first frame with a lower contrast aligns with cd_start
                                        f_cd            = c_onset_support:c_offset; % total frames with cd change + pre_onset_time support
                                        t_pad           = length(f_cd) - length(params.stim.cd.t_cmodfun);
                                        run_frames.contrast(f_cd,cd_cued_side) = cat(2,params.stim.cd.t_cmodfun, params.stim.cd.t_cmodfun(end)*ones(1,t_pad))';
                                        run_frames.button_response_cd(c_onset) = 1;
                                    else
                                        run_frames.button_response_cd(curr_frames) = 2;
                                    end
                                    
                                end
                            end % if catch
                        end % if cd trial
                    end % isempty(nsides)
                end % nr of stim events
                
                % Log fixation changes as button presses in
                % subject's time_table_master
                fix_events = find(strcmp(this_run.task_class_name,'fix'));
                if  ~isempty(fix_events)
                    % Given fixed soa and sampling without replacement,
                    % fix change can only happen every 1.4 s (42 frames) 
                    % or shorter in case we happen to run into an IBI.
                    assert(all(diff(find(abs(diff(run_frames.fix_abs_lum(run_frames.frame_event_nr>90 & run_frames.frame_event_nr<99)))>0)) <= params.stim.fix.dotmeanchange))
                    
                    fix_start_idx  = this_run.event_start(fix_events,:) + 1;    % when do trials in fixation block start 
                    fix_end_idx    = this_run.event_end(fix_events,:);          % when do trials in fixation block end (no additions because +1 for frame indexing and -1 to get the start of the last event frame cancels eachother out).
                    fix_update_idx = run_frames.fix_correct_response>0;   % when would the subject press a button to indicate luminance change          
                   
                    % find those frames where the fixation circle updated
                    % and tell the run table how many fixation changes we
                    % expect per trial row
                    fix_update_sub = find(fix_update_idx);
                    for ff = 1:length(fix_update_sub)
                        t_fix = fix_update_sub(ff);
                        t_tbl = find((fix_start_idx <= t_fix) & (t_fix <= fix_end_idx));
                        if ~isempty(t_tbl)
                            % if this is the first count, then redefine NaN
                            % as 1..
                            if isnan(this_run.nr_fix_changes(fix_events(t_tbl)))
                                this_run.nr_fix_changes(fix_events(t_tbl)) = 1;
                            else % is this is not the first count, then add one to the existing count
                                this_run.nr_fix_changes(fix_events(t_tbl)) = this_run.nr_fix_changes(fix_events(t_tbl)) + 1; % add another fixation change event to existing counter 
                            end
                        end
                    end
                    
                end

                % set non-fixation periods to NaN
                fix_block_nrs = unique(this_run.block_nr(fix_events));
                for fb = 1:length(fix_block_nrs)
                    non_fix_periods = setdiff([1:size(this_run,1)],min(find(this_run.block_nr==fix_block_nrs(fb))):max(find(this_run.block_nr==fix_block_nrs(fb))));
                    this_run.nr_fix_changes(non_fix_periods) = NaN;
                end
                % Add run_images and alpha_masks to larger cell array
                all_run_frames = cat(1,all_run_frames,run_frames);
                
                time_table_master2 = cat(1,time_table_master2,this_run);
                
            end % runs
        end % if nan session type
    end % session types
end % sessions

time_table_master = time_table_master2;


end