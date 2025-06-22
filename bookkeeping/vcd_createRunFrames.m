function [time_table_master,all_run_frames] = vcd_createRunFrames(params, time_table_master, env_type, varargin)
% [WRITE ME]
% 
%    [time_table_master,all_run_frames] = vcd_createRunFrames(params, time_table_master, env_type)
%
%
%
%

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

% define time frames from t=0 until t=22559.
% Note that the time_table_master will have the event_end of time frame 
% t=22559 as 22560, given that this is the end of last time frame (and
% technically the start of the next time frame, if there was any..)

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
                run_dur = this_run.event_end(end);
                
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
                blank_onset       = this_run.event_start(this_run.block_nr==999 | this_run.event_id==0 | this_run.event_id==90 | this_run.event_id==99); 
                blank_offset      = this_run.event_end(this_run.block_nr==999 | this_run.event_id==0 | this_run.event_id==90 | this_run.event_id==99); % includes postblank period 
                
                % Get fixation sequence
                % fix_matrix is a matrix with dims: time frames x 4, where 
                % Column 1: time points (in frames)
                % Column 2: absolute luminance values (between 0 and 255)
                % Column 3: relative change in luminance values compared to the previous time point
                % Column 4: correct button press associated with the relative luminance change (1=brighter, 2=dimmer).
                fix_matrix = vcd_createFixationSequence(params,params.stim.fix.fixsoafun, run_dur, blank_onset, blank_offset); 
                
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
                for jj = 1:length(this_run.event_id)
                    if ~isnan(this_run.event_dur(jj)) && this_run.event_dur(jj)~=0
                        dur = this_run.event_dur(jj);
                        all_events      = cat(1, all_events,      repmat(this_run.event_id(jj),        dur,1));
                        all_cued        = cat(1, all_cued,        repmat(this_run.is_cued(jj),         dur,1));
                        all_catch       = cat(1, all_catch,       repmat(this_run.is_catch(jj),        dur,1));
                        all_objectcatch = cat(1, all_objectcatch, repmat(this_run.is_objectcatch(jj),  dur,1));
                        all_crossings   = cat(1, all_crossings,   repmat(this_run.crossing_nr(jj),     dur,1));
                    end
                end
                % If numeric vector: Change zero's to NaN
                all_cued(all_cued==0) = NaN;
                
                % If logical vector: Change  NaN's to zeros
                all_catch(isnan(all_catch))             = 0;
                all_objectcatch(isnan(all_objectcatch)) = 0;
                
                % Add events to run_frames struct
                run_frames.frame_event_nr = single(all_events);       clear all_events
                run_frames.is_cued        = single(all_cued);         clear all_cued;
                run_frames.is_catch       = logical(all_catch);       clear all_catch;
                run_frames.is_objectcatch = logical(all_objectcatch); clear all_objectcatch;
                run_frames.crossingIDs    = single(all_crossings);    clear all_crossings
                
                % find stimulus and event frames
                stim_idx    = ismember(this_run.event_id, [params.exp.block.stim_epoch1_ID, params.exp.block.stim_epoch2_ID, ...
                                params.exp.block.eye_gaze_fix_ID, params.exp.block.eye_gaze_sac_target_ID, ...
                                params.exp.block.eye_gaze_pupil_black_ID, params.exp.block.eye_gaze_pupil_white_ID]);
                stim_events = this_run.event_id(stim_idx);
                stim_row    = find(stim_idx);
                
                % Preallocate space for nr of fixation changes, frame nrs,
                % contrast levels, button_response related to contrast dip/
                % we use single class to save memory
                this_run.nr_of_fix_changes    = zeros(size(this_run,1),1,'single');
                run_frames.frame_im_nr        = zeros(run_dur,2,'single');
                run_frames.contrast           = ones(run_dur,2,'single');
                run_frames.button_response_cd = zeros(run_dur,1,'single');
                
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
                                unique_im = 0;
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
                                    pre_onset_time_frames = sum(params.stim.cd.t_cmodfun==1); % time frames at initial contrast level (next one will have dip)
                                    % 20% change we will actually apply the contrast
                                    % decrement change to the cued stimulus
                                    % (uncued stimulus will never change
                                    % contrast).
                                    c_onset      = this_run.cd_start(stim_row(ii));
                                    if ~isnan(c_onset) && c_onset~=0
                                        cd_cued_side    = mod(this_run.is_cued(stim_row(ii))-1,2)+1;
                                        c_onset_support = c_onset - pre_onset_time_frames; % we shift the support function to an earlier time point such that the disp occurs at c_onset
                                        f_cd            = c_onset_support:this_run.event_end(stim_row(ii));
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
                    % given fixed interval and sampling without
                    % replacement, fix change can only happen every 1.4 s
                    % (42 frames) or 2.8 s (84 frames) in case we happen to
                    % sample the same luminance twice when restarting the
                    % sampling process of 5 lum values.
                    % assert(isequal(unique(diff(find(diff(run.fix_correct_response)>0)))', [params.stim.fix.dotmeanchange, 2*params.stim.fix.dotmeanchange]));
                    
                    fix_block_nrs  = unique(this_run.block_nr(fix_events,:))';        % should be 1 or 2 or 3 blocks per run
                    fix_update_idx = (run_frames.fix_correct_response>0);             % 1 x 22560 --> 31 or 43 fixation changes per block
                    [~,fix_block_frames] = ismember(this_run.block_nr,fix_block_nrs); % 40 trial events per block
                    fix_block_change_direction = run_frames.fix_correct_response(fix_update_idx); % 1=brighter, 2=dimmer
                    fix_block_abs_lum = run_frames.fix_abs_lum(fix_update_idx); %
                   
                   
                    % find those frames where the fixation circle updated
                    time_frames_fix_updated = t_frame(fix_update_idx);
                    fix_update_sub = find(fix_update_idx);
                    for ff = 1:length(fix_update_sub)
                        t_fix = time_frames_fix_updated(ff);
                        t_tbl = find((this_run.event_start <= t_fix) & (this_run.event_end >= t_fix));
                        this_run.nr_of_fix_changes(t_tbl) = this_run.nr_of_fix_changes(t_tbl) + 1;
                    end
                    
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