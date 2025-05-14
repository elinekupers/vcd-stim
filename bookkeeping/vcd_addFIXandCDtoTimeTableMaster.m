function [time_table_master,all_run_frames] = vcd_addFIXandCDtoTimeTableMaster(params, time_table_master, session_type)

%% ANON FUNCTIONS TO GET ONSET OF FIX AND CD

% Fixation order and fixation
fixsoafun = @() round(params.stim.fix.dotmeanchange);

% check session type
if strcmp(session_type,'MRI')
    all_sessions     = cat(3, params.exp.session.wide.ses_blocks, params.exp.session.deep.ses_blocks);
    runs_per_session = cat(1, size(params.exp.session.wide.ses_blocks,3),size(params.exp.session.deep.ses_blocks,3));
    session_types    = cat(1, params.exp.session.mri.wide.session_types,params.exp.session.mri.deep.session_types);
    nr_session_types = 2;
    
elseif strcmp(session_type,'BEHAVIOR')
    all_sessions         = params.exp.session.behavior.ses_blocks;
    runs_per_session     = params.exp.session.behavior.n_runs_per_session;
    session_types        = params.exp.session.behavior.session_types;
    nr_session_types     = 1;
end

%% Preallocate space for generated subject run frames and updated time table master
session_nrs  = unique(time_table_master.session_nr);
run_nrs      = unique(time_table_master.run_nr);

all_run_frames     = [];
time_table_master2 = [];

% loop over sessions
for ses = 1:length(session_nrs)
    
    for st = 1:nr_session_types
        
        if ~isnan(session_types(ses,st))
            % loop over runs
            for rr = 1:length(run_nrs)
                
                % grab run trials from time_table_master
                this_run = time_table_master(time_table_master.session_nr==session_nrs(ses) & ...
                    time_table_master.session_type==session_types(ses,st) & ...
                    time_table_master.run_nr==run_nrs(rr),:);
                
                % check if the last event is a post-blank
                assert(strcmp(this_run.event_name(end),'post-blank'))
                
                % get run duration (in frames)
                run_dur = this_run.event_end(end)+1;
                
                % ensure a run is of reasonable length, and is not longer than 10 min
                assert((run_dur*params.stim.presentationrate_hz)/3600 < 600)
                
                % %%%%% GENERATE FIXATION SEQUENCE %%%%%
                % every subject will get the same order.
                blank_onset  = this_run.event_start(this_run.block_nr==0);
                blank_onset(1) = 1; % include eyetracking block
                blank_offset = this_run.event_end(this_run.block_nr==0); % includes postblank period 
                blank_offset(end) = blank_offset(end) + 1;
                fix_matrix = vcd_createFixationSequence(params,fixsoafun, run_dur, blank_onset,blank_offset);
                
                % remove button responses from frozen fixation periods;
                for ii = 1:length(blank_onset)
                    fix_matrix(blank_onset(ii):blank_offset(ii),4) = 0;
                end
                
                
                run_frames = table();
                run_frames.session_nr           = repmat(ses,run_dur,1);
                run_frames.session_type         = repmat(session_types(ses,st),run_dur,1);
                run_frames.run_nr               = repmat(rr,run_dur,1);
                run_frames.fix_abs_lum          = fix_matrix(:,2);
                run_frames.fix_correct_response = fix_matrix(:,4);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % visualize fixation sequence
                if params.verbose
                    makeprettyfigures
                    figure(101); clf; set(gcf,'Position',[137,952,2424,600])
                    sgtitle(sprintf('Fixation sequence: session %02d, run %02d', session_nrs(ses),run_nrs(rr)))
                    subplot(211);
                    plot([0:1:length(run_frames.fix_abs_lum)-1].*params.stim.presentationrate_hz,run_frames.fix_abs_lum,'ko-');
                    xlabel('Time (s)'); ylabel('dot luminance');
                    ylim([0 255]); xlim([0, (length(run_frames.fix_abs_lum).*params.stim.presentationrate_hz)])
                    title('fixation dot luminance sequence')
                    
                    subplot(212);
                    plot([0:1:length(run_frames.fix_correct_response)-1].*params.stim.presentationrate_hz,  run_frames.fix_correct_response ,'ko-');
                    title('fixation dot luminance rel diff')
                    set(gca,'YTick', [0,1,2], 'YTickLabel', {'No Change','Brighter','Dimmer'})
                    xlabel('Time (s)');
                    xlim([0, (length(run_frames.fix_correct_response).*params.stim.presentationrate_hz)])
                    
                    if params.store_imgs
                        saveFigsFolder = fullfile(vcd_rootPath,'figs');
                        filename = sprintf('vcd_session%02d_run%02d_fix_sequence.png', session_nrs(ses),run_nrs(rr));
                        print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Expand subject run time table
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                all_events = []; all_cued = []; all_catch = []; all_crossings = [];
                for jj = 1:length(this_run.event_id)
                    if ~isnan(this_run.event_dur(jj)) && this_run.event_dur(jj)~=0
                        all_events = cat(1, all_events, repmat(this_run.event_id(jj), round(this_run.event_dur(jj)),1));
                        all_cued  = cat(1, all_cued,  repmat(this_run.is_cued(jj),   round(this_run.event_dur(jj)),1));
                        all_catch = cat(1, all_catch, repmat(this_run.is_catch(jj), round(this_run.event_dur(jj)),1));
                        all_crossings = cat(1, all_crossings, repmat(this_run.crossing_nr(jj), round(this_run.event_dur(jj)),1));
                    end
                end
                run_frames.frame_event_nr = single(all_events); clear all_events
                run_frames.is_cued  = single(all_cued); clear all_cued;
                run_frames.is_catch = logical(all_catch); clear all_catch;
                run_frames.crossingIDs = single(all_crossings); clear all_crossings
                
                stim_idx    = ismember(this_run.event_id, [params.exp.block.stim_epoch1_ID, params.exp.block.stim_epoch2_ID, ...
                    params.exp.block.eye_gaze_fix_ID, params.exp.block.eye_gaze_sac_target_ID, ...
                    params.exp.block.eye_gaze_pupil_black_ID, params.exp.block.eye_gaze_pupil_white_ID]);
                stim_events = this_run.event_id(stim_idx);
                stim_row    = find(stim_idx);
                
                %% Index stim cell vector with monotonic counter
                this_run.nr_of_fix_changes = zeros(size(this_run,1),1,'single');

                run_frames.frame_im_nr    = zeros(run_dur,2,'single');
                run_frames.contrast       = ones(run_dur,2,'single');
                run_frames.button_response_cd = zeros(run_dur,1,'single');
                
                for ii = 1:length(stim_row)
                    
                    % Get frame nrs
                    nframes       = round(this_run.event_dur(stim_row(ii)));
                    frame_counter = this_run.event_start(stim_row(ii))+1;
                    curr_frames   = frame_counter:(frame_counter+nframes-1);
                    
                    % Get unique stimulus image nrs fpr every frame 
                    nsides = find(~cellfun(@isempty, this_run.stim_class_name(stim_row(ii),:)));
                    
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
                    end

                    % ADD CD in frames
                    % only for STIMULI:  94 = Stim interval 1, 95 = Stim interval 2,
                    if (stim_events(ii) == params.exp.block.stim_epoch1_ID ||stim_events(ii) == params.exp.block.stim_epoch2_ID)
                        
                        % NOT A CATCH TRIAL
                        if ~this_run.is_catch(stim_row(ii))
                            % IF CONTRAST DECREMENT TASK BLOCK
                            if strcmp(this_run.task_class_name(stim_row(ii)),'cd')
                                for side = nsides
                                    % 50% change we will actually apply the contrast
                                    % decrement change to stimulus
                                    c_onset = this_run.cd_start(stim_row(ii),side);

                                    if ~isnan(c_onset) && c_onset~=0
                                        f_cd = c_onset:(c_onset+length(params.stim.cd.t_gauss)-1);
                                        run_frames.contrast(f_cd,side) = params.stim.cd.t_gauss;
                                        run_frames.button_response_cd(f_cd) = 1;
                                    else
                                        run_frames.button_response_cd(curr_frames) = 2;
                                    end
                                end
                            end
                        end % if catch
                    end
                end % events
                
                
                % Log fixation changes as button presses in
                % subject's time_table_master
                all_frames = 0:run_dur-1; % subtract one frame because table starts at t=0;
                fix_events = find(strcmp(this_run.task_class_name,'fix'));
                if  ~isempty(fix_events)
                    % given fixed interval and sampling without
                    % replacement, fix change can only happen every 1.4 s
                    % (42 frames) or 2.8 s (84 frames) in case we happen to
                    % sample the same luminance twice when restarting the
                    % sampling process of 5 lum values.
                    %                 assert(isequal(unique(diff(find(diff(run.fix_correct_response)>0)))', [params.stim.fix.dotmeanchange, 2*params.stim.fix.dotmeanchange]));
                    
                    fix_block_nrs = unique(this_run.block_nr(fix_events,:))'; % should be 1 or 2 or 3 blocks per rum
                    fix_update_idx = (run_frames.fix_correct_response>0); % 1 x 20592 --> 234 fixation changes per run
                    [~,fix_block_frames] = ismember(this_run.block_nr,fix_block_nrs); % 40 trial events per block
                    fix_block_change_direction = run_frames.fix_correct_response(fix_update_idx); % 1=brighter, 2=dimmer
                    fix_block_abs_lum = run_frames.fix_abs_lum(fix_update_idx); %
                   
                   
                    % find those frames where the fixation circle updated
                    time_frames_fix_updated = all_frames(fix_update_idx);
                    fix_update_sub = find(fix_update_idx);
                    for ff = 1:length(fix_update_sub)
                        t_fix = time_frames_fix_updated(ff);
                        t_tbl = find((this_run.event_start <= t_fix) & (this_run.event_end >= t_fix));
                        this_run.nr_of_fix_changes(t_tbl) = this_run.nr_of_fix_changes(t_tbl) + 1;
                    end
                    
                end
                
                % add timing
                run_frames.timing = single(all_frames');

                % Add run_images and alpha_masks to larger cell array
                all_run_frames = cat(1,all_run_frames,run_frames);
                
                time_table_master2 = cat(1,time_table_master2,this_run);
                
            end % runs
        end % if nan session type
    end % session types
end % sessions

time_table_master = time_table_master2;

% Store structs locally, if requested
if params.store_params
    fprintf('[%s]: Storing expanded time table for all subjects..\n',mfilename)
    saveDir = fullfile(vcd_rootPath,'workspaces','info');
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir, ...
        sprintf('time_table_master_complete_%s_%s.mat', ...
        params.disp.name, datestr(now,30))), ...
        'time_table_master','all_run_frames','-v7.3')
end



end