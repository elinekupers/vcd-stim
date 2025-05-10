function [time_table_master,all_subj_run_frames] = vcd_addFIXandCDtoTimeTableMaster(params, time_table_master)

%% ANON FUNCTIONS TO GET ONSET OF FIX AND CD

% Fixation order and fixation
fixsoafun = @() round(params.stim.fix.dotmeanchange);

% Contrast decrement gaussian time window onset
cdsoafun = @() round(params.stim.cd.meanchange + params.stim.cd.changeplusminus*(2*(rand-.5)));

%% Preallocate space for generated subject run frames and updated time table master
subject_nrs  = unique(time_table_master.subj_nr);
session_nrs  = unique(time_table_master.session_nr);
run_nrs      = unique(time_table_master.run_nr);

all_subj_run_frames = cell(length(subject_nrs),length(session_nrs),length(run_nrs));

time_table_master2 = [];


% loop over subjects
for sj = 1:length(subject_nrs)
    % loop over sessions
    for ses = 1:length(session_nrs)
        % loop over runs
        for rr = 1:length(run_nrs)
            
            % grab subj run trials from time_table_master
            subj_run = time_table_master((time_table_master.subj_nr==subject_nrs(sj) & ...
                time_table_master.session_nr==session_nrs(ses) & ...
                time_table_master.run_nr==run_nrs(rr)),:);
            
            % check if the last event is a post-blank
            assert(strcmp(subj_run.event_name(end),'post-blank'))
            
            % get run duration (in frames)
            run_dur = subj_run.event_end(end);
            
            % ensure a run is of reasonable length, and is not longer than 10 min
            assert((run_dur*params.stim.presentationrate_hz)/3600 < 600)
            
            % %%%%% GENERATE FIXATION SEQUENCE %%%%%
            fix_matrix = vcd_createFixationSequence(params,fixsoafun,run_dur);
            
            run = struct();
            run.fix_abs_lum = fix_matrix(:,2);
            run.fix_correct_response = fix_matrix(:,4);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % visualize fixation sequence
            if params.verbose
                makeprettyfigures
                figure(101); clf; set(gcf,'Position',[137,952,2424,600])
                sgtitle(sprintf('Fixation sequence: subject %03d, session %02d, run %02d', subject_nrs(sj),session_nrs(ses),run_nrs(rr)))
                subplot(211);
                plot([0:1:length(run.fix_abs_lum)-1].*params.stim.presentationrate_hz,run.fix_abs_lum,'ko-');
                xlabel('Time (s)'); ylabel('dot luminance');
                ylim([0 255]); xlim([0, (length(run.fix_abs_lum).*params.stim.presentationrate_hz)])
                title('fixation dot luminance sequence')
                
                subplot(212);
                plot([0:1:length(run.fix_correct_response)-1].*params.stim.presentationrate_hz,  run.fix_correct_response ,'ko-');
                title('fixation dot luminance rel diff')
                set(gca,'YTick', [0,1,2], 'YTickLabel', {'No Change','Brighter','Dimmer'})
                xlabel('Time (s)');
                xlim([0, (length(run.fix_correct_response).*params.stim.presentationrate_hz)])
                
                if params.store_imgs
                    saveFigsFolder = fullfile(vcd_rootPath,'figs');
                    filename = sprintf('vcd_subject%03d_session%02d_run%02d_fix_sequence.png', subject_nrs(sj),session_nrs(ses),run_nrs(rr));
                    print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Expand subject run time table
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            all_events = []; all_cued = []; all_catch = [];
            for jj = 1:length(subj_run.event_id)
                if ~isnan(subj_run.event_dur(jj)) && subj_run.event_dur(jj)~=0
                    all_events = cat(1, all_events, repmat(subj_run.event_id(jj),1,round(subj_run.event_dur(jj)))'); 
                    
                    all_cued  = cat(1, all_cued,  repmat(subj_run.is_cued(jj), 1, round(subj_run.event_dur(jj)))'); 
                    all_catch = cat(1, all_catch, repmat(subj_run.is_catch(jj),1, round(subj_run.event_dur(jj)))'); 
                end
            end
            run.frame_event_nr = all_events; clear all_events
            run.is_cued = all_cued; clear all_cued;
            run.is_catch = all_catch; clear all_catch;
            
            stim_idx    = ismember(subj_run.event_id, [91, 92,990:996]);
            stim_events = subj_run.event_id(stim_idx);
            stim_row    = find(stim_idx);
            
            subj_run.fix_lum   = cell(size(subj_run,1),1);
            subj_run.fix_start = cell(size(subj_run,1),1);
            total_nr_frames    = subj_run.event_end(end);
            run.frame_im_nr    = zeros(total_nr_frames,2);
            run.contrast       = ones(total_nr_frames,2);
            run.button_response_cd = zeros(total_nr_frames,2);
            run.is_cued        = zeros(total_nr_frames,2);
            
            % convert button response to cell vector so we can add
            % multiple responses per row.
            tmp_br = subj_run.correct_response;
            subj_run.correct_response = cell(size(subj_run,1),1);
            subj_run.correct_response = mat2cell(tmp_br, ones(size(tmp_br)));
            subj_run.cd_start         = NaN(size(subj_run,1),2);
            
            for ii = 1:length(stim_row)
                
                % get frame nrs
                nframes     = round(subj_run.event_dur(stim_row(ii)));
                frame_counter = subj_run.event_start(stim_row(ii))+1;
                            
                curr_frames = frame_counter:(frame_counter+nframes-1);
                
                % Get unique stimulus image nrs
                nsides = find(~cellfun(@isempty, subj_run.stim_class_name(stim_row(ii),:)));
                
                if ~isempty(nsides)
                    if isnan(subj_run.stim_nr_right(stim_row(ii)))
                        nsides = 1;
                    end
                
                    for side = nsides
                        if strcmp(subj_run.event_name(stim_row(ii),1),'stim1')
                            if side == 1
                                unique_im = subj_run.stim_nr_left(stim_row(ii));
                            elseif side == 2
                                unique_im = subj_run.stim_nr_right(stim_row(ii));
                            end
                        elseif strcmp(subj_run.event_name(stim_row(ii),1),'stim2')
                            unique_im = subj_run.stim2_im_nr(stim_row(ii),side);
                        else
                            unique_im = 0;
                        end

                        % fill in unique_im nr
                        run.frame_im_nr(curr_frames,side) = repmat(unique_im, length(curr_frames),1);
                    end
                end
               

    
                % STIMULI:  91 = Stim interval 1, 92 = Stim interval 2,
                if (stim_events(ii) == 91 ||stim_events(ii) == 92)
                    
                    % NOT A CATCH TRIAL
                    if ~subj_run.is_catch(stim_row(ii))
                        
                        % IF CONTRAST DECREMENT TASK BLOCK
                        if strcmp(subj_run.task_class_name(stim_row(ii)),'cd')

                            for side = nsides
                            
                                % reset cd_change to false
                                cd_change = false;
                            
                                % 50% change we will actually apply the contrast
                                % decrement change to stimulus
                                if rand(1)>params.stim.cd.prob
                                    cd_change = true;

                                    % Get onset of contrast decrement within the
                                    % stimulus period
                                    c_onset = feval(cdsoafun);

                                   
                                    f_cd = c_onset:(c_onset+length(params.stim.cd.t_gauss)-1);
                                    run.contrast(curr_frames(f_cd),side) = params.stim.cd.t_gauss;
                                    
                                    
                                    % log decrement in contrast column
                                    subj_run.cd_start(stim_row(ii),side) = c_onset;
                                
                                    % if this happens to be the cued side,
                                    % then log the correct response
                                    if isequal(subj_run.is_cued(stim_row(ii)),side)
                                        if cd_change
                                            % Add button response to start of decrement onset
                                            subj_run.correct_response{stim_row(ii)} = 1; % button 1: yes there was a change
                                        else
                                            % Add button response: no
                                            subj_run.correct_response{stim_row(ii)} = 2; % button 2: no there was a change (apply to all time points in the block)
                                        end
                                    end
                                end
                            
                            end
                        end

                    else % catch blocks
                        subj_run.correct_response{stim_row(ii)} = 0;
                    end % if catch
                    
                    % if response event is next, inherent response from stim event.
                    if subj_run.event_id(stim_row(ii)+1) == params.exp.block.response_ID  && isnan(subj_run.correct_response{stim_row(ii)+1})
                        subj_run.correct_response{stim_row(ii)+1} = subj_run.correct_response{stim_row(ii)};
                    end
                end
                % update frame counter
                frame_counter = curr_frames(end)+1;
            end % events
            
            
            
            % Log fixation changes as button presses in
            % subject's time_table_master
            fix_events = find(strcmp(subj_run.task_class_name,'fix'));
            if  ~isempty(fix_events)
                % given fixed interval and sampling without
                % replacement, fix change can only happen every 1.4 s
                % (42 frames) or 2.8 s (84 frames) in case we happen to
                % sample the same luminance twice when restarting the
                % sampling process of 5 lum values.
                assert(isequal(unique(diff(find(diff(run.fix_correct_response)>0)))', [params.stim.fix.dotmeanchange, 2*params.stim.fix.dotmeanchange]));
                
                fix_block_nrs = unique(subj_run.block_nr(fix_events,:))'; % should be 1 or 2 or 3 blocks per rum
                fix_update_idx = (run.fix_correct_response>0); % 1 x 20592 --> 234 fixation changes per run
                [~,fix_block_frames] = ismember(subj_run.block_nr,fix_block_nrs); % 40 trial events per block
                fix_block_change_direction = run.fix_correct_response(fix_update_idx); % 1=brighter, 2=dimmer
                fix_block_abs_lum = run.fix_abs_lum(fix_update_idx); %
                
                % nr of block frames should be exactly 1 or more single-stim presentation
                % block duration
                %                     assert(isequal(mod((max(subj_run.event_end(fix_block_frames>0))-min(subj_run.event_start(fix_block_frames>0)))+1,params.exp.block.total_single_epoch_dur),0))
                % get all the frames for the entire run
                all_frames = 0:total_nr_frames-1; % subtract one frame because table starts at t=0;
                % find those frames where the fixation circle updated
                time_frames_fix_updated = all_frames(fix_update_idx);
                fix_update_sub = find(fix_update_idx);
                for ff = 1:length(fix_update_sub)
                    t_fix = time_frames_fix_updated(ff); 
                    t_tbl = find((subj_run.event_start <= t_fix) & (subj_run.event_end >= t_fix));
                    
                    subj_run.fix_lum{t_tbl}   = cat(2, subj_run.fix_lum{t_tbl}, fix_block_abs_lum(ff));
                    subj_run.fix_start{t_tbl} = cat(2, subj_run.fix_start{t_tbl}, t_fix);
                end
                
                for ff_block = 1:length(fix_block_nrs)
                    fix_events2 = (fix_block_frames==ff_block);
                    % get start and end time frames for fixation block from time table master
                    fix_block_start_end = [min(subj_run.event_start(fix_events2)),  max(subj_run.event_end(fix_events2))];
                    % see what fix block frames overlap with frames where the fixation circle changed
                    fix_block_changes_idx = ((time_frames_fix_updated >= fix_block_start_end(1)) & (time_frames_fix_updated <= fix_block_start_end(2)));
                    fix_block_changes_times = time_frames_fix_updated(fix_block_changes_idx);
                    fix_block_changes_correct_response = fix_block_change_direction(fix_block_changes_idx);
                    
                    for ff = 1:length(fix_block_changes_correct_response)
                        t_fix = fix_block_changes_times(ff);
                        t_response = fix_block_changes_correct_response(ff);
                        
                        t_tbl = find((subj_run.event_start <= t_fix) & (subj_run.event_end >= t_fix));
                        if ~isempty(t_tbl)
                            if cellfun(@isnan, subj_run.correct_response(t_tbl))
                                subj_run.correct_response{t_tbl} = t_response; % replace nan with correct response
                            else
                                subj_run.correct_response{t_tbl} = cat(2, subj_run.correct_response(t_tbl), t_response); % or we add the second response
                            end
                        end
                    end
                end
            end
            
            % add timing
            run.timing = [0:1:length(run.frame_im_nr)-1]';

            % Add run_images and alpha_masks to larger cell array
            all_subj_run_frames{subject_nrs(sj),session_nrs(ses),run_nrs(rr)} = run;
            
            time_table_master2 = cat(1,time_table_master2,subj_run);
            
        end % runs
    end % sessions
end % subjects

time_table_master = time_table_master2;

% Store structs locally, if requested
if params.store_params
    fprintf('[%s]: Storing expanded time table for all subjects..\n',mfilename)
    saveDir = fullfile(vcd_rootPath,'workspaces','info');
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir, ...
        sprintf('time_table_master2_%s_%s.mat', ...
        params.disp.name, datestr(now,30))), ...
        'time_table_master','all_subj_run_frames','-v7.3')
end



end