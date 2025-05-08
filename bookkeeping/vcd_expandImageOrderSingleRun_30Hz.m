function all_subj_run_frames = vcd_expandImageOrderSingleRun_30Hz(params, ...
    time_table_master,subject_nrs,session_nrs,run_nrs,varargin)
% Functions to expand sequence of events from individual images into 30Hz
% framelocked presentation rate.
%
%  subj_run = vcd_expandImageOrderSingleRun_30Hz(params, ...
%    time_table_master, subject_nrs, session_nrs, run_nrs, ...
%                       ['all_run_images'     , <all_run_images>], ...
%                       ['all_run_alpha_masks', <all_run_alpha_masks>], ...
%                       ['fixsoafun'      , <fixsoafun>], ...
%                       ['cdsoafun'       , <cdsoafun>], ...
%                       ['load_params'    , <load_params>], ...
%                       ['store_params'   , <store_params>])

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'                , @isstruct);
p0.addRequired('time_table'            , @istable);
p0.addRequired('subject_nrs'           , @isnumeric);
p0.addRequired('session_nrs'           , @isnumeric);
p0.addRequired('run_nrs'               , @isnumeric);
p0.addParameter('all_run_images'       , {}  , @iscell);
p0.addParameter('all_run_alpha_masks'  , {}  , @iscell);
p0.addParameter('all_run_im_nrs'  , {}  , @iscell);
p0.addParameter('fixsoafun'            , []  , @(x) isa(x,'function_handle'));
p0.addParameter('cdsoafun'             , []  , @(x) isa(x,'function_handle'));
p0.addParameter('load_params'          , true, @islogical);
p0.addParameter('store_params'         , true, @islogical);
p0.addParameter('session_type', 'MRI', @(x) ismember(x,{'MRI','BEHAVIOR'}));

% Parse inputs
p0.parse(params, time_table_master, subject_nrs, session_nrs, run_nrs, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

%% Define empty variables
if isempty(fixsoafun)
    % Fixation order and fixation
    fixsoafun = @() round(params.stim.fix.dotmeanchange + (params.stim.fix.dotchangeplusminus*(2*(rand-.5))));
end

if isempty(cdsoafun)
    % Contrast decrement gaussian time window onset
    cdsoafun = @() round(params.stim.cd.meanchange + params.stim.cd.changeplusminus*(2*(rand-.5)));
end

%% Check subject, session and run nrs to process

% Check if ses nr is member of all sessions
if strcmp(session_type, 'MRI')
    ses_ok = ismember([1:params.exp.session.n_total_sessions],session_nrs);
    run_ok = ismember([1:params.exp.session.deep.n_runs_per_session],run_nrs);
    
elseif strcmp(session_type, 'BEHAVIOR')
    ses_ok = ismember([1:params.exp.session.n_behavioral_sessions],session_nrs);
end

if isempty(ses_ok)
    error('[%s]: Session number is outside the range!\n',mfilename);
end

% If left undefined, we will load images for first session
if isempty(run_ok)
    error('[%s]: Run number is outside the range!\n',mfilename);
end

%% Load params if requested and we can find the file
if load_params
    % if we load multiple subjects, sessions or runs then we preallocate a
    % cell to concatenate them
    if length(subject_nrs)>1 || length(session_nrs)>1 || length(run_nrs)>1
        all_subj_run_frames      = cell(length(subject_nrs),length(session_nrs),length(run_nrs));
    end
    
    for sj = 1:length(subject_nrs)
        subjDir = fullfile(vcd_rootPath,'workspaces','info',sprintf('subj%03d', subject_nrs(sj)));
        
        for ses = 1:length(session_nrs)
            for rr = 1:length(run_nrs)
                fname = sprintf('subj%03d_ses%02d_run%02d_frames_%s*.mat',...
                    params.disp.name,subject_nrs(sj), session_nrs(ses), run_nrs(rr));
                
                d = dir(fullfile(subjDir,fname));
                fprintf('[%s]: Found %d subj run frames .mat file(s)\n',mfilename,length(d));
                if ~isempty(d)
                    if length(d) > 1
                        warning('[%s]: Multiple .mat files! Will pick the most recent one.\n', mfilename);
                    end
                    fprintf('[%s]: Loading subj run frames .mat file: %s\n', mfilename, d(end).name);
                    subj0 = load(fullfile(d(end).folder,d(end).name),'subj_run');
                else
                    error('[%s]: Can''t find subj run frames file!\n', mfilename)
                end
                
                all_subj_run_frames{subject_nrs(sj),session_nrs(ses),run_nrs(rr)} = subj0.subj_run;
                clear subj0
            end
        end
    end
    
else
    % Preallocate space for generated subject run frames
    all_subj_run_frames = cell(length(subject_nrs),length(session_nrs),length(run_nrs));
    
    for sj = 1:length(subject_nrs)
        for ses = 1:length(session_nrs)
            for rr = 1:length(run_nrs)
                
                if isempty(all_run_images{subject_nrs(sj),session_nrs(ses),run_nrs(rr)}) || isempty(all_run_alpha_masks{subject_nrs(sj),session_nrs(ses),run_nrs(rr)})
                    
                    % load subject run images
                    subjDataDir = fullfile(vcd_rootPath,'workspaces','info',sprintf('subj%03d',subject_nrs(sj)));
                    filename    = sprintf('subj%03d_ses%02d_run%02d_images_%s*.mat', ...
                        subject_nrs(sj), session_nrs(ses), run_nrs(rr), params.disp.name);
                    
                    d = dir(fullfile(subjDataDir,filename));
                    fprintf('[%s]: Found %d  subject run image .mat file(s)\n',mfilename,length(d));
                    if ~isempty(d)
                        if length(d) > 1
                            warning('[%s]: Multiple subject run image .mat files! Will pick the most recent one.\n', mfilename);
                        end
                        fprintf('[%s]: Loading  subject run .mat file: %s\n', mfilename, d(end).name);
                        a = load(fullfile(d(end).folder,d(end).name),'run_images', 'run_alpha_masks','run_im_nr');
                        all_run_images = a.run_images;
                        all_run_alpha_masks = a.run_alpha_masks;
                        all_run_im_nrs = a.run_im_nrs;
                    else
                        error('[%s]: Can''t find subject run image file!\n', mfilename)
                    end
                end
                
                
                % grab subj run trials from time_table_master
                subj_run = time_table_master((time_table_master.subj_nr==subject_nrs(sj) & ...
                    time_table_master.session_nr==session_nrs(ses) & ...
                    time_table_master.run_nr==run_nrs(rr)),:);
                
                run_images      = all_run_images{subject_nrs(sj),session_nrs(ses),run_nrs(rr)};
                run_alpha_masks = all_run_alpha_masks{subject_nrs(sj),session_nrs(ses),run_nrs(rr)};
                run_ims         = all_run_im_nrs{subject_nrs(sj),session_nrs(ses),run_nrs(rr)};
                
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
                run.button_response_fix = fix_matrix(:,4);
                
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
                    plot([0:1:length(run.fix_button_response)-1].*params.stim.presentationrate_hz,  run.fix_button_response ,'ko-');
                    title('fixation dot luminance rel diff')
                    set(gca,'YTick', [0,1,2], 'YTickLabel', {'No Change','Brighter','Dimmer'})
                    xlabel('Time (s)');
                    xlim([0, (length(run.fix_button_response).*params.stim.presentationrate_hz)])
                    
                    if params.store_imgs
                        saveFigsFolder = fullfile(vcd_rootPath,'figs');
                        filename = sprintf('vcd_subject%03d_session%02d_run%02d_fix_sequence.png', subject_nrs(sj),session_nrs(ses),run_nrs(rr));
                        print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Expand subject run time table
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % available run images should match the stimulus events
                assert(isequal(sum(ismember(subj_run.event_id(~subj_run.is_catch), [91, 92])), size(run_images,1)))
                
                % create a larger cell where each row is a frame
                total_nr_frames    = subj_run.event_end(end);
                run.images         = cell(total_nr_frames,2);
                run.masks          = cell(total_nr_frames,2);
                run.frame_event_nr = zeros(total_nr_frames,2);
                run.frame_im_nr    = zeros(total_nr_frames,2);
                
                
                jj = 1;
                frame_counter = 2;
                
                stim_events = subj_run.event_id(ismember(subj_run.event_id, [91, 92,990:996]));
                stim_row = find(ismember(subj_run.event_id, [91, 92,990:996]));
                
                for ii = 1:length(stim_events)
                    
                    % get frame nrs
                    nframes     = round(subj_run.event_dur(stim_row(ii)));
                    curr_frames = frame_counter:(frame_counter+nframes-1);
                    run.frame_event_nr(curr_frames,:) = stim_events(ii);
                    
                    % EYETRACKING TARGETS: 990 = center circle, 991:995 = saccade targets, 996 = pupil
                    if ismember(stim_events(ii),[990:996])
                        
                        %% TODO
                        % ADD EYETRACKING IMAGES
                        run.images(frame_counter,1) = {NaN};
                        run.masks(frame_counter,1)  = {NaN};
                        run.frame_im_nr(curr_frames,1) = NaN;
                        
                        % STIMULI:  91 = Stim interval 1, 92 = Stim interval 2,
                    elseif (stim_events(ii) == 91 ||stim_events(ii) == 92)
                        
                        if ~subj_run.is_catch(stim_row(ii))
                            
                            % insert pixel image (uint8)
                            if any(strcmp(subj_run.stim_class_name(stim_row(ii),:),{'rdk'})) && ((ndims(run_images{jj,1}) == 4) || (ndims(run_images{jj,2}) == 4)) % we are dealing with rdk, which has time dim
                                
                                nsides = find(strcmp(subj_run.stim_class_name(stim_row(ii),:),{'rdk'})==true);
                                clear rdk_images
                                for side = nsides
                                    rdk_images(:,side) = squeeze(mat2cell(run_images{jj,side}, size(run_images{jj,side},1), ...
                                        size(run_images{jj,side},2), size(run_images{jj,side},3), ones(1,size(run_images{jj,side},4))));
                                end
                                
                                rdk_images2 = cell(nframes,2);
                                if nsides == [1,2]
                                    for oob = 1:size(rdk_images,1)
                                        rdk_images2(oob,:) = [rdk_images(oob,1),rdk_images(oob,2)];
                                    end
                                elseif nsides == 1
                                    for oob = 1:size(rdk_images,1)
                                        rdk_images2(oob,:) = [rdk_images(oob,1),run_images(jj,2)];
                                    end
                                elseif nsides == 2
                                    for oob = 1:size(rdk_images,1)
                                        rdk_images2(oob,:) = [run_images(jj,1),rdk_images(oob,2)];
                                    end
                                end
                                
                                if size(rdk_images2,1) < size(rdk_images2,2)
                                    rdk_images2 = rdk_images2';
                                end
                                
                                f_im    = rdk_images2;
                                f_masks = repmat(run_alpha_masks(jj,:), nframes, 1);
                                
                            else
                                f_im     = run_images(jj,:);
                                f_masks  = run_alpha_masks(jj,:);
                            end
                            
                            
                            % APPLY CONTRAST DECREMENT
                            if strcmp(subj_run.task_class_name(stim_row(ii)),'cd')
                                % 50% change we will actually apply the contrast
                                % decrement change to stimulus
                                cd_change = false;
                                
                                if rand(1)>params.stim.cd.prob
                                    cd_change = true;
                                    
                                    [f_im_cd, c_onset] = vcd_applyContrastDecrement(params, cdsoafun, subj_run.stim_class_name{stim_row(ii),1}, f_im);
                                    run.images(curr_frames,:) = f_im_cd;
                                    
                                    f_cd = c_onset:(c_onset+length(params.stim.cd.t_gauss)-1);
                                    
                                    % log decrement in contrast column
                                    subj_run.contrast(curr_frames(f_cd)) = params.stim.cd.t_gauss;
                                    
                                    % Add button response to start of decrement onset
                                    subj_run.button_response(curr_frames(c_onset-1)) = 2; % button 2: no there was a change (yet)
                                    subj_run.button_response(curr_frames(c_onset):curr_frames(end)) = 1; % button 1: yes there was a change
                                    subj_run.cd_start(curr_frames(c_onset)) = c_onset; % log precise onset
                                    
                                else
                                    % Add button response: no
                                    subj_run.button_response(curr_frames) = 2; % button 2: no there was a change (apply to all time points in the block)
                                end
                                
                                
                            else % no contrast decrement
                                if size(f_im,1) > 1
                                    % add multiple images
                                    run.images(curr_frames,:) = f_im;
                                    run.masks(curr_frames,:)  = f_masks;
                                    run.frame_im_nr(curr_frames,:) = repmat(run_ims(jj,:), length(curr_frames),1);
                                else
                                    % add single image at the start
                                    run.images(curr_frames(1),:) = f_im;
                                    run.masks(curr_frames(1),:)  = f_masks;
                                    run.frame_im_nr(curr_frames,:) = repmat(run_ims(jj,:), length(curr_frames),1);
                                end
                                
                                subj_run.button_response = 
                                
                            end
                            
                        else % catch blocks
                            run.images(curr_frames,:) = repmat({0},length(curr_frames),2);
                            run.masks(curr_frames,:)  = repmat({0},length(curr_frames),2);
                            run.frame_im_nr(curr_frames,:) = repmat(run_ims(jj,:), length(curr_frames),1);
                        end % if catch
                        
                        % count one image nr
                        jj = jj+1;

                    else % non stim-events
                        % do nothing
                    end
                    
                    % update frame counte
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
                    assert(isequal(unique(diff(find(diff(run.button_response_fix)>0)))', [params.stim.fix.dotmeanchange, 2*params.stim.fix.dotmeanchange]));
                    
                    fix_block_nrs = unique(subj_run.block_nr(fix_events,:)); % should be 1 or 2 or 3 blocks per rum
                    fix_update_idx = (run.button_response_fix>0); % 1 x 10298 --> 238 fixation changes per run
                    [~,fix_block_frames] = ismember(subj_run.block_nr,fix_block_nrs); % 40 trial events
                    fix_block_change_direction = run.button_response_fix(fix_update_idx); % 1=brighter, 2=dimmer
                    fix_block_abs_lum = run.fix_abs_lum(fix_update_idx); %
                    
                    % nr of block frames should be exactly 1 or more single-stim presentation
                    % block duration
                    assert(isequal(mod((max(subj_run.event_end(fix_block_frames>0))-min(subj_run.event_start(fix_block_frames>0)))+1,params.exp.block.total_single_epoch_dur),0))
                    % get all the frames for the entire run
                    all_frames = 0:total_nr_frames;
                    % find those frames where the fixation circle updated
                    time_frames_fix_updated = all_frames(fix_update_idx);
                    % get start and end time frames for fixation block from time table master
                    fix_block_start_end = [min(subj_run.event_start(fix_events)),  max(subj_run.event_end(fix_events))];
                    % see what fix block frames overlap with frames where the fixation circle changed
                    fix_block_changes_idx = ((time_frames_fix_updated >= fix_block_start_end(1)) & (time_frames_fix_updated <= fix_block_start_end(2)));
                    fix_block_changes_times = time_frames_fix_updated(fix_block_changes_idx);
                    fix_block_abs_lum = fix_block_abs_lum(fix_block_changes_idx);
                    fix_block_changes_correct_response = fix_block_change_direction(fix_block_changes_idx);
                    
                    for ff = 1:length(fix_block_changes_correct_response)
                        t_fix = fix_block_changes_times(ff);
                        t_response = fix_block_changes_correct_response(ff);
                        
                        t_tbl = find((subj_run.event_start <= t_fix) & (subj_run.event_end >= t_fix));
                        if ~isempty(t_tbl)
                            subj_run.button_response(t_tbl) = t_response;
                            subj_run.fix_lum(t_tbl) = fix_block_abs_lum(ff);
                            subj_run.fix_start(t_tbl) = fix_block_changes_times(ff);
                        end
                    end
                end
                
                % add timing
                run.timing = [0:1:length(run.frame_im_nr)-1]';
                
                % Store structs locally, if requested
                if params.store_params
                    fprintf('[%s]: Storing expanded time table for subject %003d..\n',mfilename, subject_nrs(sj))
                    saveDir = fullfile(vcd_rootPath,'workspaces','info',sprintf('subj%03d',subject_nrs(sj)));
                    if ~exist(saveDir,'dir'), mkdir(saveDir); end
                    save(fullfile(saveDir, ...
                        sprintf('subj%03d_ses%02d_run%02d_frames_%s_%s.mat', ...
                        subject_nrs(sj), session_nrs(ses), run_nrs(rr), params.disp.name, datestr(now,30))), ...
                        'subj_run','run','-v7.3')
                end
                
                
                % Add run_images and alpha_masks to larger cell array
                all_subj_run_frames{subject_nrs(sj),session_nrs(ses),run_nrs(rr)} = run;
                

                time_table_master.button ((time_table_master.subj_nr==subject_nrs(sj) & ...
                    time_table_master.session_nr==session_nrs(ses) & ...
                    time_table_master.run_nr==run_nrs(rr)),:) = subj_run;
                
            end % runs
        end % sessions
    end % subjects
end % load params or not

