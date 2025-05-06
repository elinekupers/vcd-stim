function all_subj_run_frames = vcd_expandImageOrderSingleRun_30Hz(params, ...
    time_table_master,subject_nrs,session_nrs,run_nrs,varargin)
% Functions to expand sequence of events from individual images into 30Hz
% framelocked presentation rate.
%
%  subj_run_frames = vcd_expandImageOrderSingleRun_30Hz(params, ...
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
                    subj0 = load(fullfile(d(end).folder,d(end).name),'subj_run_frames');
                else
                    error('[%s]: Can''t find subj run frames file!\n', mfilename)
                end
                
                all_subj_run_frames{subject_nrs(sj),session_nrs(ses),run_nrs(rr)} = subj0.subj_run_frames;
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
                        load(fullfile(d(end).folder,d(end).name),'run_images', 'run_alpha_masks');
                    else
                        error('[%s]: Can''t find subject run image file!\n', mfilename)
                    end
                end
                
                
                % grab subj run
                subj_run = time_table_master((time_table_master.subj_nr==subject_nrs(sj) & ...
                    time_table_master.session_nr==session_nrs(ses) & ...
                    time_table_master.run_nr==run_nrs(rr)),:);
                
                run_images = all_run_images{subject_nrs(sj),session_nrs(ses),run_nrs(rr)};
                run_alpha_masks = all_run_alpha_masks{subject_nrs(sj),session_nrs(ses),run_nrs(rr)};
                
                % check if the last event is a post-blank
                assert(strcmp(subj_run.event_name(end),'post-blank'))
                
                % get run duration (in frames)
                run_dur = subj_run.event_end(end)+2;
                
                assert((run_dur*params.stim.presentationrate_hz)/3600 < 600) % ensure a run is of reasonable length, and is not longer than 10 min
                
                % copy condition order table headers and scrub content
                sz = [run_dur size(subj_run(1,:),2)];
                subj_run_frames = vcd_preallocateNaNTable(sz(1), sz(2), subj_run(1,:), []);
                
                max_merge_cols = find(strcmp(subj_run.Properties.VariableNames,'islure'));
                
                % expand table with info about frame timing
                subj_run_frames.frame_nr   = [0:(run_dur-1)]';
                subj_run_frames.images     = cell(size(subj_run_frames,1),2);
                subj_run_frames.masks      = cell(size(subj_run_frames,1),2);
                subj_run_frames.fix_abs_lum = NaN(size(subj_run_frames,1),1);
                subj_run_frames.fix_button_response = NaN(size(subj_run_frames,1),1);
                subj_run_frames.button_response = NaN(size(subj_run_frames,1),1);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% %%%%% GENERATE FIXATION SEQUENCE %%%%%
                
                [~,fix_abs_lum,~,button_response_fix] = vcd_createFixationSequence(params,fixsoafun,run_dur);
                
                % Add columns to subject time table
                subj_run_frames.fix_abs_lum = fix_abs_lum;
                subj_run_frames.fix_button_response = button_response_fix;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % visualize fixation sequence
                if params.verbose
                    makeprettyfigures
                    figure(101); clf; set(gcf,'Position',[137,952,2424,600])
                    sgtitle(sprintf('Fixation sequence: subject %03d, session %02d, run %02d', subject_nrs(sj),session_nrs(ses),run_nrs(rr)))
                    subplot(211);
                    plot([0:1:length(subj_run_frames.fix_abs_lum)-1].*params.stim.presentationrate_hz,subj_run_frames.fix_abs_lum,'ko-');
                    xlabel('Time (s)'); ylabel('dot luminance');
                    ylim([0 255]); xlim([0, (length(subj_run_frames.fix_abs_lum).*params.stim.presentationrate_hz)])
                    title('fixation dot luminance sequence')
                    
                    subplot(212);
                    plot([0:1:length(subj_run_frames.fix_button_response)-1].*params.stim.presentationrate_hz,  subj_run_frames.fix_button_response ,'ko-');
                    title('fixation dot luminance rel diff')
                    set(gca,'YTick', [0,1,2], 'YTickLabel', {'No Change','Brighter','Dimmer'})
                    xlabel('Time (s)');
                    xlim([0, (length(subj_run_frames.fix_button_response).*params.stim.presentationrate_hz)])
                    
                    if params.store_imgs
                        saveFigsFolder = fullfile(vcd_rootPath,'figs');
                        filename = sprintf('vcd_subject%03d_session%02d_run%02d_fix_sequence.png', subject_nrs(sj),session_nrs(ses),run_nrs(rr));
                        print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                    end
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Expand subject run time table
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                im_nr = 1;
                frame_counter = 2;
                
                for ii = 1:length(subj_run.event_id)
                    
                    % EYETRACKING TARGETS: 990 = center circle, 991:995 = saccade targets, 996 = pupil
                    if ismember(subj_run.event_id(ii),[990:996])
                        
                        % get frame nrs
                        nframes = round(subj_run.event_dur(ii));
                        curr_frames = frame_counter:(frame_counter+nframes-1);
                        
                        % repeat row in table and add to framelocked table
                        subj_run_frames(curr_frames,1:size(subj_run,2)) = repmat(subj_run(ii,:), nframes,1);
                        
                        %% TODO 
                        %ADD EYETRACKING IMAGES
                        subj_run_frames.images(curr_frames,:) = repmat({NaN}, length(curr_frames),2);
                        subj_run_frames.masks(curr_frames,:)  = repmat({NaN}, length(curr_frames),2);
                        
                    % STIMULI:  91 = Stim interval 1, 92 = Stim interval 2,
                    elseif (subj_run.event_id(ii) == 91 || subj_run.event_id(ii) == 92)
                        
                        % get frame nrs
                        nframes = round(subj_run.event_dur(ii));
                        curr_frames = frame_counter:(frame_counter+nframes-1);
                        
                        % repeat row in table and add to framelocked table
                        subj_run_frames(curr_frames,1:size(subj_run,2)) = repmat(subj_run(ii,:), nframes,1);
                        

                        % insert pixel image (uint8)
                        if any(strcmp(subj_run.stim_class_name(ii,:),{'rdk'})) && ((ndims(run_images{im_nr,1}) == 4) || (ndims(run_images{im_nr,2}) == 4)) % we are dealing with rdk, which has time dim
                            
                            nsides = find(strcmp(subj_run.stim_class_name(ii,:),{'rdk'})==true);
                            clear rdk_images
                            for side = 1:length(nsides)
                            
                                rdk_images(:,side) = squeeze(mat2cell(run_images{im_nr,side}, size(run_images{im_nr,side},1), ...
                                                          size(run_images{im_nr,side},2), size(run_images{im_nr,side},3), ones(1,size(run_images{im_nr,side},4))));
                            end
                            
                            rdk_images2 = cell(nframes,2);
                            if nsides == [1,2]
                                for oob = 1:size(rdk_images,1)
                                    rdk_images2(oob,:) = [rdk_images(oob,1),rdk_images(oob,2)];
                                end
                            elseif nsides == 1
                                for oob = 1:size(rdk_images,1)
                                    rdk_images2(oob,:) = [rdk_images(oob,1),run_images(im_nr,2)];
                                end
                            elseif nsides == 2
                                for oob = 1:size(rdk_images,1)
                                    rdk_images2(oob,:) = [run_images(im_nr,1),rdk_images(oob,1)];
                                end
                            end
                            
                            if size(rdk_images2,1) < size(rdk_images2,2)
                                rdk_images2 = rdk_images2';
                            end
                            f_im    = rdk_images2;
                            f_masks = repmat(run_alpha_masks(im_nr,:), nframes, 1);
                            
                        else
                            f_im     = repmat(run_images(im_nr,:),nframes, 1);
                            f_masks  = repmat(run_alpha_masks(im_nr,:), nframes, 1);
                        end
                        
                        im_nr = im_nr+1;
                        
                        % add images to subject frame table
                        subj_run_frames.images(curr_frames,:) = f_im;
                        subj_run_frames.masks(curr_frames,:)  = f_masks;
                        
                        % APPLY CONTRAST DECREMENT
                        if strcmp(subj_run.task_class_name(ii),'cd')
                            % 50% change we will actually apply the contrast
                            % decrement change to stimulus
                            cd_change = false;
                            
                            if rand(1)>params.stim.cd.prob
                                cd_change = true;
                                
                                [f_im_cd, c_onset] = vcd_applyContrastDecrement(params, cdsoafun, subj_run.stim_class_name{ii,1}, f_im);
                                subj_run_frames.images(curr_frames,:) = f_im_cd;
                                
                                f_cd = c_onset:(c_onset+length(params.stim.cd.t_gauss)-1);
                                
                                % log decrement in contrast column
                                subj_run_frames.contrast(curr_frames(f_cd)) = params.stim.cd.t_gauss;
                                
                                % Add button response to start of decrement onset
                                subj_run_frames.button_response(curr_frames(c_onset)) = 1; % yes there was a change
                                
                                
                            end
                        end % if cd
                        
                        
                    else % NON-STIM EVENTS
                        
                        % get frame nrs
                        nframes = round(subj_run.event_dur(ii));
                        curr_frames = frame_counter:(frame_counter+nframes-1);
                        
                        % repeat row in table and add to framelocked table
                        subj_run_frames(curr_frames,1:size(subj_run,2)) = repmat(subj_run(ii,:), nframes,1);
                        
                        % update timing
                        subj_run_frames.event_dur(curr_frames) = 1;
                        if ii == 1 && strcmp(subj_run_frames.event_name(curr_frames(1)),'pre-blank')
                            subj_run_frames.event_start(curr_frames(1)) = 1;
                            subj_run_frames.event_end(curr_frames(1)) = subj_run_frames.event_start(curr_frames(1))+ subj_run_frames.event_dur(curr_frames(1));
                        else
                            subj_run_frames.event_start(curr_frames(1)) = subj_run_frames.event_end(curr_frames(1)-1);
                            subj_run_frames.event_end(curr_frames(1)) = subj_run_frames.event_start(curr_frames(1))+ subj_run_frames.event_dur(curr_frames(1));
                        end
                        
                        for curr_fr = curr_frames(2:end)
                            subj_run_frames.event_start(curr_fr) = subj_run_frames.event_end(curr_fr-1);
                            subj_run_frames.event_end(curr_fr) = subj_run_frames.event_start(curr_fr)+ subj_run_frames.event_dur(curr_fr);
                        end
                        
                        
                        if subj_run_frames.event_id(curr_frames(1)) == params.exp.block.response_ID % response window
                            if strcmp(subj_run_frames.task_class_name(curr_frames(1)-1),'cd')
                                if ~cd_change
                                    % Add button response: no
                                    subj_run_frames.button_response(curr_frames) = 2; % button 2: no there was a change
                                elseif cd_change
                                    % Add button response: yes
                                    subj_run_frames.button_response(curr_frames) = 1; % button 1: yes there was a change
                                end
                            end
                        end
                        
                    end
                    
                    frame_counter = curr_frames(end)+1;
                end % events
                
                fix_block = find(strcmp(subj_run_frames.stim_class_name,'fix'));
                if  ~isempty(fix_block)
                    fix_block_nrs = unique(subj_run_frames.block_nr(fix_block,:));
                    
                    for ff = 1:length(fix_block_nrs)
                        fix_block_frames = subj_run_frames.block_nr==fix_block_nrs(ff);
                        subj_run_frames.button_response(fix_block_frames) = subj_run_frames.fix_button_response(fix_block_frames);
                    end
                end
                
                
                % Store structs locally, if requested
                if params.store_params
                    fprintf('[%s]: Storing expanded time table for subject %003d..\n',mfilename, subject_nrs(sj))
                    saveDir = fullfile(vcd_rootPath,'workspaces','info',sprintf('subj%03d',subject_nrs(sj)));
                    if ~exist(saveDir,'dir'), mkdir(saveDir); end
                    save(fullfile(saveDir, ...
                        sprintf('subj%03d_ses%02d_run%02d_frames_%s_%s.mat', ...
                        subject_nrs(sj), session_nrs(ses), run_nrs(rr), params.disp.name, datestr(now,30))), ...
                        'subj_run_frames','-v7.3')
                end
                
                
                % Add run_images and alpha_masks to larger cell array
                all_subj_run_frames{subject_nrs(sj),session_nrs(ses),run_nrs(rr)} = subj_run_frames;
                
            end % runs
        end % sessions
    end % subjects
end % load params or not

