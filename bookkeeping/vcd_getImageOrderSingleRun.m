function [run_images, run_alpha_masks, all_images] = vcd_getImageOrderSingleRun(params, ...
    time_table_master, all_run_frames, subj_nr, ses_nr, ses_type, run_nr, env_type, varargin)
% VCD function to load uint8 stimuli, alpha transparency masks, and
% corresponding image nrs for each trial
%
%   [run_images, run_alpha_masks, all_images] =
%   vcd_getImageOrderSingleRun( ...
%       params, time_table_master, subj_nr, ses_nr, ses_type, run_nr, env_type, ...
%       ['all_images', <all_images>],['verbose', <verbose>],...
%       ['store_run_images',<store_run_images>],['env_type',<env_type>])
%
% Note image files are stored locally (workspaces > stimuli), and need to
% be created prior to running this function (see s_createStim.m) for every
% monitor display setup, such that it has the <params.disp.name> in the
% stimulus file name.


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'             , @isstruct);
p0.addRequired('time_table_master'  , @istable);
p0.addRequired('all_run_frames'     , @istable);
p0.addRequired('subj_nr'            , @isnumeric);
p0.addRequired('ses_nr'             , @isnumeric);
p0.addRequired('ses_type'           , @isnumeric);
p0.addRequired('run_nr'             , @isnumeric);
p0.addRequired('env_type'           , @(x) ismember(x,{'MRI','BEHAVIOR'}));
p0.addParameter('all_images'        , struct(), @isstruct);
p0.addParameter('verbose'           , false, @islogical);
p0.addParameter('savestim'          , false, @islogical); % store the actual stimuli we load for stimulus presentation

% Parse inputs
p0.parse(params, time_table_master, all_run_frames, subj_nr, ses_nr, ses_type, run_nr, env_type, varargin{:});

% Rename variables into general params struct
sfun = @(x) sprintf('%s = p0.Results.%s;',x,x);
cellfun(@eval, cellfun(sfun, fieldnames(p0.Results), 'UniformOutput', false))
clear sfun

%% Check if experimental parameters are already defined and load it if needed
if ~isfield(params,'exp') || isempty(params.exp)
    warning('[%s]: Can''t find exp field in params! Will run vcd_getSessionParams.m', mfilename);
    params.exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false, 'disp_name',params.disp.name);
end

% check if stimulus parameters are already defined and load it if needed
if ~isfield(params,'stim') || isempty(params.stim)
    warning('[%s]: Can''t find stim field in params! Will run vcd_getStimParams.m', mfilename);
    params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose', false, 'disp_name',params.disp.name);
end

%% Check subject, session and run nrs to process

% Check if subject nr is member of total subjects
subj_ok = ismember([0:999],subj_nr);
if isempty(subj_ok)
    error('\n[%s]: Subject number is outside the range!\n',mfilename);
end

% Check if ses nr is member of all sessions
if strcmp(env_type, 'MRI')
    if ses_nr == 1
        ses_ok = ismember([params.exp.session.n_wide_sessions],ses_nr);
        run_ok = ismember([1:params.exp.session.mri.wide.n_runs_per_session(ses_nr,ses_type)],run_nr); 
    else
        ses_ok = ismember([params.exp.session.n_deep_sessions],ses_nr);
        run_ok = ismember([1:params.exp.session.mri.deep.n_runs_per_session(ses_nr,ses_type)],run_nr); 
    end
elseif strcmp(env_type, 'BEHAVIOR')
    if params.is_demo
        ses_ok = ismember([1:params.exp.session.n_demo_sessions],ses_nr);
        run_ok = ismember([1:params.exp.session.demo.n_runs_per_session(ses_nr,ses_type)],run_nr);
    else
        ses_ok = ismember([1:params.exp.session.n_behavioral_sessions],ses_nr);
        run_ok = ismember([1:params.exp.session.behavior.n_runs_per_session],run_nr);
    end
end

if isempty(ses_ok)
    error('\n[%s]: Session number is outside the range!\n',mfilename);
end

% If left undefined, we will load images for first session
if isempty(run_ok)
    error('\n[%s]: Run number is outside the range!\n',mfilename);
end


%% Load params if requested and we can find the file
tic
fprintf('\n[%s]: Loading stimuli..\n',mfilename);
    

% Check if we need to load more images
% If we did not provide images.. then create the field names of the all_images struct
if isempty(all_images) || isempty(fieldnames(all_images))
    all_images = struct('gabor',[],'rdk',[],'dot',[],'obj',[],'ns',[],...
        'fix',[], 'info',[], 'eye',[],'alpha',[]);
end

% Check if we defined folder names
if ~isfield(params, 'instrfolder') || isempty(params.instrfolder)
    params.instrfolder = fullfile(vcd_rootPath,'workspaces','instructions');
end

if ~isfield(params, 'stimfolder') || isempty(params.stimfolder)
    params.stimfolder = fullfile(vcd_rootPath,'workspaces','stimuli');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% INSTRUCTIONS IM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load all stored instruction cues if they aren't already in the all_images
% struct
if ~isfield(all_images,'instr') || isempty(all_images.instr)
    if verbose; fprintf('[%s]: Loading instruction images (text + icons)..\n',mfilename); end

    d = dir(fullfile(params.instrfolder,'instruction_images', '*.png'));
    
    %  instr_im: 4D array: [x,y, 3, nr_crossings]ß
    all_images.instr = uint8(zeros(700,700,3,length(params.exp.crossings)+1));
    
    for aa = 1:length(d)
        a1 = imread(fullfile(d(aa).folder,d(aa).name));
        if ismatrix(a1)
            a1 = repmat(a1,[1 1 3]);
        end
        all_images.instr(:,:,:,aa) = a1;
        all_images.info.instr(aa)  = {d(aa).name};
        all_images.alpha.instr(aa) = NaN;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIX IM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load stored fixation dot images if needed
if ~isfield(all_images,'fix') || isempty(all_images.fix)
    if verbose; fprintf('[%s]: Loading fixation dot images..\n',mfilename); end

    %  FIX CIRCLE: 5D array: [x,y, 3, lum, rim types]
    d = dir(fullfile(params.stimfolder,params.disp.name, sprintf('fix_%s*.mat',params.disp.name)));
    a1 = load(fullfile(d(end).folder,d(end).name), 'fix_im','mask','info');

    all_images.fix       = a1.fix_im; 
    all_images.info.fix  = a1.info; 
    all_images.alpha.fix = a1.mask;
end

%%%%%%%%%%%%%%%%%%%%%%%%% EYETRACKING BLOCK STIM %%%%%%%%%%%%%%%%%%%%%%%%%%
% Load stored eyetracking block target images if needed
if ~isfield(all_images,'eye') || isempty(all_images.eye)
     if verbose; fprintf('[%s]: Loading eyetracking target images..\n',mfilename); end
    
    % ET TARGETS: 4D array: [x,y, 3, type]
    d = dir(sprintf('%s*.mat', params.stim.el.stimfile));
    a = load(fullfile(d(end).folder,d(end).name), 'sac_im','pupil_im_white','pupil_im_black');
    all_images.eye.sac_im    = a.sac_im; 
    all_images.eye.pupil_im_white = a.pupil_im_white;
    all_images.eye.pupil_im_black = a.pupil_im_black;
    clear a d;
end


% Get subject time table from master table
run_table = time_table_master((time_table_master.session_nr==ses_nr ...
    & time_table_master.session_type==ses_type ...
    & time_table_master.run_nr==run_nr),:);

run_frames = all_run_frames((all_run_frames.session_nr==ses_nr ...
    & all_run_frames.session_type==ses_type ...
    & all_run_frames.run_nr==run_nr),:);

%% Predefine cell arrays for stimuli
run_images      = cell(size(run_frames.frame_im_nr)); % time frames x stim loc (1:left and 2:right side stimuli)
run_alpha_masks = run_images;

% Get stimulus rows in time table                
stim_row    = find(ismember(run_table.event_id, [params.exp.block.stim_epoch1_ID, params.exp.block.stim_epoch2_ID, ...
                    params.exp.block.eye_gaze_fix_ID, params.exp.block.eye_gaze_sac_target_ID, ...
                    params.exp.block.eye_gaze_pupil_black_ID, params.exp.block.eye_gaze_pupil_white_ID]));


% Loop over all stimulus rows
if verbose; fprintf('[%s]: Insert stimuli for each trial into "run_frames"..\n',mfilename); end
 
for ii = 1:length(stim_row)
    
    fprintf('.')
    
    frame_counter = run_table.event_start(stim_row(ii))+1; % t=1 is 0, but we can't use 0 as index
    assert(run_frames.timing(frame_counter)==run_table.event_start(stim_row(ii)));
    
    nframes     = round(run_table.event_dur(stim_row(ii)));
    curr_frames = frame_counter:(frame_counter+nframes-1);
    
    % Get stimulus image nrs
    if isnan(run_table.stim_nr_right(stim_row(ii)))
        nsides = 1;
    else
        nsides = 2;
    end
    
    for side = 1:nsides
        
        if ~isempty(run_table.stim_class_name{stim_row(ii),side})
            % Get stimulus class name
            stimClass = run_table.stim_class_name{stim_row(ii),side};
            
            % fill in unique_im nr
            unique_im = run_frames.frame_im_nr(curr_frames(1),side);
            
            switch stimClass
                case 'gabor'
                     if verbose; fprintf('[%s]: Loading gabor images..\n',mfilename); end
                    if isempty(all_images.gabor)
                        % GABORS: 6D array: [x,y,8 orient, 4 phase,3 contrast, og + 4 delta]
                        d = dir(sprintf('%s*.mat', params.stim.gabor.stimfile));
                        a = load(fullfile(d(end).folder,d(end).name), 'gabors','info'); % also contains 'masks', but we don't need those..
                        
                        all_images.gabor = a.gabors; 
                        all_images.info.gabor  = a.info; 
%                         all_images.alpha.gabor = a.masks; 
                        clear a d;
                    end
                    
                    if strcmp(run_table.event_name(stim_row(ii)),'stim1') && run_table.is_catch(stim_row(ii)) == 0
                        
                        % GABORS: 6D array: [x,y,orient,contrast,phase,delta]
                        idx0 = find(all_images.info.gabor.unique_im==unique_im);
                        gbr_ori      = all_images.info.gabor.orient_deg(idx0);
                        gbr_contrast = all_images.info.gabor.contrast(idx0);
                        gbr_phase    = all_images.info.gabor.phase_deg(idx0);
                        ori_idx      = all_images.info.gabor.orient_i(idx0);
                        con_idx      = all_images.info.gabor.contrast_i(idx0);
                        
                        % check if stim params match
                        assert(isequal( run_table.orient_dir(stim_row(ii),side) , gbr_ori));
                        assert(isequal( run_table.contrast(stim_row(ii),side), gbr_contrast));
                        assert(isequal( run_table.gbr_phase(stim_row(ii),side), gbr_phase));
                        
                        run_images{curr_frames(1),side} = all_images.gabor(:,:,:,ori_idx,con_idx,1);
                        run_alpha_masks{curr_frames(1),side} = []; %all_images.alpha.gabor(:,:,ori_idx,con_idx,1);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'wm') ...
                            && run_table.is_catch(stim_row(ii)) == 0
                        delta_deg = run_table.stim2_delta(stim_row(ii),side);
                        
                        idx0 = find( (all_images.info.gabor.orient_deg == (run_table.orient_dir(stim_row(ii),side))) & ...
                            (all_images.info.gabor.phase_deg==run_table.gbr_phase(stim_row(ii),side)) & ...
                            (all_images.info.gabor.contrast==max(params.stim.gabor.contrast)) & ...
                            (all_images.info.gabor.delta_deg==delta_deg) );
                        
                        ori_idx      = all_images.info.gabor.orient_i(idx0);
                        con_idx      = max(all_images.info.gabor.contrast_i); % we only use highest contrast for gabor wm test im
                        test_im      = all_images.info.gabor.unique_im(idx0);
                        assert(isequal(test_im,run_table.stim2_im_nr(stim_row(ii),side)))
                        
                        % check if stim params match
                        assert(isequal( run_table.orient_dir(stim_row(ii),side) , all_images.info.gabor.orient_deg(idx0)));
                        assert(isequal( run_table.gbr_phase(stim_row(ii),side),   all_images.info.gabor.phase_deg(idx0)));
                        assert(isequal( delta_deg, all_images.info.gabor.delta_deg(idx0)));
                        
                        delta_idx = 1+ find(delta_deg == params.stim.gabor.delta_from_ref); % add one because 6th dimension in all_images.gabor index delta = [0,1,2,3,4]
                        
                        run_images{curr_frames(1),side}      = all_images.gabor(:,:,:,ori_idx,con_idx, delta_idx);
                        run_alpha_masks{curr_frames(1),side} = []; %all_images.alpha.gabor(:,:,ori_idx,con_idx, delta_idx);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'ltm') ...
                            && run_table.is_catch(stim_row(ii)) == 0
                        
                        
                        if run_table.islure(stim_row(ii)) % same stim class
                            run_images{curr_frames(1),side} = [];
                            run_alpha_masks{curr_frames(1),side} = [];
                            
                        else % other stim class
                            run_images{curr_frames(1),side} = [];
                            run_alpha_masks{curr_frames(1),side} = [];
                        end
                        
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'img') ...
                            && run_table.is_catch(stim_row(ii)) == 0
                        
                        if run_table.response(stim_row(ii),side) % yes -- dots overlap
                            run_images{curr_frames(1),side}      = [];
                            run_alpha_masks{curr_frames(1),side} = [];
                            
                        else % no -- dots overlap
                            run_images{curr_frames(1),side}      = [];
                            run_alpha_masks{curr_frames(1),side} = [];
                        end
                        
                        test_im = run_table.stim2_im_nr(stim_row(ii),side);
                        
                    elseif run_table.is_catch(stim_row(ii)) == 1
                        % catch trial
                        run_images{curr_frames,side}      = uint8(zeros(1,1,1));
                        run_alpha_masks{curr_frames,side} = uint8(zeros(1,1,1));
                    end
                    
                    
                    
                    
                case 'rdk'
                    if verbose; fprintf('[%s]: Loading rdk images..\n',mfilename); end
                    
                    if strcmp(run_table.event_name(stim_row(ii)),'stim1') && run_table.is_catch(stim_row(ii)) == 0
                        
                        % attempt to store rdk movies
%                         tbl_idx = find(ismember([params.stim.rdk.unique_im_nrs_core,params.stim.rdk.unique_im_nrs_wm_test],unique_im));
%                         
%                         if ~isfield(all_images.info, 'rdk')
%                             all_images.info.rdk = [];
%                         end
%                             
%                         if ~isempty(all_images.rdk(:,:,:,:,tbl_idx)) && ...
%                            ~isempty(all_images.info.rdk) && ismember(unique_im,all_images.info.rdk.unique_im)
%                             
%                        
%                            
%                             all_images.rdk(:,:,:,:,tbl_idx) = frames;
%                             
%                             
%                         else

                            % RDKs: 130 mat separate files: 8 directions x 3 coherence levels x 5 deltas (0 + 4 deltas)
                            stimDir = dir(fullfile(sprintf('%s*',params.stim.rdk.stimfile)));
                            filename = sprintf('%04d_vcd_rdk_ori%02d_coh%02d_delta%02d',...
                                unique_im, ...
                                find(run_table.orient_dir(stim_row(ii),side) == params.stim.rdk.dots_direction), ...
                                find(run_table.rdk_coherence(stim_row(ii),side) == params.stim.rdk.dots_coherence), ...
                                0);
                            
                            stimfile = fullfile(stimDir(1).folder,stimDir(1).name,sprintf('%s.mat', filename));
                            if exist(stimfile,'file')
                                load(stimfile, 'frames','rdk_info'); % also contains 'mask', but we don't need that..
                            else
                                error('[%s]: Can''t find RDK stim file!')
                            end
                            
                            % ensure stimulus duration
                            frames =  frames(:,:,:,1:params.stim.rdk.duration);
                            
                            %                         all_images.info.rdk = cat(1,all_images.info.rdk, rdk_info);
                            %                         [~,reorder_idx] = sort(all_images.info.rdk.unique_im,'ascend');
                            %                         all_images.info.rdk = all_images.info.rdk(reorder_idx,:);
                            
%                         end
                        
                        % check if stim description matches
                        idx0 = find(rdk_info.unique_im == unique_im);
                        dot_motdir = rdk_info.dot_motdir_deg(idx0);
                        dot_coh    = rdk_info.dot_coh(idx0);
                        assert(isequal(run_table.rdk_coherence(stim_row(ii),side),dot_coh));
                        assert(isequal(run_table.orient_dir(stim_row(ii),side),dot_motdir));


                        % expand rdk movies into frames
                        rdk_images = squeeze(mat2cell(frames, size(frames,1), ...
                            size(frames,2), size(frames,3), ones(1,size(frames,4))));
                        
                        if size(rdk_images,1) < size(rdk_images,2)
                            rdk_images = rdk_images';
                        end
                        
                        % rdk_masks = repmat({mask}, size(rdk_images,1),1); % we are ignoring masks for rdks for now
                        
                        run_images(curr_frames,side) = rdk_images;
                        run_alpha_masks(curr_frames,side) = cell(size(rdk_images,1), 1);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'wm') ...
                            && run_table.is_catch(stim_row(ii)) == 0
                        
                        delta_idx = find(run_table.stim2_delta(stim_row(ii),side) == params.stim.rdk.delta_from_ref); % here we do not +1 because we search for the file name with "delta01", "delta02", etc.
                        
                        
                        if side == 1
                            corresponding_unique_im = run_table.stim_nr_left(stim_row(ii));
                            test_im = run_table.stim2_im_nr(stim_row(ii),side);
                        elseif side == 2
                            corresponding_unique_im = run_table.stim_nr_right(stim_row(ii));
                            test_im = run_table.stim2_im_nr(stim_row(ii),side);
                        end
                        
                        % get stim features
%                         corresponding_unique_im_idx = find(rdk_info.unique_im==corresponding_unique_im);
                        updated_ori = run_table.stim2_orient_dir(stim_row(ii),side);
                        delta_test = run_table.stim2_delta(stim_row(ii),side);
                        og_ori     = mod(updated_ori-delta_test,360);

                        % RDKs: 130 mat files: 8 directions x 3 coherence levels x 5 deltas (0 + 4 deltas)
                        stimDir = dir(fullfile(sprintf('%s*',params.stim.rdk.stimfile)));
                        filename = sprintf('%04d_vcd_rdk_ori%02d_coh%02d_delta%02d',...
                            test_im, ...
                            find(og_ori == params.stim.rdk.dots_direction), ...
                            find(max(params.stim.rdk.dots_coherence) == params.stim.rdk.dots_coherence), ...
                            delta_idx);
                        
                        stimfile = fullfile(stimDir(1).folder,stimDir(1).name,sprintf('%s.mat', filename));
                        if exist(stimfile,'file')
                            load(stimfile, 'frames','rdk_info'); % also contains 'mask', but we don't need that..
                        else
                            error('[%s]: Can''t find RDK stim file!')
                        end
                        
                        
                        % check if stim description matches
                        assert(isequal(rdk_info.dot_motdir_deg, updated_ori));
                        
                        % ensure duration
                        frames =  frames(:,:,:,1:params.stim.rdk.duration);

                        
                        % expand rdk movies into frames
                        rdk_images = squeeze(mat2cell(frames, size(frames,1), ...
                            size(frames,2), size(frames,3), ones(1,size(frames,4))));
                        
                        if size(rdk_images,1) < size(rdk_images,2)
                            rdk_images = rdk_images';
                        end
                        
%                         rdk_masks = repmat({mask}, size(rdk_images,1), 1);
                        
                        run_images(curr_frames,side) = rdk_images;
                        run_alpha_masks(curr_frames,side) = cell(size(rdk_images,1), 1);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'ltm') ...
                            && run_table.is_catch(stim_row(ii)) == 0
                        
                        if run_table.islure(stim_row(ii),side) % same stim class
                            pair_im = [];
                            
                            %                                         % expand rdk movies into frames
                            %                                         rdk_images = squeeze(mat2cell(frames, size(frames,1), ...
                            %                                             size(frames,2), size(frames,3), ones(1,size(frames,4))));
                            %
                            %                                         if size(rdk_images,1) < size(rdk_images,2)
                            %                                             rdk_images = rdk_images';
                            %                                         end
                            %
                            %                                         rdk_masks = repmat(mask, size(rdk_images,1), 1);
                            run_images(curr_frames,side) = rdk_images;
                            run_alpha_masks(curr_frames,side) = cell(size(rdk_images,1), 1);
                            
                        else % other stim class
                            pair_im = [];
                            
                            %                                         % expand rdk movies into frames
                            %                                         rdk_images = squeeze(mat2cell(frames, size(frames,1), ...
                            %                                             size(frames,2), size(frames,3), ones(1,size(frames,4))));
                            %
                            %                                         if size(rdk_images,1) < size(rdk_images,2)
                            %                                             rdk_images = rdk_images';
                            %                                         end
                            %
                            %                                         rdk_masks = repmat(mask, size(rdk_images,1), 1);
                            run_images(curr_frames,side) = rdk_images;
                            run_alpha_masks(curr_frames,side) = cell(size(rdk_images,1), 1);
                            
                        end
                        
                        
                    elseif run_table.is_catch(stim_row(ii)) == 1
                        % catch trial
                        run_images(curr_frames,side)      = uint8(zeros(1,1,1));
                        run_alpha_masks(curr_frames,side) = uint8(zeros(1,1));
                    end
                    
                case 'dot'
                    if verbose; fprintf('[%s]: Loading dot images..\n',mfilename); end
                    if isempty(all_images.dot)
                        % Simple dot: 2D array: [x,y]
                        d = dir(sprintf('%s*.mat', params.stim.dot.stimfile));
                        a=load(fullfile(d(end).folder,d(end).name), 'dot','mask','info');
                        all_images.dot = a.dot; 
                        all_images.alpha.dot = a.mask; 
                        all_images.info.dot = a.info; 
                        clear a d;
                    end
                    
                    
                    if run_table.is_catch(stim_row(ii)) == 0
                        run_images{curr_frames(1),side} = all_images.dot;
                        run_alpha_masks{curr_frames(1),side} = all_images.alpha.dot;
                        
                        if strcmp(run_table.event_name(stim_row(ii)),'stim2')
                            test_im = run_table.stim2_im_nr(stim_row(ii),side);
                        end
                    elseif run_table.is_catch(stim_row(ii)) == 1
                        % catch trial
                        run_images{curr_frames,side}      = uint8(zeros(1,1,1));
                        run_alpha_masks{curr_frames,side} = uint8(zeros(1,1));
                    end
                    
                case 'obj'
                    if verbose; fprintf('[%s]: Loading object images..\n',mfilename); end
                    if isempty(all_images.obj)
                        % Complex objects: 4D array: [x,y,16 object, og + 10 rotation]
                        d = dir(sprintf('%s*.mat', params.stim.obj.stimfile));
                        a = load(fullfile(d(end).folder,d(end).name), 'objects','masks','objects_catch','masks_catch','info');
                        all_images.obj            = a.objects; 
                        all_images.objcatch       = a.objects_catch;
                        all_images.alpha.obj      = a.masks; 
                        all_images.alpha.objcatch = a.masks_catch;
                        all_images.info.obj       = a.info; 
                        clear a d;
                    end
                    
                    if strcmp(run_table.event_name(stim_row(ii)),'stim1') && run_table.is_catch(stim_row(ii)) == 0 ...
                            && (run_table.is_objectcatch(stim_row(ii)) == 0 || isnan(run_table.is_objectcatch(stim_row(ii))))
                       
                        if any(strcmp(all_images.info.obj.Properties.VariableNames, 'img_quiz_dot_specialcore_stim_nr'))
                            idx = find( (all_images.info.obj.super_cat == run_table.super_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.basic_cat == run_table.basic_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.sub_cat == run_table.sub_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.abs_rot == run_table.orient_dir(stim_row(ii),side)) & ...
                                isnan(all_images.info.obj.img_quiz_dot_specialcore_stim_nr) );
                        else
                            idx = find( (all_images.info.obj.super_cat == run_table.super_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.basic_cat == run_table.basic_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.sub_cat == run_table.sub_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.abs_rot == run_table.orient_dir(stim_row(ii),side)) );
                        end
                        
                        
                        obj_super = all_images.info.obj.super_cat_name(idx);
                        obj_basic = all_images.info.obj.basic_cat_name(idx);
                        obj_sub   = all_images.info.obj.sub_cat_name(idx);
                        obj_rotation = all_images.info.obj.abs_rot(idx);
                        
                        % check if stim idx matches
                        assert(isequal(all_images.info.obj.unique_im(idx),     unique_im));
                        assert(isequal(obj_super,   run_table.super_cat_name(stim_row(ii),side)));
                        assert(isequal(obj_basic,   run_table.basic_cat_name(stim_row(ii),side)));
                        assert(isequal(obj_sub,     run_table.sub_cat_name(stim_row(ii),side)));
                        assert(isequal(obj_rotation,run_table.orient_dir(stim_row(ii),side)));
                        
                        run_images{curr_frames(1),side}      = all_images.obj(:,:,:,unique_im==params.stim.obj.unique_im_nrs_core,1);
                        run_alpha_masks{curr_frames(1),side} = all_images.alpha.obj(:,:,unique_im==params.stim.obj.unique_im_nrs_core,1);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim1') && run_table.is_catch(stim_row(ii)) == 0 ...
                            && run_table.is_objectcatch(stim_row(ii)) == 1
                            
                            if any(strcmp(all_images.info.obj.Properties.VariableNames, 'img_quiz_dot_specialcore_stim_nr'))
                                idx = find( (all_images.info.obj.super_cat == run_table.super_cat(stim_row(ii),side)) & ...
                                    (all_images.info.obj.basic_cat == run_table.basic_cat(stim_row(ii),side)) & ...
                                    (all_images.info.obj.sub_cat == run_table.sub_cat(stim_row(ii),side)) & ...
                                    (all_images.info.obj.abs_rot == run_table.orient_dir(stim_row(ii),side)) & ...
                                    (all_images.info.obj.rel_rot == 0) & ...
                                    isnan(all_images.info.obj.img_quiz_dot_specialcore_stim_nr) );
                            else
                                idx = find( (all_images.info.obj.super_cat == run_table.super_cat(stim_row(ii),side)) & ...
                                    (all_images.info.obj.basic_cat == run_table.basic_cat(stim_row(ii),side)) & ...
                                    (all_images.info.obj.sub_cat == run_table.sub_cat(stim_row(ii),side)) & ...
                                    (all_images.info.obj.abs_rot == run_table.orient_dir(stim_row(ii),side)) & ...
                                    (all_images.info.obj.rel_rot == 0));
                            end
                            obj_super = all_images.info.obj.super_cat_name(idx);
                            obj_basic = all_images.info.obj.basic_cat_name(idx);
                            obj_sub   = all_images.info.obj.sub_cat_name(idx);
                            obj_rotation = all_images.info.obj.abs_rot(idx);

                            % check if stim idx matches
                            assert(isequal(all_images.info.obj.unique_im(idx),     unique_im));
                            assert(isequal(obj_super,   run_table.super_cat_name(stim_row(ii),side)));
                            assert(isequal(obj_basic,   run_table.basic_cat_name(stim_row(ii),side)));
                            assert(isequal(obj_sub,     run_table.sub_cat_name(stim_row(ii),side)));
                            assert(isequal(obj_rotation,run_table.orient_dir(stim_row(ii),side)));
                      
                            if ismember(unique_im,params.stim.obj.unique_im_nrs_core)
                                run_images{curr_frames(1),side}      = all_images.obj(:,:,:,unique_im==params.stim.obj.unique_im_nrs_core, 1);
                                run_alpha_masks{curr_frames(1),side} = all_images.alpha.obj(:,:,unique_im==params.stim.obj.unique_im_nrs_core, 1);
                                
                                % or the case when we deal with the objectcatch
                            elseif ismember(unique_im,params.stim.obj.unique_im_nrs_objcatch)
                                
                                assert(all_images.info.obj.is_objectcatch(idx)==1)
                                
                                tmp_im = strsplit(run_table.condition_name{stim_row(ii),side},'-');
                                tmp_im = str2num(tmp_im{2});
                                assert(strcmp(all_images.info.obj.sub_cat_name(all_images.info.obj.unique_im==tmp_im),obj_sub))
                                corresponding_im = find(all_images.info.obj.unique_im==tmp_im);
                                assert(isequal(corresponding_im, find(tmp_im==params.stim.obj.unique_im_nrs_core)));
                                obj_catch_mat = reshape(params.stim.obj.unique_im_nrs_objcatch,[],16)';
                                obj_catch_im_nrs = obj_catch_mat(corresponding_im,:);
                                obj_catch_rot_idx = find(unique_im==obj_catch_im_nrs);
                                assert(isequal(params.stim.obj.catch_rotation(corresponding_im,find(unique_im==obj_catch_im_nrs)),obj_rotation))
                                run_images{curr_frames(1),side}      = all_images.objcatch(:,:,:,corresponding_im, obj_catch_rot_idx);
                                run_alpha_masks{curr_frames(1),side} = all_images.alpha.objcatch(:,:,corresponding_im, obj_catch_rot_idx);
                            end
                            
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'wm') ...
                            && run_table.is_catch(stim_row(ii)) == 0 && (run_table.is_objectcatch(stim_row(ii)) == 0 || isnan(run_table.is_objectcatch(stim_row(ii))))
                        
                        delta_idx  = run_table.stim2_delta(stim_row(ii),side);
                        ref_dir    = run_table.orient_dir(stim_row(ii),side);
                        test_dir   = run_table.stim2_orient_dir(stim_row(ii),side);
                        delta      = test_dir - ref_dir;
                        assert(isequal(delta_idx,delta));
                        % grab stim 1 unique im nr
                        if side == 1
                            corresponding_unique_im = run_table.stim_nr_left(stim_row(ii))==params.stim.obj.unique_im_nrs_core;
                        elseif side == 2
                            corresponding_unique_im = run_table.stim_nr_right(stim_row(ii))==params.stim.obj.unique_im_nrs_core;
                        end
                        % check if stim idx matches
                        if any(strcmp(all_images.info.obj.Properties.VariableNames, 'img_quiz_dot_specialcore_stim_nr'))
                            idx = find( (all_images.info.obj.super_cat == run_table.super_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.basic_cat == run_table.basic_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.sub_cat == run_table.sub_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.abs_rot == ref_dir) & ...
                                isnan(all_images.info.obj.img_quiz_dot_specialcore_stim_nr) );
                        else
                            idx = find( (all_images.info.obj.super_cat == run_table.super_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.basic_cat == run_table.basic_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.sub_cat == run_table.sub_cat(stim_row(ii),side)) & ...
                                (all_images.info.obj.abs_rot == ref_dir) );
                        end
                        test_im         = run_table.stim2_im_nr(stim_row(ii),side);
                        delta_idx0      = find(delta==params.stim.obj.delta_from_ref) +1; % add one because 5th dimension index delta = [0,1,2,3,4] 
                        assert(strcmp(run_table.sub_cat_name(stim_row(ii),side), all_images.info.obj.sub_cat_name(idx)));
                        
                        tmp_test_im     = reshape(params.stim.obj.unique_im_nrs_wm_test,4,[]);
                        assert(isequal(test_im,tmp_test_im(delta_idx0 -1 ,corresponding_unique_im))) % subtract one here because tmp_test_im does not contain delta00
                        
                        run_images{curr_frames(1),side}      = all_images.obj(:,:,:,corresponding_unique_im,delta_idx0);
                        run_alpha_masks{curr_frames(1),side} = all_images.alpha.obj(:,:,corresponding_unique_im,delta_idx0);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'ltm') ...
                            && run_table.is_catch(stim_row(ii)) == 0
                        
                        if run_table.islure(stim_row(ii),side) % same stim class
                            pair_im = [];
                            run_images{curr_frames(1),side} = [];
                            run_alpha_masks{curr_frames(1),side} = [];
                            
                        else % other stim class
                            pair_im = [];
                            run_images{curr_frames(1),side} = [];
                            run_alpha_masks{curr_frames(1),side} = [];
                        end
                        
                        
                    elseif run_table.is_catch(stim_row(ii)) == 1
                        run_images{curr_frames,side}      = uint8(zeros(1,1,1));
                        run_alpha_masks{curr_frames,side} = uint8(zeros(1,1));
                    end
                    
                    
                case 'ns'
                    if verbose; fprintf('[%s]: Loading scene images..\n',mfilename); end
                    if isempty(all_images.ns)
                        % NS: 6D array: [x,y,3, 5 superordinate cat, 2 ns_loc, 3 obj_loc]
                        % CBlind: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 change images];
                        % Lures: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 lure images];
                        d = dir(sprintf('%s*.mat', params.stim.ns.stimfile));
                        a = load(fullfile(d(end).folder,d(end).name), 'scenes','wm_im','info'); % also contains 'ltm_lures' but we don't use those
                        all_images.ns.scenes = a.scenes; 
%                         all_images.ns.lures = a.ltm_lures; 
                        all_images.ns.wm_im = a.wm_im; 
                        all_images.info.ns = a.info; 
                        clear a d;
                    end
                    
                    % get indices for fourth, fifth, and sixth dimension of image
                    % array
                    i4 = run_table.super_cat(stim_row(ii));
                    i5 = run_table.basic_cat(stim_row(ii));
                    i6 = run_table.sub_cat(stim_row(ii));
                    
                    if strcmp(run_table.event_name(stim_row(ii)),'stim1') && run_table.is_catch(stim_row(ii)) == 0
                        
                        if any(strcmp(all_images.info.ns.Properties.VariableNames, 'img_quiz_dot_specialcore_stim_nr'))
                            idx0 = (all_images.info.ns.super_cat == run_table.super_cat(stim_row(ii),1) & ...
                                (all_images.info.ns.basic_cat == run_table.basic_cat(stim_row(ii),1) ) & ...
                                (all_images.info.ns.sub_cat == run_table.sub_cat(stim_row(ii),1)) & ...
                                (ismember(all_images.info.ns.unique_im, params.stim.ns.unique_im_nrs_core)) & ...
                                isnan(all_images.info.ns.img_quiz_dot_specialcore_stim_nr) );
                        else
                            idx0 = (all_images.info.ns.super_cat == run_table.super_cat(stim_row(ii),1) & ...
                                (all_images.info.ns.basic_cat == run_table.basic_cat(stim_row(ii),1) ) & ...
                                (all_images.info.ns.sub_cat == run_table.sub_cat(stim_row(ii),1)) & ...
                                (ismember(all_images.info.ns.unique_im, params.stim.ns.unique_im_nrs_core)));
                        end
                        
                        obj_super = all_images.info.ns.super_cat(idx0);
                        obj_basic = all_images.info.ns.basic_cat(idx0);
                        obj_sub   = all_images.info.ns.sub_cat(idx0);
                        
                        % check if stim idx matches
                        assert(isequal(all_images.info.ns.unique_im(idx0),unique_im));
                        assert(isequal(obj_super,run_table.super_cat(stim_row(ii),1)))
                        assert(isequal(obj_basic,run_table.basic_cat(stim_row(ii),1)))
                        assert(isequal(obj_sub,run_table.sub_cat(stim_row(ii),1)))
                        
                        % third dim has image and alpha mask
                        run_images{curr_frames(1),side} = all_images.ns.scenes(:,:,:,i4,i5,i6);
                        run_alpha_masks{curr_frames(1),side} = [];
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'wm') ...
                            && run_table.is_catch(stim_row(ii)) == 0
                        
                        wm_change = run_table.stim2_delta(stim_row(ii),1);
                        wm_im_nr  = run_table.stim2_im_nr(stim_row(ii), 1); % {'easy_added','hard_added','easy_removed','hard_removed'}
                        im_nr = (all_images.info.ns.unique_im == wm_im_nr);
                        info_name = params.stim.ns.change_im_name{wm_change == params.stim.ns.change_im};
                        info_name2 = strsplit(all_images.info.ns.filename{im_nr},'/');
                        
                        assert(strcmp(info_name2{2},sprintf('%s.png',info_name)));
                        assert(isequal(i4, run_table.super_cat(stim_row(ii),1)));
                        assert(isequal(i5, run_table.basic_cat(stim_row(ii),1)));
                        assert(isequal(i6, run_table.sub_cat(stim_row(ii),1)));
                        
                        run_images{curr_frames(1),side} = all_images.ns.wm_im(:,:,:,i4,i5,i6,wm_change == params.stim.ns.change_im);
                        run_alpha_masks{curr_frames(1),side} = [];
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'ltm') ...
                            && run_table.is_catch(stim_row(ii)) == 0
                        
                        if run_table.islure(stim_row(ii),1)
                            
                            lure_im = run_table.stim2_delta{stim_row(ii),1};
                            im_idx = find(strcmp(pair_im,params.stim.ns.ltm_lure_im)); % {'lure01','lure02', 'lure03', 'lure04'};
                            info_name = all_images.info.ns.(sprintf('lure%02d',im_idx));
                            pair_im = run_table.stim2_im_nr(stim_row(ii),1);
                            
                            assert(isequal(sprintf('%s.png',lure_im),info_name{im_idx}));
                            assert(isequal(i4, run_table.super_cat(stim_row(ii),1)));
                            assert(isequal(i5, run_table.basic_cat(stim_row(ii),1)));
                            assert(isequal(i6, run_table.sub_cat(stim_row(ii),1)));
                            
                            run_images{curr_frames(1),side} = all_images.ns.lures(:,:,:,i4,i5,i6,lure_im);
                            
                        else
                            pair_im = run_table.ltm_pair(stim_row(ii),1);
                            
                            run_images{curr_frames(1),side} = [];
                            run_alpha_masks{curr_frames(1),side} = [];
                        end
                        
                        
                    elseif run_table.is_catch(stim_row(ii)) == 1
                        % catch trial
                        run_images{curr_frames,side}      = uint8(zeros(1,1,1));
                        run_alpha_masks{curr_frames,side} = [];
                    end
            end % stim class
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%% APPLY CONTRAST DECREMENT %%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(run_table.task_class_name(stim_row(ii)),'cd')
                % 20% change we will actually apply the contrast
                % decrement change to stimulus (this probability is already
                % determined by vcd_addFIXandCDtoTimeTableMaster.m).
                if ~isnan(run_table.cd_start(stim_row(ii)))
                    if run_table.cd_start(stim_row(ii))~=0 && (1+mod(run_table.is_cued(stim_row(ii))-1,2))==side
                        stmclass   = run_table.stim_class_name{stim_row(ii),side};
                        rel_onset  = run_table.cd_start(stim_row(ii)) - run_table.event_start(stim_row(ii)) + 1;
                        pre_onset_frames = sum(params.stim.cd.t_cmodfun==params.stim.cd.t_cmodfun(1));
                        t_cmod_pad = length([ ((rel_onset + length(params.stim.cd.t_cmodfun)-pre_onset_frames - 1)+1) : params.exp.trial.stim_array_dur]);
                        if strcmp(stmclass,'rdk')
                            [f_im_cd, f_mask_cd] = vcd_applyContrastDecrement(params, ...
                                rel_onset, stmclass, run_images(curr_frames,side), ...
                                'input_mask', run_alpha_masks(curr_frames,side), ...
                                't_cmod_pad',t_cmod_pad); % give all frames
                        else
                            [f_im_cd, f_mask_cd] = vcd_applyContrastDecrement(params,...
                                rel_onset, stmclass, run_images(curr_frames(1),side), ...
                                'input_mask', run_alpha_masks(curr_frames(1),side), ...
                                't_cmod_pad',t_cmod_pad); % give first (and only frame)
                        end
                        run_images(curr_frames,side) = f_im_cd;
                        run_alpha_masks(curr_frames,side) = f_mask_cd;
                        clear f_im_cd f_mask_cd;
                       
                    end
                end
            end
            
        end % isempty
    end % side
end % stim idx


% Store stimuli and masks if requested
if savestim
    fprintf('[%s]: Storing run images..\n',mfilename)
    saveDir = fullfile(vcd_rootPath,'data',env_type,sprintf('vcd_subj%03d',subj_nr), sprintf('vcd_subj%03d_ses%02d',subj_nr, ses_nr));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    % grab eye, fixation and instructions images as well..
    eye_im         = all_images.eye; 
    fix_im.im      = all_images.fix;
    fix_im.alpha   = all_images.alpha.fix;
    fix_im.info    = all_images.info.fix;
    instr_im.im    = all_images.instr;
    instr_im.alpha = all_images.alpha.instr;
    instr_im.info  = all_images.info.instr;
    
    save(fullfile(saveDir, ...
        sprintf('subj%03d_ses%02d_%s_run%02d_images_%s.mat', ...
        subj_nr, ses_nr, choose(ses_type==1,'A','B'),run_nr, params.disp.name)), ...
        'run_images','run_alpha_masks','run_table','run_frames','eye_im','instr_im','fix_im','-v7.3')
end

% tell the user we are done and report time it took to load images
fprintf('Done! ');  toc;
fprintf('\n');

return

