function [run_images, run_alpha_masks] = vcd_getImageOrderSingleRun(params, ...
    time_table_master, all_run_frames, subj_nr, ses_nr, ses_type, run_nr, varargin)
% VCD function to load uint8 stimuli, alpha transparency masks, and
% corresponding image nrs for each trial
%
%   [all_run_images, all_run_alpha_masks,all_run_im_nr] =
%   vcd_getImageOrderSingleRun( ...
%       params, time_table_master, subj_nr, ses_nr, ses_type, run_nr, ...
%       ['images', <images>], ['load_params',<load_params>], ...
%       ['store_params',<store_params>],['session_env',<session_env>])
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
p0.addParameter('images'      , []  , @isstruct);
p0.addParameter('store_params', true, @islogical);
p0.addParameter('session_env', 'MRI', @(x) ismember(x,{'MRI','BEHAVIOR'}));

% Parse inputs
p0.parse(params, time_table_master, all_run_frames, subj_nr, ses_nr, ses_type, run_nr, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

%% Check subject, session and run nrs to process

% Check if subject nr is member of total subjects
subj_ok = find(ismember([1:params.exp.total_subjects],subj_nr));
if isempty(subj_ok)
    error('[%s]: Subject number is outside the range!\n',mfilename);
end

% Check if ses nr is member of all sessions
if strcmp(session_env, 'MRI')
    ses_ok = ismember([1:params.exp.session.n_total_sessions],ses_nr);
    run_ok = ismember([1:params.exp.session.deep.n_runs_per_session],run_nr);
    
elseif strcmp(session_env, 'BEHAVIOR')
    ses_ok = ismember([1:params.exp.session.n_behavioral_sessions],ses_nr);
    run_ok = ismember([1:params.exp.session.behavior.n_runs_per_session],run_nr);
end

if isempty(ses_ok)
    error('[%s]: Session number is outside the range!\n',mfilename);
end

% If left undefined, we will load images for first session
if isempty(run_ok)
    error('[%s]: Run number is outside the range!\n',mfilename);
end


%% Load params if requested and we can find the file

% Check if we need to load images
% If we did not provide images.. and params.images is empty. then
% create the field names of the images struct
if isempty(images) || isempty(fieldnames(images))
    images = struct('gabor',[],'rdk',[],'dot',[],'obj',[],'ns',[],...
        'fix',[], 'info',[], 'eye',[],'alpha',[]);
else
    % If we provided images, then add it to params struct
    % !!! WARNING that this will override the existing params.images field!)
    if isfield(params, 'images') &&  ~isempty(fieldnames(params.images))
        images = params.images;  params = rmfield(params, 'images');
    end
end

% Load stored fixation dot images if needed
if ~isfield(images,'fix') || isempty(images.fix)
    fprintf('[%s]: Loading fixation dot images..\n',mfilename);
    
    % FIX: 5D array: [x,y, 3, 5 lum, 2 widths]
    d = dir(sprintf('%s*.mat', params.stim.fix.stimfile));
    load(fullfile(d(end).folder,d(end).name), 'fix_im','mask','info');
    images.fix       = fix_im; clear fix_im;
    images.info.fix  = info; clear info;
    images.alpha.fix = mask; clear mask;
end

% Load stored fixation dot images if needed
if ~isfield(images,'eye') || isempty(images.eye)
    fprintf('[%s]: Loading eyetracking target images..\n',mfilename);
    
    % FIX: 5D array: [x,y, 3, 5 lum, 2 widths]
    d = dir(sprintf('%s*.mat', params.stim.el.stimfile));
    a = load(fullfile(d(end).folder,d(end).name), 'sac_im','pupil_im_white','pupil_im_black');
    images.eye.sac_im    = a.sac_im; 
    images.eye.pupil_im_white  = a.pupil_im_white;
    images.eye.pupil_im_black = a.pupil_im_black;
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
run_images      = cell(size(run_frames.frame_im_nr)); % second dim represent left and right side stimuli
run_alpha_masks = run_images;

% Get stimulus rows in time table
stim_events = run_table.event_id(ismember(run_table.event_id, [params.exp.block.stim_epoch1_ID, params.exp.block.stim_epoch2_ID, ...
                    params.exp.block.eye_gaze_fix_ID, params.exp.block.eye_gaze_sac_target_ID, ...
                    params.exp.block.eye_gaze_pupil_black_ID, params.exp.block.eye_gaze_pupil_white_ID]));
stim_row = find(ismember(run_table.event_id, [params.exp.block.stim_epoch1_ID, params.exp.block.stim_epoch2_ID, ...
                    params.exp.block.eye_gaze_fix_ID, params.exp.block.eye_gaze_sac_target_ID, ...
                    params.exp.block.eye_gaze_pupil_black_ID, params.exp.block.eye_gaze_pupil_white_ID]));


% Loop over all stimulus rows
fprintf('[%s]: Load stimuli for each trial..\n',mfilename); tic;
for ii = 1:length(stim_row)
    
    frame_counter = run_table.event_start(stim_row(ii))+1; % t=1 is 0, but we can't use 0 as index
    
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
                    fprintf('[%s]: Loading gabor images..\n',mfilename);
                    if isempty(images.gabor)
                        % GABORS: 6D array: [x,y,8 orient, 4 phase,3 contrast, og + 4 delta]
                        d = dir(sprintf('%s*.mat', params.stim.gabor.stimfile));
                        a = load(fullfile(d(end).folder,d(end).name), 'gabors','masks','info');
                        
                        images.gabor = a.gabors; 
                        images.info.gabor  = a.info; 
                        images.alpha.gabor = a.masks; 
                        clear a d;
                    end
                    
                    if strcmp(run_table.event_name(stim_row(ii)),'stim1') && run_table.is_catch(stim_row(ii)) == 0
                        
                        % GABORS: 6D array: [x,y,orient,contrast,phase,delta]
                        idx0 = find(images.info.gabor.unique_im==unique_im);
                        gbr_ori      = images.info.gabor.orient_deg(idx0);
                        gbr_contrast = images.info.gabor.contrast(idx0);
                        gbr_phase    = images.info.gabor.phase_deg(idx0);
                        ori_idx      = images.info.gabor.orient_i(idx0);
                        con_idx      = images.info.gabor.contrast_i(idx0);
                        % check if stim params match
                        assert(isequal( run_table.orient_dir(stim_row(ii),side) , gbr_ori));
                        assert(isequal( run_table.contrast(stim_row(ii),side), gbr_contrast));
                        assert(isequal( run_table.gbr_phase(stim_row(ii),side), gbr_phase));
                        
                        run_images{curr_frames(1),side} = images.gabor(:,:,:,ori_idx,con_idx,1);
                        run_alpha_masks{curr_frames(1),side} = images.alpha.gabor(:,:,ori_idx,con_idx,1);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'wm') && run_table.is_catch(stim_row(ii)) == 0
                        delta_deg = run_table.stim2_delta(stim_row(ii),side);
                        
                        idx0 = find( (images.info.gabor.orient_deg == (run_table.orient_dir(stim_row(ii),side))) & ...
                            (images.info.gabor.phase_deg==run_table.gbr_phase(stim_row(ii),side)) & ...
                            (images.info.gabor.contrast==run_table.contrast(stim_row(ii),side)) &...
                            (images.info.gabor.delta_deg==delta_deg) );
                        
                        ori_idx      = images.info.gabor.orient_i(idx0);
                        con_idx      = images.info.gabor.contrast_i(idx0);
                        test_im      = images.info.gabor.unique_im(idx0);
                        assert(isequal(test_im,run_table.stim2_im_nr(stim_row(ii),side)))
                        
                        % check if stim params match
                        assert(isequal( run_table.orient_dir(stim_row(ii),side) , images.info.gabor.orient_deg(idx0)));
                        assert(isequal( run_table.contrast(stim_row(ii),side),    images.info.gabor.contrast(idx0)));
                        assert(isequal( run_table.gbr_phase(stim_row(ii),side),   images.info.gabor.phase_deg(idx0)));
                        assert(isequal( delta_deg, images.info.gabor.delta_deg(idx0)));
                        
                        delta_idx = 1+ find(delta_deg == params.stim.gabor.delta_from_ref);
                        
                        run_images{curr_frames(1),side}      = images.gabor(:,:,:,ori_idx,con_idx, delta_idx);
                        run_alpha_masks{curr_frames(1),side} = images.alpha.gabor(:,:,ori_idx,con_idx, delta_idx);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'ltm') && run_table.is_catch(stim_row(ii)) == 0
                        
                        
                        if run_table.islure(stim_row(ii)) % same stim class
                            run_images{curr_frames(1),side} = [];
                            run_alpha_masks{curr_frames(1),side} = [];
                            
                        else % other stim class
                            run_images{curr_frames(1),side} = [];
                            run_alpha_masks{curr_frames(1),side} = [];
                        end
                        
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'img') && run_table.is_catch(stim_row(ii)) == 0
                        
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
                        run_alpha_masks{curr_frames,side}      = uint8(zeros(1,1,1));
                    end
                    
                    
                    
                    
                case 'rdk'
                    fprintf('[%s]: Loading rdk images..\n',mfilename);
                    
                    if strcmp(run_table.event_name(stim_row(ii)),'stim1') && run_table.is_catch(stim_row(ii)) == 0
                        
                        % RDKs: 130 mat files: 8 directions x 3 coherence levels x 5 deltas (0 + 4 deltas)
                        %                                     stimDir = dir(fullfile(sprintf('%s*',params.stim.rdk.stimfile)));
                        %                                     filename = sprintf('%04d_vcd_rdk_ori%02d_coh%02d_delta%02d',...
                        %                                         unique_im, ...
                        %                                         find(run_table_table.orient_dir(stim_row(ii),side) == params.stim.rdk.dots_direction), ...
                        %                                         find(run_table_table.rdk_coherence(stim_row(ii),side) == params.stim.rdk.dots_coherence), ...
                        %                                         0);
                        stimDir = dir(fullfile(sprintf('%s*',params.stim.rdk.stimfile)));
                        filename = sprintf('%04d_vcd_rdk_ori%02d_coh%02d_delta%02d',...
                            unique_im, ...
                            find(run_table.orient_dir(stim_row(ii),side) == params.stim.rdk.dots_direction), ...
                            find(run_table.rdk_coherence(stim_row(ii),side) == params.stim.rdk.dots_coherence), ...
                            0);
                        
                        stimfile = fullfile(stimDir(1).folder,stimDir(1).name,sprintf('%s.mat', filename));
                        if exist(stimfile,'file')
                            load(stimfile, 'frames','mask', 'rdk_info');
                        else
                            error('[%s]: Can''t find RDK stim file!')
                        end
                        
                        
                        %                                     if ~isfield(images,'info') || ~isfield(images.info,'rdk')
                        %                                         infofile = dir(fullfile(sprintf('%s*',params.stim.rdk.infofile)));
                        %                                         if ~isempty(infofile)
                        %                                             images.info.rdk = readtable(fullfile(infofile(end).folder,infofile(end).name));
                        %                                         else
                        %                                             error('[%s]: Can''t find RDK info file!')
                        %                                         end
                        %                                     end
                        
                        
                        % check if stim description matches
%                         idx0 = find(rdk_info.unique_im == unique_im);
%                         dot_motdir = rdk_info.dot_motdir_deg(idx0);
%                         dot_coh    = rdk_info.dot_coh(idx0);
%                         assert(isequal(run_table.rdk_coherence(stim_row(ii),side),dot_coh));
%                         assert(isequal(run_table.orient_dir(stim_row(ii),side),dot_motdir));
                        
                        frames =  frames(:,:,:,1:params.stim.rdk.duration);

                        % expand rdk movies into frames
                        rdk_images = squeeze(mat2cell(frames, size(frames,1), ...
                            size(frames,2), size(frames,3), ones(1,size(frames,4))));
                        
                        if size(rdk_images,1) < size(rdk_images,2)
                            rdk_images = rdk_images';
                        end
                        
                        rdk_masks = repmat({mask}, size(rdk_images,1), 1);
                        
                        run_images(curr_frames,side) = rdk_images;
                        run_alpha_masks(curr_frames,side) = rdk_masks;
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'wm') && run_table.is_catch(stim_row(ii)) == 0
                        
                        delta_idx = find(run_table.stim2_delta(stim_row(ii),side) == params.stim.rdk.delta_from_ref);
                        
                        
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
                        og_ori     = updated_ori-delta_test;
                        
                        % check if stim description matches
%                         assert(isequal(rdk_info.dot_coh(corresponding_unique_im_idx), run_table.rdk_coherence(stim_row(ii),side)));
%                         assert(isequal(rdk_info.dot_motdir_deg(corresponding_unique_im_idx), og_ori));
                        
                        % RDKs: 130 mat files: 8 directions x 3 coherence levels x 5 deltas (0 + 4 deltas)
                        stimDir = dir(fullfile(sprintf('%s*',params.stim.rdk.stimfile)));
                        filename = sprintf('%04d_vcd_rdk_ori%02d_coh%02d_delta%02d',...
                            test_im, ...
                            find(og_ori == params.stim.rdk.dots_direction), ...
                            find(run_table.rdk_coherence(stim_row(ii),side) == params.stim.rdk.dots_coherence), ...
                            delta_idx);
                        
                        stimfile = fullfile(stimDir(1).folder,stimDir(1).name,sprintf('%s.mat', filename));
                        if exist(stimfile,'file')
                            load(stimfile, 'frames','mask');
                        else
                            error('[%s]: Can''t find RDK stim file!')
                        end
                        
                        % ensure duration
                        frames =  frames(:,:,:,1:params.stim.rdk.duration);

                        
                        % expand rdk movies into frames
                        rdk_images = squeeze(mat2cell(frames, size(frames,1), ...
                            size(frames,2), size(frames,3), ones(1,size(frames,4))));
                        
                        if size(rdk_images,1) < size(rdk_images,2)
                            rdk_images = rdk_images';
                        end
                        
                        rdk_masks = repmat({mask}, size(rdk_images,1), 1);
                        
                        run_images(curr_frames,side) = rdk_images;
                        run_alpha_masks(curr_frames,side) = rdk_masks;
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'ltm') && run_table.is_catch(stim_row(ii)) == 0
                        
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
                            run_images(curr_frames,side) = []; %rdk_images;
                            run_alpha_masks(curr_frames,side) = []; %rdk_masks;
                            
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
                            run_images(curr_frames,side) = []; %rdk_images;
                            run_alpha_masks(curr_frames,side) = []; %rdk_masks;
                            
                        end
                        
                        
                    elseif run_table.is_catch(stim_row(ii)) == 1
                        % catch trial
                        run_images(curr_frames,side)      = uint8(zeros(1,1,1));
                        run_alpha_masks(curr_frames,side) = uint8(zeros(1,1));
                    end
                    
                case 'dot'
                    fprintf('[%s]: Loading dot images..\n',mfilename);
                    if isempty(images.dot)
                        % Simple dot: 2D array: [x,y]
                        d = dir(sprintf('%s*.mat', params.stim.dot.stimfile));
                        a=load(fullfile(d(end).folder,d(end).name), 'dot','mask','info');
                        images.dot = a.dot; 
                        images.alpha.dot = a.mask; 
                        images.info.dot = a.info; 
                        clear a d;
                    end
                    
                    
                    if run_table.is_catch(stim_row(ii)) == 0
                        run_images{curr_frames(1),side} = images.dot;
                        run_alpha_masks{curr_frames(1),side} = images.alpha.dot;
                        
                        if strcmp(run_table.event_name(stim_row(ii)),'stim2')
                            test_im = run_table.stim2_im_nr(stim_row(ii),side);
                        end
                    elseif run_table.is_catch(stim_row(ii)) == 1
                        % catch trial
                        run_images{curr_frames,side}      = uint8(zeros(1,1,1));
                        run_alpha_masks{curr_frames,side} = uint8(zeros(1,1));
                    end
                    
                case 'obj'
                    fprintf('[%s]: Loading object images..\n',mfilename);
                    if isempty(images.obj)
                        % Complex objects: 4D array: [x,y,16 object, og + 10 rotation]
                        d = dir(sprintf('%s*.mat', params.stim.obj.stimfile));
                        a = load(fullfile(d(end).folder,d(end).name), 'objects','masks','info');
                        images.obj = a.objects; 
                        images.alpha.obj = a.masks; 
                        images.info.obj = a.info; 
                        clear a d;
                    end
                    
                    if strcmp(run_table.event_name(stim_row(ii)),'stim1') && run_table.is_catch(stim_row(ii)) == 0
                        
                        idx = find( (images.info.obj.super_cat == run_table.super_cat(stim_row(ii),side)) & ...
                            (images.info.obj.basic_cat == run_table.basic_cat(stim_row(ii),side)) & ...
                            (images.info.obj.sub_cat == run_table.sub_cat(stim_row(ii),side)) & ...
                            (images.info.obj.abs_rot == run_table.orient_dir(stim_row(ii),side)) );
                        
                        obj_super = images.info.obj.super_cat_name(idx);
                        obj_basic = images.info.obj.basic_cat_name(idx);
                        obj_sub   = images.info.obj.sub_cat_name(idx);
                        obj_rotation = images.info.obj.abs_rot(idx);
                        
                        % check if stim idx matches
                        assert(isequal(images.info.obj.unique_im(idx),     unique_im));
                        assert(isequal(obj_super,   run_table.super_cat_name(stim_row(ii),side)));
                        assert(isequal(obj_basic,   run_table.basic_cat_name(stim_row(ii),side)));
                        assert(isequal(obj_sub,     run_table.sub_cat_name(stim_row(ii),side)));
                        assert(isequal(obj_rotation,run_table.orient_dir(stim_row(ii),side)));
                        
                        run_images{curr_frames(1),side}      = images.obj(:,:,:,unique_im==params.stim.obj.unique_im_nrs_core,1);
                        run_alpha_masks{curr_frames(1),side} = images.alpha.obj(:,:,unique_im==params.stim.obj.unique_im_nrs_core,1);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'wm') && run_table.is_catch(stim_row(ii)) == 0
                        
                        delta_idx  = run_table.stim2_delta(stim_row(ii),side);
                        ref_dir    = run_table.orient_dir(stim_row(ii),side);
                        test_dir   = run_table.stim2_orient_dir(stim_row(ii),side);
                        delta      = test_dir - ref_dir;
                        
                        % dangerous.. grab stim 1 unique im nr
                        if side == 1
                            corresponding_unique_im = run_table.stim_nr_left(stim_row(ii))==params.stim.obj.unique_im_nrs_core;
                        elseif side == 2
                            corresponding_unique_im = run_table.stim_nr_right(stim_row(ii))==params.stim.obj.unique_im_nrs_core;
                        end
                        % check if stim idx matches
                        idx = find( (images.info.obj.super_cat == run_table.super_cat(stim_row(ii),side)) & ...
                            (images.info.obj.basic_cat == run_table.basic_cat(stim_row(ii),side)) & ...
                            (images.info.obj.sub_cat == run_table.sub_cat(stim_row(ii),side)) & ...
                            (images.info.obj.abs_rot == ref_dir) );
                        
                        test_im = run_table.stim2_im_nr(stim_row(ii),side);
                        delta_idx0 = find(delta==params.stim.obj.delta_from_ref);
                        assert(strcmp(run_table.sub_cat_name(stim_row(ii),side), images.info.obj.sub_cat_name(idx)));
                        tmp_test_im = reshape(params.stim.obj.unique_im_nrs_wm_test,4,[]);
                        assert(isequal(test_im,tmp_test_im(delta_idx0,corresponding_unique_im)))
                        run_images{curr_frames(1),side} = images.obj(:,:,:,corresponding_unique_im,delta_idx0);
                        run_alpha_masks{curr_frames(1),side} = images.alpha.obj(:,:,corresponding_unique_im,delta_idx0);
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'ltm') && run_table.is_catch(stim_row(ii)) == 0
                        
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
                    fprintf('[%s]: Loading scene images..\n',mfilename);
                    if isempty(images.ns)
                        % NS: 6D array: [x,y,3, 5 superordinate cat, 2 ns_loc, 3 obj_loc]
                        % CBlind: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 change images];
                        % Lures: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 lure images];
                        d = dir(sprintf('%s*.mat', params.stim.ns.stimfile));
                        a = load(fullfile(d(end).folder,d(end).name), 'scenes','ltm_lures','wm_im','info');
                        images.ns.scenes = a.scenes; 
                        images.ns.lures = a.ltm_lures; 
                        images.ns.wm_im = a.wm_im; 
                        images.info.ns = a.info; 
                        clear a d;
                    end
                    
                    % get indices for fourth, fifth, and sixth dimension of image
                    % array
                    i4 = run_table.super_cat(stim_row(ii));
                    i5 = run_table.basic_cat(stim_row(ii));
                    i6 = run_table.sub_cat(stim_row(ii));
                    
                    if strcmp(run_table.event_name(stim_row(ii)),'stim1') && run_table.is_catch(stim_row(ii)) == 0
                        idx0 = (images.info.ns.super_cat == run_table.super_cat(stim_row(ii),1) & ...
                            (images.info.ns.basic_cat == run_table.basic_cat(stim_row(ii),1) ) & ...
                            (images.info.ns.sub_cat == run_table.sub_cat(stim_row(ii),1)) & ...
                            (ismember(images.info.ns.unique_im, params.stim.ns.unique_im_nrs_core)));
                        
                        obj_super = images.info.ns.super_cat(idx0);
                        obj_basic = images.info.ns.basic_cat(idx0);
                        obj_sub   = images.info.ns.sub_cat(idx0);
                        
                        % check if stim idx matches
                        assert(isequal(images.info.ns.unique_im(idx0),unique_im));
                        assert(isequal(obj_super,run_table.super_cat(stim_row(ii),1)))
                        assert(isequal(obj_basic,run_table.basic_cat(stim_row(ii),1)))
                        assert(isequal(obj_sub,run_table.sub_cat(stim_row(ii),1)))
                        
                        % third dim has image and alpha mask
                        run_images{curr_frames(1),side} = images.ns.scenes(:,:,:,i4,i5,i6);
                        run_alpha_masks{curr_frames(1),side} = [];
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'wm') && run_table.is_catch(stim_row(ii)) == 0
                        
                        wm_change = run_table.stim2_delta(stim_row(ii),1);
                        wm_im_nr = run_table.stim2_im_nr(stim_row(ii), 1); % {'easy_added','hard_added','easy_removed','hard_removed'}
                        im_nr = (images.info.ns.unique_im== wm_im_nr);
                        info_name = params.stim.ns.change_im_name{wm_change};
                        info_name2 = strsplit(images.info.ns.filename{im_nr},'/');
                        
                        assert(strcmp(info_name2{2},sprintf('%s.png',info_name)));
                        assert(isequal(i4, run_table.super_cat(stim_row(ii),1)));
                        assert(isequal(i5, run_table.basic_cat(stim_row(ii),1)));
                        assert(isequal(i6, run_table.sub_cat(stim_row(ii),1)));
                        
                        run_images{curr_frames(1),side} = images.ns.wm_im(:,:,:,i4,i5,i6,wm_change);
                        run_alpha_masks{curr_frames(1),side} = [];
                        
                    elseif strcmp(run_table.event_name(stim_row(ii)),'stim2') && strcmp(run_table.task_class_name(stim_row(ii)),'ltm') && run_table.is_catch(stim_row(ii)) == 0
                        
                        if run_table.islure(stim_row(ii),1)
                            
                            lure_im = run_table.stim2_delta{stim_row(ii),1};
                            im_idx = find(strcmp(pair_im,params.stim.ns.ltm_lure_im)); % {'lure01','lure02', 'lure03', 'lure04'};
                            info_name = images.info.ns.(sprintf('lure%02d',im_idx));
                            pair_im = run_table.stim2_im_nr(stim_row(ii),1);
                            
                            assert(isequal(sprintf('%s.png',lure_im),info_name{im_idx}));
                            assert(isequal(i4, run_table.super_cat(stim_row(ii),1)));
                            assert(isequal(i5, run_table.basic_cat(stim_row(ii),1)));
                            assert(isequal(i6, run_table.sub_cat(stim_row(ii),1)));
                            
                            run_images{curr_frames(1),side} = images.ns.lures(:,:,:,i4,i5,i6,lure_im);
                            
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
            
            % APPLY CONTRAST DECREMENT
            if strcmp(run_table.task_class_name(stim_row(ii)),'cd')
                % 50% change we will actually apply the contrast
                % decrement change to stimulus
  
                if run_table.cd_start(stim_row(ii),side)~=0
                    stmclass = run_table.stim_class_name{stim_row(ii),side};
                    rel_onset = run_table.cd_start(stim_row(ii),side) - run_table.event_start(stim_row(ii));
                    
                    if strcmp(stmclass,'rdk')
                        f_im_cd = vcd_applyContrastDecrement(params, rel_onset, stmclass, run_images(curr_frames,side)); % give all frames
                        run_images(curr_frames,side) = f_im_cd;
                    else
                        f_im_cd = vcd_applyContrastDecrement(params, rel_onset, stmclass, run_images(curr_frames(1),side)); % give first (and only frame)
                        run_images(curr_frames(1),side) = f_im_cd;
                    end
                    clear f_im_cd;
                end
            end
            
        end % isempty
    end % side
end % stim idx

% delete empty rows
%             empty_rows = cellfun(@isempty, run_images(:,1),'UniformOutput',false);
%             delete_me = logical(cell2mat(empty_rows));
%             run_images(delete_me,:) = [];
%             run_alpha_masks(delete_me,:) = [];

% Store structs if requested
if store_params
    fprintf('[%s]: Storing trial data..\n',mfilename)
    saveDir = fullfile(vcd_rootPath,'workspaces','info',sprintf('subj%03d',subj_nr));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir, ...
        sprintf('subj%03d_ses%02d_%s_run%02d_images_%s.mat', ...
        subj_nr, ses_nr, choose(ses_type==1,'A','B'),run_nr, params.disp.name)), ...
        'run_images','run_alpha_masks','-v7.3')
end

fprintf('Done! ');  toc;
fprintf('\n');
% end % load params


return

