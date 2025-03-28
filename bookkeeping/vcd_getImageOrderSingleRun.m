function [run_images, run_alpha_masks] = vcd_getImageOrderSingleRun(params, subj_time_table, subj_nr, session_nr, run_nr, varargin)
% VCD function to load in the RGB images for each trial
%
%   [run_images, run_alpha_masks] = vcd_getImageOrderSingleRun( ...
%       params, subj_time_table, subj_nr, session_nr, run_nr, 'images',[images])

% Note image files are stored locally (workspaces > stimuli), and need to
% be created prior to running this function (see s_createStim.m) for every
% monitor display setup, such that it has the <params.disp.name> in the
% stimulus file name.


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'             , @isstruct);
p0.addRequired('subj_time_table'    , @istable);
p0.addRequired('subj_nr'            , @isnumeric);
p0.addRequired('session_nr'         , @isnumeric);
p0.addRequired('run_nr'             , @isnumeric);
p0.addParameter('images'      , []  , @isstruct);
p0.addParameter('load_params' , true, @islogical);
p0.addParameter('store_params', true, @islogical);

% Parse inputs
p0.parse(params,subj_time_table,subj_nr, session_nr, run_nr, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0


%% Load params if requested and we can find the file
if load_params
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('image_order_%s*.mat',params.disp.name)));
    fprintf('[%s]: Found %d image_order .mat file(s)\n',mfilename,length(d));
    if ~isempty(d)
        if length(d) > 1
            warning('[%s]: Multiple .mat files! Will pick the most recent one', mfilename);
        end
        fprintf('[%s]: Loading image_order .mat file: %s\n', mfilename, d.name);
        load(fullfile(d(end).folder,d(end).name));
    else
        error('[%s]: Can''t find image_order file!', mfilename)
    end
else
    
    
    %% Predefine cell arrays for stimuli
    run_images = cell(size(subj_time_table,1),2); % second dim represent left and right side stimuli
    run_alpha_masks = run_images;
    
    % Check if we need to load images
    if isfield(params, 'images')
        images = params.images;  params = rmfield(params, 'images');
    end
    
    if isempty(images) || isempty(fieldnames(images))
%         images = vcd_loadAllUniqueImages(params);
        images = struct('gabor',[],'rdk',[],'dot',[],'obj',[],'ns',[],...
            'fix',[], 'info',[], 'image_order',[],'alpha',[]);
    end
    
    
    % Load stored fixation dot images
    if ~isfield(images,'fix') || isempty(images.fix)
        fprintf('[%s]: Loading fixation dot images..\n',mfilename);
        
        % FIX: 5D array: [x,y, 3, 5 lum, 2 widths]
        d = dir(sprintf('%s*.mat', params.stim.fix.stimfile));
        load(fullfile(d(end).folder,d(end).name), 'fix_im','mask','info');
        images.fix       = fix_im; clear fix_im;
        images.info.fix  = info; clear info;
        images.alpha.fix = mask; clear mask;
    end

    
    % Get stimulus rows in time table
    stim_idx = (strcmp(subj_time_table.event_name,'stim1') | strcmp(subj_time_table.event_name,'stim2'));
    stim_sub = find(stim_idx);
    trial_nr = 1;
    
    % Loop over all stimulus rows
    fprintf('[%s]: Load stimuli for each trial..\n',mfilename); tic;
    for ii = 1:length(stim_sub)

        % Get stimulus class name
        stimClass = subj_time_table.stim_class_name{stim_sub(ii)};
        
        % Get stimulus location
        if subj_time_table.stimloc(stim_sub(ii))==1  % left
            side = 1;
        elseif subj_time_table.stimloc(stim_sub(ii))==2 % right
            side = 2;
        else % center (for scenes)
            side = 1;
        end

        switch stimClass
            case 'gabor'
                
                if isempty(images.gabor)
                    fprintf('[%s]: Loading gabor images..\n',mfilename);
                    % GABORS: 6D array: [x,y,8 orient, 4 phase,3 contrast, og + 4 delta]
                    d = dir(sprintf('%s*.mat', params.stim.gabor.stimfile));
                    load(fullfile(d(end).folder,d(end).name), 'gabors','masks','info');
                    
                    images.gabor = gabors; clear gabors;
                    images.info.gabor  = info; clear info;
                    images.alpha.gabor = masks; clear masks;
                end
                
                if strcmp(subj_time_table.event_name(stim_sub(ii)),'stim1')
                    
                    % GABORS: 6D array: [x,y,orient,contrast,phase,delta]
                    unique_im = subj_time_table.unique_im_nr(stim_sub(ii));
                    
                    idx0 = find(images.info.gabor.unique_im==unique_im);
                    gbr_ori      = images.info.gabor.ori_deg(idx0);
                    gbr_contrast = images.info.gabor.contrast(idx0);
                    gbr_phase    = images.info.gabor.phase_deg(idx0);
                    
                    % check if stim params match
                    assert(isequal( subj_time_table.orient_dir(stim_sub(ii)) , gbr_ori));
                    assert(isequal( subj_time_table.contrast(stim_sub(ii)), gbr_contrast));
                    assert(isequal( subj_time_table.gbr_phase(stim_sub(ii)), gbr_phase));
                    
                    
                    run_images{trial_nr,side} = images.gabor(:,:,:,unique_im,1);
                    run_alpha_masks{trial_nr,side} = images.alpha.gabor(:,:,unique_im,1);
                    
                elseif strcmp(subj_time_table.event_name(stim_sub(ii)),'stim2') && strcmp(subj_time_table.task_class_name(stim_sub(ii)),'wm')
                    delta_deg = subj_time_table.stim2_delta(stim_sub(ii));
                    
                    idx0 = find( (images.info.gabor.ori_deg == (subj_time_table.orient_dir(stim_sub(ii)) - delta_deg)) & ...
                        (images.info.gabor.phase_deg==subj_time_table.gbr_phase(stim_sub(ii))) & ...
                        (images.info.gabor.contrast==subj_time_table.contrast(stim_sub(ii))) &...
                        (images.info.gabor.delta_deg==delta_deg) );
                    
                    % check if stim params match
                    assert(isequal( subj_time_table.orient_dir(stim_sub(ii)) , images.info.gabor.ori_deg(idx0) + images.info.gabor.delta_deg(idx0)));
                    assert(isequal( subj_time_table.contrast(stim_sub(ii)),    images.info.gabor.contrast(idx0)));
                    assert(isequal( subj_time_table.gbr_phase(stim_sub(ii)),   images.info.gabor.phase_deg(idx0)));
                    assert(isequal( delta_deg, images.info.gabor.delta_deg(idx0)));                    
                    
                    delta_idx = 1+ find(delta_deg == params.stim.gabor.delta_from_ref);
                    run_images{trial_nr,side} = images.gabor(:,:,:,unique_im, delta_idx);
                    run_alpha_masks{trial_nr,side} = images.alpha.gabor(:,:,unique_im, delta_idx);
                    
                elseif strcmp(subj_time_table.event_name(stim_sub(ii)),'stim2') && strcmp(subj_time_table.task_class_name(stim_sub(ii)),'ltm')
                    
                    if subj_time_table.islure(stim_sub(ii)) % same stim class
                        pair_im = [];
                        run_images{trial_nr,side} = [];
                        run_alpha_masks{trial_nr,side} = [];
                        
                    else % other stim class
                        pair_im = [];
                        run_images{trial_nr,side} = [];
                        run_alpha_masks{trial_nr,side} = [];
                    end
                end
                
            case 'rdk'
                fprintf('[%s]: Loading rdk images..\n',mfilename);

                if strcmp(subj_time_table.event_name(stim_sub(ii)),'stim1')
                    
                    if ~isfield(images,'info') || ~isfield(images.info,'rdk')
                        infofile = dir(fullfile(sprintf('%s*',params.stim.rdk.infofile)));
                        if ~isempty(infofile)
                            images.info.rdk = readtable(fullfile(infofile.folder,infofile.name));
                        else
                            error('[%s]: Can''t find RDK info file!')
                        end
                    end
                    % check if stim description matches
                    idx0 = find(images.info.rdk.unique_im == subj_time_table.unique_im_nr(stim_sub(ii)));
                    dot_motdir = images.info.rdk.dot_dir(idx0);
                    dot_coh    = images.info.rdk.dot_coh(idx0);
                    assert(isequal(subj_time_table.rdk_coherence(stim_sub(ii)),dot_coh));
                    assert(isequal(subj_time_table.orient_dir(stim_sub(ii)),dot_motdir));
                    
                    % RDKs: 130 mat files: 8 directions x 3 coherence levels x 5 deltas (0 + 4 deltas)
                    stimDir = dir(fullfile(sprintf('%s*',params.stim.rdk.stimfile)));
                    filename = sprintf('%02d_rdk_ori%02d_coh%02d_delta%02d',...
                        subj_time_table.unique_im_nr(stim_sub(ii)), ...
                        find(subj_time_table.orient_dir(stim_sub(ii)) == params.stim.rdk.dots_direction), ...
                        find(subj_time_table.rdk_coherence(stim_sub(ii)) == params.stim.rdk.dots_coherence), ...
                        0);
                    
                    stimfile = fullfile(stimDir(1).folder,stimDir(1).name,sprintf('%s.mat', filename));
                    if exist(stimfile,'file')
                        load(stimfile, 'frames','mask');
                    else
                        error('[%s]: Can''t find RDK stim file!')
                    end
                    
                    % each file contains a 4D array: [x,y,3,frames]
                    run_images{trial_nr,side} = frames;
                    run_alpha_masks{trial_nr,side} = mask;
                    
                elseif strcmp(subj_time_table.event_name(stim_sub(ii)),'stim2') && strcmp(subj_time_table.task_class_name(stim_sub(ii)),'wm')
                    
                    delta_idx = subj_time_table.stim2_delta(stim_sub(ii));
                    
                    % check if stim description matches
                    assert(equal(delta_idx, find(subj_time_table.stim2_delta(stim_sub(ii)) == params.stim.rdk.delta_from_ref)));
                    
                    % RDKs: 130 mat files: 8 directions x 3 coherence levels x 5 deltas (0 + 4 deltas)
                    stimDir = dir(fullfile(sprintf('%s*',params.stim.rdk.stimfile)));
                    filename = sprintf('%02d_rdk_ori%02d_coh%02d_delta%02d',...
                        subj_time_table.unique_im_nr(stim_sub(ii)), ...
                        find(subj_time_table.orient_dir(stim_sub(ii)) == params.stim.rdk.dots_direction), ...
                        find(subj_time_table.rdk_coherence(stim_sub(ii)) == params.stim.rdk.dots_coherence), ...
                        delta_idx);
                    
                    stimfile = fullfile(stimDir(1).folder,stimDir(1).name,sprintf('%s.mat', filename));
                    if exist(stimfile,'file')
                        load(stimfile, 'frames','mask');
                    else
                        error('[%s]: Can''t find RDK stim file!')
                    end
                    
                    % each file contains a 4D array: [x,y,3,frames]
                    run_images{trial_nr,side} = frames;
                    run_alpha_masks{trial_nr,side} = mask;
                    
                elseif strcmp(subj_time_table.event_name(stim_sub(ii)),'stim2') && strcmp(subj_time_table.task_class_name(stim_sub(ii)),'ltm')
                    
                    if subj_time_table.islure(stim_sub(ii)) % same stim class
                        pair_im = [];
                        run_images{trial_nr,side} = [];
                        run_alpha_masks{trial_nr,side} = [];
                        
                    else % other stim class
                        pair_im = [];
                        run_images{trial_nr,side} = [];
                        run_alpha_masks{trial_nr,side} = [];
                    end
                    
                end
                
                
                
            case 'dot'
                
                if isempty(images.dot)
                    % Simple dot: 2D array: [x,y]
                    d = dir(sprintf('%s*.mat', params.stim.dot.stimfile));
                    load(fullfile(d(end).folder,d(end).name), 'simple_dot','mask','info');
                    images.dot = simple_dot; clear simple_dot;
                    images.alpha.dot = mask; clear mask;
                    images.info.dot = info; clear info;
                end
                
                run_images{trial_nr,side} = images.dot;
                run_alpha_masks{trial_nr,side} = images.alpha.dot;
                
                
            case 'obj'
                
                if isempty(images.obj)
                    fprintf('[%s]: Loading object images..\n',mfilename);
                    % Complex objects: 4D array: [x,y,16 object, og + 10 rotation]
                    d = dir(sprintf('%s*.mat', params.stim.obj.stimfile));
                    load(fullfile(d(end).folder,d(end).name), 'objects','masks','info','im_order');
                    images.obj = objects; clear objects;
                    images.alpha.obj = masks; clear masks;
                    images.info.obj = info; clear info;
                end
                
                unique_im = subj_time_table.unique_im_nr(stim_sub(ii));
                
                if strcmp(subj_time_table.event_name(stim_sub(ii)),'stim1')
                    
                    
                    idx = find( (image_info.obj.superordinate_i == subj_time_table.super_cat(stim_sub(ii))) & ...
                        (image_info.obj.basic_i == subj_time_table.basic_cat(stim_sub(ii))) & ...
                        (image_info.obj.subordinate_i == subj_time_table.sub_cat(stim_sub(ii))) & ...
                        (image_info.obj.rot_abs == subj_time_table.orient_dir(stim_sub(ii))) );
                    
                    obj_super = images.info.obj.superordinate(idx);
                    obj_basic = images.info.obj.basic(idx);
                    obj_sub   = images.info.obj.subordinate(idx);
                    obj_rotation = images.info.obj.rot_abs(idx);
                    
                    % check if stim idx matches
                    assert(isequal(idx,         unique_im));
                    assert(isequal(obj_super,   subj_time_table.super_cat_name(stim_sub(ii))));
                    assert(isequal(obj_basic,   subj_time_table.basic_cat_name(stim_sub(ii))));
                    assert(isequal(obj_sub,     subj_time_table.sub_cat_name(stim_sub(ii))));
                    assert(isequal(obj_rotation,subj_time_table.orient_dir(stim_sub(ii))));
                    
                    run_images{trial_nr,side} = images.obj(:,:,:,unique_im,1);
                    run_alpha_masks{trial_nr,side} = images.alpha.obj(:,:,unique_im,1);
                    
                elseif strcmp(subj_time_table.event_name(stim_sub(ii)),'stim2') && strcmp(subj_time_table.task_class_name(stim_sub(ii)),'wm')
                    
                    delta_idx = subj_time_table.stim2_delta(stim_sub(ii));
                    delta     = subj_time_table.stim2_delta(stim_sub(ii));
                    
                    % check if stim idx matches
                    idx = find( (image_info.obj.superordinate_i == subj_time_table.super_cat(stim_sub(ii))) & ...
                        (image_info.obj.basic_i == subj_time_table.basic_cat(stim_sub(ii))) & ...
                        (image_info.obj.subordinate_i == subj_time_table.sub_cat(stim_sub(ii))) & ...
                        (image_info.obj.rot_abs == subj_time_table.orient_dir(stim_sub(ii))) );
                    
                    delta_idx0 = find(delta,params.stim.obj.delta_from_ref);
                    %                 delta_idx1 = find(delta, [-10:2:10]);
                    assert(isequal(delta, images.info.obj(idx,delta_idx0)));
                    
                    run_images{trial_nr,side} = images.obj(:,:,:,unique_im,delta_idx);
                    run_alpha_masks{trial_nr,side} = images.alpha.obj(:,:,unique_im,delta_idx);
                    
                elseif strcmp(subj_time_table.event_name(stim_sub(ii)),'stim2') && strcmp(subj_time_table.task_class_name(stim_sub(ii)),'ltm')
                    
                    if subj_time_table.islure(stim_sub(ii)) % same stim class
                        pair_im = [];
                        run_images{trial_nr,side} = [];
                        run_alpha_masks{trial_nr,side} = [];
                        
                    else % other stim class
                        pair_im = [];
                        run_images{trial_nr,side} = [];
                        run_alpha_masks{trial_nr,side} = [];
                    end
                    
                end
                
                
            case 'ns'
                
                if isempty(images.ns)
                    fprintf('[%s]: Loading scene images..\n',mfilename);
                    % NS: 6D array: [x,y,3, 5 superordinate cat, 2 ns_loc, 3 obj_loc]
                    % CBlind: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 change images];
                    % Lures: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 lure images];
                    d = dir(sprintf('%s*.mat', params.stim.ns.stimfile));
                    load(fullfile(d(end).folder,d(end).name), 'scenes','lures','cblind','info','im_order');
                    images.ns = scenes; clear scenes;
                    images.lures = lures; clear lures;
                    images.cblind = cblind; clear cblind;
                    images.info.ns = info; clear info;
                end
                
                unique_im = subj_time_table.unique_im_nr(stim_sub(ii));
                
                % get indices for fourth, fifth, and sixth dimension of image
                % array
                i4 = subj_time_table.super_cat(stim_sub(ii));
                i5 = subj_time_table.basic_cat(stim_sub(ii));
                i6 = subj_time_table.sub_cat(stim_sub(ii));
                
                if strcmp(subj_time_table.event_name(stim_sub(ii)),'stim1')
                    
                    idx0 = (images.info.ns.superordinate_i == subj_time_table.super_cat(stim_sub(ii)) & ...
                        (images.info.ns.ns_loc_i == subj_time_table.basic_cat(stim_sub(ii)) ) & ...
                        (images.info.ns.obj_loc_i == subj_time_table.sub_cat(stim_sub(ii))) );
                    
                    obj_super = images.info.ns.superordinate_i(idx0);
                    obj_basic = images.info.ns.ns_loc_i(idx0);
                    obj_sub   = images.info.ns.obj_loc_i(idx0);
                    
                    % check if stim idx matches
                    assert(isequal(find(idx0),subj_time_table.unique_im_nr(stim_sub(ii))));
                    assert(isequal(obj_super,subj_time_table.super_cat(stim_sub(ii))))
                    assert(isequal(obj_basic,subj_time_table.basic_cat(stim_sub(ii))))
                    assert(isequal(obj_sub,subj_time_table.sub_cat(stim_sub(ii))))
                    
                    % third dim has image and alpha mask
                    run_images{trial_nr,side} = images.ns(:,:,:,i4,i5,i6);
                    run_alpha_masks{trial_nr,side} = [];
                    
                elseif strcmp(subj_time_table.event_name(stim_sub(ii)),'stim2') && strcmp(subj_time_table.task_class_name(stim_sub(ii)),'wm')
                    
                    change_im = subj_time_table.stim2_delta(stim_sub(ii));
                    change_im_name = subj_time_table.stim2{stim_sub(ii)};

                    im_idx    = find(strcmp(change_im_name,params.stim.ns.change_im)); % {'easy_added','hard_added','easy_removed','hard_removed'}
                    info_name = images.info.ns.(sprintf('change_img%02d',im_idx));
                    assert(isequal(sprintf('%s.png',change_im_name),info_name{im_idx}));
                    
                    run_images{trial_nr,side} = images.cblind(:,:,:,i4,i5,i6,im_idx);
                    run_alpha_masks{trial_nr,side} = [];
                    
                elseif strcmp(subj_time_table.event_name(stim_sub(ii)),'stim2') && strcmp(subj_time_table.task_class_name(stim_sub(ii)),'ltm')
                    
                    if subj_time_table.islure(stim_sub(ii))
                        lure_im = subj_time_table.stim2_delta{stim_sub(ii)};
                        im_idx = find(strcmp(pair_im,params.stim.ns.lure_im)); % {'lure01','lure02', 'lure03', 'lure04'};
                        info_name = images.info.ns.(sprintf('lure%02d',im_idx));
                        assert(isequal(sprintf('%s.png',lure_im),info_name{im_idx}));
                        
                        run_images{trial_nr,side} = images.lures(:,:,:,i4,i5,i6,im_idx);
                        
                    else
                        pair_im = subj_time_table.ltm_pair(stim_sub(ii));
                        
                        run_images{trial_nr,side} = [];
                        run_alpha_masks{trial_nr,side} = [];
                    end
                end
        end
        
        if side == 2 
            trial_nr = trial_nr + 1;
        elseif subj_time_table.stimloc(stim_sub(ii))==3 && ...
                subj_time_table.block_local_trial_nr(stim_sub(ii))~=subj_time_table.block_local_trial_nr(stim_sub(ii)+1)
            trial_nr = trial_nr + 1;
        end
    end
    
    % delete empty rows
    empty_rows = cellfun(@isempty, run_images(:,1),'UniformOutput',false);
    delete_me = logical(cell2mat(empty_rows));
    run_images(delete_me,:) = [];
    run_alpha_masks(delete_me,:) = [];
    
    
    % Store structs if requested
    if params.store_params
        fprintf('[%s]: Storing trial data..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info',sprintf('subj%03d',subj_nr));
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir, ...
            sprintf('subj%03d_ses%02d_run%02d_images_%s_%s.mat', ...
            subj_nr, session_nr, run_nr, params.disp.name, datestr(now,30))), ...
            'run_images','run_alpha_masks')
    end
    
    
end

fprintf('Done! ');  toc;
fprintf('\n');

return

