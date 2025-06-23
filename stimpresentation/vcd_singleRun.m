function [data, all_images] = vcd_singleRun(subj_nr, ses_nr, ses_type, run_nr, dispName, varargin)
% VCD stimulus presentation function to prepare running one VCD core 
% stimulus run:
% 
%    [data, all_images] = vcd_singleRun(subj_nr, ses_nr, ses_type, run_nr, dispName, varargin)
% 
% INPUTS (mandatory, for optional inputs, see input parser below):
%  subj_nr      : subject number (integral number between 1-999)
%  ses_nr       : session number (integral number between 1-27)
%  ses_type     : session type (either 1 for version A, 2 for version B)
%  run_nr       : run number (integral number between 1-15)
%  dispName     : display name to load display params. Choose from: 
%                 '7TAS_BOLDSCREEN32', 'KKOFFICE_AOCQ3277',
%                 'PPROOM_EIZOFLEXSCAN', 'EKHOME_ASUSVE247','CCNYU_VIEWPIXX3D'
%
% OUTPUTS:
%   data        : struct with behavioral button presses and monitor
%                  refresh rate timing, as well as other parameters.
%   all_images  : struct with all VCD images in uint8 pixels, to
%                  avoid additional loading time when executing multiple runs.
%
% Examples:
%  [data, all_images] = vcd_singleRun(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'env_type','BEHAVIOR', 'wantdatabypass',true)
%
% AUTHOR:
%  Written by Eline Kupers @ UMN (kupers@umn.edu) (2024,2025)
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
% MANDATORY INPUTS
p.addRequired('subj_nr' , @isnumeric); % subject number (integral number between 1-999)
p.addRequired('ses_nr'  , @isnumeric); % session number (integral number between 1-27)
p.addRequired('ses_type', @isnumeric); % session type (either 1 for version A, 2 for version B)
p.addRequired('run_nr'  , @isnumeric); % run number (integral number between 1-15)
p.addRequired('dispName', @(x) ismember(x,{'7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','PPROOM_EIZOFLEXSCAN','EKHOME_ASUSVE247','CCNYU_VIEWPIXX3D'})) % display name to get the right display params

% OPTIONAL INPUTS
p.addParameter('savedatafolder'     , ''        , @ischar);                      % place to store data with today's date
p.addParameter('behaviorfile'       , ''        , @ischar);                      % filename used for stored matlab file with behavioral data (name will add today's date)
p.addParameter('eyelinkfile'        , ''        , @ischar);                      % filename used for stored eyelink edf file with eyetracking data (name will add today's date)
p.addParameter('stim'               , []        , @isstruct);                    % if you don't want to reload stimuli, you need stim.im to containing uint8 images (time frames x 2 cell array). 
p.addParameter('loadparams'         , true      , @islogical)                    % (boolean) whether load stim/condition params or regenerate
p.addParameter('storeparams'        , true      , @islogical)                    % whether to store stimulus params
p.addParameter('laptopkey'          , -3        , @isnumeric);                   % listen to all keyboards/boxes (is this similar to k=-3;?)
p.addParameter('wanteyetracking'    , false     , @islogical);                   % whether to try to hook up to the eyetracker
p.addParameter('wantdatabypass'     , false     , @islogical);                   % whether to skip the experiment and just save dummy .mat file
p.addParameter('deviceNr'           , -3        , @isnumeric);                   % kbWait/kbCheck input device number to listen to. Default = -3, listen to all devices. Previously: vcd_checkDevices(params.deviceNr, params.device_check);
p.addParameter('device_check'       , 'both'    , @char);                        % what type of devices do we want to check for button presses: 'external','internal', or 'both'
p.addParameter('triggerkey'         , {'5','t'}, @(x) iscell(x) || isstring(x))  % key(s) that starts the experiment
p.addParameter('triggerkeyname'     , '''5'' or ''t''', @isstring)               % for display only
p.addParameter('userkeys'           , {'1','2','3','4'}, @(x) iscell(x) || isstring(x)) % key(s) that participants are expected to push
p.addParameter('offsetpix'          , [0 0]     , @isnumeric);                   % offset of screen in pixels [10 20] means move 10-px right, 20-px down
p.addParameter('movieflip'          , [0 0]     , @isnumeric)                    % whether to flip up-down, whether to flip left-right
p.addParameter('wantsynctest'       , true      , @islogical)                    % whether we want to run the PTB sync test or not
p.addParameter('savestim'           , false     , @islogical)                    % whether we want to store matlab file with stimuli and timing
p.addParameter('loadstimfromrunfile', false     , @islogical)                    % whether we want to load stim from run file
p.addParameter('ptbMaxVBLstd'       , 0.0009    , @isnumeric)                    % what standard deviation for screen flip duration do we allow?
p.addParameter('env_type'           , []        , @(x) ismember(x, {'MRI','BEHAVIOR'})); % are we running the 'BEHAVIOR' (PProom) or 'MRI' (7TAS) version of the VCD core experiment?
p.addParameter('timetable_file'     , ''        , @ischar);                      % what randomization file are we loading? file should exist in   
p.addParameter('all_images'         , struct()  , @isstruct);                    % preloaded all_images in a single struct (to save time)
p.addParameter('verbose'            , true      , @islogical)                    % (boolean) whether to print out text in command window. Default = true. 
p.addParameter('store_imgs'         , false     , @islogical)                    % (boolean) whether to save figures locally. Default = false.    
p.addParameter('infofolder'         , fullfile(vcd_rootPath,'workspaces','info')        , @ischar); % where are the *_info.csv file(s)?
p.addParameter('stimfolder'         , fullfile(vcd_rootPath,'workspaces','stimuli')     , @ischar); % where are the mat-files with store stimuli?
p.addParameter('instrfolder'        , fullfile(vcd_rootPath,'workspaces','instructions'), @ischar); % where are the txt and png files with task instructions?

% Parse inputs
p.parse(subj_nr, ses_nr, ses_type, run_nr, dispName, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    eval(sprintf('params.%s = p.Results.%s;', rename_me{ff},rename_me{ff}));
end
clear rename_me ff p

% deal with movieflip
if params.movieflip(1) && params.movieflip(2) %#ok<NODEF>
    flipfun = @(x) flipdim(flipdim(x,1),2); %#ok<*DFLIPDIM>
elseif params.movieflip(1)
    flipfun = @(x) flipdim(x,1);
elseif params.movieflip(2)
    flipfun = @(x) flipdim(x,2);
else
    flipfun = @(x) x;
end

% release images
all_images = params.all_images; params = rmfield(params, 'all_images');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERIPHERALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(params, 'disp') || isempty(params.disp)
    params.disp = vcd_getDisplayParams(params.dispName);
end

if params.wantsynctest
    skipsync = 0; % if wantsynctest = true, we do the psychtoolbox synctest
else
    skipsync = 1; % if wantsynctest = false, we skip the psychtoolbox synctest
end

% pton input arguments are:
% 1: [width, height, framerate, bitdepth]
% 2: winsize (fraction: default is full extent)
% 3: clutfile -- 0 for linear CLUT (-2 for squaring CLUT for BOLDSCREEN to simulate normal monitors --> NOTE: we do this manually!)
% 4: skipsync (bool: 0 is false -- do not skip the text, 1 is true -- skip the test)
% 5: wantstereo (bool: default is false)
if strcmp(params.disp.name, 'PPROOM_EIZOFLEXSCAN')
    % apparently PP room monitor native refresh rate shows up as 0 (but is
    % actually 60 Hz), due to the Bits# stimulus processor.
    % PProom EIZOFLEXScan screen ptonparams are expected to be {[1920 1200 0 24],[], 0, 0}
    ptonparams = {[params.disp.w_pix params.disp.h_pix 0 24],[],params.disp.clut, skipsync};
else % '7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','EKHOME_ASUSVE247', 'CCNYU_VIEWPIXX3D'
    % Nova1x32 coil with BOLDscreen and big eye mirrors ptonparams are expected to be {[1920 1080 120 24],[], 0, 0}
    ptonparams = {[params.disp.w_pix params.disp.h_pix params.disp.refresh_hz 24],[],params.disp.clut, skipsync};
end

% %%%%%%%%% SETUP RNG %%%%%%%%%
rand('seed', sum(100*clock)); %#ok<RAND>
randn('seed', sum(100*clock)); %#ok<RAND>
params.rng.rand  = rand;
params.rng.randn = randn;

% Very hacky but alais.. If we want to run the vcd-core experiment in a 
%  different environment without generating the stimuli for that particular 
%  display setup, we need to change disp name to load the stimulusfiles
%  that are default on VCD Google drive / Github.
%  For example, if we want to run the PPROOM_EIZOFLEXSCAN stimuli on 
%  CCNYU_VIEWPIXX3D..
if strcmp(params.env_type,'BEHAVIOR') && strcmp(params.disp.name,'CCNYU_VIEWPIXX3D')
    params.disp.name = 'PPROOM_EIZOFLEXSCAN';
end

% Infer env_type from display name if empty
if isempty(params.env_type) 
    if strcmp(params.disp.name,{'CCNYU_VIEWPIXX3D','PPROOM_EIZOFLEXSCAN','KKOFFICE_AOCQ3277','EKHOME_ASUSVE247'})
        params.env_type = 'BEHAVIOR';  
    elseif strcmp(params.disp.name,'7TAS_BOLDSCREEN32')
        params.env_type = 'MRI';
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% LOAD STIM & EXP PARAMS %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get stimulus params
if ~isfield(params, 'stim') || isempty(params.stim)
    if params.loadparams
        d = dir(fullfile(params.infofolder,sprintf('stim_%s*.mat',params.disp.name)));
        if  isempty(d)
            warning('[%s]: Can''t find stim params for %s, will reload params without overwriting params',mfilename,params.disp.name);
            params.stim = vcd_getStimParams('disp_name', params.disp.name, ...
                'load_params',false, ...
                'store_params', params.storeparams);
        else
            load(fullfile(d(end).folder,d(end).name),'stim');
            params.stim = stim; clear stim; %#ok<NODEF>
        end
    else
        params.stim  = vcd_getStimParams('disp_name', params.disp.name, ...
            'load_params',params.loadparams, ....
            'store_params', params.storeparams);
    end
end

% Get experimental session params
if ~isfield(params, 'exp') ||  isempty(params.exp)
    if params.loadparams
        d = dir(fullfile(params.infofolder,'exp*.mat'));
        if ~isempty(d)
            load(fullfile(d(end).folder,d(end).name),'exp');
            params.exp = exp; clear exp;
        else
            params.exp = vcd_getSessionParams('disp_name', params.disp.name, ...
                'load_params', false,...
                'store_params', params.storeparams);
        end
    else
        params.exp = vcd_getSessionParams('disp_name', params.disp.name,...
            'load_params', false, ...
            'store_params', params.storeparams);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% LOAD/CREATE TIME TABLE MASTER %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we look user wants to load stimuli from file
if params.loadstimfromrunfile
    % Images are in the format of:
    % subj001_ses01_A_run01_images_PPROOM_EIZOFLEXSCAN_20250505T184109.mat
    if isempty(params.savedatafolder)
        params.savedatafolder = fullfile(vcd_rootPath,'data',params.env_type, ...
            sprintf('vcd_subj%03d',params.subj_nr), sprintf('vcd_subj%03d_ses%02d',params.subj_nr, params.ses_nr));
    end
    d = dir(fullfile(params.savedatafolder,  ...
        sprintf('subj%03d_ses%02d_%s_run%02d_images_%s.mat', ...
        params.subj_nr, params.ses_nr, choose(params.ses_type==1,'A','B'), params.run_nr, params.disp.name)));
    if ~isempty(d)
        a1 = load(fullfile(d(end).folder,d(end).name));
        stim.im    = a1.run_images;
        stim.masks = a1.run_alpha_masks;
        stim.eye   = cat(4, a1.eye_im.sac_im, a1.eye_im.pupil_im_black, a1.eye_im.pupil_im_white);
        all_images.fix        = a1.fix_im.im;
        all_images.alpha.fix  = a1.fix_im.alpha;
        all_images.instr      = a1.instr_im.im;
        all_images.info.instr = a1.instr_im.info;
        run_frames = a1.run_frames;
        run_table  = a1.run_table;
    else
        error('[%s]: Can''t find file with run_images!!',mfilename)
    end
    clear d a1;
    
else % if not, then we look user pointed to a timetable_file

    if ~isempty(params.timetable_file)
        %% %%%%%%%%%%%%%% LOAD EXISTING TIME_TABLE_MASTER %%%%%%%%%%%%%%%%
        d = dir(fullfile(params.timetable_file));
        if ~isempty(d)
            load(fullfile(d(end).folder,d(end).name),'time_table_master','all_run_frames');
        else
            error('[%s]: Can''t find time table master!!',mfilename)
        end
        
    else % if not, then we create one on the fly 
        
        %% %%%%%%%%% CREATE SUBJECT TIME_TABLE_MASTER %%%%%%%%%%%%%%%
        % create folder to store time tables once we've created them
        timetablefiledir = strsplit(params.savedatafolder,'_ses');
        params.timetablefiledir = timetablefiledir{1};
        clear timetablefiledir;

        % see if there is a subject condition_master_shuffled, which is the precursor of time_table_master
        % and create the time_table_master from there..
        d = dir(fullfile(params.timetablefiledir, '%s_condition_master*.mat',sprintf('vcd_subj%03d',params.subj_nr)));
        if ~isempty(d)
            load(fullfile(d(end).folder,d(end).name),'condition_master_shuffled');

            [time_table_master,all_run_frames] = ...
                vcd_createRunTimeTables(params, ...
                'load_params',false, ...
                'store_params',true, ...
                'condition_master',condition_master_shuffled,...
                'env_type',params.env_type, ...
                'saveDir',params.timetablefiledir, ...
                'subj_id',sprintf('vcd_subj%03d',params.subj_nr));

        else
            %% %% CREATE CONDITION_MASTER AND TIME_TABLE_MASTER %%%%%%%%%
            % if there is no condition_master_shuffled for this subject, then
            % we create both condition_master_shuffled and time_table_master on the spot.
            d = dir(fullfile(vcd_rootPath,'workspaces','info', sprintf('condition_master*%s*.mat', params.disp.name)));
            a1 = load(fullfile(d(end).folder,d(end).name),'condition_master');

            % make randomization file!
            [~,~, time_table_master, all_run_frames] = vcd_createSessions(params,...
                'load_params',false, ...
                'store_params',true, ...
                'condition_master',a1.condition_master,...
                'env_type',params.env_type, ...
                'saveDir',params.timetablefiledir, ...
                'subj_id',sprintf('vcd_subj%03d',params.subj_nr));

            % clear up
            clear a1 d
        end
    end

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% TRUNCATE TIME_TABLE_MASTER %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Subselect the run frames and run table from bigger tables
    % (all_run_frames and time_table_master) with the entire experiment
    run_frames = all_run_frames(all_run_frames.session_nr == params.ses_nr & ...
        all_run_frames.session_type == params.ses_type & ...
        all_run_frames.run_nr ==params.run_nr,:);
    
    run_table = time_table_master(time_table_master.session_nr==params.ses_nr & ...
        time_table_master.session_type==params.ses_type &...
        time_table_master.run_nr==params.run_nr,:);
    
    clear time_table_master all_run_frames
 
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD RUN STIMULI    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If you give the function a struct called "stim" with "im" field, we will
    % not load the images..
    if ~exist('stim','var') || ~isfield(scan,'im') || isempty(stim.im) || isempty(stim.eye)
        
        % we want to load stimuli on the fly (takes 15-60 seconds depending on the run)
        if ~exist('all_images','var') && isempty(all_images)
            all_images = struct();
        end
        
        [images, masks, all_images] = vcd_getImageOrderSingleRun(params, ...
            run_table, run_frames, params.subj_nr, params.ses_nr, params.ses_type, params.run_nr, params.env_type, ...
            'all_images', all_images, 'savestim', params.savestim);
        
        
        % Insert images and mask into stim struct
        stim.im    = images;
        stim.masks = masks;
        stim.eye   = cat(4, all_images.eye.sac_im, all_images.eye.pupil_im_black, all_images.eye.pupil_im_white);
        
        clear images masks
    end
    
end



% Get fixation stim rect
fix_thin_rect0   = CenterRect([0 0 size(all_images.fix,1) size(all_images.fix,2)], [0 0 params.disp.w_pix params.disp.h_pix]);
fix_thick_rect0  = CenterRect([0 0 size(all_images.fix,1) size(all_images.fix,2)], [0 0 params.disp.w_pix params.disp.h_pix]);
fix_thin_rect    = {fix_thin_rect0 + repmat(params.offsetpix,1,2)}; clear fix_thin_rect0
fix_thick_rect   = {fix_thick_rect0 + repmat(params.offsetpix,1,2)}; clear fix_thick_rect0

% merge fixation dot image and mask
fix_thin_full         = cell(1,size(all_images.fix,4));
fix_thick_full_white  = cell(1,size(all_images.fix,4));
fix_thick_left        = cell(1,size(all_images.fix,4));
fix_thick_right       = cell(1,size(all_images.fix,4));
fix_thick_both        = cell(1,size(all_images.fix,4));

for ll = 1:size(all_images.fix,4) % loop over luminance values
    fix_thin_full{ll}        = feval(flipfun,  cat(3, all_images.fix(:,:,:,ll,1), all_images.alpha.fix(:,:,1)));
    fix_thick_full_white{ll} = feval(flipfun,  cat(3, all_images.fix(:,:,:,ll,2), all_images.alpha.fix(:,:,2)));
    fix_thick_left{ll}       = feval(flipfun,  cat(3, all_images.fix(:,:,:,ll,3), all_images.alpha.fix(:,:,3)));
    fix_thick_right{ll}      = feval(flipfun,  cat(3, all_images.fix(:,:,:,ll,4), all_images.alpha.fix(:,:,4)));
    fix_thick_both{ll}       = feval(flipfun,  cat(3, all_images.fix(:,:,:,ll,5), all_images.alpha.fix(:,:,5)));
end

% only mean gray luminance // no transparency
fix_thick_full_black{1} = feval(flipfun,  all_images.fix(:,:,:,ll,6));

% add fix images and rects to struct
fix_im = struct();
fix_im.fix_thin_full         = fix_thin_full;         clear fix_thin_full
fix_im.fix_thick_full_white  = fix_thick_full_white;  clear fix_thick_full
fix_im.fix_thick_left        = fix_thick_left;        clear fix_thick_left
fix_im.fix_thick_right       = fix_thick_right;       clear fix_thick_right
fix_im.fix_thick_both        = fix_thick_both;        clear fix_thick_both
fix_im.fix_thick_full_black  = fix_thick_full_black;  clear fix_thick_full_black

fix_im.fix_thin_rect         = fix_thin_rect;         clear fix_thin_rect
fix_im.fix_thick_rect        = fix_thick_rect;        clear fix_thick_rect




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% IMAGE XY CENTER, SIZE, OFFSET %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('scan','var') || ~isfield(scan, 'rects') || isempty(scan.rects)
    % recenter x,y-center coordinates if needed
    if any(params.offsetpix~=0) || isempty(params.offsetpix)
        xc = (ptonparams{1}(1)/2) + params.offsetpix(1);
        yc = (ptonparams{1}(2)/2) + params.offsetpix(2);
        % store offset
        params.stim.xc = xc;
        params.stim.yc = yc;
    else
        params.stim.xc = (ptonparams{1}(1)/2);
        params.stim.yc = (ptonparams{1}(2)/2);
    end
    
    % To calculate stim rects, we need stim [x,y]-center in pixels, taking
    % display size, stimulus aperture size, and fixation offset into account.
    centers = cell(size(stim.im,1),2);
    apsize  = cell(size(stim.im,1),2);
    
    % To avoid making multiple copies of the same static image, we create a
    % matrix called "im_IDs" which will refer to the static images stored 
    % in stim.im for every run frames. (We will do the same for centers and
    % aperture sizes..
    im_IDs = NaN(size(stim.im,1),2);

    % Check which stimulus sides (left/right) are empty or not for each run
    % time frame.
    stim_frames = (~cellfun(@isempty, stim.im));
    
    % Now loop over stimuli
    for nn = 1:size(stim.im,1)

        % We only do get centers/sizes for stimulus events or images in the
        % eye tracking block
        if ismember(run_frames.frame_event_nr(nn), ...
                [params.exp.block.stim_epoch1_ID, params.exp.block.stim_epoch2_ID, ...
                params.exp.block.eye_gaze_fix_ID,params.exp.block.eye_gaze_pupil_white_ID, ...
                params.exp.block.eye_gaze_pupil_black_ID, params.exp.block.eye_gaze_sac_target_ID])
            
            % check if we have left and right stimulus locations or only
            % one central stimulus location..
            numSides = find(stim_frames(nn,:));
            if ~isempty(numSides)
                for side = numSides % 1:left/2:right
                    
                    % If it is a gabor..
                    if ismember(run_frames.frame_im_nr(nn,side), [params.stim.gabor.unique_im_nrs_core,params.stim.gabor.unique_im_nrs_wm_test])
                        centers{nn,side} = [params.stim.gabor.x0_pix(side) + params.stim.xc, ... % x-coord (pixels)
                            params.stim.gabor.y0_pix(side) + params.stim.yc]; % y-coord (pixels)
                        
                    % If it is an RDK..
                    elseif ismember(run_frames.frame_im_nr(nn,side), [params.stim.rdk.unique_im_nrs_core,params.stim.rdk.unique_im_nrs_wm_test])
                        centers{nn,side} = [params.stim.rdk.x0_pix(side) + params.stim.xc, ...
                            params.stim.rdk.y0_pix(side) + params.stim.yc];
                        
                    % If it is an single dot..
                    elseif ismember(run_frames.frame_im_nr(nn,side), [params.stim.dot.unique_im_nrs_core,params.stim.dot.unique_im_nrs_wm_test])
                        % deal with dot pol2cart
                        if ismember(run_frames.frame_im_nr(nn,side), params.stim.dot.unique_im_nrs_wm_test)
                            test_im2D  = reshape(params.stim.dot.unique_im_nrs_wm_test,4,[]);
                            [~,test_im_idx] = ismember(run_frames.frame_im_nr(nn,side), test_im2D);
                            [x,y] = ind2sub([size(test_im2D,1),size(test_im2D,2)],test_im_idx);
                            dot_x = params.stim.dot.x0_pix_delta(x,y);
                            dot_y = params.stim.dot.y0_pix_delta(x,y);
                        elseif ismember(run_frames.frame_im_nr(nn,side), params.stim.dot.unique_im_nrs_core)
                            [~,core_im_idx] = ismember(run_frames.frame_im_nr(nn,side), params.stim.dot.unique_im_nrs_core);
                            dot_x = params.stim.dot.x0_pix(core_im_idx);
                            dot_y = params.stim.dot.y0_pix(core_im_idx);
                        end
                        centers{nn,side} = [dot_x,dot_y];
                        
                    % If it is an object..
                    elseif ismember(run_frames.frame_im_nr(nn,side), [params.stim.obj.unique_im_nrs_core,params.stim.obj.unique_im_nrs_wm_test])
                        centers{nn,side} = [params.stim.obj.x0_pix(side) + params.stim.xc, ...
                            params.stim.obj.y0_pix(side) + params.stim.yc];
                        
                    % If it is a natural scene..
                    elseif ismember(run_frames.frame_im_nr(nn,side), ...
                            [params.stim.ns.unique_im_nrs_core,params.stim.ns.unique_im_nrs_wm_test,params.stim.ns.unique_im_nrs_wm_test,params.stim.ns.unique_im_nrs_ltm_lures])
                        centers{nn,side} = [params.stim.ns.x0_pix + params.stim.xc, params.stim.ns.y0_pix + params.stim.yc];
                        
                    % If it is a catch trial
                    elseif isnan(run_frames.frame_im_nr(nn,side)) || (run_frames.frame_im_nr(nn,side)==0)
                        centers{nn,side} = [NaN, NaN];
                    end
                    
                    % ADD STIMULUS SIZE
                    if isnan(run_frames.frame_im_nr(nn,side))
                        % If it is a catch trial
                        apsize{nn,side}  = [NaN, NaN];
                        
                    else % ADD IMAGE ID, copy size/im_ID/center if we deal with static stim.

                        % If we have RDKs, then we have a new image every time frame so we can just use the counter
                        if ismember(run_frames.frame_im_nr(nn,side), [params.stim.rdk.unique_im_nrs_core,params.stim.rdk.unique_im_nrs_wm_test])
                            
                            % Add image ID
                            im_IDs(nn,side) = nn;
                            
                            % apsize 1: image width -- second dim (pixels), apsize 2: image height -- first dim (pixels)
                            apsize{nn,side} = [size(stim.im{nn,side},2), size(stim.im{nn,side},1)];
                            
                        % if we deal with CD tasks..
                        elseif ismember(run_frames.crossingIDs(nn),find(~cellfun(@isempty, regexp(params.exp.crossingnames,'cd-*'))))
                            if ~isempty(stim.im(nn,side))
                                
                                % Add image ID
                                im_IDs(nn,side) = nn;
                            
                                % Add aperture size
                                apsize{nn,side} = [size(stim.im{nn,side},2), size(stim.im{nn,side},1)];
                            
                               if isempty(stim.im{nn+1,side})
                                   im_IDs(nn:(nn+params.stim.stimdur_frames-1),side)  = nn;
                                   centers(nn:(nn+params.stim.stimdur_frames-1),side) = mat2cell(centers{im_IDs(nn,side),side},1,2);
                                   apsize(nn:(nn+params.stim.stimdur_frames-1),side)  = mat2cell([size(stim.im{im_IDs(nn,side),side},2), size(stim.im{im_IDs(nn,side),side},1)],1,2);
                               end
                            else
                                im_IDs(nn:(nn+params.stim.stimdur_frames-1),side)  = nn;
                                centers(nn:(nn+params.stim.stimdur_frames-1),side) = mat2cell(centers{im_IDs(nn,side),side},1,2);
                                apsize(nn:(nn+params.stim.stimdur_frames-1),side)  = mat2cell([size(stim.im{im_IDs(nn,side),side},2), size(stim.im{im_IDs(nn,side),side},1)],1,2);
                            end
                            
                        % If we only have a static image, repeat im_ID,
                        % centers and aperture size for upcoming frames
                        elseif isnan(im_IDs(nn-1,side))
                            im_IDs(nn:(nn+params.stim.stimdur_frames-1),side)  = nn;
                            centers(nn:(nn+params.stim.stimdur_frames-1),side) = mat2cell(centers{im_IDs(nn,side),side},1,2);
                            apsize(nn:(nn+params.stim.stimdur_frames-1),side)  = mat2cell([size(stim.im{im_IDs(nn,side),side},2), size(stim.im{im_IDs(nn,side),side},1)],1,2);
                        end
                    end
                    
                    % COMBINE IMAGES AND MASKS, FLIP IMAGES IF REQUESTED
                    if isempty(stim.masks{nn,side})
                        stim.im{nn,side} = feval(flipfun, stim.im{nn,side});
                    else
                        stim.im{nn,side} = feval(flipfun, cat(3, stim.im{nn,side}, stim.masks{nn,side}));
                        stim.masks{nn,side} = [];
                    end
                end % numSides
            end % isempty numSides
        end % is stim event
    end % nn time frame loop
end % if scan struct exists


% Add info to run frames and stim struct
run_frames.im_IDs = im_IDs;  clear im_ID
stim.centers      = centers; clear centers;
stim.apsize       = apsize;  clear apsize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% CENTER RECT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stim.rects = cell(size(stim.centers));
nSides = unique(run_frames.is_cued(~isnan(run_frames.is_cued)));
if length(nSides)==1 && nSides==3
    nSides = 1; % if we happen to only have NS blocks..
elseif length(nSides)>1 && any(nSides==3)
    nSides = [1,2];
else
    nSides = nSides';
end
for side = nSides
    % Find the non-empty center and size cells for each stimulus side
    nonemptycenters   = ~cellfun(@isempty, stim.centers(:,side));
    nonemptysizes     = ~cellfun(@isempty, stim.apsize(:,side));
    % We expect centers and sizes to have the same nr of empty cells..
    assert(isequal(nonemptycenters,nonemptysizes)); clear nonemptysizes
    % Now only select the nonempty centers and sizes
    centers_shortlist = stim.centers(nonemptycenters,side);
    apsize_shortlist  = stim.apsize(nonemptycenters,side);
    centers_mat       = cell2mat(centers_shortlist);
    size_mat          = cell2mat(apsize_shortlist);
    destination_mat = [centers_mat(:,1) - size_mat(:,1)./2, ... x1 top-left-x
                       centers_mat(:,2) - size_mat(:,2)./2, ... y1 top-left-y
                       centers_mat(:,1) + size_mat(:,1)./2, ... x2 bottom-right-x
                       centers_mat(:,2) + size_mat(:,2)./2];  % y2 bottom-right-y
                            
    % Create stimulus "rects" for PTB. 
    % CenterRect cannot handle empty cells or NaNs
    % stimsize        = [width, height] in pixels
    % destinationRect = [top-left-x, top-left-y, bottom-right-x, bottom-right-y] pixel location of the stimulus
    rects_mat       = CenterRect([zeros(size(size_mat,1),1) zeros(size(size_mat,1),1), size_mat(:,1) size_mat(:,2)], ...
                        [destination_mat(:,1), destination_mat(:,2), ...
                         destination_mat(:,3), destination_mat(:,4)]);
    rects_shortlist = mat2cell(rects_mat,ones(size(rects_mat,1),1));
    
    % Insert the "rects" into the struct
    stim.rects(nonemptycenters,side) = rects_shortlist;
end
% Clear some memory
clear centers_shortlist apsize_shortlist rects_shortlist ...
        centers_mat size_mat destination_mat rects_mat nSides side



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK INSTRUCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
introscript.im   = all_images.instr(:,:,:,1);
introscript.rect = CenterRect([0 0 size(introscript.im,1) size(introscript.im,2)], ...
                              [0 0 params.disp.w_pix,params.disp.h_pix]);
 
% Find the crossing_nrs for each stimulus block
taskIDs   = unique(run_table.crossing_nr);
taskIDs   = taskIDs(~isnan(taskIDs));
taskIDs   = taskIDs(taskIDs~=0); % black periods
taskIDs   = taskIDs(taskIDs~=999); % eyetracking block events
taskNames = params.exp.crossingnames(taskIDs);

% Now load the task instructions from file.
taskscript = struct();
taskscript.im   = cell(1,length(taskIDs));
taskscript.rect = cell(1,length(taskIDs));
for nn = 1:length(taskIDs)
    
    % load in instruction image
    taskscript.im{nn}   = all_images.instr(:,:,:,taskIDs(nn)+1);  % +1 because introtext is first image
    taskscript.rect{nn} = CenterRect([0 0 size(taskscript.im{nn},1) size(taskscript.im{nn},2)], [0 0 params.disp.w_pix,params.disp.h_pix]);
    
    % Check we got the right file
    taskName = strrep(taskNames{nn},'-','_');
    assert(isequal(all_images.info.instr{taskIDs(nn)+1}, sprintf('%02d_runvcdcore_%s.png', taskIDs(nn),taskName))); % +1 because introtext is first image
end

% Clear some memory
clear taskIDs taskNames d taskName

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  DEFINE FILENAMES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get time of stim call
timeofshowstimcall = datestr(now,30);

% Create subj###_ses## folder where behavioral, VBL timing, and eyetracking
% data will be stored
if isempty(params.savedatafolder)
    params.savedatafolder = fullfile(vcd_rootPath,'data',params.env_type,...
        sprintf('vcd_subj%03d',params.subj_nr), ...
        sprintf('vcd_subj%03d_ses%02d',params.subj_nr, params.ses_nr));
end

% Create subj###_ses## folder if it doesn't exist
if ~exist(params.savedatafolder,'dir'), mkdir(params.savedatafolder); end

% Create behavioral matlab file name if user didn't define it yet (we do
% this at the end because we want to use the "timeofshowstimcall")
if ~isfield(params,'behaviorfile') || isempty(params.behaviorfile)
    params.behaviorfile = sprintf('behavior_%s_vcd_subj%03d_ses%02d_%s_run%02d.mat',...
        timeofshowstimcall, params.subj_nr, params.ses_nr,choose(params.ses_type==1,'A','B'),params.run_nr);
end

% Create eyelink file if user didn't define one yet.
if params.wanteyetracking
    if ~isfield(params,'eyelinkfile') || isempty(params.eyelinkfile)
        params.eyelinkfile = sprintf('eye_%s_vcd_subj%03d_ses%02d_%s_run%02d.edf', ...
            timeofshowstimcall,params.subj_nr,params.ses_nr,choose(params.ses_type==1,'A','B'),params.run_nr);
    end
else
    params.eyelinkfile = [];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START EXPERIMENT! %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tell the user how many frames we expect for this run
fprintf(' *** Expected duration of this run is ***\n')
fprintf(' %d time frames (%d Hz) \n',size(run_frames,1), params.stim.presentationrate_hz);
fprintf(' %3.2f seconds (or %d minutes and %02d seconds) \n',size(run_frames,1)/60, floor(size(run_frames,1)/3600),rem(size(run_frames,1),3600)/60);


if ~params.wantdatabypass
    
  % If we want to run the actual experiment
  [data, ~] = vcd_showStimulus(...
      params, ptonparams, ...
      fix_im, ...
      stim, ...
      run_frames, ...
      run_table, ...
      introscript, ...
      taskscript, ...
      params.deviceNr);

else % or just generate the files + dummy data

  % calc
  numtimeframes = sum(run_table.event_dur);                       % number of timeframes expected
  idealexpdur   = numtimeframes/params.stim.presentationrate_hz;  % ideal exp dur in seconds
  mfi           = 1/params.stim.presentationrate_hz;              % idealized flip interval in seconds
  
  % create dummy data
  data = struct();
  data.wantframefiles    = 0;
  data.detectinput       = 1;
  data.forceglitch       = 0;
  data.timeKeys          = {739780.00   {'absolutetimefor0'}
                            -0.016      {'trigger'}
                            -0.014      {'t'}
                            idealexpdur {'DONE'}};
  data.timing.mfi        = mfi;
  data.timing.glitchcn   = 0;
  data.timing.timeframes = linspacefixeddiff(0,mfi,numtimeframes);
  data.timing.starttime  = 12162.00;
  save(fullfile(params.savedatafolder,params.behaviorfile),'params','run_table','run_frames','data');

end

return


