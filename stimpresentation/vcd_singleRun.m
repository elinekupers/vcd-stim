function [] = vcd_singleRun(subj_nr, ses_nr, run_nr, varargin)


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p = inputParser;
p.addRequired('subj_nr'         , @isnumeric); % subject number
p.addRequired('ses_nr'          , @isnumeric); % session number
p.addRequired('ses_type'        , @isnumeric); % session type
p.addRequired('run_nr'          , @isnumeric); % nun number
p.addParameter('scan'           , struct()  , @isstruct);                    % struct with exp im rects, centers,  and sizes for ptb
p.addParameter('savedatadir'    , []        , @ischar);                      % place to store data with today's date
p.addParameter('behaviorfile'   , []        , @ischar);                      % filename to store behavioral data with today's date
p.addParameter('eyelinkfile'    , []        , @ischar);                      % where the eyelink edf file can be obtained
p.addParameter('loadparams'     , true      , @islogical)                    % whether load stim/condition params or regenerate
p.addParameter('storeparams'    , true      , @islogical)                    % whether to store stimulus params
p.addParameter('infofolder'     , fullfile(vcd_rootPath,'workspaces','info'), @ischar);     % where the *_info.csv file is
p.addParameter('stimfolder'     , fullfile(vcd_rootPath,'workspaces','stimuli'), @ischar);  % where the images can be obtained
p.addParameter('instrtextdir'   , fullfile(vcd_rootPath,'workspaces','instructions'), @ischar); % where the task instructions can be obtained
p.addParameter('laptopkey'      , -3        , @isnumeric);                   % listen to all keyboards/boxes (is this similar to k=-3;?)
p.addParameter('wanteyetracking', false     , @islogical);                   % whether to try to hook up to the eyetracker
p.addParameter('deviceNr'       , []        , @isnumeric);                   % kbWait/Check input device number
p.addParameter('device_check'   , 'both'    , @char);                        % what type of devices do we want to check for button presses: 'external','internal', or 'both'
p.addParameter('triggerkey'     , {'5%','t'}, @(x) iscell(x) || isstring(x)) % key that starts the experiment
p.addParameter('triggerkeyname' , '''5'' or ''t''', @isstring)               % for display only
p.addParameter('offsetpix'      , [0 0]     , @isnumeric);                   % offset of screen in pixels [10 20] means move 10-px right, 20-px down
p.addParameter('movieflip'      , [0 0]     , @isnumeric)                    % whether to flip up-down, whether to flip left-right
p.addParameter('debugmode'      , false     , @islogical)                    % whether to use debug mode (no BOLDscreen, no eyelink)
p.addParameter('dispName'       , '7TAS_BOLDSCREEN32' , @ischar)             % display params: 7TAS_BOLDSCREEN32, KKOFFICE_AOCQ3277, PPROOM_EIZOFLEXSCAN, 'EKHOME_ASUSVE247'
p.addParameter('savestim'       , false     , @islogical)                    % whether we want to store temp file with stimuli and timing
p.addParameter('loadstimfromrunfile', false     , @islogical)                % whether we want to load stim from run file

% Parse inputs
p.parse(subj_nr, ses_nr, run_nr, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    %     eval([sprintf('%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
    eval([sprintf('params.%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p

% release scan and timing var
scan = params.scan; rmfield(params,'scan');

% deal with movieflip
if params.movieflip(1) && params.movieflip(2)
    flipfun = @(x) flipdim(flipdim(x,1),2); %#ok<*DFLIPDIM>
elseif params.movieflip(1)
    flipfun = @(x) flipdim(x,1);
elseif params.movieflip(2)
    flipfun = @(x) flipdim(x,2);
else
    flipfun = @(x) x;
end


%% %%%%%%%%%%%%% PERIPHERALS %%%%%%%%%%%%%

if ~isfield(params, 'disp') || isempty(params.disp)
    params.disp = vcd_getDisplayParams(params.dispName); % BOLDSCREEN is default
end


% EK HACK START ---
if params.debugmode
    params.disp.name = 'PPROOM_EIZOFLEXSCAN';
end
% EK HACK END ---

if params.debugmode % skip synctest
    skipsync = 1;
else
    skipsync = 0;
end

% Get device nr for KbCheck
deviceNr = vcd_checkDevices(params.deviceNr, params.device_check);


% Nova1x32 coil with BOLDscreen and big eye mirrors
% expected to be {[1920 1080 120 24],[], 0, 0}
% 1: [width, height, framerate, bitdepth]
% 2: winsize (fraction: default is full extent)
% 3: clutfile -- 0 for linear CLUT (-2 for squaring CLUT for BOLDSCREEN to simulate normal monitors --> NB: we do this manually!)
% 4: skipsync (bool: 0 is false, 1 is true)
% 5: wantstereo (bool: default is false)
if strcmp(params.disp.name, 'PPROOM_EIZOFLEXSCAN')
    % apparently PP room monitor native refresh rate show up as 0 (but is 60 Hz)
    ptonparams = {[params.disp.w_pix params.disp.h_pix 0 24],[],params.disp.clut, skipsync};
else
    ptonparams = {[params.disp.w_pix params.disp.h_pix params.disp.refresh_hz 24],[],params.disp.clut, skipsync};
end

% %%%%%%%%% SETUP RNG %%%%%%%%%
rand('seed', sum(100*clock));
randn('seed', sum(100*clock));
params.rng.rand = rand;
params.rng.randn = randn;

%% %%%%%%%%%%%%% STIM PARAMS %%%%%%%%%%%%%

% Infer session type
if strcmp(params.disp.name, '7TAS_BOLDSCREEN32')
    session_env = 'MRI';
    
    % Buttonbox / keyboard
    params.ignorekeys = KbName({params.triggerkey});  % dont record TR triggers as subject response
    
elseif strcmp(params.disp.name,'PPROOM_EIZOFLEXSCAN')
    session_env = 'BEHAVIOR';
else
    session_env = 'MRI';
end


% Stimulus params
if ~isfield(params, 'stim') || isempty(params.stim)
    if params.loadparams
        d = dir(fullfile(params.infofolder,sprintf('stim_%s*.mat',params.disp.name)));
        if  isempty(d)
            warning('[%s]: Can''t find stim params for %s, will reload params without overwriting params',mfilename,params.disp.name);
            params.stim   = vcd_getStimParams('disp_name', params.disp.name, ...
                'load_params',false, ...
                'store_params', params.storeparams);
            
        else
            load(fullfile(d(end).folder,d(end).name),'stim');
            params.stim = stim; clear stim;
        end
    else
        params.stim   = vcd_getStimParams('disp_name', params.disp.name, ...
            'load_params',params.loadparams, ....
            'store_params', params.storeparams);
    end
end

% Session params
if ~isfield(params, 'exp') ||  isempty(params.exp)
    if params.loadparams
        d = dir(fullfile(params.infofolder,'exp*.mat'));
        if ~isempty(d)
            load(fullfile(d(end).folder,d(end).name),'exp');
            params.exp = exp; clear exp;
        else
            params.exp = vcd_getSessionParams('disp_name', params.disp.name,'load_params',false, 'store_params', params.storeparams);
        end
    else
        params.exp = vcd_getSessionParams('disp_name', params.disp.name,'load_params',false, 'store_params', params.storeparams);
    end
end

% TIME TABLE MASTER
if ~exist('time_table_master','var') ||  isempty(params.exp)
    if params.loadparams
        d = dir(fullfile(params.infofolder,'time_table_master_complete*.mat'));
        if ~isempty(d)
            load(fullfile(d(end).folder,d(end).name),'time_table_master','all_run_frames');
        end
    end
    
else
    error('[%s]: Can''t find time table master!!',mfilename)
end

%% %%%%%%%%% LOAD IMAGES & ORDER

if ~exist('scan','var') || ~isfield(scan,'exp_im') || isempty(scan.exp_im)
    
    if params.loadstimfromrunfile
        
        % Images are in the format of:
        % subj001_ses01_A_run01_images_PPROOM_EIZOFLEXSCAN_20250505T184109.mat
        d = dir(fullfile(params.infofolder,sprintf('subj%03d',params.subj_nr),  ...
            sprintf('subj%03d_ses%02d_%s_run%02d_images_%s_*.mat', ...
            params.subj_nr, params.ses_nr, choose(params.ses_type==1,'A','B'), params.run_nr, params.disp.name)));
        a = load(fullfile(d(end).folder,d(end).name));
        images = a.images;
        masks = a.masks;
        clear a d;
        
    else
        
        [images, masks] = vcd_getImageOrderSingleRun(params, ...
            time_table_master, all_run_frames, params.subj_nr, params.ses_nr, params.ses_type, params.run_nr, ...
            'store_params', false,'session_env', session_env);
    end
    
    run_frames = all_run_frames(all_run_frames.session_nr == params.ses_nr & ...
                                all_run_frames.session_type == params.ses_type & ...
                            all_run_frames.run_nr ==params.run_nr,:);
    run_table = time_table_master(time_table_master.session_nr==params.ses_nr ... 
        & time_table_master.session_type==params.ses_type ... 
        & time_table_master.run_nr==params.run_nr,:);
    
    stim.im    = images;
    stim.masks = masks;
    clear images masks
    
end

%% %%%%%%%%% IMAGE XY CENTER OFFSET

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
    
    % Get [x,y]-center in pixels of peripheral stimuli given display size, and
    % offset. Get stimulus aperture size..
    centers = cell(size(stim.im,1),2);
    apsize  = cell(size(stim.im,1),2);
    
    
    for nn = 1:size(stim.im,1)
        
        if ismember(run_frames.frame_event_nr(nn), [params.exp.block.stim_epoch1_ID, params.exp.block.stim_epoch2_ID, ...
                params.exp.block.eye_gaze_fix_ID,params.exp.block.eye_gaze_pupil_white_ID, params.exp.block.eye_gaze_pupil_black_ID, params.exp.block.eye_gaze_sac_target_ID])
            
            numSides = find(~cellfun(@isempty, stim.im(nn,:)));
            
            if ~isempty(numSides)
                
                for side = numSides
                    
                    if ismember(run_frames.frame_im_nr(nn,side), [params.stim.gabor.unique_im_nrs_core,params.stim.gabor.unique_im_nrs_wm_test])
                        centers{nn,side} = [params.stim.gabor.x0_pix(side) + params.stim.xc, ... % x-coord (pixels)
                            params.stim.gabor.y0_pix(side) + params.stim.yc]; % y-coord (pixels)
                        
                    elseif ismember(run_frames.frame_im_nr(nn,side), [params.stim.rdk.unique_im_nrs_core,params.stim.rdk.unique_im_nrs_wm_test])
                        centers{nn,side} = [params.stim.rdk.x0_pix(side) + params.stim.xc, ...
                            params.stim.rdk.y0_pix(side) + params.stim.yc];
                        
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
                        
                        centers{nn,side} = [dot_x + params.stim.xc, dot_y + params.stim.yc];
                        
                    elseif ismember(run_frames.frame_im_nr(nn,side), [params.stim.obj.unique_im_nrs_core,params.stim.obj.unique_im_nrs_wm_test])
                        centers{nn,side} = [params.stim.obj.x0_pix(side) + params.stim.xc, ...
                            params.stim.obj.y0_pix(side) + params.stim.yc];
                        
                    elseif ismember(run_frames.frame_im_nr(nn,side), ...
                            [params.stim.ns.unique_im_nrs_core,params.stim.ns.unique_im_nrs_wm_test,params.stim.ns.unique_im_nrs_wm_test,params.stim.ns.unique_im_nrs_ltm_lures])
                        centers{nn,side} = [params.stim.ns.x0_pix + params.stim.xc, params.stim.ns.y0_pix + params.stim.yc];
                    
                    elseif isnan(run_frames.frame_im_nr(nn,side)) || (run_frames.frame_im_nr(nn,side)==0)
                        centers{nn,side} = [NaN, NaN];
                    end
                    
                    if isnan(run_frames.frame_im_nr(nn,side))
                        apsize{nn,side}  = [NaN, NaN];
                    else
                        % apsize 1: image width (pixels), apsize 2: image height (pixels)
                        apsize{nn,side}  = [size(stim.im{nn,side},2), size(stim.im{nn,side},1)];
                    end
                    
                    % COMBINE IMAGE AND MASKS, FLIP IM IF REQUESTED
                    if isempty(stim.masks{nn,side})
                        stim.im{nn,side} = feval(flipfun, stim.im{nn,side});
                    else
                        stim.im{nn,side} = feval(flipfun, cat(3, stim.im{nn,side}, stim.masks{nn,side}));
                        stim.masks{nn,side} = [];
                    end
                    
                end
            end
        end
    end
end


stim.centers = centers; clear centers;
stim.apsize  = apsize; clear apsize;

%%%%%%% CENTER RECT
nonemptycenters = ~cellfun(@isempty, stim.centers(:,1));
nonemptysizes = ~cellfun(@isempty, stim.apsize(:,1));
assert(isequal(nonemptycenters,nonemptysizes)); clear nonemptysizes
centers_shortlist = stim.centers(nonemptycenters,:);
apsize_shortlist  = stim.apsize(nonemptycenters,:);

% insert NaNs for second column when using single square stimulus (nat scene)
if size(centers_shortlist,2)==2
    centers_shortlist(cellfun(@isempty,centers_shortlist(:,2)),2) = ...
        repmat({[NaN,NaN]},size(find(cellfun(@isempty,centers_shortlist(:,2))),1),1);
end
if size(apsize_shortlist,2)==2
    apsize_shortlist(cellfun(@isempty,apsize_shortlist(:,2)),2) = ...
        repmat({[NaN,NaN]},size(find(cellfun(@isempty,apsize_shortlist(:,2))),1),1);
end
rects_shortlist = cell(size(apsize_shortlist));

for side = [1,2]
    rects_shortlist(:,side) = ...
        cellfun(@(stimsize,stimcenter) CenterRectOnPoint([0 0 stimsize(1) stimsize(2)], stimcenter(1), stimcenter(2)), ...
        apsize_shortlist(:,side),centers_shortlist(:,side), 'UniformOutput', false);
end

stim.rects = cell(size(stim.centers));
stim.rects(nonemptycenters,:) = rects_shortlist;

%% Accumulate stimulus idx
im_IDs     = NaN(length(run_frames.frame_im_nr),1);
empty_rows = cellfun(@isempty, stim.im(:,1));
non_empty_rows = ~empty_rows;
empty_rows = find(empty_rows);
non_empty_rows = find(non_empty_rows);

im_counter = 1;
for mm = 1:length(non_empty_rows)
    tmp = run_frames.frame_im_nr(non_empty_rows(mm),:)~=[0,0];
    if any(tmp)
        im_IDs(non_empty_rows(mm):(non_empty_rows(mm)+params.stim.gabor.duration-1)) = non_empty_rows(mm);
    end
    im_counter = im_counter+1;
end

run_frames.im_IDs = im_IDs; clear im_IDs im_cnt empty_rows

%% %%%%%%%%%%%%% BACKGROUND IM %%%%%%%%%%%%%
%  BACKGROUND: 4D array: [x,y, 3, num images]
% input 2: 'puzzle'  (square central image + peripheral apertures,
%          'dotring' (simple dot iso eccen donut),
%          'comb'    (puzzle + simple dot iso-eccen ring overlayed)
% input 3: 'skinny'  (no space between stim and edge of background) or
%          'fat'     (+2 deg from stim)
% input 4: number of unique noise background images
% input 5: offset of center in pixels [x,y]
% create background image that adjusts for shifted background if needed

if ~exist('bckground','var') || isempty(bckground)
    d = dir(fullfile(params.stimfolder,params.disp.name, sprintf('bckgrnd_%s*.mat',params.disp.name)));
    a = load(fullfile(d(end).folder, d(end).name),'bckgrnd_im');
    bckground = a.bckgrnd_im(:,:,:,params.run_nr);
    clear a d
end

if any(params.offsetpix~=[0,0])
    bckground = vcd_pinknoisebackground(params, 'comb', 'fat', 1, params.offsetpix);
end

bckground = feval(flipfun,bckground);

%% %%%%%%%%%%%%% ET IM %%%%%%%%%%%%%
%  ET TARGETS: 4D array: [x,y, 3, type]
% Load stored eyetracking block target images if needed
if ~exist('eye_im','var') || isempty(eye_im)
    fprintf('[%s]: Loading eyetracking target images..\n',mfilename);
    
    % FIX: 5D array: [x,y, 3, 5 lum, 2 widths]
    d = dir(sprintf('%s*.mat', params.stim.el.stimfile));
    a = load(fullfile(d(end).folder,d(end).name), 'sac_im','pupil_im_white','pupil_im_black');
    eye_im.sac_im          = a.sac_im; 
    eye_im.pupil_im_white  = a.pupil_im_white;
    eye_im.pupil_im_black  = a.pupil_im_black;
    clear a d;
end


%% %%%%%%%%%%%%% FIX IM %%%%%%%%%%%%%
%  FIX CIRCLE: 5D array: [x,y, 3, lum, type]

if ~exist('fix','var') || isempty(fix) || ~isfield(fix,'im')
    d = dir(fullfile(params.stimfolder,params.disp.name, sprintf('fix_%s*.mat',params.disp.name)));
    a = load(fullfile(d(end).folder, d(end).name),'fix_im','mask');
    fix_im = a.fix_im;
    fix_mask = a.mask;
    clear a d
end


% make fixation dot texture
fix_thin_full   = cell(1,size(fix_im,4));
fix_thick_full  = cell(1,size(fix_im,4));
fix_thick_left  = cell(1,size(fix_im,4));
fix_thick_right = cell(1,size(fix_im,4));
fix_thick_both  = cell(1,size(fix_im,4));

fix_thin_rect0   = CenterRect([0 0 round(size(fix_im,1)) round(size(fix_im,2))], [0 0 params.disp.w_pix params.disp.h_pix]);
fix_thick_rect0  = CenterRect([0 0 round(size(fix_im,1)) round(size(fix_im,2))], [0 0 params.disp.w_pix params.disp.h_pix]);
fix_thin_rect   = {fix_thin_rect0 + repmat(params.offsetpix,1,2)};
fix_thick_rect  = {fix_thick_rect0 + repmat(params.offsetpix,1,2)};

for ll = 1:size(fix_im,4) % loop over luminance values
    fix_thin_full{ll}   = feval(flipfun,  cat(3, fix_im(:,:,:,ll,1), fix_mask(:,:,1)));
    fix_thick_full{ll}  = feval(flipfun,  cat(3, fix_im(:,:,:,ll,2), fix_mask(:,:,2)));
    fix_thick_left{ll}  = feval(flipfun,  cat(3, fix_im(:,:,:,ll,3), fix_mask(:,:,3)));
    fix_thick_right{ll} = feval(flipfun,  cat(3, fix_im(:,:,:,ll,4), fix_mask(:,:,4)));
    fix_thick_both{ll}  = feval(flipfun,  cat(3, fix_im(:,:,:,ll,5), fix_mask(:,:,5)));
end

clear fix_im fix_mask;
fix_im = struct();
fix_im.fix_thin_full   = fix_thin_full; clear fix_thin_full
fix_im.fix_thick_full  = fix_thick_full; clear fix_thick_full
fix_im.fix_thick_left  = fix_thick_left; clear fix_thick_left
fix_im.fix_thick_right = fix_thick_right; clear fix_thick_right
fix_im.fix_thick_both  = fix_thick_both; clear fix_thick_both
fix_im.fix_thin_rect   = fix_thin_rect; clear fix_thin_rect
fix_im.fix_thick_rect  = fix_thick_rect; clear fix_thick_rect


%% %%%%%%%%%%%%% TASK INSTRUCTIONS %%%%%%%%%%%%%

introscript = fullfile(params.instrtextdir,'00_runvcdcore_subjectinstructions.txt');
if ~exist(introscript,'file')
    error('[%s]: Can''t find instructions text file!',mfilename')
end

taskIDs = unique(run_table.crossing_nr);
taskIDs = taskIDs(~isnan(taskIDs));
taskIDs = taskIDs(taskIDs~=0);
taskIDs = taskIDs(taskIDs~=999);
taskNames = params.exp.crossingnames(taskIDs);

for nn = 1:length(taskIDs)
    taskName = strrep(taskNames{nn},'-','_');
    d = dir(fullfile(params.instrtextdir,sprintf('%02d_runvcdcore_%s.txt', taskIDs(nn),taskName)));
    taskscript{nn} = fullfile(d.folder,d.name);
end


%% %%%%%%%%%%%%% INIT SCREEN %%%%%%%%%%%%%
Screen('Preference', 'SyncTestSettings', 0.0004); %.0004
oldCLUT = pton(ptonparams{:});

win  = firstel(Screen('Windows'));
oldPriority = Priority(MaxPriority(win));

rect = Screen('Rect',win); % what is the total rect
% rect = CenterRect(round([0 0 rect(3)*winsize rect(4)*winsize]),rect);



%% %%%%%%%%%%%%% Eyelink stuff %%%%%%%%%%%%%
% ANON EYE FUN for SYNC TIME
if params.wanteyetracking
    tfunEYE     = @() Eyelink('Message','SYNCTIME');
    
    if ~isfield(params,'eyelinkfile') || isempty(params.eyelinkfile)
        params.eyelinkfile = fullfile(params.savedatadir,sprintf('eye_%s_vcd_subj-%s_run-%d.edf',datestr(now,30),params.subj_nr,params.run_nr));
    end
else
    tfunEYE = @() fprintf('EXP STARTS.\n');
end



% INIT/CALIBRATION/VALIDATION
if params.wanteyetracking && ~isempty(params.eyelinkfile)
    
    % initialize
    assert(EyelinkInit()==1);
    
    %     % Check TCP/IP address
    %     if strcmp(dispName,'7TAS_BOLDSCREEN32')
    %         Eyelink('SetAddress','100.1.1.1') %% copied from MGS exp
    %     elseif strcmp(dispName,'PPROOM_EIZOFLEXSCAN')
    %         Eyelink('SetAddress','100.1.1.1') %% <--- CHECK THIS
    %     end
    
    % Get eyelink default params
    el = EyelinkInitDefaults(win);
    if ~isempty(el.callback)
        PsychEyelinkDispatchCallback(el);
    end
    
    % Update default EYELINK params with VCD needs
    el = vcd_setEyelinkParams(el);
    EyelinkUpdateDefaults(el);
    
    % Get window size
    wwidth = rect(3); wheight = rect(4);  % returns in pixels
    
    % Ensure window size is what we think it is
    assert(isequal(wwidth,params.disp.w_pix))
    assert(isequal(wheight,params.disp.h_pix))
    assert(isequal(round(wwidth/2),params.disp.xc))
    assert(isequal(round(wheight/2),params.disp.yc))
    
    % Tell the experimentor
    fprintf('Pixel size of window is width: %d, height: %d.\n',wwidth,wheight);
    fprintf('Pixel center offset of window is [x,y]=[%d,%d].\n',params.offsetpix(1),params.offsetpix(2));
    
    % recenter EL coordinates if needed (
    % EK: should we just use updated params.stim.xc/yc??
    if any(params.offsetpix~=[0,0]) || isempty(params.offsetpix)
        xc_off = round(wwidth/2) + params.offsetpix(1);
        yc_off = round(wheight/2) + params.offsetpix(2);
    else
        xc_off = round(wwidth/2);
        yc_off = round(wheight/2);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Customize calibration points (include pixel shift, change size and
    % color)
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld',0,0,wwidth-1,wheight-1); % X,Y coordinates left/top/right/bottom of display area
    Eyelink('message','DISPLAY_COORDS %ld %ld %ld %ld',0,0,wwidth-1,wheight-1);
    % IF WE DON"T WANT CUSTOM DOT POSITIONS: Set number of calibration/validation dots and spread: horizontal-only(H) or horizontal-vertical(HV) as H3, HV3, HV5, HV9 or HV13
    Eyelink('command','calibration_type = HV5'); % horizontal-vertical 5-points.
    Eyelink('command','generate_default_targets = NO');
    Eyelink('command','calibration_samples  = 5');
    Eyelink('command','calibration_sequence = 0,1,2,3,4');
    Eyelink('command','calibration_targets  = %d,%d %d,%d %d,%d %d,%d %d,%d',...
        xc_off,yc_off,  ... center x,y
        xc_off + params.stim.el.point2point_distance_pix, yc_off, ... horz shift right
        xc_off - params.stim.el.point2point_distance_pix, yc_off, ... horz shift left
        xc_off, yc_off + params.stim.el.point2point_distance_pix, ... vert shift down
        xc_off, yc_off - params.stim.el.point2point_distance_pix); %  vert shift up
    Eyelink('command','validation_samples = 5');
    Eyelink('command','validation_sequence = 0,1,2,3,4');
    Eyelink('command','validation_targets  = %d,%d %d,%d %d,%d %d,%d %d,%d',...
        xc_off,yc_off,  ... center x,y
        xc_off + params.stim.el.point2point_distance_pix, yc_off, ... horz shift right
        xc_off - params.stim.el.point2point_distance_pix, yc_off, ... horz shift left
        xc_off, yc_off + params.stim.el.point2point_distance_pix, ... vert shift down
        xc_off, yc_off - params.stim.el.point2point_distance_pix); %  vert shift up
    
    Eyelink('command','active_eye = LEFT');
    Eyelink('command','binocular_enabled','NO')
    Eyelink('command','enable_automatic_calibration','NO'); % force manual calibration sequencing, if yes, provide Eyelink('command','automatic_calibration_pacing=1500');
    Eyelink('command','recording_parse_type = GAZE'); %from manual (default)
    Eyelink('command','sample_rate = %d', 1000); % hz
    Eyelink('driftcorrect_cr_disable','YES'); % yes to disable drift correction -- we don't want that!
    %  EyelinkDoDriftCorrection(el); % No drift correction.
    % other ways of drawing things in EL:     Eyelink('Command','draw_box %d %d %d %d %d',xPos-boundary,yPos-boundary,xPos+boundary,yPos+boundary,colorFrm);
    
    % what events (columns) are recorded in EDF:
    Eyelink('command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    % what samples (columns) are recorded in EDF:
    Eyelink('command','file_sample_data = LEFT,RIGHT,GAZE,GAZERES,PUPIL,AREA,STATUS');
    % events available for real time:
    Eyelink('command','link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    % samples available for real time:
    Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,GAZERES,PUPIL,AREA,STATUS');
    
    % make temp name and open
    eyetempfile = sprintf('%s.edf', datestr(now, 'HHMMSS')); %less than 8 digits!
    fprintf('Saving eyetracking data to %s.\n',eyetempfile);
    Eyelink('Openfile',eyetempfile);  % NOTE THIS TEMPORARY FILENAME. REMEMBER THAT EYELINK REQUIRES SHORT FILENAME!
    
    % Send preamble to EDF header
    preamble = sprintf('add_file_preamble_text ''VCD experiment @ CMRR %s, subject %03d, session %03d, run: %02d.''', ...
        params.disp.name, params.subj_nr, params.ses_nr, params.run_nr);
    Eyelink('command', preamble);
    
    % Do calibration & validation!
    checkcalib = input('Do you want to do a calibration (0=no, 1=yes)? ','s');
    if isequal(checkcalib,'1')
        fprintf('Please perform calibration. When done, press the output/record button.\n');
        EyelinkDoTrackerSetup(el);
    end
    
    
    % Start recording
    Eyelink('StartRecording');
    
end


%% START EXPERIMENT

% call ptviewmovie
timeofshowstimcall = datestr(now,30);

% GO!
[data,getoutearly] = vcd_showStimulus(...
    win, rect, params, ...
    scan, ...
    bckground, ...
    fix_im, ...
    eye_im, ...
    stim, ...
    run_frames, ...
    run_table, ...
    introscript, ...
    taskscript, ...
    tfunEYE, ...
    deviceNr, ...
    oldCLUT,...
    oldPriority);

%% CLEAN UP AND SAVE

fprintf('RUN ENDED: %4.4f.\n',GetSecs);

% Close out eyelink
if params.wanteyetracking
    if ~isempty(eyetempfile)
        
        % Close eyelink and record end
        Eyelink('StopRecording');
        Eyelink('message', sprintf('EXP END %d',GetSecs));
        Eyelink('CloseFile');
        status = Eyelink('ReceiveFile',eyetempfile, params.savedatadir, 1);
        fprintf('ReceiveFile status %d\n', status);
        
        % RENAME DOWNLOADED FILE TO THE FINAL FILENAME
        mycmd=['mv ' params.savedatadir '/' eyetempfile ' ' params.eyelinkfile];
        system(mycmd);
        if status <= 0, fprintf('\n\nProblem receiving edf file\n\n');
        else
            fprintf('Data file ''%s'' can be found in ''%s''\n', params.eyelinkfile, pwd);
        end
        Eyelink('ShutDown'); % Do we need to shut down???
        
    end
end

% Create folder where data will be stored
if isempty(params.savedatadir)
    params.savedatadir = fullfile(vcd_rootPath,'data', ...
        sprintf('%s_vcd_subj%d_ses%02d',...
        timeofshowstimcall,params.subj_nr, params.ses_nr));
end

if ~exist('params.savedatadir','dir'), mkdir(params.savedatadir); end

% Create filename where behavioral and timing data will be stored
if ~isfield(params,'behaviorfile') || isempty(params.behaviorfile)
    params.behaviorfile = sprintf('%s_vcd_subj%d_ses%02d_run%02d.mat', ...
        timeofshowstimcall,params.subj_nr,params.ses_nr,params.run_nr);
end

% figure out names of all variables except 'stim', 'fix', and
% 'bckground'
vars = whos;
vars = {vars.name};
vars = vars(cellfun(@(x) ~isequal(x,'fix_im','bckground','stim'),vars));


% Save data (button presses, params, etc)
save(fullfile(params.savedatadir,params.behaviorfile),vars{:}, '-v7.3');


% check the timing
if getoutearly == 0 %if we completed the experiment
    fprintf('Experiment duration was %4.3f.\n',data.timing.endtime);
    slack = [-0.0004, 0.0004]; % (seconds) as much as we allow PTB SyncTime to vary
    expectedduration = (length(data.timing.timeframes)/params.disp.refresh_hz) + slack;
    if data.timing.endtime > expectedduration(1) && data.timing.endtime < expectedduration(2)
        fprintf('Timing was ok and within [%d - %d] sec slack\n',slack(1),slack(2));
    else
        fprintf('ERROR !!! Timing was OFF!!! Difference between expected and recorded is %3.2f seconds\n', expectedduration-data.timing.endtime);
    end
    
    % Visualize timing
    vcd_checkMonitorTiming(data)
   
    if exist(fullfile(params.savedatadir,params.behaviorfile),'file')
        d = dir(fullfile(vcd_rootPath,'tmp_data*.mat'));
        if ~isempty(d)
            delete(fullfile(d(end).folder,d(end).name))
            fprintf('Deleted tmp file because we stored the complete behavioral file.');
        end
    end
    
end





