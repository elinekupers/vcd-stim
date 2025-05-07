function [] = vcd_singleRun(subjID, sesID, runnum, varargin)


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p = inputParser;
p.addRequired('subjID'          , @isnumeric); % subject number
p.addRequired('sesID'           , @isnumeric); % session number
p.addRequired('runnum'          , @isnumeric); % nun number
p.addParameter('scan'           , struct()  , @isstruct);                    % struct with optimized images for ptb
p.addParameter('timing'         , struct()  , @isstruct);                    % struct with frame timing for ptb
p.addParameter('images'         , struct()  , @isstruct);                    % struct with images (we preload before calling vcd_singleRun to avoid delays).
p.addParameter('savedatadir'    , []        , @ischar);                      % place to store data with today's date
p.addParameter('behaviorfile'   , []        , @ischar);                      % filename to store behavioral data with today's date
p.addParameter('eyelinkfile'    , []        , @ischar);                      % where the eyelink edf file can be obtained
p.addParameter('loadparams'     , true      , @islogical)                    % whether load stim/condition params or regenerate
p.addParameter('storeparams'    , true      , @islogical)                    % whether to store stimulus params
p.addParameter('overwrite_randomized_params',false, @islogical);             % whether to overwrite randomization of certain stimulus parameters
p.addParameter('infofolder'     , fullfile(vcd_rootPath,'workspaces','info'), @ischar);     % where the *_info.csv file is
p.addParameter('stimfolder'     , fullfile(vcd_rootPath,'workspaces','stimuli'), @ischar);  % where the images can be obtained
p.addParameter('instrtextdir'   , fullfile(vcd_rootPath,'workspaces','instructions'), @ischar); % where the task instructions can be obtained
p.addParameter('laptopkey'      , -3        , @isnumeric);                   % listen to all keyboards/boxes (is this similar to k=-3;?)
p.addParameter('wanteyetracking', false     , @islogical);                   % whether to try to hook up to the eyetracker
p.addParameter('deviceNr'       , []        , @isnumeric);                   % kbWait/Check input device number
p.addParameter('triggerkey'     , {'5%','t'}, @(x) iscell(x) || isstring(x)) % key that starts the experiment
p.addParameter('triggerkeyname' , '''5'' or ''t''', @isstring)               % for display only
p.addParameter('offsetpix'      , [0 0]     , @isnumeric);                   % offset of screen in pixels [10 20] means move 10-px right, 20-px down
p.addParameter('movieflip'      , [0 0]     , @isnumeric)                    % whether to flip up-down, whether to flip left-right
p.addParameter('debugmode'      , false     , @islogical)                    % whether to use debug mode (no BOLDscreen, no eyelink)
p.addParameter('dispName'       , '7TAS_BOLDSCREEN32' , @ischar)             % display params: 7TAS_BOLDSCREEN32, KKOFFICE_AOCQ3277, PPROOM_EIZOFLEXSCAN, 'EKHOME_ASUSVE247'
p.addParameter('savestimtiming' , false     , @islogical)                    % whether we want to store tempfile with stimulus timing (no stimuli)
p.addParameter('savestim'       , false     , @islogical)                    % whether we want to store temp file with stimuli and timing
p.addParameter('loadstimtiming' , false     , @islogical)                    % whether we want to store temp stim timing file

% Parse inputs
p.parse(subjID, sesID, runnum, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    %     eval([sprintf('%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
    eval([sprintf('params.%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p

% release scan and timing var
scan = params.scan; rmfield(params,'scan');
timing = params.timing; rmfield(params,'timing');

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

if params.debugmode
    params.disp.name = 'PPROOM_EIZOFLEXSCAN';
end

devices = PsychHID('Devices');
for vv = 1:length(devices)
    if strcmp(devices(vv).usageName,'Keyboard') && ~isempty(regexp(devices(vv).product,'\w*Internal Keyboard\w*'))
        deviceNr_tmp.internal = vv;
    elseif strcmp(devices(vv).usageName,'Keyboard') && ~isempty(regexp(devices(vv).product,'\w*USB\w*'))
        deviceNr_tmp.external = vv;
    end
end

if params.debugmode % skip synctest
    skipsync = 1;
    if isempty(params.deviceNr) || isempty(params.deviceNr_tmp)
        deviceNr = -3; % listen to all
        warning('[%s]: No specified device nrs, will listen to all devices!\n',mfilename)
    elseif ~isfield(deviceNr_tmp,'internal') && ~isfield(deviceNr_tmp,'external')
        deviceNr = -3; % listen to all
        warning('[%s]: No internal or external device nrs found, will listen to all devices!\n',mfilename)
    else
        if isfield(deviceNr_tmp,'external') && ~isempty(deviceNr_tmp.external) && ...
                (~isfield(deviceNr_tmp,'internal') || isempty(deviceNr_tmp.internal))
            deviceNr = deviceNr_tmp.external;
            fprintf('[%s]: Using external device number(s): %d, %s %s\n',mfilename,deviceNr,devices(deviceNr).product, devices(deviceNr).manufacturer)
        elseif isfield(deviceNr_tmp,'internal') && ~isempty(deviceNr_tmp.internal) && ...
                (~isfield(deviceNr_tmp,'external') || isempty(deviceNr_tmp.external))
            deviceNr = deviceNr_tmp.internal;
            fprintf('[%s]: Using internal device number(s): %d %s %s\n',mfilename,deviceNr,devices(deviceNr).product, devices(deviceNr).manufacturer)
        elseif isfield(deviceNr_tmp,'external') && isfield(deviceNr_tmp,'internal') && ...
                ~isempty(deviceNr_tmp.external) && ~isempty(deviceNr_tmp.internal)
            deviceNr = [deviceNr_tmp.internal, deviceNr_tmp.external];
            fprintf('[%s]: Using external and internal device number(s): %d, %s %s\n',mfilename,deviceNr,devices(deviceNr).product, devices(deviceNr).manufacturer)
        end
    end
else
    skipsync = 0;
    if isempty(params.deviceNr) && ~isempty(deviceNr_tmp)
        deviceNr = deviceNr_tmp.external;
        fprintf('[%s]: Using external device number: %d, %s %s\n',mfilename,deviceNr,devices(deviceNr).product, devices(deviceNr).manufacturer)
    elseif  isempty(params.deviceNr) && isempty(deviceNr_tmp.external)
        error('[%s]: Cannot find external device number!',mfilename)
    end
end

% Nova1x32 coil with BOLDscreen and big eye mirrors
% expected to be {[1920 1080 120 24],[], 0, 0}
% 1: [width, height, framerate, bitdepth]
% 2: winsize (fraction: default is full extent)
% 3: clutfile -- 0 for linear CLUT (-2 for squaring CLUT for BOLDSCREEN to simulate normal monitors --> NB: we do this manually!)
% 4: skipsync (bool: 0 is false, 1 is true)
% 5: wantstereo (bool: default is false)
ptonparams = {[params.disp.w_pix params.disp.h_pix params.disp.refresh_hz 24],[],params.disp.clut, skipsync};

% Buttonbox / keyboard
params.ignorekeys = KbName({params.triggerkey});  % dont record TR triggers as subject responses


%% %%%%%%%%% SETUP RNG %%%%%%%%%

rand('seed', sum(100*clock));
randn('seed', sum(100*clock));
params.rng.rand = rand;
params.rng.randn = randn;

%% %%%%%%%%%%%%% STIM PARAMS %%%%%%%%%%%%%

% Infer session type
if strcmp(params.disp.name, '7TAS_BOLDSCREEN32')
    session_type = 'MRI';
    
elseif strcmp(params.disp.name,'PPROOM_EIZOFLEXSCAN')
    session_type = 'BEHAVIORAL';
else
    session_type = 'MRI';
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
        load(fullfile(d(end).folder,d(end).name),'exp');
        params.exp = exp; clear exp;
    else
        params.exp = vcd_getSessionParams(params,params.loadparams,params.storeparams);
    end
end

%% %%%%%%%%% LOAD IMAGES & ORDER

if ~exist('scan','var') || ~isfield(scan,'exp_im') || isempty(scan.exp_im)
    % Images are in the format of:
    % subj001_ses01_run01_images_PPROOM_EIZOFLEXSCAN_20250505T184109.mat
    d = dir(fullfile(params.infofolder,sprintf('subj%03d',params.subjID),  ...
        sprintf('subj%03d_ses%02d_run%02d_frames_%s_*.mat', ...
        params.subjID, params.sesID, params.runnum, params.disp.name)));
    a = load(fullfile(d(end).folder,d(end).name));
    scan = a.subj_run_frames; clear a;
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
    centers = cell(size(scan.frame_nr)); %cell(size(timing.seq_stim));
    apsize  = cell(size(scan.frame_nr)); %cell(size(timing.seq_stim));
    
    im_w_mask = uint8.empty(length(scan.frame_nr),0);
    im_w_mask = repmat(mat2cell(im_w_mask,ones(10298,1)),1,2);
    
    for nn = 1:length(scan.frame_nr)
        
        if strcmp(scan.event_name(nn),'stim1') || strcmp(scan.event_name(nn),'stim2')
            
            numSides = find(~cellfun(@isempty, scan.images(nn,:)));
            
            for side = numSides
                
                if strcmp(scan.stim_class_name{nn,side},'gabor')
                    centers{nn,side} = [params.stim.gabor.x0_pix(side) + params.stim.xc, ... % x-coord (pixels)
                        params.stim.gabor.y0_pix(side) + params.stim.yc]; % y-coord (pixels)
                    
                elseif strcmp(scan.stim_class_name{nn,side},'rdk')
                    centers{nn,side} = [params.stim.rdk.x0_pix(side) + params.stim.xc, ...
                        params.stim.rdk.y0_pix(side) + params.stim.yc];
                    
                elseif strcmp(scan.stim_class_name{nn,side},'dot')
                    
                    % deal with dot pol2cart
                    if strcmp(scan.event_name(nn),'stim2')
                        delta_idx = find(scan.stim2_delta(nn,side)==params.stim.dot.delta_from_ref);
                        
                        tmp             = reshape(params.stim.dot.unique_im_nrs_wm_test,4,[]);
                        [~,test_im_idx] = ismember(scan.stim2_im_nr(nn,side),params.stim.dot.unique_im_nrs_wm_test);
                        test_im_sub = ind2sub(test_im_idx, size(tmp,1),size(tmp,2));
                        
                        checkme = [params.stim.dot.x0_pix_delta(test_im_sub(2),delta_idx), params.stim.dot.y0_pix_delta(test_im_sub(2),delta_idx)];
                        dot_angle = scan.stim2_orient_dir(nn,side);
                        [dot_x,dot_y] = pol2cart(deg2rad(dot_angle),params.stim.dot.iso_eccen);
                        
                    elseif strcmp(scan.event_name(nn),'stim1')
                        dot_angle = scan.orient_dir(nn,side);
                        xy_idx = find(dot_angle == params.stim.dot.ang_deg);
                        dot_x = params.stim.dot.x0_pix(xy_idx);
                        dot_y = params.stim.dot.y0_pix(xy_idx);
                    end
                    
                    centers{nn,side} = [dot_x + params.stim.xc, dot_y + params.stim.yc];
                    
                elseif strcmp(scan.stim_class_name{nn,side},'obj')
                    centers{nn,side} = [params.stim.obj.x0_pix(side) + params.stim.xc, ...
                        params.stim.obj.y0_pix(side) + params.stim.yc];
                    
                elseif strcmp(scan.stim_class_name{nn,side},'ns')
                    centers{nn,side} = [params.stim.ns.x0_pix + params.stim.xc, params.stim.ns.y0_pix + params.stim.yc];
                elseif isnan(scan.stim_class_name{nn,side})
                    centers{nn,side} = [NaN, NaN];
                end
                
                if isnan(scan.stim_class_name{nn,side})
                    apsize{nn,side}  = [NaN, NaN];
                else
                    % apsize 1: image width (pixels), apsize 2: image height (pixels)
                    apsize{nn,side}  = [size(scan.images{nn,side},2), size(scan.images{nn,side},1)];
                end
                
                % COMBINE IMAGE AND MASKS, FLIP IM IF REQUESTED
                if isempty(scan.masks{nn,side})
                    scan.images{nn,side} = feval(flipfun, scan.images{nn,side});
                else
                    scan.images{nn,side} = feval(flipfun, cat(3, scan.images{nn,side}, scan.masks{nn,side}));
                    scan.masks{nn,side} = [];
                end
                
            end
        end
    end
end


scan.centers = centers; clear centers;
scan.apsize  = apsize; clear apsize;

%%%%%%% CENTER RECT
nonemptycells = ~cellfun(@isempty, scan.centers(:,1));
emptycells2 = ~cellfun(@isempty, scan.apsize(:,1));
assert(isequal(nonemptycells,emptycells2)); clear emptycells2
centers_shortlist = scan.centers(nonemptycells,:);
apsize_shortlist  = scan.apsize(nonemptycells,:);

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

scan.rects = cell(size(scan.centers));
scan.rects(nonemptycells,:) = rects_shortlist;




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
    bckground = a.bckgrnd_im(:,:,:,params.runnum);
end

if any(params.offsetpix~=[0,0])
    bckground = vcd_pinknoisebackground(params, 'comb', 'fat', 1, params.offsetpix);
end

bckground = feval(flipfun,bckground);

%% %%%%%%%%%%%%% FIX IM %%%%%%%%%%%%%
%  FIX CIRCLE: 5D array: [x,y, 3, lum, type]

if ~exist('fix_im','var') || isempty(fix_im)
    d = dir(fullfile(params.stimfolder,params.disp.name, sprintf('fix_%s*.mat',params.disp.name)));
    a = load(fullfile(d(end).folder, d(end).name),'fix_im','mask');
    fix_im = a.fix_im;
    fix_mask = a.mask;
end


% make fixation dot texture
fix_thin_full   = cell(1,size(fix_im,4));
fix_thick_full  = cell(1,size(fix_im,4));
fix_thick_left  = cell(1,size(fix_im,4));
fix_thick_right = cell(1,size(fix_im,4));
fix_thick_both  = cell(1,size(fix_im,4));
fix_thin_rect   = CenterRect([0 0 round(size(fix_im,1)) round(size(fix_im,2))], [0 0 params.disp.w_pix params.disp.h_pix]);
fix_thick_rect  = CenterRect([0 0 round(size(fix_im,1)) round(size(fix_im,2))], [0 0 params.disp.w_pix params.disp.h_pix]);
fix_thin_rect   = fix_thin_rect + repmat(params.offsetpix,1,2);
fix_thick_rect  = fix_thick_rect + repmat(params.offsetpix,1,2);

for ll = 1:size(fix_im,4) % loop over luminance values
    fix_thin_full{ll}   = feval(flipfun,  cat(3, fix_im(:,:,:,ll,1), fix_mask(:,:,1)));
    fix_thick_full{ll}  = feval(flipfun,  cat(3, fix_im(:,:,:,ll,2), fix_mask(:,:,2)));
    fix_thick_left{ll}  = feval(flipfun,  cat(3, fix_im(:,:,:,ll,3), fix_mask(:,:,3)));
    fix_thick_right{ll} = feval(flipfun,  cat(3, fix_im(:,:,:,ll,4), fix_mask(:,:,4)));
    fix_thick_both{ll}  = feval(flipfun,  cat(3, fix_im(:,:,:,ll,5), fix_mask(:,:,5)));
end

clear fix_im;
fix_im = struct('fix_thin_full',[], 'fix_thick_full', [], ...
    'fix_thick_left', [], 'fix_thick_right', [], ...
    'fix_thick_both', [], 'fix_thin_rect', [], 'fix_thick_rect', []);
fix_im.fix_thin_full   = fix_thin_full;
fix_im.fix_thick_full  = fix_thick_full;
fix_im.fix_thick_left  = fix_thick_left;
fix_im.fix_thick_right = fix_thick_right;
fix_im.fix_thick_both  = fix_thick_both;
fix_im.fix_thin_rect   = fix_thin_rect;
fix_im.fix_thick_rect  = fix_thick_rect;

%% Save images and/or timing just in case we run into an error?
if params.savestimtiming
    fprintf('[%s]: Saving stimulus timing..',mfilename)
    tic
    timing_data = struct('scan',scan, 'params',params,'timing',timing);
    timing_data.scan.exp_im = []; timing_data.scan.exp_im_masks = [];
    timing_data.scan.fix_im = []; timing_data.scan.fix_im_masks = [];
    timing_data.scan.bckground_im = [];
    save(fullfile(params.savedatadir,sprintf('tmp_expim_%s.mat',datestr(now,30))),'timing_data','-v7.3');
    clear tmp_data
    fprintf('done!\n'); toc
end

% if params.savestim
%     fprintf('[%s]: Saving stimuli and timing.. may take a minute!',mfilename)
%     tic
%     timing_data = struct('scan',scan, 'params',params,'timing',timing);
%     save(fullfile(params.savedatadir,sprintf('tmp_expim_%s.mat',datestr(now,30))),'timing_data','-v7.3');
%     clear tmp_data
%     fprintf('done!\n'); toc
% end

%% %%%%%%%%%%%%% TASK INSTRUCTIONS %%%%%%%%%%%%%

introscript = fullfile(params.instrtextdir,'00_runvcdcore_subjectinstructions.txt');
if ~exist(introscript,'file')
    error('[%s]: Can''t find instructions text file!',mfilename')
end

taskIDs = unique(scan.block_ID);
taskIDs = taskIDs(~isnan(taskIDs));
taskIDs = taskIDs(taskIDs~=0);
taskNames = params.exp.crossingnames(taskIDs);

for nn = 1:length(taskIDs)
    taskName = strrep(taskNames{nn},'-','_');
    d = dir(fullfile(params.instrtextdir,sprintf('%02d_runvcdcore_%s.txt', taskIDs(nn),taskName)));
    taskscript{nn} = fullfile(d.folder,d.name);
end


%% %%%%%%%%%%%%% INIT SCREEN %%%%%%%%%%%%%
Screen('Preference', 'SyncTestSettings', .0004);
oldCLUT = pton(ptonparams{:});

win  = firstel(Screen('Windows'));
oldPriority = Priority(MaxPriority(win));

rect = Screen('Rect',win); % what is the total rect
% rect = CenterRect(round([0 0 rect(3)*winsize rect(4)*winsize]),rect);
[win, rect] = Screen('OpenWindow',max(Screen('Screens')),params.stim.bckgrnd_grayval,rect);





%% %%%%%%%%%%%%% Eyelink stuff %%%%%%%%%%%%%
% ANON EYE FUN for SYNC TIME
if params.wanteyetracking
    tfunEYE     = @() cat(2,fprintf('EXP STARTS.\n'),Eyelink('Message','SYNCTIME'));
    
    if ~isfield(params,'eyelinkfile') || isempty(params.eyelinkfile)
        params.eyelinkfile = fullfile(params.savedatadir,sprintf('eye_%s_vcd_subj-%s_run-%d.edf',datestr(now,30),params.subjID,params.runnum));
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
    
    % Update default param struct with VCD needs
    el = vcd_setEyelinkParams(el);
    EyelinkUpdateDefaults(el);
    
    % Get window size
    wwidth = rect(3); wheight = rect(4);  % returns in pixels
    
    % Ensure window size is what we think it is
    assert(isequal(wwidth,disp.w_pix))
    assert(isequal(wheight,disp.h_pix))
    assert(isequal(round(wwidth/2),disp.xc))
    assert(isequal(round(wheight/2),disp.yc))
    
    % Tell the experimentor
    fprintf('Pixel size of window is width: %d, height: %d.\n',wwidth,wheight);
    fprintf('Pixel center offset of window is [x,y]=[%d,%d].\n',params.offsetpix(1),params.offsetpix(2));
    
    % recenter EL coordinates if needed (
    % EK: should we just use updated params.stim.xc/yc??
    if ~iszero(params.offsetpix) || isempty(params.offsetpix)
        xc_off = round(wwidth/2) + params.offsetpix(1);
        yc_off = round(wheight/2) + params.offsetpix(2);
    else
        xc_off = round(wwidth/2);
        yc_off = round(wheight/2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Customize calibration points (include pixel shift, change size and
    % color)
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
    
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld',0,0,wwidth-1,wheight-1); % X,Y coordinates left/top/right/bottom of display area
    Eyelink('message','DISPLAY_COORDS %ld %ld %ld %ld',0,0,wwidth-1,wheight-1);
    % IF WE DON"T WANT CUSTOM DOT POSITIONS: Set number of calibration/validation dots and spread: horizontal-only(H) or horizontal-vertical(HV) as H3, HV3, HV5, HV9 or HV13
    %     Eyelink('command','calibration_type = HV5'); % horizontal-vertical 5-points.
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
    introscript, ...
    taskscript, ...
    tfunEYE, ...
    deviceNr, ...
    oldCLUT,...
    oldPriority);

%% CLEAN UP AND SAVE

% Close out eyelink
if params.wanteyetracking
    if ~isempty(eyetempfile)
        
        % before we close out the eyelink, we send one more syntime message
        Eyelink('Message',eval(tfunEND));
        
        % Close eyelink and record end
        Eyelink('StopRecording');
        Eyelink('message', sprintf('END %d',tGetSecs));
        Eyelink('CloseFile');
        status = Eyelink('ReceiveFile',eyetempfile, params.savedatadir, 1);
        fprintf('ReceiveFile status %d\n', status);
        fprintf('RUN ENDED: %4.4f.\n',GetSecs);
        
        % RENAME DOWNLOADED FILE TO THE FINAL FILENAME
        mycmd=['mv ' params.savedatadir '/' eyetempfile ' ' params.savedatadir '/' params.eyelinkfile];
        system(mycmd);
        if status <= 0, fprintf('\n\nProblem receiving edf file\n\n');
        else
            fprintf('Data file ''%s'' can be found in ''%s''\n', params.eyelinkfile, pwd);
        end
        Eyelink('ShutDown'); % Do we need to shut down???
        
        data.totalTime = ts-startTime;
        data.timeKeys = [data.timeKeys; {ts 'end'}];
    end
else
    fprintf('RUN ENDED: %4.4f.\n',GetSecs);
end

% Create folder where data will be stored
if isempty(params.savedatadir)
    params.savedatadir = fullefile(vcd_rootPath,'data', ...
        sprintf('%s_vcd_subj%d_ses%02d',...
        timeofshowstimcall,params.subjID, params.sesID));
end
if ~exist('params.savedatadir','dir'), mkdir(params.savedatadir); end

% figure out names of all variables except 'images'
vars = whos;
vars = {vars.name};
vars = vars(cellfun(@(x) ~isequal(x,'images'),vars));

% Create filename where behavioral and timing data will be stored
if ~isfield(params,'behaviorfile') || isempty(params.behaviorfile)
    params.behaviorfile = sprintf('%s_vcd_subj%d_ses%02d_run%02d.mat', ...
        timeofshowstimcall,params.subjID,params.sesID,params.runnum);
end

% Save data (button presses, params, etc)
save(fullfile(params.savedatadir,params.behaviorfile),vars{:});

% Clean up
ShowCursor;
Screen('CloseAll');
ptoff(oldCLUT);

% check the timing
if getoutearly == 0 %if we completed the experiment
    fprintf('Experiment duration was %4.3f.\n',data.timing.endtime);
    slack = [-5 5].*params.stim.framedur_s;
    expectedduration = (timing.trig_timing(end)+params.stim.framedur_s) + slack;
    
    if data.timing.endtime > expectedduration(1) && data.timing.endtime < expectedduration(2)
        fprintf('Timing was ok and within [%d - %d] sec slack\n',slack(1),slack(2));
    else
        fprintf('ERROR !!! Timing was OFF!!! Difference between expected and recorded is %3.2f s\n', expectedduration-data.timing.endtime);
    end
end
