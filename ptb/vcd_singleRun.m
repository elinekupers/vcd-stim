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
p.addParameter('savetempstimuli', false     , @islogical)                    % whether we want to store temp stim file
p.addParameter('loadtempstimuli', false     , @islogical)                    % whether we want to store temp stim file

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

%% %%%%%%%%%%%%% PERIPHERALS %%%%%%%%%%%%%

if ~isfield(params, 'disp') || isempty(params.disp)
    params.disp = vcd_getDisplayParams(params.dispName); % BOLDSCREEN is default
end

devices = PsychHID('Devices');
for vv = 1:length(devices)
    if strcmp(devices(vv).usageName,'Keyboard') && ~isempty(regexp(devices(vv).product,'\w*Internal Keyboard\w*'))
        deviceNr.internal = vv;
    elseif strcmp(devices(vv).usageName,'Keyboard') && ~isempty(regexp(devices(vv).product,'\w*USB\w*'))
        deviceNr.external = vv;
    end
end

if params.debugmode % skip synctest
    skipsync = 1;
    if isempty(params.deviceNr)
        params.deviceNr = -3; % listen to all
    elseif ~isfield(deviceNr,'internal') && ~isfield(deviceNr,'external')
        params.deviceNr = -3; % listen to all
    else
        if isfield(deviceNr,'external') && ~isempty(deviceNr.external) && ...
                (~isfield(deviceNr,'internal') || isempty(deviceNr.internal))
            params.deviceNr = deviceNr.external;
        elseif isfield(deviceNr,'internal') && ~isempty(deviceNr.internal) && ...
                (~isfield(deviceNr,'external') || isempty(deviceNr.external))
            params.deviceNr = deviceNr.internal;
        elseif isfield(deviceNr,'external') && isfield(deviceNr,'internal') && ...
                ~isempty(deviceNr.external) && ~isempty(deviceNr.internal)
            params.deviceNr = [deviceNr.internal, deviceNr.external];
        end
    end
else
    skipsync = 0;
    if isempty(params.deviceNr) && ~isempty(deviceNr.external)
        params.deviceNr = deviceNr.external;
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

% Stimulus params
if ~isfield(params, 'stim') || isempty(params.stim)
    if params.loadparams
        d = dir(fullfile(params.infofolder,'stim*.mat'));
        load(fullfile(d(end).folder,d(end).name),'stim');
        params.stim = stim; clear stim;
    else
        params.stim   = vcd_getStimParams('all',disp.name,params.loadparams,params.storeparams);
    end
end

% Trial params
% **** Block / trial order ****
% We need to define all possible combinations. Depending on the session,
% run, and subject, we need to define the combinations we want to show,
% converted into miniblocks, where trials are sampling the manipulations in
% a pseudo-random way. In each run, we have manipulations that we prioritize
% to fully sample, otherwise it is difficult to compare conditions (e.g.,
% we want to sample all contrast levels within the run).

if ~isfield(params, 'trials') || isempty(params.trials)
    if params.loadparams
        d = dir(fullfile(params.infofolder,'trials*.mat'));
        load(fullfile(d(end).folder,d(end).name),'all_trials');
        params.trials = all_trials; clear all_trials;
    else
        params.trials = vcd_makeTrials(params,params.loadparams,params.storeparams);
    end
end

% Session params
if ~isfield(params, 'exp') ||  isempty(params.exp)
    if params.loadparams
        d = dir(fullfile(params.infofolder,'exp_session*.mat'));
        load(fullfile(d(end).folder,d(end).name),'exp_session');
        params.exp = exp_session; clear exp_session;
    else
        params.exp = vcd_getSessionParams(params,params.loadparams,params.storeparams);
    end
end

%% %%%%%%%%%%%%% GET MINIBLOCK ORDER %%%%%%%%%%%%%
% We have master trial sequence, and a sequence where blocks are shuffled
% given the subject.

if params.loadparams
    d = dir(fullfile(params.infofolder,'subject_sessions_*.mat'));
    load(fullfile(params.infofolder,d(end).name),'subject_sessions'); % pick the last one we saved
    % subject_sessions is a struct
    subj_run = subject_sessions(params.sesID,params.subjID).run(params.runnum);
else
    subject_sessions = vcd_getSessionParams(params,params.loadparams,params.storeparams);
    subj_run = subject_sessions(params.sesID,params.subjID).run(params.runnum);
end

for jj = 1:length(params.exp.stimClassLabels)
    d = dir(sprintf('%s*.csv',params.stim.(params.exp.stimClassLabels{jj}).infofile));
    image_info.(params.exp.stimClassLabels{jj}) = readtable(fullfile(d(end).folder,d(end).name));
end



%% %%%%%%%%% INDEX IMAGES THAT NEED TO BE LOADED
[run_image_order,im_seq_order] = vcd_getImageOrder(subj_run.block, image_info, params);

if ~exist('scan','var') || ~isfield(scan,'exp_im') || isempty(scan.exp_im)
    % LOAD AND RESIZE IMAGES etc
    
    % exp_im is a cell with dims:
    % blocks x trials x locations (1:l, 2:r) x stim epoch (first or second)
    [exp_im, alpha_masks, images] = vcd_loadRunImages(run_image_order, subj_run.block, params);
    
    scan.exp_im = exp_im; clear exp_im
    scan.exp_im_masks = alpha_masks; clear alpha masks;
end

%% %%%%%%%%%%%%% FIXATION IM/PARAMS %%%%%%%%%%%%%

if ~exist('scan','var') || ~isfield(scan,'fix_im') || isempty(scan.fix_im)
    
    % Load stored fixation dot images
    if isempty(images.fix)
        fprintf('[%s]: Loading fixation dot images..\n',mfilename);
        % FIX: 5D array: [x,y, 3, 5 lum, 2 widths]
        d = dir(sprintf('%s*.mat', params.stim.fix.stimfile));
        load(fullfile(d(end).folder,d(end).name), 'fix_im','mask','info');
        images.fix = fix_im; clear fix_im;
        images.info.fix = info; clear info;
        images.alpha.fix = mask; clear mask;
    end
    
    % rescale fix dot if needed:
    if ~isempty(params.stim.fix.dres) && ...
            length(params.stim.fix.dres)==2 && ...
            sum(params.stim.fix.dres)~=0
        
        fprintf('[%s]: Resampling fixation dot images..\n',mfilename);
        old_fix_sz = size(images.fix);
        p.stim.(stimClass).iscolor
        tmp_im = reshape(images.fix, ...
            old_fix_sz(1),... x
            old_fix_sz(2),... y
            3, ... % rgb
            []); % 5 luminance and 2 rim widths
        
        numIn = size(tmp_im,4);
        tmp_im = double(tmp_im);
        for kk = 1:numIn
            tmp0 = imresize(tmp_im(:,:,:,kk),p.stim.fix.dres); %% DEFAULT IS BICUBIC
            
            % square for CLUT
            tmp0 = double(tmp0);
            pix_range = [min(tmp0),max(tmp0)];
            
            % if pix lum range is 0-1
            if pix_range(2)<=1
                img  = floor((tmp0*255)+params.stim.bckgrnd_grayval);
                img = img.^2;
                
                % if pix lum range is 1-255
            elseif pix_range(2)<=255
                tmp0 = tmp0.^2;
            end
            im_rz(:,:,:,kk) = uint8(tmp0);
        end
        
        images.fix = reshape(im_rz,old_fix_sz(1),old_fix_sz(2),old_fix_sz(3),old_fix_sz(4),old_fix_sz(5));
        
        tmp_alpha = double(images.alpha.fix);
        tmp_alpha = imresize(tmp_alpha,p.stim.fix.dres); %% DEFAULT IS BICUBIC
        images.alpha.fix = uint8(tmp_alpha);
    end
    
    if ~isfield(scan,'fix_im') || ~isempty(scan.fix_im)
        scan.fix_im = images.fix;
        scan.fix_alpha_mask = images.alpha.fix;
    end
    
end

%% %%%%%%%%% TIMING

if ~exist('timing','var') || ~isfield(timing,'seq_stim') || isempty(timing.seq_stim)
    timing = vcd_getImageTiming_framelocked30Hz(params, subj_run, im_seq_order, scan.exp_im, scan.fix_im, scan.exp_im_masks);
    %     timing = vcd_getImageTiming_systemclocklocked(params, subj_run, im_seq_order,  scan.exp_im, scan.exp_im_masks, scan.fix_im, scan.fix_im_masks);
    
    cellblock = struct2cell(subj_run.block);
    cellblock = squeeze(cellblock);
    
    celltiming = (timing.seq_stim);
    
    st_ID = cell2mat(cellblock(2,:));
    st_start = []; st_end = [];
    for ii = 1:length(st_ID)
        st_start = cat(1,st_start, cellblock{6,ii}.onset_time(1));
        st_end = cat(1,st_end, cellblock{6,ii}.run_time(end));
    end
    
    timing.block = [st_ID', st_start,st_end];
    
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
    
    % 1-30 = all non-NS stim-task crossings
    % 31-39 = NS stim-task crossings
    gb_IDs   = find(~cellfun(@isempty, (regexp(params.exp.stimTaskLabels,'-gabor','ONCE'))))';
    rdk_IDs  = find(~cellfun(@isempty, (regexp(params.exp.stimTaskLabels,'-rdk','ONCE'))))';
    dot_IDs  = find(~cellfun(@isempty, (regexp(params.exp.stimTaskLabels,'-dot','ONCE'))))';
    cobj_IDs = find(~cellfun(@isempty, (regexp(params.exp.stimTaskLabels,'-cobj','ONCE'))))';
    ns_IDs   = find(~cellfun(@isempty, (regexp(params.exp.stimTaskLabels,'-ns','ONCE'))))';
    
    % Get [x,y]-center in pixels of peripheral stimuli given display size, and
    % offset. Get stimulus aperture size..
    centers = cell(size(timing.seq_stim));
    apsize  = cell(size(timing.seq_stim));
    
    unique_blocks =  [unique(timing.trig_block,'stable')];
    
    wm_blocks = find(~cellfun(@isempty, regexp(params.exp.stimTaskLabels,'wm')));
    ltm_blocks = find(~cellfun(@isempty, regexp(params.exp.stimTaskLabels,'ltm')));
    img_blocks = find(~cellfun(@isempty, regexp(params.exp.stimTaskLabels,'img')));
    
    trial_counter = 1;
    block_counter = 1;
    
    for nn = 1:length(timing.trig_block)
        block_nr = timing.trig_block(nn);
        
        if (block_nr > 0) && (block_nr < 90)
            prev_block_counter = block_counter;
            block_counter = find(unique_blocks==block_nr);
            
            if prev_block_counter ~= block_counter % new block!
                trial_counter = 0;
            end
            
            if  (timing.trig_stim(nn,1)<90)
                if timing.trig_stim(nn,1) ~= timing.trig_stim(nn-1,1)
                    trial_counter = trial_counter + 1;
                    if timing.trig_stim(nn-1,1)==96
                        delay_period = 1;
                    else
                        delay_period = 0;
                    end
                end
                
                
                numSides = length(unique(timing.trig_stim(nn,:)));
                
                for side = 1:numSides
                    
                    if ~isempty(intersect(block_nr,gb_IDs))
                        centers{nn,side}(1) = params.stim.gabor.x0_pix(side) + params.stim.xc; % x-coord (pixels)
                        centers{nn,side}(2) = params.stim.gabor.y0_pix(side) + params.stim.yc; % y-coord (pixels)
                        
                    elseif ~isempty(intersect(block_nr,rdk_IDs))
                        centers{nn,side}(1) = params.stim.rdk.x0_pix(side) + params.stim.xc;
                        centers{nn,side}(2) = params.stim.rdk.y0_pix(side) + params.stim.yc;
                        
                    elseif ~isempty(intersect(block_nr,dot_IDs))
                        
                        if ~isempty(intersect(block_nr,[wm_blocks,ltm_blocks,img_blocks]))
                            trial_counter0 = round(trial_counter./2);
                        end
                        
                        % deal with dot pol2cart
                        if delay_period && ~isnan(cellblock{4,block_counter}(trial_counter0,side).ref_delta)
                            dot_angle = cellblock{4,block_counter}(trial_counter0,side).ref_delta;
                            [dot_x,dot_y] = pol2cart(deg2rad(dot_angle),params.stim.dot.iso_eccen);
                        elseif isfield(cellblock{4,block_counter}(trial_counter,side),'loc_deg')
                            dot_angle = cellblock{4,block_counter}(trial_counter,side).loc_deg;
                            xy_idx = find(dot_angle == params.stim.dot.loc_deg);
                            dot_x = params.stim.dot.x0_pix(xy_idx);
                            dot_y = params.stim.dot.y0_pix(xy_idx);
                        end
                        
                        centers{nn,side}(1) = dot_x + params.stim.xc;
                        centers{nn,side}(2) = dot_y + params.stim.yc;
                        
                    elseif ~isempty(intersect(block_nr,cobj_IDs))
                        centers{nn,side}(1) = params.stim.cobj.x0_pix(side) + params.stim.xc;
                        centers{nn,side}(2) = params.stim.cobj.y0_pix(side) + params.stim.yc;
                        
                    elseif ~isempty(intersect(block_nr,ns_IDs))
                        centers{nn,side}(1) = params.stim.ns.x0_pix(side) + params.stim.xc;
                        centers{nn,side}(2) = params.stim.ns.y0_pix(side) + params.stim.yc;
                    end
                    
                    apsize{nn,side}(1)  = size(timing.trig_seq_exp_im_w_cd{nn}{1},2); % apsize 1: image width (pixels)
                    apsize{nn,side}(2)  = size(timing.trig_seq_exp_im_w_cd{nn}{1},1); % apsize 2: height (pixels)
                    %                     disp(apsize{nn,side})
                end
            end
        end
    end
end


scan.centers = centers; clear centers trial_counter block_counter;
scan.apsize  = apsize; clear apsize;

%%%%%%% CENTER RECT
nonemptycells = ~cellfun(@isempty, scan.centers(:,1));
emptycells2 = ~cellfun(@isempty, scan.apsize(:,1));
assert(isequal(nonemptycells,emptycells2)); clear emptycells2
scan.centers_shortlist = scan.centers(nonemptycells,:);
scan.apsize_shortlist = scan.apsize(nonemptycells,:);

% insert NaNs for second column when using single square stimulus (nat scene)
scan.centers_shortlist(cellfun(@isempty,scan.centers_shortlist(:,2)),2) = ...
    repmat({[NaN,NaN]},size(find(cellfun(@isempty,scan.centers_shortlist(:,2))),1),1);

scan.apsize_shortlist(cellfun(@isempty,scan.apsize_shortlist(:,2)),2) = ...
    repmat({[NaN,NaN]},size(find(cellfun(@isempty,scan.apsize_shortlist(:,2))),1),1);

scan.rects_shortlist = cell(size(scan.apsize_shortlist));

for side = [1,2]
    scan.rects_shortlist(:,side) = cellfun(@(stimsize,stimcenter) CenterRectOnPoint([0 0 stimsize(1) stimsize(2)], stimcenter(1), stimcenter(2)), ...
        scan.apsize_shortlist(:,side),scan.centers_shortlist(:,side), 'UniformOutput', false);
end

scan.rects = cell(size(scan.centers));
scan.rects(nonemptycells,:) = scan.rects_shortlist;




%% %%%%%%%%%%%%% BACKGROUND IM %%%%%%%%%%%%%
%  BACKGROUND: 3D array: [x,y, num images]
% input 2: 'puzzle'  (square central image + peripheral apertures,
%          'dotring' (simple dot iso eccen donut),
%          'comb'    (puzzle + dotring overlayed)
% input 3: 'skinny'  (no space between stim and edge of background) or
%          'fat'     (+2 deg from stim)
% input 4: number of unique noise background images
% input 5: offset of center in pixels [x,y]
% create background image that adjusts for shifted background if needed
if ~exist('scan','var') || ~isfield(scan, 'bckground') || isempty(scan.bckground)
    scan.bckground = vcd_pinknoisebackground(params, 'comb', 'fat', 1, params.offsetpix);
end

%% Save images and timing just in case we run into an error?
if params.savetempstimuli && ~params.loadtempstimuli % assume we don't to resave the same stimuli
    save(fullfile(params.savedatadir,sprintf('tmp_expim_%s.mat',datestr(now,30))),'scan','timing','params','-v7.3')
end

%% %%%%%%%%%%%%% TASK INSTRUCTIONS %%%%%%%%%%%%%

% inputs for showinstructionscreen: (fileIn,txtOffset=250,instTextWrap=75,textsize=25,background=128,offset=[0,0])
% introscript  = sprintf('showinstructionscreen(''%s'',''00_runvcdcore_subjectinstructions.txt'',250,75,25,%d,[%d %d])',params.instrtextdir,params.stim.bckgrnd_grayval,params.offsetpix(1),params.offsetpix(2));
introscript = fullfile(params.instrtextdir,'00_runvcdcore_subjectinstructions.txt');
if ~exist(introscript,'file')
    error('[%s]: Can''t find instructions text file!',mfilename')
end

tasks_to_run = squeeze(struct2cell(subj_run.block));

taskNames = tasks_to_run(1,:);
tasksIDs = cell2mat(tasks_to_run(2,:));

task_idx = find(cellfun(@isempty, regexp(taskNames,'blank')));
for nn = 1:length(task_idx)
    d = dir(fullfile(params.instrtextdir,sprintf('%02d_runvcdcore*.txt', tasksIDs(task_idx(nn)))));
    taskscript{nn} = fullfile(d.folder,d.name);
end


%% %%%%%%%%%%%%% INIT SCREEN %%%%%%%%%%%%%
Screen('Preference', 'SyncTestSettings', .0004);
oldclut = pton(ptonparams{:});

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
    [wwidth,wheight] = [rect(3), rect(4)];  % returns in pixels
    
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
        xc_off + params.stim.el.point2point_distance, yc_off, ... horz shift right
        xc_off - params.stim.el.point2point_distance, yc_off, ... horz shift left
        xc_off, yc_off + params.stim.el.point2point_distance, ... vert shift down
        xc_off, yc_off - params.stim.el.point2point_distance); %  vert shift up
    Eyelink('command','validation_samples = 5');
    Eyelink('command','validation_sequence = 0,1,2,3,4');
    Eyelink('command','validation_targets  = %d,%d %d,%d %d,%d %d,%d %d,%d',...
        xc_off,yc_off,  ... center x,y
        xc_off + params.stim.el.point2point_distance, yc_off, ... horz shift right
        xc_off - params.stim.el.point2point_distance, yc_off, ... horz shift left
        xc_off, yc_off + params.stim.el.point2point_distance, ... vert shift down
        xc_off, yc_off - params.stim.el.point2point_distance); %  vert shift up
    
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
[data,getoutearly] = vcd_showStimulus(win, rect,...
    params, ...
    scan, ...
    timing, ...
    introscript, ...
    taskscript, ...
    tfunEYE, ...
    params.deviceNr, ...
    oldclut,  oldPriority);


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
ptoff(oldclut);

% check the timing
if getoutearly == 0 %if we completed the experiment
    fprintf('Experiment duration was %4.3f.\n',data.timing.endtime);
    slack = [-5 5].*params.stim.fps;
    expectedduration = (timing.trig_timing(end)+params.stim.fps) + slack;
    
    if data.timing.endtime > expectedduration(1) && data.timing.endtime < expectedduration(2)
        fprintf('Timing was ok and within [%d - %d] sec slack\n',slack(1),slack(2));
    else
        fprintf('ERROR !!! Timing was OFF!!! Difference between expected and recorded is %3.2f s\n', expectedduration-data.timing.endtime);
    end
end
