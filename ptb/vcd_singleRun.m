function [status,images,maskimages] = vcd_singleRun(subjID, sesID, runnum, varargin)


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p = inputParser;
p.addRequired('subjID'          , @isnumeric); % subject number 
p.addRequired('sesID'           , @isnumeric); % session number 
p.addRequired('runnum'          , @isnumeric); % nun number
p.addParameter('images'         , struct()  , @isstruct);  % struct with images (we preload before calling vcd_singleRun to avoid delays).
p.addParameter('savedatadir'    , []        , @ischar);       % place to store data with today's date
p.addParameter('behaviorfile'   , []        , @ischar);       % filename to store behavioral data with today's date
p.addParameter('eyelinkfile'    , []        , @ischar);       % where the eyelink edf file can be obtained
p.addParameter('loadparams'     , true      , @islogical)       % whether load stim/condition params or regenerate
p.addParameter('infofolder'     , fullfile(vcd_rootPath,'workspaces','info'), @ischar); % where the *_info.csv file is
p.addParameter('stimfolder'     , fullfile(vcd_rootPath,'workspaces','stimuli'), @ischar); % where the images can be obtained
p.addParameter('instrtextdir'   , fullfile(vcd_rootPath,'workspaces','instructions'), @ischar); % where the task instructions can be obtained
p.addParameter('laptopkey'      , -3        , @isnumeric);      % listen to all keyboards/boxes (is this similar to k=-3;?)
p.addParameter('wanteyetracking', false     , @islogical);      % whether to try to hook up to the eyetracker
p.addParameter('triggerkey'     , {'5%','t'}, @(x) iscell(x) || isstring(x)) % key that starts the experiment
p.addParameter('triggerkeyname' , '''5'' or ''t''', @isstring)  % for display only
p.addParameter('offsetpix'      , [0 0]     , @isnumeric);      % offset of screen in pixels [10 20] means move 10-px right, 20-px down
p.addParameter('movieflip'      , [0 0]     , @isnumeric)       % whether to flip up-down, whether to flip left-right
p.addParameter('debugmode'      , false     , @islogical)       % whether to use debug mode (no BOLDscreen, no eyelink)
p.addParameter('storeparams'    , true      , @islogical)       % whether to store stimulus params
p.addParameter('dispName'       , '7TAS_BOLDSCREEN32' , @ischar) % display params: 7TAS_BOLDSCREEN32, KKOFFICE_AOSQ3277, PPROOM_EIZOFLEXSCAN

% Parse inputs
p.parse(subjID, sesID, runnum, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
%     eval([sprintf('%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
    eval([sprintf('params.%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p


%% %%%%%%%%%%%%% PERIPHERALS %%%%%%%%%%%%%

if ~isfield(params, 'disp') || isempty(params.disp)
     params.disp = vcd_getDisplayParams(params.dispName); % BOLDSCREEN is default
end

if params.debugmode % skip syntest
    skipsync = true;
else
    skipsync = false;
end

% Nova1x32 with BOLDscreen and big eye mirrors
% expected to be {[1920 1080 120 24],[], 0, 0} 
% 1: [width, height, framerate, bitdepth]
% 2: winsize (fraction: default is full extent)
% 3: clutfile -- 0 for linear CLUT (choose -2 for squaring CLUT because we need to simulate normal monitors)
% 4: skipsync (bool: 0 is false, 1 is true)
% 5: wantstereo (bool: default is false)
ptonparams = {[params.disp.w_pix params.disp.h_pix params.disp.refresh_hz 24],[],0, skipsync}; 

% Buttonbox / keyboard
params.ignorekeys = KbName({params.triggerkey});  % dont record TR triggers as subject responses

%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SETUP RNG %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('twister', sum(100*clock));
params.rand = rand;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% STIM PARAMS %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **** Block / trial order ****
% We need to define all possible combinations 
% Depending on the session, run, and subject, we need to define the 
% combinations we want to show, converted into miniblocks, where trials are
% sampling the manipulations in a pseudo-random way. In each run, 
% we have manipulations that we prioritize to fully sample, otherwise it is
% difficult to compare conditions (e.g., we want to sample all contrast
% levels within the run).
% 
% We want a master trial function and a function thats in the trial
% specifics given the subject and session nr

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

%% INDEX IMAGES THAT NEED TO BE LOADED
run_image_order = vcd_getImageOrder(subj_run.block, image_info, params);

%% LOAD AND RESIZE IMAGES etc
% exp_im is a cell with dims:
% blocks x trials x locations (1:l, 2:r) x stim epoch (first or second)
[exp_im, images] = vcd_loadRunImages(run_image_order, subj_run.block, params);

%% %%%%%%%%%%%%% FIXATION IM/PARAMS %%%%%%%%%%%%%
% Fixation order and fixation
fixsoafun = @() round(params.stim.fix.dotmeanchange*params.stim.fps + params.stim.fix.dotchangeplusminus*(2*(rand-.5))*params.stim.fps);

% Load stored fixation dot images
if isempty(images.fix)
    fprintf('[%s]: Loading fixation dot images..\n',mfilename);
    % FIX: 5D array: [x,y, 3, 5 lum, 2 widths]
    d = dir(sprintf('%s*.mat', params.stim.fix.stimfile));
    load(fullfile(d(end).folder,d(end).name), 'fix_im','info');
    images.fix = fix_im; clear fix_im;
    images.info.fix = info; clear info;
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
end
fix_im = images.fix;


%% TIMING
seq = {}; % task_cue_ID = 97; ITI_ID = 98; % IBI_ID = 99;
seq_timing = []; 
spatial_cue = []; 

% 6 fields (name, ID, within_session_repeat, trial, trial_type, timing)
% 8 blocks: 1:run, 2:block, 3:stimtaskID, 4:unique_im, 5:spatial_cue, 6:onset_time, 7:event_dur, 8:run_time
cellblock = squeeze(struct2cell(subj_run.block));

for bb = 1:size(cellblock,2)
    tmp_timing = cellblock{6,bb};
    im_nr = tmp_timing.unique_im;
    
    if ~iscell(im_nr)
        % convert to cell
        im_nr = num2cell(im_nr);
    end
    
    
    % check for nans
    for xi = 1:length(im_nr)
        tmp_im = im_nr{xi};
        if isnan(tmp_im)
            tmp_im(isnan(tmp_im))=0;
            im_nr{xi} = tmp_im;
        end
        clear tmp_im
    end  
    seq = cat(1, seq, im_nr);

    seq_timing = cat(1,seq_timing, tmp_timing.run_time);
    spatial_cue = cat(1,spatial_cue, tmp_timing.spatial_cue);
end


%% START HERE (decide on fps)

trig_timing = [0:params.stim.fps:(seq_timing(end)/params.stim.fps)-(params.stim.fps)]'; % seconds
trig = zeros(size(trig_timing,1),2);
for tt = 1:length(seq_timing)
    event_time = seq_timing(tt);
    t_idx = find(trig_timing == event_time);
    if isempty(t_idx)
        error('[%s]: Event_time doesn''t match monitor refresh rate', mfilename)
    end
    if length(seq{tt})==2
        trig(t_idx,:) = seq{tt};
    else
        trig(t_idx,:) = repmat(seq{tt},1,2);
    end
end


for tt = 1:length(trig)
    fix_seq(tt) = fixsoafun;
end


% contrast change seq
cd_seq = [];

% repmat(Expand(params.taskColor,1,2), length(params.seq)/4,1);



%% IMAGE XY CENTER OFFSET
% recenter x,y-center coordinates if needed
if ~iszero(params.offsetpix) || isempty(params.offsetpix)
    xc = (ptonparams{1}(1)/2) + params.offsetpix(1);
    yc = (ptonparams{1}(2)/2) + params.offsetpix(2);
    
    % store offset
    params.stim.xc = xc;
    params.stim.yc = yc;

else
    params.stim.xc = (ptonparams{1}(1)/2);
    params.stim.yc = (ptonparams{1}(2)/2);
end


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
images.bckground = vcd_pinknoisebackground(p, 'comb', 'fat', 1, params.offsetpix); 


%% EK: Do we still want this??
% Get a vector with number of images in each class
% for fn = 1:length(params.exp.stimClassLabels)
%     sz = size(images.(params.exp.stimClassLabels{fn}));
%     if params.stim.(params.exp.stimClassLabels{fn}).iscolor
%         numinclass(fn) = prod(sz(4:end));
%     else
%         numinclass(fn) = prod(sz(3:end));
%     end
% end


%% %%%%%%%%%%%%% EYELINK PARAMS %%%%%%%%%%%%%

% EYE FUN for SYNC TIME
if wanteyetracking
    tfunEYE     = @() cat(2,fprintf('EXP STARTS.\n'),Eyelink('Message','SYNCTIME'));

    if ~isfield(params,'eyelinkfile') || isempty(params.eyelinkfile) 
        params.eyelinkfile = fullfile(params.savedatadir,sprintf('eye_%s_vcd_subj-%s_run-%d.edf',datestr(now,30),params.subjID,params.runnum));
    end
else
    tfunEYE = @() fprintf('EXP STARTS.\n');
end





%% %%%%%%%%%%%%% TASK INSTRUCTIONS %%%%%%%%%%%%%

%%%%%%%%%%%%% Q: add folder ?? --> instrtextfolder
introscript  = 'showinstructionscreen(''runvcdcore_presubjectinstructions.txt'',250,75,25,78)';  % inputs are (fileIn,txtOffset,instTextWrap,textsize,background)

tasks_to_run = struct2cell(subj_run.run);

tasksIDs = cell2mat(tasks_to_run(:,:,2));

for nn = 1:length(tasksIDs)
    d = dir(fullfile(instrtextdir,sprintf('%02d_runvcdcore*.txt', tasksIDs(nn))));
    
    taskscript{nn} = sprintf('showinstructionscreen(''%s'',250,75,25,78)',fullfile(d.folder,d.name));
end





%% %%%%%%%%%%%%% INIT SCREEN %%%%%%%%%%%%%
oldclut = pton(ptonparams{:});


% Flag incase we want to quit early
getoutearly = 0; 

%% %%%%%%%%%%%%% EYELINK: initialize, setup, calibrate, and start %%%%%%%%%%%%%
if wanteyetracking && ~isempty(params.eyelinkfile)
    
    assert(EyelinkInit()==1);
    win = firstel(Screen('Windows'));
    
    if strcmp(dispName,'7TAS_BOLDSCREEN32')
        Eyelink('SetAddress','100.1.1.1') %% copied from MGS exp
    elseif strcmp(dispName,'PPROOM_EIZOFLEXSCAN')
        Eyelink('SetAddress','100.1.1.1') %% <--- CHECK THIS
    end
    
    %eyelink defaults
    el = EyelinkInitDefaults(win);
    if ~isempty(el.callback)
        PsychEyelinkDispatchCallback(el);
    end
    
    el = vcd_setEyelinkParams;
    
    EyelinkUpdateDefaults(el);
    [wwidth,wheight] = Screen('WindowSize',win);  % returns in pixels
    
    assert(isequal(wwidth,disp.w_pix))
    assert(isequal(wheight,disp.h_pix))
    assert(isequal(round(wwidth/2),disp.xc))
    assert(isequal(round(wheight/2),disp.yc))
    
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
    Eyelink('command','file_sample_data = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS');
    % events available for real time:
    Eyelink('command','link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    % samples available for real time:
    Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');

    % make temp name and open
    eyetempfile = sprintf('%s.edf', datestr(now, 'HHMMSS')); %less than 8 digits!
    fprintf('Saving eyetracking data to %s.\n',eyetempfile);
    Eyelink('Openfile',eyetempfile);  % NOTE THIS TEMPORARY FILENAME. REMEMBER THAT EYELINK REQUIRES SHORT FILENAME!
    
    checkcalib = input('Do you want to do a calibration (0=no, 1=yes)? ','s');
    if isequal(checkcalib,'1')
        fprintf('Please perform calibration. When done, press the output/record button.\n');
        EyelinkDoTrackerSetup(el);
    end
    
    Eyelink('StartRecording');
    
end

%% START EXPERIMENT

% call ptviewmovie
timeofshowstimcall = datestr(now);

% GO!
[timeframes, timekeys, digitrecord, trialoffsets] = ...
    vcd_showStimulus(...
        params, ...
        exp_im, ...
        bckground_im, ...
        fix_im, ...
        seq, ...
        seq_timing, ...
        fix_seq, ...
        trig, ...
        introscript, ...
        taskscript, ...
        fixsoafun, ...
        movieflip, ...
        tfunEYE);


%% CLEAN UP AND SAVE

if isempty(params.savedatadir)
    params.savedatadir = fullefile(vcd_rootPath,'data',sprintf('%s_vcd_subj%d_ses%02d',timeofshowstimcall,params.subjID, params.sesID));
end

if ~exist('params.savedatadir','dir'), mkdir(params.savedatadir); end

tENDfun   = @() fprintf('RUN ENDED.\n');

% Close out eyelink
if ~isempty(eyetempfile)
    
    % before we close out the eyelink, we send one more syntime message
    Eyelink('Message',eval(tENDfun));
    
    % Close eyelink and record end
    Eyelink('StopRecording');
    ts = GetSecs;
    Eyelink('message', sprintf('END %d',ts));
    Eyelink('CloseFile');
    status = Eyelink('ReceiveFile',eyetempfile, params.savedatadir, 1);
    fprintf('ReceiveFile status %d\n', status);
    
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
else
    eval(tENDfun);
end


% figure out names of all variables except 'images' and 'maskimages' and others   [MAKE THIS INTO A FUNCTION?]
vars = whos;
vars = {vars.name};
vars = vars(cellfun(@(x) ~isequal(x,'images'),vars));

if ~isfield(params,'behaviorfile') || isempty(params.behaviorfile)
    params.behaviorfile = sprintf('%s_vcd_subj%d_ses%02d_run%02d.mat',timeofshowstimcall,params.subjID,params.sesID,params.runnum);
end

% Save data (button presses, params, etc)
save(fullfile(params.savedatadir,params.behaviorfile),vars{:});

ShowCursor;
Screen('CloseAll');

% unsetup PT
ptoff(oldclut);

% check the timing
if getoutearly == 0 %if we completed the experiment
%     fprintf('Experiment duration was %.3f.\n',data.totalTime);
% %     expectedduration = [-0.2, 0.2] + (length(X)/10); %[67.0 67.4];  % expected duration for 5 point calibration 4 trials
%     
%     if data.totalTime > expectedduration(1) && data.totalTime < expectedduration(2)
%       fprintf('Timing was fine!\n');
%     else
%       fprintf('ERROR! Timing was off!!!!');
%     end
end
