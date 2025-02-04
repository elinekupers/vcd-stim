function [images,maskimages] = vcd_singleRun(subjnum, runnum, varargin)


%% Read inputs
p = inputParser;
p.addRequired ('subjnum'        , @isnumeric);
p.addRequired ('runnum'         , @isnumeric);
p.addParameter('root'           , vcd_rootPath, @isstring);      % root folder
p.addParameter('folder'         , datestr(now,30), @isstring);  % today's date
p.addParameter('infofile'       , fullfile(vcd_rootPath,'stimuli','workspaces','classinfo.mat'), @isstring); % where the classinfo.mat file is
p.addParameter('stimfile'       , fullfile(vcd_rootPath,'stimuli','workspaces','stimfile.mat'), @isstring); % where the 'images' variables can be obtained
p.addParameter('eyelinkfile'    , []        , @isstring);       % where the eyelink edf file can be obtained
p.addParameter('laptopKey'      , -3        , @isnumeric);      % listen to all keyboards/boxes (is this similar to k=-3;?)
p.addParameter('wanteyetracking', false     , @islogical);      % whether to try to hook up to the eyetracker
p.addParameter('triggerkey'     , {'5%','t'}, @(x) iscell(x) || isstring(x)) % key that starts the experiment
p.addParameter('triggerkeyname' , '''5'' or ''t''', @isstring)  % for display only
p.addParameter('offset_pix'     , [0 0]     , @isnumeric);      % offset of screen in pixels [10 20] means move 10-px right, 20-px down
p.addParameter('movieflip'      , [0 0]     , @isnumeric)       % whether to flip up-down, whether to flip left-right
p.addParameter('debugmode'      , false     , @islogical)       % whether to use debug mode (no BOLDscreen, no eyelink)


% Parse inputs
p.parse(subjnum, runnum, varargin{:});

% Rename variables into general params struct
params    = struct();
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    eval([sprintf('params.%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff

%% DISPLAY PARAMS

if params.debugmode
    disp = vcd_getDisplayParams('KKOFFICEQ3277');
    wanteyetracking = false;
    Screen('Preference', 'SkipSyncTests', 1)
else
    disp = vcd_getDisplayParams(); % BOLDSCREEN is default
end

% Nova1x32 with BOLDscreen and big eye mirrors:
ptonparams = {[disp.w_pix disp.h_pix disp.refresh_hz 24],[],-2};  % {[1920 1080 120 24],[],-2} BOLDSCREEN (squaring CLUT because we need to simulate normal monitors)

%% KEYPRESS PARAMS
params.ignorekeys = KbName({params.triggerkey});  % dont record TR triggers as subject responses

%% STIM PARAMS
params.stim = vcd_getStimParams('all');

% ANON FUNCTIONS
soafun = @() round(params.stim.fix.dotmeanchange*params.stim.fps + params.stim.fix.dotchangeplusminus*(2*(rand-.5))*params.stim.fps);

%% Other PARAMS
setupscript    = 'showinstructionscreen(''runvcdcore_subjectinstructions.txt'',250,75,25,78)';  % inputs are (fileIn,txtOffset,instTextWrap,textsize,background)

%% Set-up rand
rand('twister', sum(100*clock));
params.rand = rand;

%% EYELINK PARAMS

% EYE FUN for SYNC TIME
if params.wanteyetracking
    tfun = @() cat(2,fprintf('STIMULUS STARTED.\n'),Eyelink('Message','SYNCTIME'));
    if ~isfield(params, 'eyelinkfile') || isempty(params.eyelinkfile) 
        eyelinkfile = fullfile(vcd_rootPath,'data',sprintf('eye_%s_vcd_subj-%s_run-%d.edf',datestr(now,'yyyymmdd-HHMMSS'),subjnum,runnum));
    else
        eyelinkfile = params.eyelinkfile;
    end
else
    tfun = @() fprintf('STIMULUS STARTED.\n');  % this is executed when the expt starts
end


%% Block / trial order

% We need to define all possible combinations 
% Depending on the session, run, and subject, we need to define the 
% combinations we want to show, converted into miniblocks, where trials are
% sampling the manipulations in a pseudo-random way. In each run, 
% we have manipulations that we prioritize to fully sample, otherwise it is
% difficult to compare conditions (e.g., we want to sample all contrast
% levels within the run).

% We want a master trial function and a function thats in the trial specifics given the subject and session nr 
trials = vcd_makeTrials(params);

%% Fixation order and fixation color




%% Task instructons



%% Other stuff
% %% Set scale factor (ERK: this is already defined in stim params)
% if ~isempty(dres) && length(dres)==1
%     scfactor = -dres;
% else
%     scfactor = [];
% end

%% LOAD in the stimuli

images = struct('gabor',[],'rdk',[],'dot',[],'cobj',[], 'ns',[]);

% ERK: Why are GABORS in 4 bins the same??
% GABORS: 4 oriented bins (cells), each cell has [x,y,contrast,phase,angle]

images.gabor = load(params.stim.gabor.stimfile);
images.rdk   = load(params.stim.rdk.stimfile, 'rdk');
images.dot   = load(params.stim.dot.stimfile, 'dot');
images.cobj  = load(params.stim.cobj.stimfile, 'objects');
images.ns    = load(params.stim.ns.stimfile, 'scenes');


%% Initalize screen
oldclut = pton(ptonparams{:});

% [windowPtr,center,blankColor] = doScreen; % Opens screen, hides cursor, etc

% Flag incase we want to quit early
getoutearly = 0; 

%% EYELINK: initialize, setup, calibrate, and start 
if wanteyetracking && ~isempty(eyelinkfile)
    
    assert(EyelinkInit()==1);
    win = firstel(Screen('Windows'));
    
    %eyelink defaults
    el = EyelinkInitDefaults(win);
    if ~isempty(el.callback)
        PsychEyelinkDispatchCallback(el);
    end
    
    el = vcd_setEyelinkParams;
    
    EyelinkUpdateDefaults(el);
    [wwidth,wheight] = Screen('WindowSize',win);  % returns in pixels
    
%     % recenter EL coordinates if needed
%     if ~iszero(offset_pix) || isempty(offset_pix)
%         centerXY = [round(wwidth/2), round(wheight/2)] + offset_pix;
%         
%     end

    
    fprintf('Pixel size of window is width: %d, height: %d.\n',wwidth,wheight);
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld',0,0,wwidth-1,wheight-1); % X,Y coordinates left/top/right/bottom of display area
    Eyelink('message','DISPLAY_COORDS %ld %ld %ld %ld',0,0,wwidth-1,wheight-1);
    % Set number of calibration/validation dots and spread: horizontal-only(H) or horizontal-vertical(HV) as H3, HV3, HV5, HV9 or HV13
    Eyelink('command','calibration_type = HV5'); % horizontal-vertical 5-points.
    Eyelink('command','active_eye = LEFT');
    Eyelink('command','binocular_enabled','NO')
    Eyelink('command','enable_automatic_calibration','NO'); % force manual calibration sequencing, if yes, provide Eyelink('command','automatic_calibration_pacing=1500');
    Eyelink('command','recording_parse_type = GAZE'); %from manual (default)
    Eyelink('command','sample_rate = %d', 1000); % hz
    Eyelink('driftcorrect_cr_disable','YES'); % disable drift correction -- we don't want that!
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
    % note that we expect that something should probably issue the command:
    %   Eyelink('Message','SYNCTIME');
    % before we close out the eyelink.

end

%% START EXPERIMENT

% call ptviewmovie
timeofshowstimcall = datestr(now);

% GO!
[timeframes, timekeys, digitrecord, trialoffsets] = ...
    vcd_showStimulus(p, trials, images, setupscript, soafun);




%% CLEAN UP AND SAVE

% Close out eyelink
if ~isempty(tmp_eyelinkfile)
    
    % Close eyelink and record end
    Eyelink('StopRecording');
    ts = GetSecs;
    Eyelink('message', sprintf('END %d',ts));
    Eyelink('CloseFile');
    pathToSaveFileIn = fullfile(vcd_rootPath,'data',sessionFolder);
    status = Eyelink('ReceiveFile',tmp_eyelinkfile, pathToSaveFileIn,1);
    fprintf('ReceiveFile status %d\n', status);
    % RENAME DOWNLOADED FILE TO THE FINAL FILENAME
    mycmd=['mv ' pathToSaveFileIn '/' tmp_eyelinkfile ' ' pathToSaveFileIn '/' eyelinkfile]; 
    system(mycmd);
    if status <= 0, fprintf('\n\nProblem receiving edf file\n\n');
    else
        fprintf('Data file ''%s'' can be found in ''%s''\n', subject.eyefilename, pwd);
    end
    Eyelink('ShutDown'); % Do we need to shut down???
    
    subject.totalTime = ts-startTime;
    subject.timeKeys = [subject.timeKeys; {ts 'end'}];
end


% figure out names of all variables except 'images' and 'maskimages' and others   [MAKE THIS INTO A FUNCTION?]
vars = whos;
vars = {vars.name};
vars = vars(cellfun(@(x) ~isequal(x,'images') & ~isequal(x,'maskimages') ...
                       & ~isequal(x,'validlocations') & ~isequal(x,'A') ...
                       & ~isequal(x,'trialtask') & ~isequal(x,'mashgaponecyclefun'),vars));

subject.outfile = sprintf('sub-%s_run-%d_%s.mat', subject.name, subject.runnum, subject.timestamp);

% Save data (button presses, params, etc)
save(outfile,vars{:});

ShowCursor;
Screen('CloseAll');

% unsetup PT
ptoff(oldclut);


% check the timing
if getoutearly == 0 %if we completed the experiment
    fprintf('Experiment duration was %.3f.\n',subject.totalTime);
    expectedduration = [-0.2, 0.2] + (length(dot_cond)/10); %[67.0 67.4];  % expected duration for 5 point calibration 4 trials
    
    if subject.totalTime > expectedduration(1) && subject.totalTime < expectedduration(2)
      fprintf('Timing was fine!\n');
    else
      fprintf('ERROR! Timing was off!!!!');
    end
end
