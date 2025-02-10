function [status,images,maskimages] = vcd_singleRun(subjID, sesID, runnum, varargin)


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p = inputParser;
p.addRequired ('subjID'        , @isnumeric);
p.addRequired ('sesID'          , @isnumeric);
p.addRequired ('runnum'         , @isnumeric);
p.addParameter('root'           , vcd_rootPath, @isstring);      % root folder
p.addParameter('tmpFolder'      , datestr(now,30), @isstring);  % today's date
p.addParameter('loadparams'     , true      , @islogical)       % whether load stim/condition params or regenerate
p.addParameter('infofolder'     , fullfile(vcd_rootPath,'workspaces','info'), @isstring); % where the *_info.csv file is
p.addParameter('stimfolder'     , fullfile(vcd_rootPath,'workspaces','stimuli'), @isstring); % where the images can be obtained
p.addParameter('instrtextfolder', fullfile(vcd_rootPath,'workspaces','instructions'), @isstring); % where the task instructions can be obtained
p.addParameter('eyelinkfile'    , []        , @isstring);       % where the eyelink edf file can be obtained
p.addParameter('laptopkey'      , -3        , @isnumeric);      % listen to all keyboards/boxes (is this similar to k=-3;?)
p.addParameter('wanteyetracking', false     , @islogical);      % whether to try to hook up to the eyetracker
p.addParameter('triggerkey'     , {'5%','t'}, @(x) iscell(x) || isstring(x)) % key that starts the experiment
p.addParameter('triggerkeyname' , '''5'' or ''t''', @isstring)  % for display only
p.addParameter('offsetpix'      , [0 0]     , @isnumeric);      % offset of screen in pixels [10 20] means move 10-px right, 20-px down
p.addParameter('movieflip'      , [0 0]     , @isnumeric)       % whether to flip up-down, whether to flip left-right
p.addParameter('debugmode'      , false     , @islogical)       % whether to use debug mode (no BOLDscreen, no eyelink)
p.addParameter('storeparams'   , true      , @islogical)       % whether to store stimulus params

% Parse inputs
p.parse(subjID, sesID, runnum, varargin{:});

% Rename variables into general params struct
params    = struct();
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    eval([sprintf('params.%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p

%% %%%%%%%%%%%%% PERIPHERALS PARAMS %%%%%%%%%%%%%

if params.debugmode
    dispName = 'KKOFFICEQ3277';
    disp = vcd_getDisplayParams('KKOFFICEQ3277');
    wanteyetracking = false;
    Screen('Preference', 'SkipSyncTests', 1);
else
    disp = vcd_getDisplayParams(); % BOLDSCREEN is default
end

% Nova1x32 with BOLDscreen and big eye mirrors:
ptonparams = {[disp.w_pix disp.h_pix disp.refresh_hz 24],[],-2};  % {[1920 1080 120 24],[],-2} BOLDSCREEN (squaring CLUT because we need to simulate normal monitors)

% Buttonbox / keyboard
params.ignorekeys = KbName({params.triggerkey});  % dont record TR triggers as subject responses

% SETUP RNG
rand('twister', sum(100*clock));
params.rand = rand;

%% %%%%%%%%%%%%% STIM PARAMS %%%%%%%%%%%%%

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

if params.loadparams
    d = dir(fullfile(params.infofolder,'stim*.mat'));
    load(fullfile(d(end).folder,d(end).name),'stim');
    params.stim = stim; clear stim;
    
    d = dir(fullfile(params.infofolder,'trials*.mat'));
    load(fullfile(d(end).folder,d(end).name),'all_trials');
    params.trials = all_trials; clear all_trials;
    
    params.exp = vcd_getSessionParams;
    
else
    params.stim   = vcd_getStimParams('all',disp.name,true);
    params.exp    = vcd_getSessionParams;
    params.trials = vcd_makeTrials(params);
end


%% %%%%%%%%%%%%% LOAD STIMULI %%%%%%%%%%%%%

if ~exist('images','var') || isempty(images)

    images = struct('gabor',[],'rdk',[],'dot',[],'cobj',[], 'ns',[]);

    % GABORS: 6D array: [x,y,orient,contrast,phase,delta]
    load(params.stim.gabor.stimfile, 'gabors'); 
    images.gabor = gabors; clear gabors;
    
    % RDKs: 8 directions x 3 coherence levels cell array: [x,y,3, frames]
    load(params.stim.rdk.stimfile, 'rdk'); 
    images.rdk = rdk; clear rdk;
    
    % Simple dot: 2D array: [x,y]
    load(params.stim.dot.stimfile, 'simple_dot'); 
    images.dot = simple_dot; clear simple_dot;
    
    % Complex objects: 4D array: [x,y,object,rotation]
    load(params.stim.cobj.stimfile, 'objects'); 
    images.cobj = objects; clear objects;
    
    % NS: 5D array: [x,y,superordinate cat, ns_loc, obj_loc]
    load(params.stim.ns.stimfile, 'scenes'); 
    images.ns = scenes; clear scenes;

%%%%%%%%%%%%% EK TO DO %%%%%%%%%%%%%
%       % load maskimages/????
%       load(stimfile,'maskimages');
%       if ~exist('maskimages','var')
%           maskimages = {};
%       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      % resize if desired
      for fn = 1:length(params.exp.stimClassLabels)
      
          if ~isempty(params.stim.(params.exp.stimClassLabels{fn}).dres) && ...
                  length(params.stim.(params.exp.stimClassLabels{fn}).dres)==2
          tic;
          fprintf('resampling the stimuli; this may take a while');
            % GRAY
          if strcmp(params.exp.stimClassLabels{fn},'gabor') || ...
                  strcmp(params.exp.stimClassLabels{fn},'rdk') || ...
                  strcmp(params.exp.stimClassLabels{fn},'dot') || ...
                  strcmp(params.exp.stimClassLabels{fn},'cobj')

              tmp_im = reshape(params.stim.(stimClassLabels{fn}), ...
                  size(params.stim.(stimClassLabels{fn}),1),... x
                  size(params.stim.(stimClassLabels{fn}),2),... y
                  []);
              iscolor = false;
              % COLOR
          elseif strcmp(params.exp.stimClassLabels{fn},'ns')
              tmp_im = reshape(params.stim.(stimClassLabels{fn}), ...
                  size(params.stim.(stimClassLabels{fn}),1),... x
                  size(params.stim.(stimClassLabels{fn}),2),... y
                  3, ...
                  []);
              iscolor = true;
          end
          
          szIm = size(tmp_im); numIm = szIm(end);
          images.(stimClassLabels{fn}) = cell(1,numIm);
          for p = 1:numIm
              statusdots(p,numIm);
              if iscolor
                  temp = imresize(tmp_im(:,:,:,p),params.stim.(stimClassLabels{fn}).dres);
              else
                  temp = imresize(tmp_im(:,:,p),params.stim.(stimClassLabels{fn}).dres);
              end        
              images.(stimClassLabels{fn}){p} = temp;
          end
          
          %%%%%%%%%%%%% EK TO DO %%%%%%%%%%%%%
%           for p = 1:length(maskimages)
%               temp = cast([],class(maskimages{p}));
%               for q=1:size(maskimages{p},3)
%                   temp(:,:,q) = imresize(maskimages{p}(:,:,q),dres);
%               end
%               maskimages{p} = temp;
%           end
          fprintf('done!\n');
          toc
          end
      end

end

% Get a vector with number of images in each class
for fn = 1:length(params.exp.stimClassLabels)
    sz = size(images.(params.exp.stimClassLabels{fn}));
    if params.stim.(params.exp.stimClassLabels{fn}).iscolor
        numinclass(fn) = prod(sz(4:end));
    else
        numinclass(fn) = prod(sz(3:end));
    end
end

%% %%%%%%%%%%%%% DEFINE MINIBLOCK ORDER %%%%%%%%%%%%%

if params.loadparams
    d = dir(fullfile(params.infofolder,'subject_sessions_*.mat'));
    load(fullfile(params.infofolder,d(end).name),'subject_sessions'); % pick the last one we saved
    
    subj_session = subject_sessions(params.sesID,params.subjID).run(params.runnum,:);
    
else
    subject_sessions = vcd_getStimParams('all',dispName,true);
    subj_session = subject_sessions(params.sesID,params.subjID).run(params.runnum,:); 
end

% INDEX IMAGES THAT NEED TO BE LOADED
im_order = vcd_getRunImagesToLoad(subj_session,images,params);

%% %%%%%%%%%%%%% EYELINK PARAMS %%%%%%%%%%%%%

% EYE FUN for SYNC TIME
if params.wanteyetracking
    tfun = @() cat(2,fprintf('STIMULUS STARTED.\n'),Eyelink('Message','SYNCTIME'));
    if ~isfield(params, 'eyelinkfile') || isempty(params.eyelinkfile) 
        eyelinkfile = fullfile(vcd_rootPath,'data',sprintf('eye_%s_vcd_subj-%s_run-%d.edf',datestr(now,'yyyymmdd-HHMMSS'),subjID,runnum));
    else
        eyelinkfile = params.eyelinkfile;
    end
else
    tfun = @() fprintf('STIMULUS STARTED.\n');  % this is executed when the expt starts
end

%% %%%%%%%%%%%%% FIXATION PARAMS %%%%%%%%%%%%%

% Fixation order and fixation
soafun = @() round(params.stim.fix.dotmeanchange*params.stim.fps + params.stim.fix.dotchangeplusminus*(2*(rand-.5))*params.stim.fps);




%% %%%%%%%%%%%%% TASK INSTRUCTIONS %%%%%%%%%%%%%

%%%%%%%%%%%%% Q: add folder ?? --> instrtextfolder

setupscript  = 'showinstructionscreen(''runvcdcore_presubjectinstructions.txt'',250,75,25,78)';  % inputs are (fileIn,txtOffset,instTextWrap,textsize,background)

for ii = 1:length(params.exp.stimTaskLabels)
    if regexp('*scc*',params.exp.stimTaskLabels{ii})
        script_stimTask(ii) = sprintf('runvcdcore_scc.txt');
    else
        script_stimTask(ii) = sprintf('runvcdcore_%s.txt',params.exp.stimTaskLabels{ii});
    end
end





%% %%%%%%%%%%%%% INIT SCREEN %%%%%%%%%%%%%
oldclut = pton(ptonparams{:});

%%%%%%%%%%%%% Q: Do we still need this? %%%%%%%%%%%%%
% [windowPtr,center,blankColor] = doScreen; % Opens screen, hides cursor, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flag incase we want to quit early
getoutearly = 0; 

%% %%%%%%%%%%%%% EYELINK: initialize, setup, calibrate, and start %%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%% TODO: EK FIX %%%%%%%%%%%%%
    
%     % recenter EL coordinates if needed
%     if ~iszero(offset_pix) || isempty(offset_pix)
%         centerXY = [round(wwidth/2), round(wheight/2)] + offset_pix;
%         
%     end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    vcd_showStimulus(p, images, im_order, setupscript, soafun);




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
