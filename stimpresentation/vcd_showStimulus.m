function [data,getoutearly,run_frames,run_table] = vcd_showStimulus(...
    params, ...
    ptonparams, ...
    fix, ...
    stim, ...
    run_frames, ...
    run_table, ...
    introscript, ...
    taskscript, ...
    deviceNr)

%% PTB functionality notes:
% Inputs to DrawTexture: 
%   windowPointer, texturePointer(s), [sourceRect], destRects, rotAngles, filterModes, globalAlphas, modulateColors, textureShader, specialFlags, auxParameters]);
% Inputs to MakeTextures/DrawTextures (plural) instead of DrawTexture (single): 
%   * texturePointers need to be: n vector (where n is the number of textures)
%   * destinationRects need to be: 4 row x n columns (where n is the number of textures)
%   * Make sure to adjust have n x alpha and n x color inputs.
%  example: Screen('DrawTextures',win, stim_textures', [], catcell(1,stim.rects(run_frames.im_IDs(framecnt,1),:))', [0;0],[], [1;1], 255*ones(2,3)');
% Inputs to DrawFormattedText inputs are winptr, tstring, [sx], [sy], [color], [wrapat], [flipHorizontal], [flipVertical], [vSpacing], [righttoleft], [winRect]
%
%  <framecolor> can also be size(<frameorder>,2) x 1 with values in [0,1] indicating an alpha change.
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%  SET CONSTANTS / FLAGS / COUNTERS  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% internal constants
fliplead = 10/1000;  % min amount of time to allocate prior to flip

% Set counters
glitchcnt         = 0;
when              = 0;

% Set logical flags
getoutearly       = false;  % we are not getting out early (yet). Will change to true when experimenter presses ESCAPE button
allowforceglitch  = false;  % 0 means do nothing special. [1 D] means allow keyboard input 'p' to force a glitch of duration D secs.
forceglitch       = false;  % useful for testing timing
wantframefiles    = false;  % do we want to print every single frame flipped on the screen (to create a video of the experiment).
detectinput       = true;   % do we want detect the button presses prior to the experiment starts (to test if button box is working)

% Preallocate space for key presses and timestamps
timekeys          = {};
frameorder        = 1:size(stim.im,1);
timeframes        = NaN(1, floor(size(frameorder,2)-1)+1);
stim_texture      = [];
stim_textures     = [];

% make tmpdir
if wantframefiles
    tmpDir     = '~/Desktop/tmp/'; %#ok<UNRCH>
    framefiles = {'~/Desktop/tmp/frame%05d.png', []};
    if ~exist(tmpDir,'dir'), mkdir(tmpDir); end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   INIT SCREEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Screen('Preference', 'SyncTestSettings', params.ptbMaxVBLstd); % params.ptbMaxVBLstd defines what deviation from empirical monitor refresh rate do we allow before calling it a missed flip (and throwing an error)
oldCLUT     = pton(ptonparams{:});
win         = firstel(Screen('Windows'));
oldPriority = Priority(MaxPriority(win));
rect        = Screen('Rect',win); % get total screen rect   % alternatively: rect = CenterRect(round([0 0 rect(3)*winsize rect(4)*winsize]),rect);

% OTHER PTB stuff
mfi            = Screen('GetFlipInterval',win);     % empirical refresh rate determined when we "opened" the ptb window pointer
frameduration  = round(params.stim.framedur_s/mfi); % 60 Hz presentation, 1 time frame for office/psph monitors (60 Hz), 2 time frames for BOLDscreen (120 Hz);

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);
Screen('TextSize', win, params.disp.fontsize);
Screen('TextStyle', win, 0);

% Copied from PTB: Priority value corresponds to the proportion of the
% video frame period  during which MATLAB is guranteed 100 percent of CPU
% time: GuaranteedCPUtimePerFramePeriod = priority/10 * frame period.
% If your computer has multiple video displays running at different frame 
% rates then Priority will choose the shortest frame period.
Priority(9); % 9 = max value possible for Mac OS X. 
HideCursor; % now hide cursors for all screen windows

% IMPORTANT to ensure proper functioning (flush caches, etc.)
clear PsychHID;
clear KbCheck;
clear KbWait;

% Run functions as first time running them always takes more time
GetSecs;
now;
ceil(1);
fprintf('');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% CREATE TEXTURES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Create fixation textures prior to exp onset (as we need them throughout the experiment)
% We make fixation dot textures for each luminance value and rim thickness.
fix_texture.thin_full         = cell(1,size(fix.fix_thin_full,2));
fix_texture.thick_full_white  = cell(1,size(fix.fix_thin_full,2));
fix_texture.thick_left        = cell(1,size(fix.fix_thin_full,2));
fix_texture.thick_right       = cell(1,size(fix.fix_thin_full,2));
fix_texture.thick_both        = cell(1,size(fix.fix_thin_full,2));
for ll = 1:size(fix.fix_thin_full,2) % loop over luminance values
    fix_texture.thin_full{ll}         = Screen('MakeTexture',win,fix.fix_thin_full{ll});
    fix_texture.thick_full_white{ll}  = Screen('MakeTexture',win,fix.fix_thick_full_white{ll});
    fix_texture.thick_left{ll}        = Screen('MakeTexture',win,fix.fix_thick_left{ll});
    fix_texture.thick_right{ll}       = Screen('MakeTexture',win,fix.fix_thick_right{ll});
    fix_texture.thick_both{ll}        = Screen('MakeTexture',win,fix.fix_thick_both{ll});
end

fix_texture.thick_full_black{1}  = Screen('MakeTexture',win,fix.fix_thick_full_black{1});


% Prepare fixation texture vector outside the flip loop
fix_tex    = cell(length(stim.im),1);
fix_rect   = fix_tex;
im_tex     = fix_tex;
im_rect    = fix_tex;
framecolor = fix_tex;
task_tex   = fix_tex;
task_rect  = fix_tex;

unique_crossingIDs = unique(run_frames.crossingIDs(~isnan(run_frames.crossingIDs)));
unique_crossingIDs(ismember(unique_crossingIDs,[999,0])) = [];


% Make textures for task instruction images once
alltasktex = {};
for p = 1:length(unique_crossingIDs)
    alltasktex{p} = Screen('MakeTexture', win, taskscript.im{p});
end

% Make textures for eyetracking block images once
allimtex = {};
for p = 1:size(stim.eye,4)
    allimtex{p} = Screen('MakeTexture', win, stim.eye(:,:,:,p));
end

clear p

% Loop over time frames
for nn = 1:size(run_frames.frame_event_nr,1)
    
    % Get event number for this time frame
    eventID = run_frames.frame_event_nr(nn);
    if isnan(eventID)
        eventID = 0;
    end
    
    % set up fixation dot textures
    lum_idx = find(run_frames.fix_abs_lum(nn)==params.stim.fix.dotlum);
    
    switch eventID
    case 92
        if run_frames.is_cued(nn)==1
            fix_tex(nn)  = fix_texture.thick_left(lum_idx);
            fix_rect(nn) = fix.fix_thick_rect;
        elseif run_frames.is_cued(nn)==2
            fix_tex(nn)  = fix_texture.thick_right(lum_idx);
            fix_rect(nn) = fix.fix_thick_rect;
        elseif run_frames.is_cued(nn)==3
            fix_tex(nn)  = fix_texture.thick_both(lum_idx);
            fix_rect(nn) = fix.fix_thick_rect;
        end
    case 90 % TASK CUE -- NOTE: no fixation circle!
%         fix_tex(nn)  = []; % NOTE: if left empty, row nn will shrink fix_tex length by 1. WE DON'T WANT THAT as we want the length to match the run duration (hence it's commented out).
%         fix_rect(nn) = [];
%         run_frames.fix_abs_lum(nn) = NaN; % remove 128 as absolute luminance because there is no fixation target twinkle
    case {91 93 94 95 96 97 98} % ALL STIMULUS EVENTS + ITI (thick fixation circle rim)
        fix_tex(nn)  = fix_texture.thick_full_white(lum_idx);
        fix_rect(nn) = fix.fix_thick_rect;
    case 99 % IBI (thin fixation circle rim)
        fix_tex(nn)  = fix_texture.thin_full(lum_idx);
        fix_rect(nn) = fix.fix_thin_rect;
    otherwise
        if eventID >= 990 % eyetracking targets / pupil black/white displays
%           fix_tex(nn) = []; % NOTE: if left empty, row nn will shrink fix_tex length by 1. WE DON'T WANT THAT as we want the length to match the run duration (hence it's commented out).
%           fix_rect(nn) = []; 
%           run_frames.fix_abs_lum(nn) = NaN; % remove 128 as absolute luminance because there is no fixation target twinkle
        elseif eventID == 0 % pre/post blank rest period (thin fixation circle rim)
            fix_tex(nn)  = fix_texture.thin_full(lum_idx);
            fix_rect(nn) = fix.fix_thin_rect;
        end
    end
    
    switch eventID
    % task_cue_ID           = 90; % Text on display to instruct subject
    % post_task_cue_ITI_ID  = 91; % time between task cue and first trial of the block (thick white fixation circle rim)
    % spatial_cue_ID        = 92; % fixation dot turning red on either L/R/both sides  
    % pre_stim_blank_ID     = 93; % blank period between spatial cue and stimulus onset
    % stim_epoch1_ID        = 94; % stim onset (1st interval)
    % stim_epoch2_ID        = 95; % stim onset (2nd interval after delay)
    % delay_ID              = 96; % delay period between two stimulus epochs
    % response_ID           = 97; % time for subject to respond
    % ITI_ID                = 98; % inter-trial interval
    % IBI_ID                = 99; % inter-block interval

        % Draw background + fix dot on top
        case {0, 91, 92, 93, 94, 95, 96, 97, 98, 99}
            im_tex{nn}     = fix_tex{nn}; 
            im_rect{nn}    = fix_rect{nn};
            framecolor{nn} = 255*ones(1,3); 

        case 90 % task_cue_ID
            % Get instructions from png file
            task_tex{nn}   = alltasktex{run_frames.crossingIDs(nn)==unique_crossingIDs};
            task_rect{nn}  = taskscript.rect{run_frames.crossingIDs(nn)==unique_crossingIDs};
            framecolor{nn} = 255*ones(1,3);
        
        % Draw background with eyetracking target
        case {990,991} % eye_gaze_fix_ID = 990,991; % central fixation "rest" and "target"
            im_tex{nn}      = allimtex{1};
            im_rect{nn}     = rect;
            framecolor{nn}  = 255*ones(1,3);
            
        case 992 % eye_gaze_sac_target_ID = left
            im_tex{nn}      = allimtex{2};
            im_rect{nn}     = rect;
            framecolor{nn}  = 255*ones(1,3);
            
        case 993 % eye_gaze_sac_target_ID = right
            im_tex{nn}      = allimtex{3};
            im_rect{nn}     = rect;
            framecolor{nn}  = 255*ones(1,3);
            
        case 994 % eye_gaze_sac_target_ID = up
            im_tex{nn}     = allimtex{4};
            im_rect{nn}    = rect;
            framecolor{nn} = 255*ones(1,3);
            
        case 995 % eye_gaze_sac_target_ID = down
            im_tex{nn}     = allimtex{5};
            im_rect{nn}    = rect;
            framecolor{nn} = 255*ones(1,3);
            
        case 996 % eye_gaze_pupil_ID is black
            im_tex{nn}     = allimtex{6};
            im_rect{nn}    = rect;
            framecolor{nn} = 255*ones(1,3);
                         
        case 997 % eye_gaze_pupil_ID is white
            im_tex{nn}     = allimtex{7};
            im_rect{nn}    = rect;
            framecolor{nn} = 255*ones(1,3);
    end
end

% Use thick rim fixation circle for pre-run instruction screen.
fix_tex_preruninstr  = fix_texture.thick_full_black{1}; % mid gray luminance (128)
fix_rect_preruninstr = fix.fix_thick_rect{1};


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  EYELINK SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~params.wanteyetracking
    % No need for EYE fun
    tfunEYE            = @() fprintf('\n');
else
    % ANON EYE FUN for SYNC TIME
    tfunEYE     = @() Eyelink('Message','SYNCTIME');

    % Initialize Eyelink
    et_ok = EyelinkInit;
    assert(et_ok==1);

    % Get Eyelink default params
    el = EyelinkInitDefaults(win);
    if ~isempty(el.callback)
        PsychEyelinkDispatchCallback(el);
    end
    
    % Update default Eyelink params with VCD needs
    el = vcd_setEyelinkParams(el);
    EyelinkUpdateDefaults(el);
    
    % Get window size
    wwidth = rect(3); wheight = rect(4);  % returns in pixels
    
    % Ensure window size is what we think it is
    assert(isequal(wwidth,params.disp.w_pix))
    assert(isequal(wheight,params.disp.h_pix))
    assert(isequal(round(wwidth/2),params.disp.xc))
    assert(isequal(round(wheight/2),params.disp.yc))
    
    % Tell the user
    fprintf('Pixel size of window is width: %d, height: %d.\n',wwidth,wheight);
    fprintf('Pixel center offset of window is [x,y]=[%d,%d].\n',params.offsetpix(1),params.offsetpix(2));
    
    % Recenter EL coordinates if needed
    % EK: should we just use updated params.stim.xc/yc??
    if any(params.offsetpix~=[0,0]) || isempty(params.offsetpix)
        xc_off = round(wwidth/2) + params.offsetpix(1);
        yc_off = round(wheight/2) + params.offsetpix(2);
    else
        xc_off = round(wwidth/2);
        yc_off = round(wheight/2);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Customize calibration points (include pixel shift, change size and color)
    Eyelink('command',sprintf('screen_phys_coords = %3.1f, %3.1f, %3.1f, %3.1f', ...
        params.disp.el_monitor_size(1),params.disp.el_monitor_size(2),params.disp.el_monitor_size(3),params.disp.el_monitor_size(4)));  % monitor size in millimeters (center to left, top, right, and bottom). PProom: [-260.0, 162.5, 260.0, -162.5]
    Eyelink('command','screen_distance = %ld %ld', params.disp.el_screen_distance(1),params.disp.el_screen_distance(2)); % distance in millimeters from eye to top and bottom edge of the monitor. PProom: [1003 1003]
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld',0,0,wwidth-1,wheight-1); % X,Y coordinates left/top/right/bottom of display area
    Eyelink('message','DISPLAY_COORDS %ld %ld %ld %ld',0,0,wwidth-1,wheight-1);
    % IF YOU DON'T WANT CUSTOM DOT POSITIONS: Set number of calibration/validation dots and spread: horizontal-only(H) or horizontal-vertical(HV) as H3, HV3, HV5, HV9 or HV13
    Eyelink('command','calibration_type = HV5'); % horizontal-vertical 5-points.
    Eyelink('command','generate_default_targets = NO');
    Eyelink('command','calibration_samples  = 6');
    Eyelink('command','calibration_sequence = 0,1,2,3,4,5');
    Eyelink('command','calibration_targets  = %d,%d %d,%d %d,%d %d,%d %d,%d',...
        xc_off,yc_off,  ... center x,y
        xc_off + params.stim.el.point2point_distance_pix, yc_off, ... horz shift right
        xc_off - params.stim.el.point2point_distance_pix, yc_off, ... horz shift left
        xc_off, yc_off + params.stim.el.point2point_distance_pix, ... vert shift down
        xc_off, yc_off - params.stim.el.point2point_distance_pix); %  vert shift up
    Eyelink('command','validation_samples = 6');
    Eyelink('command','validation_sequence = 0,1,2,3,4,5');
    Eyelink('command','validation_targets  = %d,%d %d,%d %d,%d %d,%d %d,%d',...
        xc_off,yc_off,  ... center x,y
        xc_off + params.stim.el.point2point_distance_pix, yc_off, ... horz shift right
        xc_off - params.stim.el.point2point_distance_pix, yc_off, ... horz shift left
        xc_off, yc_off + params.stim.el.point2point_distance_pix, ... vert shift down
        xc_off, yc_off - params.stim.el.point2point_distance_pix); %  vert shift up
    
    Eyelink('command','active_eye = LEFT');
    Eyelink('command','binocular_enabled','NO');
    Eyelink('command','enable_automatic_calibration','NO'); % force manual calibration sequencing, if yes, provide Eyelink('command','automatic_calibration_pacing=1500');
    Eyelink('command','recording_parse_type = GAZE'); %from manual (default)
    Eyelink('command','sample_rate = %d', 1000); % hz
    Eyelink('command','driftcorrect_cr_disable = YES'); % yes to disable drift correction -- we don't want that!
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
    
    % make temp file name and open
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
    elseif isequal(checkcalib,'0')
        fprintf('Will SKIP calibration.\n');
    else
        fprintf('Did you press the wrong button?? Let''s try again.\n');
        checkcalib = input('Do you want to do a calibration (0=no, 1=yes)? ','s');
        if isequal(checkcalib,'1')
            fprintf('Please perform calibration. When done, press the output/record button.\n');
            EyelinkDoTrackerSetup(el);
        elseif isequal(checkcalib,'0')
            fprintf('Will SKIP calibration.\n');
        end
    end
    
    % Start recording
    Eyelink('StartRecording');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% LOAD PRE-RUN INSTRUCTION SCREEN %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Draw background (gray screen) 
Screen('FillRect',win,params.stim.bckgrnd_grayval,rect); % Previously: Screen('DrawTexture',win, bckrgound_texture,[], bckground_rect,[], 0, 1, 255*ones(1,3));

% Make and draw pre-run intructions from png file + fixation circle
intro_tex  = Screen('MakeTexture', win, introscript.im);
Screen('DrawTexture',win,intro_tex,[], introscript.rect, 0, [], 1, 255*ones(1,3)); % draw intro text
Screen('DrawTexture',win,fix_tex_preruninstr,[], fix_rect_preruninstr, 0, [], 1, 255*ones(1,3)); % draw thin fix circle

% Show the subject the intro screen (pre-trigger)
Screen('Flip',win);

% Close intro texture
Screen('Close',intro_tex);

fprintf('Instructions are on screen, waiting for trigger...\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORE EXP CODE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Please have participant press a button to confirm they are ready.\n');
ListenChar(2);

while 1
  [keyIsDown,~,keyCode,~] = KbCheck(deviceNr);
  if keyIsDown
    temp = KbName(keyCode);
    if any(strcmp(temp(1), params.userkeys))
      fprintf(' ** PARTICIPANT BUTTON SUCCESSFULLY DETECTED! ** \n')
      break;
    end
  end
end

% Once button is detected, we draw center saccade target circle on mean luminance gay background
pre_et_tex = Screen('MakeTexture',win, stim.eye(:,:,:,1));
Screen('DrawTexture',win,pre_et_tex,[], rect, 0, [], 1, 255*ones(1,3)); % draw thick black fix circle
Screen('Flip',win);
clear pre_et_tex

fprintf('Press trigger key to begin the movie. (consider turning off network, energy saver, software updates, psychopy.)\n');

while 1
    [~,keyCode,~] = KbWait(deviceNr,2); % deviceNr = -3 --> listen to all devices; outputs: secs,keyCode,deltaSecs
    temp = KbName(keyCode);
    
    if isempty(params.triggerkey) || any(strcmp(temp(1), params.triggerkey))
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% LOG THE START     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

feval(tfunEYE);
timekeys = [timekeys; {GetSecs 'trigger'}];
fprintf('EXP START.\n'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%% DRAW THE TEXTURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
framecnt = 0;
prevlinecnt = NaN;
prevreportcnt = NaN;
while 1
    
    framecnt = framecnt +1;    % NOTE: framecnt seems to track frame
    frame0 = floor(framecnt);
    
    % we have to wait until the last frame of the run sequence is done.
    if frame0 >= size(frameorder,2)+1        % Note that it is >= instead of == because it is possible that glitches push us too far
        while 1
            if GetSecs >= whendesired
                getoutearly = 1;
                % Log the end in eyelink
                if ~isempty(tfunEYE) % SEND Eyelink('Message','SYNCTIME'));
                    feval(tfunEYE);
                end
                    timekeys = [timekeys; {GetSecs 'DONE'}]; %#ok<AGROW>
                break;
            end
        end
    end
    
    % get out early?
    if getoutearly
        break;
    end
    
    % Get gray background
    Screen('FillRect',win,params.stim.bckgrnd_grayval,rect);
    
    % Here, we deal with making the stimulus texture / drawing text
    switch run_frames.frame_event_nr(framecnt)
        
        % 0  : pre/post blank
        % 90 : task_cue       
        % 91 : post_task_cue_ITI
        % 92 : spatial_cue
        % 93 : pre_stim_blank
        % 94 : stim_epoch1
        % 95 : stim_epoch2
        % 96 : delay
        % 97 : response
        % 98 : ITI
        % 99 : IBI 
        
        % Draw ET targets on grey/white/black background
        case {990, 991, 992, 993, 994, 995, 996, 997}
            Screen('DrawTexture',win,im_tex{framecnt},[],im_rect{framecnt},0,[],1,framecolor{framecnt});

        % Draw fix circle (thin or thick)
        case {0, 91, 92, 93, 96, 97, 98, 99}
            Screen('DrawTexture',win,im_tex{framecnt},[],im_rect{framecnt},0,[],1,framecolor{framecnt});
        
        % Draw task instruction text (we do not plot a fixation circle such
        % that subjects are encouraged to read the task instructions)
        case 90 
            Screen('DrawTexture',win,task_tex{framecnt},[], task_rect{framecnt}, 0, [], 1, framecolor{framecnt});
            
        case {94, 95} % stim IDs
            % stim.im is a cell with dims: frames x 2, where each cell has a uint8 image (1:l, 2:r)
            sides = length(run_frames.im_IDs(framecnt,~isnan(run_frames.im_IDs(framecnt,:))));
            if  sides == 1 % Make and draw one stimulus texture
                stim_texture = Screen('MakeTexture',win, stim.im{run_frames.im_IDs(framecnt,sides),sides});
                Screen('DrawTexture',win,stim_texture,[], stim.rects{run_frames.im_IDs(framecnt,sides),sides}, 0,[],1, 255*ones(1,3));
                Screen('Close',stim_texture); stim_texture = [];
            elseif sides == 2  % Make and draw two stimulus textures
                stim_textures(1) = Screen('MakeTexture',win, stim.im{run_frames.im_IDs(framecnt,1),1});
                stim_textures(2) = Screen('MakeTexture',win, stim.im{run_frames.im_IDs(framecnt,2),2});
                Screen('DrawTextures',win,stim_textures',[], catcell(1,stim.rects(run_frames.im_IDs(framecnt,1),:))', [0;0],[], [1;1], 255*ones(2,3)');
                Screen('Close',stim_textures); stim_textures = [];
            end

            % Draw fix dot on top
            Screen('DrawTexture',win,fix_tex{framecnt},[], fix_rect{framecnt}, 0,[],1, 255*ones(1,3));
    end

    % Give hint to PT that we're done drawing
    Screen('DrawingFinished',win);
    
    %%%%%%%%%%%%%%%%%%%%%%%% the main while loop that actually puts up stimuli and records button presses
    
    % Read input until we have to do the flip
    while 1

        % try to read input (instantaneous)
        if detectinput
            [keyIsDown,secs,keyCode,~] = KbCheck(deviceNr);  % previously -3 listen to all devices
            if floor(secs/5) ~= prevlinecnt  % every 5 seconds, fprintf a new line
              fprintf('$\n');
              prevlinecnt = floor(secs/5);
            end
            if keyIsDown
                
                % get the name of the key and record it
                kn = KbName(keyCode);
                timekeys = [timekeys; {secs kn}]; %#ok<AGROW>
                if floor(secs/0.1) ~= prevreportcnt
                    if iscell(kn)
                        for kk=1:length(kn)
                            fprintf('%s',kn{kk}(1));
                        end
                    else                   
                     fprintf('%s',kn(1));
                    end
                    prevreportcnt = floor(secs/0.1);
                end
                
                % check if ESCAPE was pressed
                if isequal(kn,'ESCAPE')
                    fprintf('Escape key detected.  Exiting prematurely.\n');
                    getoutearly = 1;
                    break;
                end
                
                % force a glitch?
                if allowforceglitch(1) && isequal(kn,'p')
                    WaitSecs(allowforceglitch(2));
                end
                
            end
        end
        
        % if we are in the initial case OR if we have hit the when time, then display the frame
        if when == 0 || GetSecs >= when
            
            % issue the flip command and record the empirical time
            [VBLTimestamp,~,~,~,~] = Screen('Flip',win,  0);
            timeframes(framecnt) = VBLTimestamp;
            
            % get matlab now for the very first stimulus frame
            if framecnt==1
                absnowtime = now;
            end
            
            % Detect glitch
            if when ~= 0 && (VBLTimestamp - whendesired) > (mfi * (1/2))
                glitchcnt = glitchcnt + 1;
                didglitch = 1;
            else
                didglitch = 0;
            end
            
            % get out of this loop
            break;

        end
            
     end   
        
     % write frame to file if desired
     if wantframefiles
         if isempty(framefiles{2}) %#ok<UNRCH>
             imwrite(Screen('GetImage',win),sprintf(framefiles{1},framestart));
         else
             imwrite(uint8(placematrix(zeros([framefiles{2} 3]),Screen('GetImage',win))),sprintf(framefiles{1},framestart));
         end
     end
    
     % calc
     fliplead0 = min(fliplead,(9/10)*mfi);
     
     % update when
     if didglitch
         % if there were glitches, proceed from our earlier when time.
         % set the when time to a little bit before the desired frame.
         % notice that the accuracy of the mfi is strongly assumed here.
         whendesired = whendesired + mfi * frameduration;
         when = whendesired - fliplead0; %#ok<*NASGU>  

         % if the current time is already past whendesired, we are doomed,
         % and so we have to drop a frame. here, we do it repeatedly until
         % we are in the clear. this gives us at least a chance of getting
         % back on track, but it is NOT guaranteed. ultimately, the user needs
         % do some checking of NaNs in timeframes to check for dropped frames.
         while GetSecs >= whendesired - fliplead0
           framecnt = framecnt + 1;
           whendesired = whendesired + mfi * frameduration;
           when = whendesired - fliplead0;
        end
         
     else
         % if there were no glitches, just proceed from the last recorded time
         % and set the when time to a little bit before the desired time.
         % notice that the accuracy of the mfi is only weakly assumed here,
         % since we keep resetting to the empirical VBLTimestamp.
         whendesired = VBLTimestamp + mfi * frameduration;
         when = whendesired - fliplead0;  % should we be less aggressive??
     end

end

fprintf('RUN ENDED.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% SAVE EYELINK DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if params.wanteyetracking
    if ~isempty(eyetempfile)
        % Close eyelink and record end
        Eyelink('StopRecording');
        Eyelink('message', sprintf('EXP END %d',GetSecs));
        Eyelink('CloseFile');
        status = Eyelink('ReceiveFile',eyetempfile, params.savedatafolder, 1);
        fprintf('ReceiveFile status %d\n', status);
        
        % Rename temporary EDF file to final file name
        mycmd=['mv ' params.savedatafolder '/' eyetempfile ' ' params.savedatafolder '/' params.eyelinkfile];
        system(mycmd);
        if status <= 0, fprintf('\n\nProblem receiving edf file\n\n');
        else
            fprintf('Data file ''%s'' can be found in ''%s''\n', params.eyelinkfile, pwd);
        end
        Eyelink('ShutDown'); % Here "ShutDown" means close the TCP/IP link, not actually shutting down the OS of the eyelink host machine
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% BUTTON LOGGING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adjust the times in timeframes and timekeys to be relative to the first time recorded.
% thus, time==0 corresponds to the showing of the first frame.
starttime = timeframes(1);
timeframes = timeframes - starttime;
if size(timekeys,1) > 0
    timekeys(:,1) = cellfun(@(x) x - starttime,timekeys(:,1),'UniformOutput',0);
end
timekeys = [{absnowtime 'absolutetimefor0'}; timekeys];

% NOT NECESSARY BECAUSE vcdbehavioralanalysis.m deals with this kind of stuff.
% % report basic timing information to stdout
% fprintf('we had %d glitches!\n',glitchcnt);
% fprintf('we had %d dropped frames!\n',sum(isnan(timeframes)));
% dur = (timeframes(end)-timeframes(1)) * (length(timeframes)/(length(timeframes)-1));
% fprintf('projected total run duration: %.10f\n',dur);
% fprintf('frames per second: %.10f\n',length(timeframes)/dur);


% Add button presses and monitor timing to data struct
data = struct();
data.wantframefiles         = wantframefiles;
data.detectinput            = detectinput;
data.forceglitch            = forceglitch;
data.timeKeys               = timekeys;
data.timing.mfi             = mfi;
data.timing.glitchcnt       = glitchcnt;
data.timing.timeframes      = timeframes;
data.timing.starttime       = starttime;
%data.timing.endtime         = dur;
%data.timing.empiricalfps    = length(timeframes)/dur;
data.timing.frameduration   = frameduration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% STORING DATA    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear memory from stuff we don't need
clear fix_tex fix_rect fix_texture alltasktex allimtex ...
      task_tex task_rect stim_textures nn lum_idx ll kn ...
      intro_tex im_tex im_rect framecolor fix_rect_preruninstr ...
      fix_tex_preruninstr sides temp

% remove all_images from params struct in case that sneaked in
params.all_images = [];
  
% figure out names of all variables except uint8 images
vars = whos;
vars = {vars.name};
vars = vars(cellfun(@(x) ~ismember(x,{'fix', ...
    'stim', 'all_images', ...
    'wantframefiles' 'detectinput' 'forceglitch' 'timekeys' 'mfi' 'glitchcnt' ...
    'timeframes' 'starttime' 'dur' 'frameduration'}),vars)); 

% Save data (button presses, params, etc)
save(fullfile(params.savedatafolder,params.behaviorfile),vars{:});  % '-v7.3'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORMANCE CHECKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get behavioral performance
performance = vcdbehavioralanalysis(fullfile(params.savedatafolder,params.behaviorfile));

% Get feedback display text
[fb_txt, fbtext_rect] = vcd_getFeedbackDisplay(params, rect, performance);

% Show performance to subject (draw text on gray background)
Screen('FillRect', win, params.stim.bckgrnd_grayval, rect);
DrawFormattedText(win, fb_txt, (fbtext_rect(3)/2)-350, (fbtext_rect(4)/2)-100,0,150,[],[],[],[],fbtext_rect); % inputs are winptr, tstring, sx, sy, color, wrapat, flipHorizontal, flipVertical, vSpacing, righttoleft, winRect)
Screen('Flip',win,0);

% Check monitor timing
ptviewmoviecheck(data.timing.timeframes,data.timeKeys,[],{'5' 't'});

% Check eyetracking data
if params.wanteyetracking
    eyeresults = vcdeyetrackingpreprocessing( ...
        fullfile(params.savedatafolder,params.eyelinkfile), ...
        fullfile(params.savedatafolder,params.behaviorfile), performance);
end

% Let the user decide when to end the experiment
fprintf(' *** PRESS SPACE BAR TO END RUN. *** \n');
while 1
    [~,keyCode,~] = KbWait(deviceNr,2); % deviceNr = -3; KbWait outputs are secs,keyCode,deltaSecs
    kn = KbName(keyCode);
    if isequal(kn,'space')
        fprintf('SPACE key detected.  Ending run.\n');
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% PT CLEANUP STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Restore priority and cursor
ListenChar(0);
ShowCursor;
ptoff(oldCLUT);
Priority(oldPriority);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



