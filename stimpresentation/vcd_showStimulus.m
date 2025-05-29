function [data,getoutearly,run_frames,run_table] = vcd_showStimulus(...
    win, rect, params, ...
    fix, ...
    eye_im, ...
    stim, ...
    run_frames, ...
    run_table, ...
    introscript, ...
    taskscript, ...
    tfunEYE, ...
    deviceNr, ...
    oldCLUT,...
    oldPriority, ...
    eyetempfile)

%% internal constants
fliplead = 5/1000;  % min amount of time to allocate prior to flip

%% Set flags and counters
getoutearly    = 0;
glitchcnt      = 0;
when           = 0;
forceglitch    = false;
wantframefiles = false;
detectinput    = true;

%% Preallocate space for key presses and timestamps
timekeys       = {};
digitrecord    = [];
digitframe     = [];
digitpolarity  = [];

%% PREPARE IMAGES
allowforceglitch  = 0; % 0 means do nothing special. [1 D] means allow keyboard input 'p' to force a glitch of duration D secs.
frameorder        = 1:size(stim.im,1);

% init variables, routines, constants
timeframes = NaN(1, floor(size(frameorder,2)-1)+1);

% make tmpdir
if wantframefiles
    tmpDir = '~/Desktop/tmp/'; %#ok<UNRCH>
    framefiles = {'~/Desktop/tmp/frame%05d.png', []};
    if ~exist(tmpDir,'dir'), mkdir(tmpDir); end
end

%% PTB stuff
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

% Run functions as first time running them always takes more time
GetSecs;
now;
ceil(1);
fprintf('');

% IMPORTANT to ensure proper functioning (flush caches, etc.)
clear PsychHID;
clear KbCheck;

%% Create background and fixation textures prior to exp onset (as we need them throughout the experiment)
bckground_rect    = rect; %CenterRect([0 0 round(size(bckground,1)) round(size(bckground,2))],rect);
% bckrgound_texture = Screen('MakeTexture', win, bckground); % no more background image

% make fixation dot texture
fix_texture.thin_full   = cell(1,size(fix.fix_thin_full,2));
fix_texture.thick_full  = cell(1,size(fix.fix_thin_full,2));
fix_texture.thick_left  = cell(1,size(fix.fix_thin_full,2));
fix_texture.thick_right = cell(1,size(fix.fix_thin_full,2));
fix_texture.thick_both  = cell(1,size(fix.fix_thin_full,2));
for ll = 1:size(fix.fix_thin_full,2) % loop over luminance values
    fix_texture.thin_full{ll}   = Screen('MakeTexture',win,fix.fix_thin_full{ll});
    fix_texture.thick_full{ll}  = Screen('MakeTexture',win,fix.fix_thick_full{ll});
    fix_texture.thick_left{ll}  = Screen('MakeTexture',win,fix.fix_thick_left{ll});
    fix_texture.thick_right{ll} = Screen('MakeTexture',win,fix.fix_thick_right{ll});
    fix_texture.thick_both{ll}  = Screen('MakeTexture',win,fix.fix_thick_both{ll});
end

%% Prepare background and fixation texture vector outside the flip loop
fix_tex    = cell(length(stim.im),1);
fix_rect   = fix_tex;
im_tex     = fix_tex;
im_rect    = fix_tex;
framecolor = fix_tex;
txt_tex    = fix_tex;
txt_rect   = fix_tex;

for nn = 1:size(run_frames.frame_event_nr,1)
    
    eventID = run_frames.frame_event_nr(nn);
    if isnan(eventID)
        eventID = 0;
    end
    
    % set up fixation dot textures
    lum_idx = find(run_frames.fix_abs_lum(nn)==params.stim.fix.dotlum);
    
    if eventID==92 % SPATIAL CUE
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
    elseif ismember(eventID,90) % TASK CUE -- NOTE: no fixation circle!
        fix_tex(nn)  = [];
        fix_rect(nn) = [];
    elseif ismember(eventID,[91,93,94,95,96,97,98]) % ALL STIMULUS EVENTS + ITI (thick fixation circle rim)
        fix_tex(nn)  = fix_texture.thick_full(lum_idx);
        fix_rect(nn) = fix.fix_thick_rect;
    elseif ismember(eventID,99) % IBI (thin fixation circle rim)
        fix_tex(nn)  = fix_texture.thin_full(lum_idx);
        fix_rect(nn) = fix.fix_thin_rect;
    elseif eventID >= 990 % eyetracking targets / pupil black/white displays
        fix_tex(nn) = [];
        fix_rect(nn) = []; 
    elseif eventID == 0 % pre/post blank rest period (thin fixation circle rim)
        fix_tex(nn)  = fix_texture.thin_full(lum_idx);
        fix_rect(nn) = fix.fix_thin_rect;
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

    % If we want to use DrawTextures (plural) (instead of DrawTexture (single)):
    % * TexturePointers need to be: n vector (where n is the number of textures)
    % * DestinationRects need to be: 4 row x n columns (where n is the number of textures)

        % Draw background + fix dot on top
        case {0, 91, 92, 93, 94, 95, 96, 97, 98, 99}
            
            im_tex{nn}     = fix_tex{nn}; %cat(1, bckrgound_texture, fix_tex(nn));
            im_rect{nn}    = fix_rect{nn}; %cat(1, bckground_rect, fix_rect{nn});
            framecolor{nn} = 255*ones(1,3); %255*ones(2,3); % <framecolor> can also be size(<frameorder>,2) x 1 with values in [0,1] indicating an alpha change.

        case 90 % task_cue_ID
            % Get instructions from text file
            script = taskscript{~cellfun(@isempty, ...
                regexp(taskscript,sprintf('%02d',run_frames.crossingIDs(nn)),'match'))};
            [task_instr, task_rect] = vcd_getInstructionText(params, script, rect);
            
%             im_tex{nn}  = bckrgound_texture;
%             im_rect{nn} = bckground_rect;
            txt_tex{nn}  = task_instr;
            txt_rect{nn} = task_rect;
            framecolor{nn} = 255*ones(1,3);
        
        % Draw background with eyetracking target
        case {990,991} % eye_gaze_fix_ID = 990,991; % central fixation "rest" and "target"
            im_tex{nn}  = Screen('MakeTexture',win,eye_im.sac_im(:,:,:,1));
            im_rect{nn} = rect;
            framecolor{nn} = 255*ones(1,3);
            
        case 992 % eye_gaze_sac_target_ID  = left
            im_tex{nn}  = Screen('MakeTexture',win,eye_im.sac_im(:,:,:,2));
            im_rect{nn} = rect;
            framecolor{nn} = 255*ones(1,3);
            
        case 993 % eye_gaze_sac_target_ID  = right
            im_tex{nn}  = Screen('MakeTexture',win,eye_im.sac_im(:,:,:,3));
            im_rect{nn} = rect;
            framecolor{nn} = 255*ones(1,3);
            
        case 994 % eye_gaze_sac_target_ID  = up
            im_tex{nn}  = Screen('MakeTexture',win,eye_im.sac_im(:,:,:,4));
            im_rect{nn} = rect;
            framecolor{nn} = 255*ones(1,3);
            
        case 995 % eye_gaze_sac_target_ID  = down
            im_tex{nn}  = Screen('MakeTexture',win,eye_im.sac_im(:,:,:,5));
            im_rect{nn} = rect;
            framecolor{nn} = 255*ones(1,3);
            
        case 996 % eye_gaze_pupil_ID is black
            im_tex{nn}  = Screen('MakeTexture',win,eye_im.pupil_im_black);
            im_rect{nn} = rect;
            framecolor{nn} = 255*ones(1,3);
                         
        case 997 % eye_gaze_pupil_ID is white
            im_tex{nn}  = Screen('MakeTexture',win,eye_im.pupil_im_white);
            im_rect{nn} = rect;
            framecolor{nn} = 255*ones(1,3);
    end
end

% Get pre-run intructions
[instrtext, prerun_text_rect] = vcd_getInstructionText(params, introscript, rect);

% Draw background (gray screen)
% Screen('DrawTexture',win, bckrgound_texture,[], bckground_rect,[], 0, 1, 255*ones(1,3)); % INPUTS TO DRAWTEXTURE: windowPointer, texturePointer(s), [sourceRect], destRects, rotAngles, filterModes, globalAlphas, modulateColors, textureShader, specialFlags, auxParameters]);
Screen('FillRect',win,params.stim.bckgrnd_grayval,rect);
  
% Draw intro text (pre-trigger)
DrawFormattedText(win, instrtext, 'center', (prerun_text_rect(4)/2)-50, 0, 75,[],[],[],[],prerun_text_rect);
Screen('Flip',win);

fprintf('Instructions are on screen, waiting for trigger...\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORE EXP CODE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Press trigger key to begin the movie. (consider turning off network, energy saver, software updates.)\n');

while 1
    [~,keyCode,~] = KbWait(deviceNr,2); % previously deviceNr = -3; outputs: secs,keyCode,deltaSecs
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

        % Draw background + fix circle (thin or thick) on top
        case {0, 91, 92, 93, 96, 97, 98, 99}
            Screen('DrawTexture',win,im_tex{framecnt},[],im_rect{framecnt},0,[],1,framecolor{framecnt});
            % draw background and dot textures
%             Screen('DrawTextures',win,cell2mat(im_tex{framecnt}),[],im_rect{framecnt}',[0;0],[],[1;1],framecolor{framecnt}');
            
        case 90 % task_cue_ID
            
            % draw background and background textures -- We do not plot a
            % fixation circle such that subjects are encouraged to read the
            % tak instructions
%             Screen('DrawTexture',win, im_tex{framecnt},[],im_rect{framecnt},0,[],1,framecolor{framecnt});
            
            % draw text
            % inputs are winptr, tstring, sx, sy, color, wrapat, flipHorizontal, flipVertical, vSpacing, righttoleft, winRect)
            DrawFormattedText(win, txt_tex{framecnt}, 'center', (txt_rect{framecnt}(4)/2)-25,0,75,[],[],[],[],txt_rect{framecnt});
            
        case {94, 95} % stim IDs
            % Draw stimulus textures
%             Screen('DrawTexture',win, bckrgound_texture,[], bckground_rect, 0, [], 1, 255*ones(1,3));...
        
            % stim.im is a cell with dims: frames x 2, where each cell has a uint8 image (1:l, 2:r)
            for side = 1:length(run_frames.im_IDs(framecnt,~isnan(run_frames.im_IDs(framecnt,:))))
                stim_texture = Screen('MakeTexture',win, stim.im{run_frames.im_IDs(framecnt,side),side});
                Screen('DrawTexture',win,stim_texture,[], stim.rects{run_frames.im_IDs(framecnt,side),side}, 0,[],1, 255*ones(1,3));
                Screen('Close',stim_texture);
            end

            % Draw fix dot on top
            Screen('DrawTexture',win,fix_tex{framecnt},[], fix_rect{framecnt}, 0,[],1, 255*ones(1,3));
    end

    % give hint to PT that we're done drawing
    Screen('DrawingFinished',win);
    
    %%%%%%%%%%%%%%%%%%%%%%%% the main while loop that actually puts up stimuli and records button presses
    
    % Read input until we have to do the flip
    while 1

        % try to read input (instantaneous)
        if detectinput
            [keyIsDown,secs,keyCode,~] = KbCheck(deviceNr);  % previously -3 listen to all devices
            if keyIsDown
                
                % get the name of the key and record it
                kn = KbName(keyCode);
                timekeys = [timekeys; {secs kn}]; %#ok<AGROW>
                
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
        status = Eyelink('ReceiveFile',eyetempfile, params.savedatadir, 1);
        fprintf('ReceiveFile status %d\n', status);
        
        % Rename temporary EDF file to final file name
        mycmd=['mv ' params.savedatadir '/' eyetempfile ' ' params.savedatadir '/' params.eyelinkfile];
        system(mycmd);
        if status <= 0, fprintf('\n\nProblem receiving edf file\n\n');
        else
            fprintf('Data file ''%s'' can be found in ''%s''\n', params.eyelinkfile, pwd);
        end
        Eyelink('ShutDown'); % Here "ShutDown" means close the TCP/IP link, not actually shutting down the OS on the eyelink host machine
        
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

% report basic timing information to stdout
fprintf('we had %d glitches!\n',glitchcnt);
fprintf('we had %d dropped frames!\n',sum(isnan(timeframes)));
dur = (timeframes(end)-timeframes(1)) * (length(timeframes)/(length(timeframes)-1));
fprintf('projected total run duration: %.10f\n',dur);
fprintf('frames per second: %.10f\n',length(timeframes)/dur);

% prepare output
digitrecord = {digitrecord digitframe digitpolarity};

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
data.timing.endtime         = dur;
data.timing.empiricalfps    = length(timeframes)/dur;
data.timing.frameduration   = frameduration;
data.digitrecord            = digitrecord;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% STORING DATA    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure out names of all variables except uint8 images
vars = whos;
vars = {vars.name};
vars = vars(cellfun(@(x) ~ismember(x,{'fix_im','bckground','stim','eye_im' ...
  'wantframefiles' 'detectinput' 'forceglitch' 'timekeys' 'mfi' 'glitchcnt' ...
  'timeframes' 'starttime' 'dur' 'frameduration' 'digitrecord'}),vars)); 

% Save data (button presses, params, etc)
save(fullfile(params.savedatadir,params.behaviorfile),vars{:});  % '-v7.3'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORMANCE CHECKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get behavioral performance
performance = vcdbehavioralanalysis(fullfile(params.savedatadir,params.behaviorfile));

% Get feedback display text
[fb_txt, fbtext_rect] = vcd_getFeedbackDisplay(params, rect, performance, taskscript);

% Show performance to subject (draw text on gray background)
Screen('FillRect', win, params.stim.bckgrnd_grayval, rect);
DrawFormattedText(win, fb_txt, (fbtext_rect(3)/2)-120, (fbtext_rect(4)/2)-25,0,75,[],[],[],[],fbtext_rect); % inputs are winptr, tstring, sx, sy, color, wrapat, flipHorizontal, flipVertical, vSpacing, righttoleft, winRect)
Screen('Flip',win,0);

% Check monitor timing
ptviewmoviecheck(data.timing.timeframes,data.timeKeys,[],{'5' 't'});

% Check eyetracking data
if params.wanteyetracking
    eyeresults = vcdeyetrackingpreprocessing( ...
        fullfile(params.savedatadir,params.eyelinkfile), ...
        fullfile(params.savedatadir,params.behaviorfile), performance);
end

% Let the user decide when to end the experiment
fprintf('PRESS SPACE BAR TO END RUN.\n');
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
ShowCursor;
ptoff(oldCLUT);
Priority(oldPriority);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



