function [data,getoutearly,run_frames,subj_run_table,scan] = vcd_showStimulus(...
    win, rect, params, ...
    scan, ...
    bckground, ...
    fix, ...
    eye_im, ...
    stim, ...
    run_frames, ...
    subj_run_table, ...
    introscript, ...
    taskscript, ...
    tfunEYE, ...
    deviceNr, ...
    oldCLUT,...
    oldPriority)

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
    if ~exist('tmpDir','dir'), mkdir(tmpDir); end
end

%% ptb stuff

mfi            = Screen('GetFlipInterval',win);
frameduration  = round(params.stim.framedur_s/mfi); % 30 Hz presentation, 2 frames for office/psph monitors (60 Hz), 4 frames for BOLDscreen (120 Hz);

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);
Screen('TextSize', win, params.disp.fontsize);
Screen('TextStyle', win, 0);

Priority(9);

if ~params.debugmode
    HideCursor;
else
    ShowCursor;
end

% run functions as first time running them always takes more time
GetSecs;
now;
ceil(1);
fprintf('');


%% Create background and fixation textures prior to exp onset (as we need them throughout the experiment)
bckground_rect    = rect; %CenterRect([0 0 round(size(bckground,1)) round(size(bckground,2))],rect);
bckrgound_texture = Screen('MakeTexture', win, bckground);

% make fixation dot texture
fix_texture_thin_full   = cell(1,size(fix.fix_thin_full,2));
fix_texture_thick_full  = cell(1,size(fix.fix_thin_full,2));
fix_texture_thick_left  = cell(1,size(fix.fix_thin_full,2));
fix_texture_thick_right = cell(1,size(fix.fix_thin_full,2));
fix_texture_thick_both  = cell(1,size(fix.fix_thin_full,2));
for ll = 1:size(fix.fix_thin_full,2) % loop over luminance values
    fix_texture_thin_full{ll}   = Screen('MakeTexture',win,fix.fix_thin_full{ll});
    fix_texture_thick_full{ll}  = Screen('MakeTexture',win,fix.fix_thick_full{ll});
    fix_texture_thick_left{ll}  = Screen('MakeTexture',win,fix.fix_thick_left{ll});
    fix_texture_thick_right{ll} = Screen('MakeTexture',win,fix.fix_thick_right{ll});
    fix_texture_thick_both{ll}  = Screen('MakeTexture',win,fix.fix_thick_both{ll});
end

%% create eyetracking targets
et_rect = rect; et_texture = {};
for sac = 1:size(eye_im.sac_im,4)
    et_texture{sac} = Screen('MakeTexture',win,eye_im.sac_im(:,:,:,sac));
end
et_texture{size(eye_im.sac_im,4)+1} = Screen('MakeTexture',win,eye_im.pupil_im_white);
et_texture{size(eye_im.sac_im,4)+2} = Screen('MakeTexture',win,eye_im.pupil_im_black);

%% Prepare background and fixation texture vector outside the flip loop

fix_tex    = cell(length(stim.im),1);
fix_rect   = fix_tex;
im_tex     = fix_tex;
im_rect    = fix_tex;
framecolor = fix_tex;
txt_tex    = fix_tex;
txt_rect   = fix_tex;

% run_frames.frame_event_nr = run_frames.frame_event_nr(1:end-1);
% run_frames.is_catch = run_frames.is_catch(1:end-1);
% run_frames.crossingIDs = run_frames.crossingIDs(1:end-1);

for nn = 1:size(run_frames.frame_event_nr,1)
    
    eventID = run_frames.frame_event_nr(nn);
    if isnan(eventID)
        eventID = 0;
    end
    
    % set up fixation dot textures
    lum_idx = find(run_frames.fix_abs_lum(nn)==params.stim.fix.dotlum);
    
    if eventID == 0 || isnan(run_frames.is_cued(nn)) || (run_frames.is_cued(nn)==0)
        fix_tex(nn)  = fix_texture_thin_full(lum_idx);
        fix_rect(nn) = fix.fix_thin_rect;   
    else
        if eventID==95
            if run_frames.is_cued(nn)==1
                fix_tex(nn)  = fix_texture_thick_left(lum_idx);
                fix_rect(nn) = fix.fix_thick_rect;
            elseif run_frames.is_cued(nn)==2
                fix_tex(nn)  = fix_texture_thick_right(lum_idx);
                fix_rect(nn) = fix.fix_thick_rect;
            elseif run_frames.is_cued(nn)==3
                fix_tex(nn)  = fix_texture_thick_both(lum_idx);
                fix_rect(nn) = fix.fix_thick_rect;
            end
        elseif ismember(eventID,[90,91,92,93,94,96,97])
            fix_tex(nn)  = fix_texture_thick_full(lum_idx);
            fix_rect(nn) = fix.fix_thick_rect;
        elseif ismember(eventID,[98,99])
            fix_tex(nn)  = fix_texture_thin_full(lum_idx);
            fix_rect(nn) = fix.fix_thin_rect;
        elseif (eventID > 0) || (eventID < 90)
            fix_tex(nn) = fix_texture_thick_full(lum_idx);
            fix_rect(nn) = fix.fix_thick_rect;
        elseif eventID > 990 % eyetracking target (TODO: implement actual targets)
            fix_tex(nn) = fix_texture_thin_full(lum_idx);
            fix_rect(nn) = fix.fix_thin_rect;
        end
    end
    
    switch eventID
        
%     exp.block.stim_epoch1_ID        = 91; % generic stim ID
%     exp.block.stim_epoch2_ID        = 92; % generic stim ID
%     exp.block.response_ID           = 93; % Time for subject to respond
%     exp.block.trial_start_ID        = 94; % Fixation dot thickening
%     exp.block.spatial_cue_ID        = 95; % Fixation dot turning black on L/R/both sides
%     exp.block.delay_ID              = 96; % Delay period between two stimulus epochs
%     exp.block.task_cue_ID           = 97; % Text on display to instruct subject
%     exp.block.ITI_ID                = 98; % Inter-trial interval
%     exp.block.IBI_ID                = 99; % Inter-block interval
        
        % Draw background + fix dot on top
        case {0, 90, 93, 94, 95, 96, 98, 99}
            
            % DrawTextures
            % * TexturePointers  need to be: n vector (where n is the number of textures)
            % * DestinationRects need to be: 4 row x n columns (where n is the number of textures)
            
            im_tex{nn}  = cat(1, bckrgound_texture, fix_tex(nn));
            im_rect{nn} = cat(1, bckground_rect, fix_rect{nn});
            framecolor{nn} = 255*ones(2,3); % <framecolor> can also be size(<frameorder>,2) x 1 with values in [0,1] indicating an alpha change.
            
        case 97 % task_cue_ID
            
            script = taskscript{~cellfun(@isempty, regexp(taskscript,sprintf('%02d',run_frames.crossingIDs(nn)),'match'))};
            [task_instr, task_rect] = vcd_getInstructionText(params, script, rect);
            
            im_tex{nn}  = cat(1, bckrgound_texture, fix_tex(nn));
            im_rect{nn} = cat(1, bckground_rect, fix_rect{nn});
            
            txt_tex{nn} = task_instr;
            txt_rect{nn} = task_rect;
            
            framecolor{nn} = 255*ones(2,3); % <framecolor> can also be size(<frameorder>,2) x 1 with values in [0,1] indicating an alpha change.
        
        % Draw background with eyetracking target
        % eye_gaze_fix_ID         = 990; % fixation target
        % eye_gaze_sac_target_ID  = 991:995; % central, left, right, up, down.
        % eye_gaze_pupil_ID       = 996; % white then black

        case {990,991}
            im_tex{nn}  = et_texture{1};
            im_rect{nn} = et_rext;
            framecolor{nn} = 255*ones(1,3);
            
        case 992
            im_tex{nn}  = et_texture{2};
            im_rect{nn} = et_rext;
            framecolor{nn} = 255*ones(1,3);
            
        case 993
            im_tex{nn}  = et_texture{3};
            im_rect{nn} = et_rext;
            framecolor{nn} = 255*ones(1,3);
            
        case 994
            im_tex{nn}  = et_texture{4};
            im_rect{nn} = et_rext;
            framecolor{nn} = 255*ones(1,3);
            
        case 995
            im_tex{nn}  = et_texture{5};
            im_rect{nn} = et_rext;
            framecolor{nn} = 255*ones(1,3);
            
        case 996
            im_tex{nn}  = et_texture{6};
            im_rect{nn} = et_rext;
            framecolor{nn} = 255*ones(1,3);
                        
        case 997
            im_tex{nn}  = et_texture{7};
            im_rect{nn} = et_rext;
            framecolor{nn} = 255*ones(1,3);
            
        case {91, 92}
            if run_frames.is_catch(nn) % treat catch trials as delays
                im_tex{nn} = cat(1, bckrgound_texture, fix_tex{nn});
                im_rect{nn} = cat(1, bckground_rect, fix_rect{nn});
                framecolor{nn} = 255*ones(2,3);
            end
      
    end
end

clear frame;

% Get pre-run intructions
[instrtext, prerun_text_rect] = vcd_getInstructionText(params, introscript, rect);

% Draw background + dot
% windowPointer, texturePointer(s), [sourceRect], destRects, rotAngles, filterModes, globalAlphas, modulateColors, textureShader, specialFlags, auxParameters]);
Screen('DrawTextures',win, im_tex{1},[], im_rect{1}',[], [0;0], [1;1], framecolor{1}');
    
% Draw text
DrawFormattedText(win, instrtext, 'center', (prerun_text_rect(4)/2)-50, 0, 75,[],[],[],[],prerun_text_rect);
Screen('Flip',win);

fprintf('Instructions are on screen, waiting for trigger...\n');


%% CORE EXP CODE!

fprintf('Press trigger key to begin the movie. (consider turning off network, energy saver, software updates.)\n');
% safemode = 0;
while 1
    [~,keyCode,~] = KbWait(deviceNr,2); % previously deviceNr = -3; outputs: secs,keyCode,deltaSecs
    temp = KbName(keyCode);
    
    if isempty(params.triggerkey) || any(strcmp(temp(1), params.triggerkey))
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% log the start!

feval(tfunEYE); %fprintf('EXP START'); % SEND Eyelink('Message','SYNCTIME'));
timekeys = [timekeys; {GetSecs 'trigger'}];

%% DRAW THE TEXTURES
framecnt = 0;
for frame = 1:size(frameorder,2)+1 % we add 1 to log end
    
    framecnt = framecnt +1;
    frame0 = floor(framecnt);
    
    % we have to wait until the last frame of the run sequence is done.
    if frame0 == size(frameorder,2)+1        
        while 1
            if GetSecs >= whendesired
                getoutearly = 1;
                % Log the end
                if ~isempty(tfunEYE)
                    feval(tfunEYE);
                    timekeys = [timekeys; {GetSecs 'DONE'}];
                end
                break;
            end
        end
    end
    
    % get out early?
    if getoutearly
        break;
    end
    
    switch run_frames.frame_event_nr(framecnt)
        
        % 0  : pre/post blank
        % 93 : exp_session.block.response_ID
        % 94 : exp_session.block.trial_start_ID
        % 95 : exp_session.block.spatial_cue_ID
        % 96 : exp_session.block.delay_ID
        % 97 : exp_session.block.task_cue_ID
        % 98 : exp_session.block.ITI_ID
        % 99 : exp_session.block.IBI_ID
        
        % Draw background + thin or thick fix dot on top
        case {0, 90, 93, 94, 95, 96, 98, 99, 990, 991, 992, 993, 994, 995, 996, 997}
            % draw background and dot textures
            Screen('DrawTextures',win,cell2mat(im_tex{frame}),[],im_rect{frame}',[0;0],[],[1;1],framecolor{frame}');
            
        case 97 % task_cue_ID
            
            % draw background and left/right cuing dot textures
            Screen('DrawTextures',win, cell2mat(im_tex{frame}),[],im_rect{frame}',[0;0],[],[1;1],framecolor{frame}');
            
            % draw text
            % inputs are winptr, tstring, sx, sy, color, wrapat, flipHorizontal, flipVertical, vSpacing, righttoleft, winRect)
            DrawFormattedText(win, txt_tex{frame}, 'center', (txt_rect{frame}(4)/2)-25,0,75,[],[],[],[],txt_rect{frame});
            
        case {91,92} % stim IDs
            % Draw stimulus textures
            Screen('DrawTexture',win, bckrgound_texture,[], bckground_rect, 0, [], 1, 255*ones(1,3));...
        
            % im_w_mask is a cell with dims: frames x 1, where each cell has 1 or 2 sides (1:l, 2:r)
            for side = 1:length(find(~cellfun(@isempty, stim.im(run_frames.im_IDs(frame),:))))
                stim_texture = Screen('MakeTexture',win, stim.im{run_frames.im_IDs(frame),side});
                Screen('DrawTexture',win,stim_texture,[], stim.rects{run_frames.im_IDs(frame),side}, 0,[],1, 255*ones(1,3));
            end

            % Draw fix dot on top
            Screen('DrawTexture',win,fix_tex{frame},[], fix_rect{frame}, 0,[],1, 255*ones(1,3));
            Screen('Close',stim_texture);
    end

    % give hint to PT that we're done drawing
    Screen('DrawingFinished',win);
    
    %%%%%%%%%%%%%%%%%%%%%%%% the main while loop that actually puts up stimuli and records button presses
    
    % here, deal with making the stimulus frame / texture / stuff
    % read input until we have to do the flip
    while 1
        
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
            
            % otherwise, try to read input
        else
            if detectinput
                [keyIsDown,secs,keyCode,~] = KbCheck(deviceNr);  % previously -3 listen to all devices
                if keyIsDown
                    
                    % get the name of the key and record it
                    kn = KbName(keyCode);
                    timekeys = [timekeys; {secs kn}];
                    
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
        end
     end   
        
     
     % write to file if desired
     if wantframefiles
         if isempty(framefiles{2}) %#ok<UNRCH>
             imwrite(Screen('GetImage',win),sprintf(framefiles{1},framestart));
         else
             imwrite(uint8(placematrix(zeros([framefiles{2} 3]),Screen('GetImage',win))),sprintf(framefiles{1},framestart));
         end
     end
     
     % update when
     if didglitch
         % if there were glitches, proceed from our earlier when time.
         % set the when time to 9/10 a frame before the desired frame.
         % notice that the accuracy of the mfi is strongly assumed here.
         whendesired = whendesired + mfi * frameduration;
         when = whendesired - mfi * (9/10); %#ok<*NASGU>    % RZ code: when = (when + mfi / 2) + mfi * frameduration - mfi / 2;
     else
         % if there were no glitches, just proceed from the last recorded time
         % and set the when time to 9/10 a frame before the desired time.
         % notice that the accuracy of the mfi is only weakly assumed here,
         % since we keep resetting to the empirical VBLTimestamp.
         whendesired = VBLTimestamp + mfi * frameduration;
         when = whendesired - mfi * (9/10);  % should we be less aggressive??    % RZ code: %     when = VBLTimestamp + mfi * frameduration - mfi / 2;  % should we be less aggressive??
         
     end
     
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUTTON LOGGING CLEAN UP STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
dur = (timeframes(end)-timeframes(1)) * (length(timeframes)/(length(timeframes)-1));
fprintf('projected total run duration: %.10f\n',dur);
fprintf('frames per second: %.10f\n',length(timeframes)/dur);

% prepare output
digitrecord = {digitrecord digitframe digitpolarity};

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

if getoutearly
    % Save a quick version of the data in case something fails...
    vars = whos;
    vars = {vars.name};
    vars = vars(cellfun(@(x) ~isequal(x,'stim'),vars));
    save(fullfile(vcd_rootPath,sprintf('tmp_data_%s.mat',datestr(now,30))),vars{:});
end
 

if ~getoutearly
   performance = vcd_getBehavioralPerformance(params, data, correct_response);
   [instrtext, txt_rect, params] = vcd_getMotivationText(params, performance.hitrate);

   % draw background and central fixation dot
   Screen('DrawTextures',win, cell2mat(im_tex{1}),[],im_rect{1}',[0;0],[],[1;1],framecolor{1}');
   
   % draw text
   % inputs are winptr, tstring, sx, sy, color, wrapat, flipHorizontal, flipVertical, vSpacing, righttoleft, winRect)
   DrawFormattedText(win, instrtext, 'center', (txt_rect{frame}(4)/2)-25,0,75,[],[],[],[],txt_rect{1});
   
   waitSecs(4);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PT CLEANUP STUFF
Screen('Close',win);
ptoff();

% restore priority and cursor
Priority(oldPriority);
ShowCursor;




return




