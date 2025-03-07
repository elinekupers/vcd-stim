function [data,getoutearly] = vcd_showStimulus(win, rect, params, ...
    scan, ...
    timing, ...
    introscript, ...
    taskscript, ...
    tfunEYE, ...
    deviceNr, ...
    oldCLUT,...
    oldPriority)

getoutearly    = 0;
glitchcnt      = 0;
when           = 0;
framecnt       = 0;
wantframefiles = false;
detectinput    = true;

%% Preallocate space for key presses and timestamps
timekeys       = {};
digitrecord    = [];
digitframe     = [];
digitpolarity  = [];

%% PREPARE IMAGES

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

allowforceglitch     = 0; % 0 means do nothing special. [1 D] means allow keyboard input 'p' to force a glitch of duration D secs.
frameorder           = 1:size(timing.trig_stim,1);

% init variables, routines, constants
timeframes = NaN(1, floor(size(frameorder,2)-1)+1);

% make tmpdir
if wantframefiles
    tmpDir = '~/Desktop/tmp/'; %#ok<UNRCH>
    framefiles = {'~/Desktop/tmp/frame%05d.png', []};
    if ~exist('tmpDir','dir'), mkdir(tmpDir); end
end

%% ptb stuff

scan.ifi       = Screen('GetFlipInterval',win);
mfi            = scan.ifi;
frameduration  = round(mfi)/params.stim.fps; % 30 Hz presentation, 2 frames for office/psph monitors (60 Hz), 4 frames for BOLDscreen (120 Hz);

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);
Screen('TextSize', win, 25);
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

%% Make textures
[im_tex, im_rect, txt_tex, txt_rect, framecolor] = ...
    vcd_makeTextures(params, win, rect, flipfun, scan, timing, taskscript);

% Get pre-run intructions 
[instrtext, prerun_text_rect] = vcd_getInstructionText(params, introscript, rect);
% Draw background + dot
Screen('DrawTextures',win, im_tex{1},[], im_rect{1}, 0, [], 1, framecolor{1});...

% Draw text
DrawFormattedText(win, instrtext, 'center', (prerun_text_rect(4)/2)-250,0,75,[],[],[],[],prerun_text_rect);
Screen('Flip',win);

fprintf('Instructions are on screen, waiting for trgger...\n');


%% CORE EXP CODE!

fprintf('press trigger key to begin the movie. (consider turning off network, energy saver, software updates.)\n');
% safemode = 0;
while 1
    [~,keyCode,~] = KbWait(deviceNr,2); % previously deviceNr = -3; outputs: secs,keyCode,deltaSecs
    temp = KbName(keyCode);

    if isempty(params.triggerkey) || any(strcmp(temp(1), params.triggerkey))
        break;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% log the start!

feval(tfunEYE) % Eyelink('Message','SYNCTIME'));
timekeys = [timekeys; {GetSecs 'trigger'}];

%% DRAW THE TEXTURES

for frame = 1:size(frameorder,2)+1
    
    framecnt = framecnt + 1;
    frame0 = floor(frame);
    reporttext = '';
    
    blockID = timing.trig_block(framecnt);
    
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
    
    % draw text
    DrawFormattedText(win, txt_tex{frame}, 'center', (txt_rect{frame}(4)/2)-50,0,75,[],[],[],[],txt_rect{frame});
    
    % draw textures
    Screen('DrawTextures',win,im_tex{frame,:},[],im_rect{frame,:},0,[],1,framecolor{frame});
    
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
    
%     % write to file if desired
%     if wantframefiles
%         if isempty(framefiles{2}) %#ok<UNRCH>
%             imwrite(Screen('GetImage',win),sprintf(framefiles{1},framecnt));
%         else
%             imwrite(uint8(placematrix(zeros([framefiles{2} 3]),Screen('GetImage',win))),sprintf(framefiles{1},framecnt));
%         end
%     end
    
    % update when
    if didglitch
        % if there were glitches, proceed from our earlier when time.
        % set the when time to 9/10 a frame before the desired frame.
        % notice that the accuracy of the mfi is strongly assumed here.
        whendesired = whendesired + mfi * frameduration;
        when = whendesired - mfi * (9/10); %#ok<*NASGU>
    else
        % if there were no glitches, just proceed from the last recorded time
        % and set the when time to 9/10 a frame before the desired time.
        % notice that the accuracy of the mfi is only weakly assumed here,
        % since we keep resetting to the empirical VBLTimestamp.
        whendesired = VBLTimestamp + mfi * frameduration;
        when = whendesired - mfi * (9/10);  % should we be less aggressive??
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PT CLEANUP STUFF
Screen('Close','all');
ptoff(oldCLUT);

% restore priority and cursor
Priority(oldPriority);
ShowCursor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUTTON LOGGING CLEAN UP STUFF
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
data.timeKeys               = timekeys;
data.timing.glitchcnt       = glitchcnt;
data.timing.timeframes      = timeframes;
data.timing.starttime       = starttime;
data.timing.endtime         = timeframes(end);
data.timing.empiricalrundur = dur;
data.timing.empiricalfps    = length(timeframes)/dur;
data.digitrecord            = digitrecord;


end




