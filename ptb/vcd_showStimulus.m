function [data,getoutearly] = vcd_showStimulus(params, ...
    scan, ...
    timing, ...
    introscript, ...
    taskscript, ...
    tfunEYE)

getoutearly    = 0;
glitchcnt      = 0;
frameduration  = params.stim.fps;
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
    flipfun = @(x) flipdim(flipdim(x,1),2);
elseif params.movieflip(1)
    flipfun = @(x) flipdim(x,1);
elseif params.movieflip(2)
    flipfun = @(x) flipdim(x,2);
else
    flipfun = @(x) x;
end

% ptviewmovie params
% <framecolor> (optional) is size(<frameorder>,2) x 3 with values in [0,255].
%   each row indicates a multiplication on color channels to apply for a
%   given frame.  default is 255*ones(size(<frameorder>,2),3) which means
%   to multiply each channel by 1 (i.e. do nothing special).  note that
%   when an entry in <frameorder> is 0, <framecolor> has no effect for that entry.
%   <framecolor> can also be size(<frameorder>,2) x 1 with values in [0,1]
%   indicating an alpha change.

params.stim.scfactor = 1; % positive number with the scaling to apply to the images
filtermode           = choose(params.stim.scfactor==1,0,1); % size(<frameorder>,2) x 3 with values in [0,255]. each row indicates a multiplication on color channels to apply for a
                                                            % given frame.
framecolor           = 255*ones(1,3); % <framecolor> can also be size(<frameorder>,2) x 1 with values in [0,1] indicating an alpha change.
allowforceglitch     = 0; % 0 means do nothing special. [1 D] means allow keyboard input 'p' to force a glitch of duration D secs.
frameorder           = 1:size(timing.trig_stim,1);

% init variables, routines, constants
timeframes = repmat(NaN,[1 floor(size(frameorder,2)-1)+1]);




%% ptb stuff
win  = firstel(Screen('Windows'));
oldPriority = Priority(MaxPriority(win));

rect = Screen('Rect',win); % what is the total rect
% rect = CenterRect(round([0 0 rect(3)*winsize rect(4)*winsize]),rect);
[win, rect] = Screen('OpenWindow',max(Screen('Screens')),params.stim.bckgrnd_grayval,rect);
scan.ifi = Screen('GetFlipInterval',win);
Screen('Preference', 'SyncTestSettings', .0004);
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);
Priority(9);
mfi = Screen('GetFlipInterval',win);  % re-use what was found upon initialization!
mfi = 1/round(1/mfi);

if ~params.debugmode
    HideCursor;
    deviceNr = params.deviceNr(1);
else
    deviceNr = params.deviceNr(2);
    ShowCursor;
end

%% Create background and fixation textures prior to exp onset (as we need them throughout the experiment)
bckground_rect = CenterRect([0 0 round(size(scan.bckground,1)) round(size(scan.bckground,2))],rect);
bckground_rect = rect;
bckrgound_texture = Screen('MakeTexture', win, feval(flipfun,scan.bckground));

% make fixation dot texture
fix_texture_thin_full = {};
fix_texture_thick_full = {};
fix_texture_thick_left = {};
fix_texture_thick_right = {};
fix_rect_thin = CenterRect([0 0 round(size(scan.fix_im,1)) round(size(scan.fix_im,2))],rect);
fix_rect_thick = CenterRect([0 0 round(size(scan.fix_im,1)) round(size(scan.fix_im,2))],rect);

for ll = 1:size(scan.fix_im,4) % loop over luminance values
    fix_texture_thin_full{ll} = Screen('MakeTexture',win,feval(flipfun,cat(3,scan.fix_im(:,:,:,ll,1),params.stim.fix.dotopacity.*ones(size(scan.fix_im,1),size(scan.fix_im,2)))));
    fix_texture_thick_full{ll} = Screen('MakeTexture',win,feval(flipfun,cat(3,scan.fix_im(:,:,:,ll,2),params.stim.fix.dotopacity.*ones(size(scan.fix_im,1),size(scan.fix_im,2)))));
    fix_texture_thick_left{ll} = Screen('MakeTexture',win,feval(flipfun,cat(3,scan.fix_im(:,:,:,ll,3),params.stim.fix.dotopacity.*ones(size(scan.fix_im,1),size(scan.fix_im,2)))));
    fix_texture_thick_right{ll} = Screen('MakeTexture',win,feval(flipfun,cat(3,scan.fix_im(:,:,:,ll,4),params.stim.fix.dotopacity.*ones(size(scan.fix_im,1),size(scan.fix_im,2)))));
end

% run functions as first time running them always takes more time
GetSecs;
now;
ceil(1);
fprintf('');

% display instruction screen
Screen('FillRect',win,params.stim.bckgrnd_grayval,rect);
Screen('Flip',win);

if ~isempty(introscript)
    if ischar(introscript)
        evalin('caller',introscript);
    else
        feval(introscript);
    end
end

fprintf('Instructions are on screen, waiting for trgger...\n');


%% SAFEMODE: Wait for a key press to start
fprintf('press trigger key to begin the movie. (consider turning off network, energy saver, software updates.)\n');
% safemode = 0;
while 1
    [secs,keyCode,deltaSecs] = KbWait(deviceNr,2); % previously deviceNr = -3;
    temp = KbName(keyCode);
%     if isequal(temp(1),'=')
%         if safemode
%             safemode = 0;
%             fprintf('SAFE MODE OFF (the scan can start now).\n');
%         else
%             safemode = 1;
%             fprintf('SAFE MODE ON (the scan will not start).\n');
%         end
%     else
%         if safemode
%         else
            if isempty(params.triggerkey) || any(strcmp(temp(1), params.triggerkey))
                break;
            end
%         end
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% log the start!

fprintf('STIMULUS STARTED.\n');
feval(tfunEYE) %Eyelink('Message','SYNCTIME'));
timekeys = [timekeys; {GetSecs 'trigger'}];


%% %%%%%%% MAKE TEXTURES

prev_blockID = timing.trig_block(1);

% show the stim
framecnt = 0;
for frame = 1:size(frameorder,2)+1
    
    framecnt = framecnt + 1;
    frame0 = floor(frame);
    reporttext = '';
    
    blockID = timing.trig_block(frame);
    
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
    
    % get fixation dot
    opacity_idx = timing.trig_fix(frame);
    if isnan(timing.trig_spatial_cue(frame))
        fix_tex0 = fix_texture_thin_full{opacity_idx};
        fix_rect = fix_rect_thick;
    else
        switch timing.trig_spatial_cue(frame)
            case 1
                fix_tex0 = fix_texture_thick_left{opacity_idx};
                fix_rect = fix_rect_thin;
            case 2
                fix_tex0 = fix_texture_thick_right{opacity_idx};
                fix_rect = fix_rect_thin;
            case 0
                fix_tex0 = fix_texture_thick_full{opacity_idx};
                fix_rect = fix_rect_thin;
        end
    end
    
    
    if frame==1 || prev_blockID~=blockID
        sendELmessage = true;
    else
        sendELmessage = false;
    end
    
    switch timing.trig_stim(frame)
        
        % 0  : pre/post blank
        % 94 : trial_start_ID 
        % 95 : spatial_cue_ID 
        % 96 : delay_ID
        % 97 : response_ID
        % 98 : ITI
        % 99 : IBI
        case {0, 94, 95, 96, 97, 98} 
            % Draw background and thin rim fix dot texture
            Screen('DrawTexture',win,bckrgound_texture,[], bckground_rect, 0, filtermode, 1, framecolor);...
            Screen('DrawTexture',win,fix_tex0,[], fix_rect, 0, filtermode, params.stim.fix.dotopacity, [framecolor,params.stim.fix.dotopacity]);...
        
        case 93 % task_cue_ID            
%             showinstructionscreen_pinknoisebackground(taskscript{timing.trig_block(frame)},250,75,25,...
%                 params.offsetpix,bckrgound_texture,fix_tex0,bckground_rect,fix_rect,filtermode,framecolor,params.stim.fix.dotopacity) 

            fd = fopen(taskscript{timing.trig_block(frame)}, 'rt');
            mytext = '';
            tl = fgets(fd);
            lcount = 0;
            while lcount < 48
                mytext = [mytext tl]; %#ok<*AGROW>
                tl = fgets(fd);
                lcount = lcount + 1;
            end
            fclose(fd);
            instr=mytext;

            rect_text = rect + [offset(1) offset(2) offset(1) offset(2)];
            
            % draw background + fix dot
            Screen('DrawTexture',win,bckrgound_texture,[], bckground_rect, 0, filtermode, 1, framecolor);...
            Screen('DrawTexture',win,fix_tex0,[], fix_rect, 0, filtermode, params.stim.fix.dotopacity, [framecolor,params.stim.fix.dotopacity]);...                
            Screen('TextSize', win, 25);
            Screen('TextStyle', win, 0);
            DrawFormattedText(win, instr, 'center', (rect_text(4)/2)-250,0,75,[],[],[],[],rect_text);
            
        case {1:39} % stimuli
            % %             wm_blocks = find(~cellfun(@isempty, regexp(params.exp.stimTaskLabels,'wm')));
            % %             ltm_blocks = find(~cellfun(@isempty, regexp(params.exp.stimTaskLabels,'ltm')));
            % %             img_blocks = find(~cellfun(@isempty, regexp(params.exp.stimTaskLabels,'img')));
            
            % Draw background + dot
            Screen('DrawTexture',win,bckrgound_texture,[], bckground_rect, 0, filtermode, 1, framecolor);...

            % Draw stim texture + % Draw thin rim fix dot texture
            % trig_seq_exp_im_w_cd is a cell with dims: frames x 1, where each cell has 1 or 2 sides (1:l, 2:r)
            for side = 1:length(timing.trig_seq_exp_im_w_cd{frame,:})
                txttemp = feval(flipfun,timing.trig_seq_exp_im_w_cd{frame,side}{1});
                stim_rect = scan.centers{frame,side};
                stim_texture = Screen('MakeTexture',win, txttemp);
                Screen('DrawTextures',win,stim_texture,[],stim_rect,0,filtermode,1,framecolor);
            end

            Screen('DrawTexture',win,fix_tex0,[], fix_rect, 0, filtermode, params.stim.fix.dotopacity, [framecolor,params.stim.fix.dotopacity]);...
            Screen('Close',stim_texture);
    end
    
    
    % give hint to PT that we're done drawing
    Screen('DrawingFinished',win);
 
    %%%%%%%%%%%%%%%%%%%%%%%% the main while loop that actually puts up stimuli and records button presses
%     if params.debugmode
%         % get matlab now for the very first stimulus frame
%         if framecnt==1
%             absnowtime = now;
%         end
%         
%         VBLTimestamp = Screen(win, 'Flip', win, 0);
% %         VBLTimestamp = Screen(win, 'Flip', win, absnowtime + timing.trig_timing(frame) - (mfi * (1/2)));
%         timeframes(framecnt) = VBLTimestamp;
%         
%         % report text to command window?
%         if ~isempty(reporttext)
%             fprintf(reporttext);
%         end
%         
%         % Send eyelink a message when new block is started
%         ts = GetSecs;
%         timekeys = [timekeys; {ts blockID}];
%         
%         if detectinput
%             [keyIsDown,secs,keyCode,~] = KbCheck(-3);  % all devices
%             if keyIsDown
%                 
%                 % get the name of the key and record it
%                 kn = KbName(keyCode);
%                 timekeys = [timekeys; {secs kn}];
%                 
%                 % check if ESCAPE was pressed
%                 if isequal(kn,'ESCAPE')
%                     fprintf('Escape key detected.  Exiting prematurely.\n');
%                     getoutearly = 1;
%                     break;
%                 end
%             end
%         end
%     
%     else
        
        when = 0;

        % here, deal with making the stimulus frame / texture / stuff
        % read input until we have to do the flip
        while 1

            % EK: Why would load gamma every frame???
            %           % load the gamma table (for a future frame)
            %           if ~isempty(specialcon)
            %               frameL = frame0 + specialcon{4};
            %               if frameL <= size(frameorder,2)
            %                   if frameorder(1,frameL)==0  % if blank frame, who cares, don't change
            %                   else
            %                       con = specialcon{2}(frameorder(1,frameL));
            %                       if lastsc ~= con
            %                           %            sound(sin(1:100),1);
            %                           Screen('LoadNormalizedGammaTable',win,specialcluts(:,:,allcons==con));  % don't use loadOnNextFlip!
            %                           lastsc = con;
            %                       end
            %                   end
            %               end
            %           end




            % if we are in the initial case OR if we have hit the when time, then display the frame
            if when == 0 | GetSecs >= when

                % EK: we don't need this????
                % right before we issue the flip command, deal with frameevents.
                % hopefully the frameevent doesn't slow us down.
                %                 if ~isempty(frameevents)
                %                     if frameevents(frame0)==0
                %                     else
                %                         temp = framefuncs{frameevents(frame0)};
                %                         if ischar(temp)
                %                             evalin('caller',temp);
                %                         else
                %                             feval(temp);
                %                         end
                %                     end
                %                 end

                % issue the flip command and record the empirical time
                [VBLTimestamp,~,~,~,~] = Screen('Flip',win,  0);
                timeframes(framecnt) = VBLTimestamp;

                % get matlab now for the very first stimulus frame
                if framecnt==1
                    absnowtime = now;
                end

                % report text to command window?
                if ~isempty(reporttext)
                    fprintf(reporttext);
                end

                % Send eyelink a message when new block is started
                ts = GetSecs;
                timekeys = [timekeys; {ts blockID}];
                if params.wanteyetracking && sendELmessage 
                    Eyelink('message', sprintf('BLOCKID %d %d',blockID,ts));
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
                    [keyIsDown,secs,keyCode,~] = KbCheck(-3);  % all devices
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
%     end

    % write to file if desired
%     if wantframefiles
%         if isempty(framefiles{2})
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
        when = whendesired - mfi * (9/10);
    else
        % if there were no glitches, just proceed from the last recorded time
        % and set the when time to 9/10 a frame before the desired time.
        % notice that the accuracy of the mfi is only weakly assumed here,
        % since we keep resetting to the empirical VBLTimestamp.
        whendesired = VBLTimestamp + mfi * frameduration;
        when = whendesired - mfi * (9/10);  % should we be less aggressive??
    end
    
%     % PUT THE STIMULUS ON THE SCREEN USING A FLIP
%     Screen('FillRect',win,params.stim.bckrgnd_grayval,rect);
%     Screen('Flip',win);
    
    % update when
    if didglitch
        % if there were glitches, proceed from our earlier when time.
        % set the when time to 9/10 a frame before the desired frame.
        % notice that the accuracy of the mfi is strongly assumed here.
        whendesired = whendesired + mfi * frameduration;
        when = whendesired - mfi * (9/10);
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
Screen('Close',bckrgound_texture);

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



% % do cleanup if necessary
% if ~isempty(cleanupscript)
%   if ischar(cleanupscript)
%     evalin('caller',cleanupscript);
%   else
%     feval(cleanupscript);
%   end
% end

% % do some checks
% if wantcheck
%   ptviewmoviecheck(timeframes,timekeys,[],'t');
% end



end



%% OLD STUFF MAIN EXPERIMENT


% OLD STUFF TR trigger to start experiment
% getKey(params.triggerkey,k);
% fprintf('*** TRIGGER DETECTING. NOW STARTING. ***\n');
% ts = GetSecs; %record when trigger happened
% Eyelink('message', sprintf('TRIGGER %d',ts)); % USE


% target= struct('target_color',[],'xy_pixels',[],'background_color',[]);

% if scan.trigger
%     Screen(win, 'TextSize', 24);
% else
%     Screen(win, 'TextSize', 18);
% end
%
% Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% rect = Screen('Rect',screennum);  % what is the total rect
% rect = CenterRect(round([0 0 rect(3)*winsize rect(4)*winsize]),rect);
% [win, rect] = Screen('OpenWindow',max(Screen('Screens')),127,rect);
%
% params.ifi = Screen('GetFlipInterval',win);
% Screen('Preference', 'SyncTestSettings', .0004);
% Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% Screen('Preference','TextRenderer',1);
% HideCursor;


%%%%%%% Load textures
% fTex = cell(1,length(scan.allFlips)-1);
% for n = 1:(length(scan.allFlips)-1)
%
%     if any(flip.IDs{n} > 0) % draw 1 or 4 squares
%         trialImages = flip.IDs{n};
%
%         for ii = 1:length(trialImages)
%             f = squares{trialImages(ii)};
%             fTex{n}(ii) = Screen('MakeTexture',win,f);
%             ok = Screen('PreloadTextures', win, fTex{n}(ii)); %#ok<NASGU>
%
%         end
%     end
% end

%% DISPLAY TASK INSTRUCTIONS PRE-BACKTICK
% [w,h] = RectSize(Screen('TextBounds',win,p.task));
% DrawFormattedText(win, params.task, xc-(w/2),yc-(h/2), params.taskColor(1,:),[], flipLR, flipUD);
% Screen(win, 'Flip', 0);
% Screen('FillRect',win,grayval,rect);
%%% initial window - wait for backtick
% startText = 'Waiting for trigger, get ready!';
% Screen('DrawText',win, startText, 10,10,params.textColor);
