function [data,getoutearly] = vcd_showStimulus(params, ...
    scan, ...
    timing, ...
    introscript, ...
    taskscript, ...
    tfunEYE, ...
    deviceNr)

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
    flipfun = @(x) flipdim(flipdim(x,1),2); %#ok<*DFLIPDIM>
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
timeframes = NaN(1, floor(size(frameorder,2)-1)+1);

% make tmpdir
if wantframefiles
    tmpDir = '~/Desktop/tmp/'; %#ok<UNRCH>
    framefiles = {'~/Desktop/tmp/frame%05d.png', []};
    if ~exist('tmpDir','dir'), mkdir(tmpDir); end
end

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
else
    ShowCursor;
end

%% Create background and fixation textures prior to exp onset (as we need them throughout the experiment)
% bckground_rect = CenterRect([0 0 round(size(scan.bckground,1)) round(size(scan.bckground,2))],rect);
bckground_rect    = rect;
bckrgound_texture = Screen('MakeTexture', win, feval(flipfun,scan.bckground));

% make fixation dot texture
fix_texture_thin_full = {};
fix_texture_thick_full = {};
fix_texture_thick_left = {};
fix_texture_thick_right = {};
fix_rect_thin = CenterRect([0 0 round(size(scan.fix_im,1)) round(size(scan.fix_im,2))],rect);
fix_rect_thick = CenterRect([0 0 round(size(scan.fix_im,1)) round(size(scan.fix_im,2))],rect);
for ll = 1:size(scan.fix_im,4) % loop over luminance values
    fix_texture_thin_full{ll} = Screen('MakeTexture',win,feval(flipfun,scan.fix_im(:,:,:,ll,1)));   %cat(3, ...,params.stim.fix.dotopacity.*ones(size(scan.fix_im,1),size(scan.fix_im,2)))));
    fix_texture_thick_full{ll} = Screen('MakeTexture',win,feval(flipfun,scan.fix_im(:,:,:,ll,2)));   %cat(3, ...,,params.stim.fix.dotopacity.*ones(size(scan.fix_im,1),size(scan.fix_im,2)))));
    fix_texture_thick_left{ll} = Screen('MakeTexture',win,feval(flipfun,scan.fix_im(:,:,:,ll,3)));  %cat(3, ...,params.stim.fix.dotopacity.*ones(size(scan.fix_im,1),size(scan.fix_im,2)))));
    fix_texture_thick_right{ll} = Screen('MakeTexture',win,feval(flipfun,scan.fix_im(:,:,:,ll,4)));  %cat(3, ...,params.stim.fix.dotopacity.*ones(size(scan.fix_im,1),size(scan.fix_im,2)))));
end

fix_texture_alpha = Screen('MakeTexture',win,feval(flipfun,scan.fix_alpha_mask));

%% create single vector for fix textures
% 
% % set up fixation dot textures
% opacity_idx = timing.trig_seq_fix(framecnt);
% 
% if isnan(timing.seq_spatial_cue(framecnt))
%     fix_tex0 = fix_texture_thin_full{opacity_idx};
%     fix_rect = fix_rect_thin;
% else
%     switch timing.trig_spatial_cue(framecnt)
%         case 1
%             fix_tex0 = fix_texture_thick_left{opacity_idx};
%             fix_rect = fix_rect_thick;
%         case 2
%             fix_tex0 = fix_texture_thick_right{opacity_idx};
%             fix_rect = fix_rect_thick;
%         case 0
%             fix_tex0 = fix_texture_thick_full{opacity_idx};
%             fix_rect = fix_rect_thick;
%     end
% end


%%




% run functions as first time running them always takes more time
GetSecs;
now;
ceil(1);
fprintf('');

% display instruction screen
fd = fopen(introscript, 'rt');
instrtext = '';
tl = fgets(fd);
lcount = 0;
while lcount < 48
    instrtext = [instrtext tl]; %#ok<*AGROW>
    tl = fgets(fd);
    lcount = lcount + 1;
end
fclose(fd);

rect_text = rect + [params.offsetpix(1) params.offsetpix(2) params.offsetpix(1) params.offsetpix(2)];

% draw background + fix dot
Screen('DrawTexture',win,bckrgound_texture,[], bckground_rect, 0, filtermode, 1, framecolor);...
  
Screen('TextSize', win, 25);
Screen('TextStyle', win, 0);
DrawFormattedText(win, instrtext, 'center', (rect_text(4)/2)-250,0,75,[],[],[],[],rect_text);
Screen('Flip',win);

fprintf('Instructions are on screen, waiting for trgger...\n');


%% SAFEMODE: Wait for a key press to start
fprintf('press trigger key to begin the movie. (consider turning off network, energy saver, software updates.)\n');
% safemode = 0;
while 1
    [~,keyCode,~] = KbWait(deviceNr,2); % previously deviceNr = -3; outputs: secs,keyCode,deltaSecs
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

% fprintf('STIMULUS STARTED.\n');
feval(tfunEYE) %Eyelink('Message','SYNCTIME'));
timekeys = [timekeys; {GetSecs 'trigger'}];


%% %%%%%%% MAKE TEXTURES
blockcnt = 0;
prev_blockID = timing.trig_block(1);

% show the stim
framecnt = 0;
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
    
    % get fixation dot
    opacity_idx = timing.trig_fix(framecnt);
    if isnan(timing.trig_spatial_cue(framecnt))
        fix_tex0 = fix_texture_thin_full{opacity_idx};
        fix_rect = fix_rect_thin;
    else
        switch timing.trig_spatial_cue(framecnt)
            case 1
                fix_tex0 = fix_texture_thick_left{opacity_idx};
                fix_rect = fix_rect_thick;
            case 2
                fix_tex0 = fix_texture_thick_right{opacity_idx};
                fix_rect = fix_rect_thick;
            case 0
                fix_tex0 = fix_texture_thick_full{opacity_idx};
                fix_rect = fix_rect_thick;
        end
    end
    
    if prev_blockID~=blockID
        sendELmessage = true;
        blockcnt = blockcnt+1;
    else
        sendELmessage = false;
    end
    
    if frame==1, sendELmessage = true; end
    
    % Draw background 
    Screen('DrawTexture',win,bckrgound_texture,[], bckground_rect, 0, filtermode, 1, framecolor);

    switch timing.trig_stim(frame,1)
        
        % 0  : pre/post blank
        % 93 : exp_session.miniblock.response_ID
        % 94 : exp_session.miniblock.trial_start_ID
        % 95 : exp_session.miniblock.spatial_cue_ID
        % 96 : exp_session.miniblock.delay_ID
        % 97 : exp_session.miniblock.task_cue_ID
        % 98 : exp_session.miniblock.ITI_ID
        % 99 : exp_session.miniblock.IBI_ID
   
        case {0, 93, 94, 95, 96, 98, 99} 
            % Draw fix dot on top
            Screen('DrawTexture',win,fix_tex0,[], fix_rect, 0, filtermode, params.stim.fix.dotopacity, []);
            Screen('DrawTexture',win,fix_texture_alpha,[], fix_rect_thick, 0, filtermode, [], framecolor);
            
        case 97 % task_cue_ID  
            
            fix_tex0 = fix_texture_thick_full{opacity_idx};
            fix_rect = fix_rect_thick;
            Screen('DrawTexture',win,fix_tex0,[], fix_rect, 0, filtermode, params.stim.fix.dotopacity, framecolor);...
            
            fd = fopen(taskscript{~cellfun(@isempty, regexp(taskscript,sprintf('%02d',blockID),'match'))}, 'rt');
            instrtext = '';
            tl = fgets(fd);
            lcount = 0;
            while lcount < 48
                instrtext = [instrtext tl]; %#ok<*AGROW>
                tl = fgets(fd);
                lcount = lcount + 1;
            end
            fclose(fd);

            rect_text = rect + [params.offsetpix(1) params.offsetpix(2) params.offsetpix(1) params.offsetpix(2)];
            DrawFormattedText(win, instrtext, 'center', (rect_text(4)/2)-50,0,75,[],[],[],[],rect_text);
        
            
                 
        % Draw stimulus textures
        % 1-30 = all 2 peripheral stimulus aperture stim-task crossings:
        %   01-xx = Gabors
        %   xx-xx = RDKs
        %   xx-xx = Simple dot
        %   xx-30 = Complex objects
        
        case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30} 
            if blockID <= 30
                % trig_seq_exp_im_w_cd is a cell with dims: frames x 1, where each cell has 1 or 2 sides (1:l, 2:r)
                for side = 1:length(find(~cellfun(@isempty, timing.trig_seq_exp_im_w_cd{frame})))
                    txttemp = feval(flipfun,timing.trig_seq_exp_im_w_cd{frame}{side});
                    stim_rect = scan.rects{frame,side};
                    stim_texture = Screen('MakeTexture',win, txttemp);
                    Screen('DrawTexture',win,stim_texture,[],stim_rect,0,filtermode,1,framecolor);

                    % draw mask if it is there
                    if ~isempty(timing.trig_seq_exp_im_mask{frame}{1})
                        if ~isnan(timing.trig_seq_exp_im_mask{frame}{1})
                            txttemp = feval(flipfun,timing.trig_seq_exp_im_mask{frame}{1});
                            mask_texture = Screen('MakeTexture',win, txttemp);
                            Screen('DrawTexture',win,stim_texture,[],stim_rect,0,filtermode,1,framecolor);
                        end
                    end
                end
            
            else % 31-39 = Natural Scene stim-task crossings ({31,32,33,34,35,36,37,38,39})
                txttemp = feval(flipfun,timing.trig_seq_exp_im_w_cd{frame}{1});
                stim_rect = scan.rects{frame,1};
                stim_texture = Screen('MakeTexture',win, txttemp);
                Screen('DrawTexture',win,stim_texture,[],stim_rect,0,filtermode,1,framecolor);
            end
            
        % Draw fix dot on top
        Screen('DrawTexture',win,fix_tex0,[], fix_rect, 0, filtermode, params.stim.fix.dotopacity, []);
        Screen('DrawTexture',win,fix_texture_alpha,[], fix_rect_thick, 0, filtermode, [], framecolor);
            
    end
    
    
    % give hint to PT that we're done drawing
    Screen('DrawingFinished',win);
 
    %%%%%%%%%%%%%%%%%%%%%%%% the main while loop that actually puts up stimuli and records button presses
        
        when = 0;

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
            imwrite(Screen('GetImage',win),sprintf(framefiles{1},framecnt));
        else
            imwrite(uint8(placematrix(zeros([framefiles{2} 3]),Screen('GetImage',win))),sprintf(framefiles{1},framecnt));
        end
    end
    
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
Screen('Close',stim_texture);
Screen('Close',fix_tex0);
ptoff;

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




