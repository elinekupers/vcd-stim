function [data,getoutearly] = vcd_showStimulus(win, rect, params, ...
    scan, ...
    timing, ...
    bckground, ...
    fix_im, ...
    fix_mask, ...
    introscript, ...
    taskscript, ...
    tfunEYE, ...
    deviceNr, ...
    oldCLUT,...
    oldPriority)

getoutearly    = 0;
glitchcnt      = 0;
when           = 0;
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
frameorder           = 1:size(scan.frame_nr,1);

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
frameduration  = round(mfi)/params.stim.framedur_s; % 30 Hz presentation, 2 frames for office/psph monitors (60 Hz), 4 frames for BOLDscreen (120 Hz);

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

%% Make textures

%% Create background and fixation textures prior to exp onset (as we need them throughout the experiment)
% bckground_rect = CenterRect([0 0 round(size(scan.bckground,1)) round(size(scan.bckground,2))],rect);
bckground_rect    = rect;
bckrgound_texture = Screen('MakeTexture', win, feval(flipfun,bckground));

% make fixation dot texture
fix_texture_thin_full   = {};
fix_texture_thick_full  = {};
fix_texture_thick_left  = {};
fix_texture_thick_right = {};
fix_texture_thick_both  = {};
fix_rect_thin = CenterRect([0 0 round(size(fix_im,1)) round(size(fix_im,2))],rect);
fix_rect_thick = CenterRect([0 0 round(size(fix_im,1)) round(size(fix_im,2))],rect);
for ll = 1:size(scan.fix_im,4) % loop over luminance values
    fix_texture_thin_full{ll} = Screen('MakeTexture',win,feval(flipfun,  cat(3, fix_im(:,:,:,ll,1), fix_mask.*params.stim.fix.dotopacity)));
    fix_texture_thick_full{ll} = Screen('MakeTexture',win,feval(flipfun,  cat(3, fix_im(:,:,:,ll,2), fix_mask.*params.stim.fix.dotopacity)));
    fix_texture_thick_left{ll} = Screen('MakeTexture',win,feval(flipfun,  cat(3, fix_im(:,:,:,ll,3), fix_mask.*params.stim.fix.dotopacity)));
    fix_texture_thick_right{ll} = Screen('MakeTexture',win,feval(flipfun,  cat(3, fix_im(:,:,:,ll,4), fix_mask.*params.stim.fix.dotopacity)));
    fix_texture_thick_both{ll} = Screen('MakeTexture',win,feval(flipfun,  cat(3, fix_im(:,:,:,ll,5), fix_mask.*params.stim.fix.dotopacity)));
end

%% Prepare background and fixation texture vector outside the flip loop

fix_tex    = cell(length(scan.frame_nr),1);
fix_rect   = fix_tex;
im_tex     = fix_tex;
im_rect    = fix_tex;
framecolor = fix_tex;
txt_tex    = fix_tex;
txt_rect   = fix_tex;

for frame = 1:length(scan.frame_nr)
    
    blockID = scan.event_ID(frame);
    
    % set up fixation dot textures
    lum_idx = find(scan.fix_abs_lum(frame)==params.stim.fix.dotlum);
    
    if blockID == 0 && (scan.is_cued(frame)==0 || isnan(scan.is_cued(frame)))
        fix_tex{frame} = fix_texture_thin_full{lum_idx};
        fix_rect{frame} = fix_rect_thin;
        
    else
        if scan.event_id(frame,1)==95
            if scan.is_cued(frame)==1
                fix_tex{frame}  = fix_texture_thick_left{lum_idx};
                fix_rect{frame} = fix_rect_thick;
            elseif scan.is_cued(frame)==2
                fix_tex{frame}  = fix_texture_thick_right{lum_idx};
                fix_rect{frame} = fix_rect_thick;
            elseif scan.is_cued(frame)==3
                fix_tex{frame}  = fix_texture_thick_both{lum_idx};
                fix_rect{frame} = fix_rect_thick;
            end
        elseif any(intersect(scan.event_id(frame),[93,94,96,97]))
            fix_tex{frame} = fix_texture_thick_full{lum_idx};
            fix_rect{frame} = fix_rect_thick;
        elseif any(intersect(scan.event_id(frame),[98,99]))
            fix_tex{frame} = fix_texture_thin_full{lum_idx};
            fix_rect{frame} = fix_rect_thin;
        elseif scan.event_id(frame) < 90
            fix_tex{frame} = fix_texture_thick_full{lum_idx};
            fix_rect{frame} = fix_rect_thick;
        end
    end
    
    switch timing.trig_stim(frame,1)
        
        % 0  : pre/post blank
        % 93 : exp_session.block.response_ID
        % 94 : exp_session.block.trial_start_ID
        % 95 : exp_session.block.spatial_cue_ID
        % 96 : exp_session.block.delay_ID
        % 97 : exp_session.block.task_cue_ID
        % 98 : exp_session.block.ITI_ID
        % 99 : exp_session.block.IBI_ID
        
        % Draw background + thin fix dot on top
        case {0, 93, 94, 95, 96, 98, 99}
            
            % DrawTextures
            % * TexturePointers  need to be: n vector (where n is the number of textures)
            % * DestinationRects need to be: 4 row x n columns (where n is the number of textures)
            im_tex{frame} = cat(1, bckrgound_texture, fix_tex{frame});
            im_rect{frame} = cat(1, bckground_rect, fix_rect{frame});
            framecolor{frame} = 255*ones(2,3); % <framecolor> can also be size(<frameorder>,2) x 1 with values in [0,1] indicating an alpha change.
            
        case 97 % task_cue_ID
            
            script = taskscript{~cellfun(@isempty, regexp(taskscript,sprintf('%02d',blockID),'match'))};
            [task_instr, task_rect] = vcd_getInstructionText(params, script, rect);
            
            im_tex{frame} = cat(1, bckrgound_texture, fix_tex{frame});
            im_rect{frame} = cat(1, bckground_rect, fix_rect{frame});
            
            txt_tex{frame} = task_instr;
            txt_rect{frame} = task_rect;
            
            framecolor{frame} = 255*ones(2,3); % <framecolor> can also be size(<frameorder>,2) x 1 with values in [0,1] indicating an alpha change.
            
    end
end


% Get pre-run intructions
[instrtext, prerun_text_rect] = vcd_getInstructionText(params, introscript, rect);

% Draw background + dot
% windowPointer, texturePointer(s), [sourceRect], destRects, rotAngles, filterModes, globalAlphas, modulateColors, textureShader, specialFlags, auxParameters]);
Screen('DrawTextures',win, im_tex{1},[], im_rect{1}',[], [0;0], [1;1], framecolor{1}');
    
% Draw text
DrawFormattedText(win, instrtext, 'center', (prerun_text_rect(4)/2)-50,0,75,[],[],[],[],prerun_text_rect);
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


framecnt = 0;
for frame = 1:size(frameorder,2)+1
    
    framecnt = framecnt + 1;
    frame0 = floor(framecnt);
    whendesired = scan.frame_nr(framecnt);
    
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
    
    switch timing.trig_stim(framecnt)
        
        % 0  : pre/post blank
        % 93 : exp_session.block.response_ID
        % 94 : exp_session.block.trial_start_ID
        % 95 : exp_session.block.spatial_cue_ID
        % 96 : exp_session.block.delay_ID
        % 97 : exp_session.block.task_cue_ID
        % 98 : exp_session.block.ITI_ID
        % 99 : exp_session.block.IBI_ID
        
        % Draw background + thin fix dot on top
        case {0, 93, 94, 95, 96, 98, 99}
            
            % draw background and dot textures
            Screen('DrawTextures',win,im_tex{framecnt},[],im_rect{framecnt}',[0;0],[],[1;1],framecolor{framecnt}');
            
        case 97 % task_cue_ID
            
            % draw background and left/right cuing dot textures
            Screen('DrawTextures',win, im_tex{framecnt},[],im_rect{framecnt}',[0;0],[],[1;1],framecolor{framecnt}');
            
            % draw text
            % inputs are winptr, tstring, sx, sy, color, wrapat, flipHorizontal, flipVertical, vSpacing, righttoleft, winRect)
            DrawFormattedText(win, txt_tex{framecnt}, 'center', (txt_rect{framecnt}(4)/2)-50,0,75,[],[],[],[],txt_rect{framecnt});
            
        case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
            % Draw stimulus textures
            % 1-30 = all 2 peripheral stimulus aperture stim-task crossings:
            %   01-xx = Gabors
            %   xx-xx = RDKs
            %   xx-xx = Simple dot
            %   xx-30 = Complex objects
            %   100+ images are WM/LTM query images (lures or change
            %   blindness)
            Screen('DrawTexture',win, bckrgound_texture,[], bckground_rect, 0, [], 1, 255*ones(1,3));...
            blockID = timing.trig_block(framecnt);
        
            if blockID <= 30
                % trig_seq_exp_im_w_cd is a cell with dims: frames x 1, where each cell has 1 or 2 sides (1:l, 2:r)
                for side = 1:length(find(~cellfun(@isempty, timing.trig_seq_exp_im_w_cd{framecnt})))

                    txttemp = feval(flipfun, timing.trig_seq_exp_im_w_cd{framecnt}{side});
                    stim_rect = scan.rects{framecnt,side};

                    stim_texture = Screen('MakeTexture',win, txttemp);
                    Screen('DrawTexture',win,stim_texture,[], stim_rect, 0,[],1, 255*ones(1,3));
                end

            else % 31-39 = Natural Scene stim-task crossings ({31,32,33,34,35,36,37,38,39})
                txttemp = feval(flipfun,timing.trig_seq_exp_im_w_cd{framecnt}{1}); % there is only one image
                stim_rect = scan.rects{framecnt,1};
                stim_texture = Screen('MakeTexture',win, txttemp);
                Screen('DrawTexture',win,stim_texture,[], stim_rect, 0,[],1, 255*ones(1,3));
            end

            % Draw fix dot on top
            Screen('DrawTexture',win,fix_tex{frame},[], fix_rect{frame}, 0,[],1, 255*ones(1,3));

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


return




