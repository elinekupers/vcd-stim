function [data, timeframes, digitrecord, trialoffsets] = ...
    vcd_showStimulus(params, ...
                        scan, ...
                        timing, ...
                        introscript, ...
                        taskscript, ...
                        tfunEYE)


%% preallocate space
data.timeKeys = {};
                 
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




%% ptb stuff
win  = firstel(Screen('Windows'));
rect = Screen('Rect',win); % what is the total rect
% rect = CenterRect(round([0 0 rect(3)*winsize rect(4)*winsize]),rect);
[win, rect] = Screen('OpenWindow',max(Screen('Screens')),127,rect);
scan.ifi = Screen('GetFlipInterval',win);
Screen('Preference', 'SyncTestSettings', .0004);
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Preference','TextRenderer',1);
HideCursor;
Priority(9);

Screen('FillRect',win,grayval,rect);

% display instruction screen

% do setup if necessary
% showinstructionscreen(fileIn,txtOffset,instTextWrap,textsize,background,offset)
if ~isempty(introscript)
  if ischar(introscript)
    evalin('caller',setupscript);
  else
    feval(introscript);
  end
end

% WaitSecs(1);
% Screen('FillRect',win,params.stim.bckgrnd_grayval);
% Screen('Flip',win);
% DrawFormattedText(windowPtr,introscript,'center','center',textColor);
% Screen('Flip',win);
% fprintf('Instructions are on screen, waiting for trgger...\n');

%%% PTONMOVIE
% wait for a key press to start
fprintf('press trigger key to begin the movie. (consider turning off network, energy saver, software updates.)\n');
safemode = 0;
while 1
  [secs,keyCode,deltaSecs] = KbWait(-3,2);
  temp = KbName(keyCode);
  if isequal(temp(1),'=')
    if safemode
      safemode = 0;
      fprintf('SAFE MODE OFF (the scan can start now).\n');
    else
      safemode = 1;
      fprintf('SAFE MODE ON (the scan will not start).\n');
    end
  else
    if safemode
    else
      if isempty(params.triggerkey) || isequal(temp(1),params.triggerkey)
        break;
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% log the start!

fprintf('STIMULUS STARTED.\n');
eval(tfunEYE) %Eyelink('Message','SYNCTIME'));
data.timekeys = [data.timekeys; {GetSecs 'trigger'}];


%% %%%%%%% MAKE TEXTURES

% Make background texture
bckrgound_texture = Screen('MakeTexture', win, scan.bckground_im);
Screen('DrawTexture',win,bckrgound_texture,[],movierect,0,filtermode,1,framecolor(frame0,:));
Screen('Close',bckrgound_texture);

% Screen('FillRect',win,params.stim.bckgrnd_grayval);  % REMOVED! this means do whole screen.    % ,movierect);


% for nn = 1:length(fix_im)
%     fix_texture{nn,1} = Screen('MakeTexture', win, fix_im(:,:,:,nn,1)); % thin
%     fix_texture{nn,2} = Screen('MakeTexture', win, fix_im(:,:,:,nn,2)); % thick
% end
% 
% stim_texture = {};
% for nn = 1:length(exp_im)
%     for mm = 1:length(exp_im{nn})
%         stim_texture{nn,mm} = Screen('MakeTexture',win, exp_im{nn});
%     end
% end

txttemp = feval(flipfun,images(:,:,:,frameorder(1,frame0)));
stim_texture = Screen('MakeTexture',win,txttemp);
Screen('DrawTexture',win,stim_texture,[],movierect,0,filtermode,1,framecolor(frame0,:));
Screen('Close',stim_texture);

stim_texture = Screen('MakeTexture',win,cat(3,fixationimage(:,:,:,-fixationorder(1+frame0)),uint8(fixationorder(end)*fixationalpha)));
Screen('DrawTexture',win,stim_texture{p},[],fixationrect{p},0,0);
Screen('Close',stim_texture{p});

% give hint to PT that we're done drawing
Screen('DrawingFinished',win);

%%%%%%%%%%%%%%%%%%%%%%%% the main while loop that actually puts up stimuli and records button presses

when = 0;

% here, deal with making the stimulus frame / texture / stuff
% read input until we have to do the flip
while 1

  % if we are in the initial case OR if we have hit the when time, then display the frame
  if when == 0 | GetSecs >= when
  
    % issue the flip command and record the empirical time
    [VBLTimestamp,StimulusOnsetTime,FlipTimestamp,Missed,Beampos] = Screen('Flip',win,  0);
    timeframes(framecnt) = VBLTimestamp;
    
    % get matlab's now for the very first stimulus frame
    if framecnt==1
      absnowtime = now;
    end
    
%     % report text to command window?
%     if ~isempty(reporttext)
%       fprintf(reporttext);
%     end

    % if we missed, report it
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
    [keyIsDown,secs,keyCode,deltaSecs] = KbCheck(-3);  % all devices
    if keyIsDown

      % get the name of the key and record it
      kn = KbName(keyCode);
      data.timekeys = [data.timekeys; {secs kn}];

      % check if ESCAPE was pressed
      if isequal(kn,'ESCAPE')
        fprintf('Escape key detected.  Exiting prematurely.\n');
        getoutearly = 1;
        break;
      end

    end
  end

end


% collect response and measure timing
%     if t > 1 %ignore the first cell
%         previous = t-1;
%         if dot_cond(t) ~= dot_cond(previous)
%             ts = GetSecs;
%             Eyelink('message', sprintf('STIM_ONSET %d %d',dot_cond(t),ts));
%             data.timeKeys = [data.timeKeys; {ts 'stim'}];
%             trialEnd = ts-startTime;
%             data.timePerTrial(t) = trialEnd; % recorded trial duration
%         end
%     end
%     
%     %record keys
%     [keys RT] = recordKeys(startTime+(t-1)*viewTime,viewTime,k,params.ignorekeys); %KJ lag fix
%     data.keys{t} = keys;
%     data.rt(t) = min(RT);
%     if isequal(data.keys{t}, 'ESCAPE')
%         fprintf('Escape key detected. Exiting prematurely. \n');
%         getoutearly = 1;
%         return;
%     end 

% PUT THE STIMULUS ON THE SCREEN USING A FLIP
Screen('FillRect',win,grayval,rect);
Screen('Flip',win);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we have to wait until the last frame is done.  so this is how we hack that in.
if DONE WITH EXPERIMENT
  while 1
    if GetSecs >= whendesired
      eval(tfunEYE); %Eyelink('Message','SYNCTIME');
      timekeys = [timekeys; {GetSecs 'DONE'}];
      break;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PT CLEANUP STUFF

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



% % main for loop drawing stimuli   
% for t = 1:numTrials
%     
%     
%     currTrial = [];
%     
%     Screen('FillRect', windowPtr, target(this_trial).color); 
%     Screen('DrawTextures', windowPtr, target(this_trial).Tex, [], target(this_trial).Rects);
%     Screen('Flip',windowPtr); 
    
    
    
    
    
    %% RDK screen stuff
    %             %draw on the
    %             if any(isnan(prod(pos,1))==0)
    %                 Screen('DrawDots', screen_struct.cur_window, pos, dots_struct.dot_size, dots_struct.dot_color', AP.spec.center);
    %
    %                 %update the loop
    %                 loopi = loopi + 1;
    %                 if loopi > dots_struct.interval
    %                     loopi = 1;
    %                 end
    %
    %
    %                 %flip the screen to make things
    %                 Screen('Flip', screen_struct.cur_window);
    %             end
    
    %             %flip the screen to clean it
    %             Screen('Flip', screen_struct.cur_window);
    
    %         end
    
    
    
    

    
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
