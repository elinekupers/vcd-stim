function [data, timeframes, digitrecord, trialoffsets] = ...
    vcd_showStimulus(params, ...
                        scan, ...
                        timing, ...
                        introscript, ...
                        taskscript, ...
                        tfunEYE)

getoutearly = 0;
glitchcnt   = 0;
when        = 0;

%% Preallocate space for key presses and timestamps
data.timeKeys       = {};
data.digitrecord    = [];
data.trialoffsets   = [];
data.digitframe     = [];
data.digitpolarity  = [];
                 
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
% <scfactor> (optional) is a positive number with the scaling to apply
%   to the images in <images>.  if supplied, we multiply the number
%   of pixels in each dimension by <scfactor> and then round.  we use
%   bilinear filtering when rendering the images.  default is 1, and in
%   this case, we use nearest neighbor filtering (which is presumably
%   faster than bilinear filtering).
% <allowforceglitch> (optional) is
%   0 means do nothing special
%   [1 D] means allow keyboard input 'p' to force a glitch of duration D secs.
%     note that this has an effect only if <detectinput> is on.
%     forcing glitches is useful for testing purposes.
%   default: 0.
filtermode = choose(params.stim.scfactor==1,0,1);
framecolor = 255*ones(1,3); 
allowforceglitch = 0;
frameorder = size(timing.trig_stim,2);

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
HideCursor;
Priority(9);
mfi = Screen('GetFlipInterval',win);  % re-use what was found upon initialization!
mfi = 1/round(1/mfi);


%% Create background and fixation textures prior to exp onset (as we need them throughout the experiment)
bckground_rect = CenterRect([0 0 round(size(scan.bckground_im,1)) round(size(scan.bckground_im,2))],rect);
bckrgound_texture = Screen('MakeTexture', win, scan.bckground_im);

% make fixation dot texture
fix_texture_thin_full = {};
fix_texture_thick_full = {};
fix_texture_thick_left = {};
fix_texture_thick_right = {};
fix_rect_thin = CenterRect([0 0 round(size(scan.fix_im,1)) round(size(scan.bckground_im,2))],rect);
fix_rect_thick = CenterRect([0 0 round(size(scan.fix_im,1)) round(size(scan.bckground_im,2))],rect);

for ll = 1:length(scan.fix_im,4) % loop over luminance values
    fix_texture_thin_full{ll} = Screen('MakeTexture',win,cat(3,scan.fix_im(:,:,:,ll,1),params.stim.fix.dotopacity));
    fix_texture_thick_full{ll} = Screen('MakeTexture',win,cat(3,scan.fix_im(:,:,:,ll,2),params.stim.fix.dotopacity));
    fix_texture_thick_left{ll} = Screen('MakeTexture',win,cat(3,scan.fix_im(:,:,:,ll,3),params.stim.fix.dotopacity));
    fix_texture_thick_right{ll} = Screen('MakeTexture',win,cat(3,scan.fix_im(:,:,:,ll,4),params.stim.fix.dotopacity));
end

% run functions as first time running them always takes more time
GetSecs;
now;
ceil(1);
fprintf('');

% display instruction screen
Screen('FillRect',win,params.stim.bckgrnd_grayval,rect);

if ~isempty(introscript)
  if ischar(introscript)
    evalin('caller',setupscript);
  else
    feval(introscript);
  end
end

WaitSecs(1);
Screen('Flip',win);
fprintf('Instructions are on screen, waiting for trgger...\n');


%% SAFEMODE: Wait for a key press to start
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

% show the movie
framecnt = 0;
for frame = 1:(frameorder+1)
  
    framecnt = framecnt + 1;
    frame0 = floor(frame);
    reporttext = '';

      % we have to wait until the last frame of the run sequence is done.  
      if frame0 == size(frameorder,2)+1
        while 1
          if GetSecs >= whendesired
            getoutearly = 1;
            % Log the end
            if ~isempty(tfunEYE)
              feval(tfunEYE);
              data.timekeys = [data.timekeys; {GetSecs 'trigger'}];
            end
            break;
          end
        end
      end

      % get out early?
      if getoutearly
        break;
      end

      opacity_idx = timing.fix_seq(frame);
      if timing.spatial_cue_seq==1
          fix_tex0 = fix_texture_thin_left{opacity_idx};
          fix_rect = fix_rect_thin;
      elseif timing.spatial_cue_seq==2
          fix_tex0 = fix_texture_thin_right{opacity_idx};
          fix_rect = fix_rect_thin;
      elseif timing.spatial_cue_seq==0
          fix_tex0 = fix_texture_thin_full{opacity_idx};
          fix_rect = fix_rect_thin;
      elseif isnan(timing.spatial_cue_seq)
          fix_tex0 = fix_texture_thick_full{opacity_idx};
          fix_rect = fix_rect_thick;
      end
      
      for side = 1:length(frameorder)
          

          
          switch frameorder{frame,side}
              
                
              case 0
                  % Draw background and thin rim fix dot texture
                  Screen('DrawTexture',win,bckrgound_texture,[],bckground_rect,0,filtermode,1,framecolor);
                  Screen('DrawTexture',win,fix_tex0,[],fix_rect,0,filtermode,1, params.stim.fix.dotopacity);
                  
              case 93 % task_cue_ID
                  task_ID = timing.seq;
                  
                  showinstructionscreen_pinknoisebackground(taskscript{task_ID},250,75,25,...
                      params.offsetpix,bckrgound_texture,fix_tex0,bckground_rect,fix_rect,filtermode,framecolor,params.stim.fix.dotopacity)
                  
%                   % Draw background and thin rim fix dot texture
%                   Screen('DrawTexture',win,bckrgound_texture,[],bckground_rect,0,filtermode,1,framecolor);
%                   Screen('DrawTexture',win,fix_texture_thick{frame},[],fix_rect_thin,0,filtermode,1, params.stim.fix.dotopacity);

                
              case 94 % trial_start_ID
                  % Draw background and thick rim fix dot texture
                  Screen('DrawTexture',win,bckrgound_texture,[],bckground_rect,0,filtermode,1,framecolor);
                  Screen('DrawTexture',win,fix_tex0,[],fix_rect,0,filtermode,1, params.stim.fix.dotopacity);

              case 95 % spatial_cue_ID
                  % Draw background and thick rim fix dot texture
                  Screen('DrawTexture',win,bckrgound_texture,[],bckground_rect,0,filtermode,1,framecolor);
                  Screen('DrawTexture',win,fix_tex0,[],fix_rect,0,filtermode,1, params.stim.fix.dotopacity);
                  
              case 96 % delay_ID
                  
              case 97 % response_ID 
                  
              case 98 % ITI_ID
                  % Draw background and thin rim fix dot texture
                  Screen('DrawTexture',win,bckrgound_texture,[],bckground_rect,0,filtermode,1,framecolor);
                  Screen('DrawTexture',win,fix_tex0,[],fix_rect,0,filtermode,1, params.stim.fix.dotopacity);
                  
              case 99 % IBI_ID
                  % Draw background and thin rim fix dot texture
                  Screen('DrawTexture',win,bckrgound_texture,[],bckground_rect,0,filtermode,1,framecolor);
                  Screen('DrawTexture',win,fix_tex0,[],fix_rect,0,filtermode,1, params.stim.fix.dotopacity);
                 
              case {5} % WM trials    
                  
              case {1:39} % stimuli
                  
                  % Draw background and thin rim fix dot texture
                  Screen('DrawTexture',win,bckrgound_texture,[],bckground_rect,0,filtermode,1,framecolor);
                  Screen('DrawTexture',win,fix_tex0,[],fix_rect,0,filtermode,1, params.stim.fix.dotopacity);
                  
                  % Draw stim texture
                  % exp_im is a cell with dims: blocks x trials x locations (1:l, 2:r) x stim epoch (first or second)
                  stim_texture = {};
                  for nn = 1:length(scan.exp_im)
                      for mm = 1:length(scan.exp_im{nn})
                          txttemp = feval(flipfun,scan.exp_im{frame,side,1});
                          
                          stim_texture{nn,mm} = Screen('MakeTexture',win, scan.exp_im{nn});
                          % Screen('DrawTexture',win,stim_texture,[],bckground_rect,0,filtermode,1,framecolor(frame0,:));
                          Screen('Close',stim_texture);
                      end
                  end
                  
                  
                  
          end
      end
     
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
Screen('FillRect',win,params.stim.bckrgnd_grayval,rect);
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

Screen('Close',bckrgound_texture);
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
