function [timeframes, timekeys, digitrecord, trialoffsets] = ...
    vcd_showStimulus(p, images, subj_session, setupscript, soafun)




% Get [x,y]-center in pixels of peripheral stimuli given display size, stim size and offest 
scan.centers = [disp.xc + p.stim.x0_pix + p.offset_pix(1), ...
                disp.yc + p.stim.y0_pix + p.offset_pix(2)];

            
% get information about the PT setup
win = firstel(Screen('Windows'));
rect = Screen('Rect',win);

%% PREPARE IMAGES

%%%%%%% FLIP UP/DOWN, LEFT/RIGHT
% if p.movieflip(1)  % flip up-down
%     scan.centers = repmat([disp.w_pix disp.h_pix],length(scan.centers),1)-scan.centers;
% end

%%%%%%% CENTER RECT
% scan.rects = CenterRectOnPoint([0 0 scan.squareSize scan.squareSize],scan.centers(:,1), scan.centers(:,2));

%%%%%%% MAKE TEXTURES

% ptb stuff
HideCursor;
Priority(9);

Screen('FillRect',win,grayval,rect);

% Make textures for each stimulus class
stimclass = {'gabor','dot','rdk','cobj','ns'};
num_stimclass = length(stimclass);

for cc = 1:length(stimclass)
    
    texture.(stimclass(cc)) = {};

    for p=1:size(images.gabor(),4)
        texture{p} = Screen('MakeTexture',win, images{:,:,:,p});
    end
    
    
end



fTex = cell(1,length(scan.allFlips)-1);

%%%%%%% Load textures
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





%%%% PTONMOVIE
% % wait for a key press to start
% fprintf('press trigger key to begin the movie. (consider turning off network, energy saver, software updates.)\n');
% safemode = 0;
% while 1
%   [secs,keyCode,deltaSecs] = KbWait(-3,2);
%   temp = KbName(keyCode);
%   if isequal(temp(1),'=')
%     if safemode
%       safemode = 0;
%       fprintf('SAFE MODE OFF (the scan can start now).\n');
%     else
%       safemode = 1;
%       fprintf('SAFE MODE ON (the scan will not start).\n');
%     end
%   else
%     if safemode
%     else
%       if isempty(triggerkey) || isequal(temp(1),triggerkey)
%         break;
%       end
%     end
%   end
% end


% %% NEWER SETUP SCREEN???
% Screen('FillRect',win,grayval,rect);
% 
% %%% initial window - wait for backtick
% startText = 'Waiting for trigger, get ready!';
% Screen('DrawText',win, startText, 10,10,params.textColor);
% 
% % display instruction screen
% WaitSecs(1);
% Screen('FillRect',windowPtr,blankColor);
% Screen('Flip',windowPtr);
% DrawFormattedText(windowPtr,instruction_str,'center','center',textColor);
% Screen('Flip',windowPtr);
% fprintf('waiting for trgger...\n');
% 
% % TR trigger to start experiment
% getKey(params.triggerkey,k);
% fprintf('*** TRIGGER DETECTING. NOW STARTING. ***\n');
% ts = GetSecs; %record when trigger happened
% Eyelink('message', sprintf('TRIGGER %d',ts)); % USE tfunEYE
% subject.timeKeys = [subject.timeKeys; {ts 'trigger'}];

%% MAIN EXPERIMENT

% target= struct('target_color',[],'xy_pixels',[],'background_color',[]);

if scan.trigger
    Screen(win, 'TextSize', 24);
else
    Screen(win, 'TextSize', 18);
end

Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);





[w,h] = RectSize(Screen('TextBounds',win,p.task.));
DrawFormattedText(win, p.task, xc-(w/2),yc-(h/2), p.taskColor(1,:),[], flipLR, flipUD);
Screen(win, 'Flip', 0);







% main for loop drawing stimuli   
for t = 1:numTrials
    
    
    this_trial = dot_cond(t);
    
    Screen('FillRect', windowPtr, target(this_trial).color); 
    Screen('DrawTextures', windowPtr, target(this_trial).Tex, [], target(this_trial).Rects);
    Screen('Flip',windowPtr); 
    
    
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
    
    
    
    

    % collect response and measure timing
    if t > 1 %ignore the first cell
        previous = t-1;
        if dot_cond(t) ~= dot_cond(previous)
            ts = GetSecs;
            Eyelink('message', sprintf('STIM_ONSET %d %d',dot_cond(t),ts));
            subject.timeKeys = [subject.timeKeys; {ts 'stim'}];
            trialEnd = ts-startTime;
            subject.timePerTrial(t) = trialEnd; % recorded trial duration
        end
    end
    
    %record keys
    [keys RT] = recordKeys(startTime+(t-1)*viewTime,viewTime,k,p.ignorekeys); %KJ lag fix
    subject.keys{t} = keys;
    subject.rt(t) = min(RT);
    if isequal(subject.keys{t}, 'ESCAPE')
        fprintf('Escape key detected. Exiting prematurely. \n');
        getoutearly = 1;
        return;
    end 
end
