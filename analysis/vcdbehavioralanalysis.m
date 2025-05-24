function results = vcdbehavioralanalysis(filename);

% function results = vcdbehavioralanalysis(filename);
%
% <filename> is the .mat file saved after running the VCD experiment for one run
%
% Check timing, deal with buttons and triggers, and analyze the behavioral data.
% We include a number of sanity checks (various asserts and warnings).
% Important information is reported to the command window for visual inspection.
% Finally, we return comprehensive results in <results>.
%
% The <results> struct contains the following:
%   <trialinfo> is a table (documented below)
%   <totaldur> is the total empirical duration of the experiment in seconds.
%     if partial data, this will be less than what is desired.
%   <glitchcnt> is the number of glitches during the presentation
%   <droppedcnt> is the number of dropped frames during the presentation
%   <matlabnowtime> is the absolute time corresponding to the first frame
%   <mristarttime> is the time in seconds (relative to the first frame) 
%     that the trigger was detected ('trigger' in timekeys)
%   <donetime> is the time in seconds (relative to the first frame) 
%     corresponding to the completion of the last frame ('done' in timekeys).
%     if partial data, <donetime> will be returned as NaN.
%   <timeframes> is a vector of times indicating the presentation time of
%     each stimulus frame. time=0 corresponds to the onset of the first frame.
%     this version (in contrast to what is in the raw .mat file) has undergone
%     linear interpolation to fix dropped frames (which are originally coded as NaN).
%   <triggertimes> is a vector of times that we detected the trigger
%   <userkeys> is a cell vector of possible user keys
%   <userkeycounts> is a vector of number of times that the keys were pressed.
%     note that this is inclusive of everything that happened in the run.
%     (it isn't restricted to the response windows.)
%   === the following are included if we have at least 5 detected triggers:
%   <expectedtriggers> is the total number of TRs (triggers) that are expected
%   <mdf> is the median diff trigger-time, indicating our best estimate of the empirical TR
%   <numcrazymdferrs> is the number of mdf errors that are crazy (non-multiples)
%   <mdferrsok> is 0/1 indicating whether all mdf errors seem to be innocuous dropped triggers
%   ===
%   <summary> is a table (documented below)
%
% <results.trialinfo> is a table with N rows where N is the number of trials in the run.
% Note that the fixation task generates many dot change events --- we treat all of these events
% as if they are trials, even though they are not trials in the conventional VCD sense.
%   session_nr - as usual
%   session_type - as usual
%   run_nr - as usual
%   block_nr - only positive integer cases
%   crossing_nr - only positive integer cases
%   trial_nr - for tasks other than FIX, this is just the usual trial_nr.
%              for the FIX task, this is set to NaN for each dot change event.
%   is_catch - as usual. for the special FIX dot change events, this is set to 0.
%   onset_start - frame time corresponding to stimulus onset. (this is either
%                 stim1 or stim2 where appropriate, and for the special FIX dot change
%                 events, stimulus onset refers to the change in the dot luminance.)
%   onset_abstime - stimulus onset time as a MATLAB serial date number (units are days).
%   correct_response - 1 through 4 indicating the correct response
%   change_mind - 0/1 indicating whether the subject pressed more than one unique button.
%                 in the case of FIX dot change events, subjects are not allowed to change
%                 their mind, and change_mind is set to NaN.
%   button_pressed - 1 through 4 indicating the official button pressed by the subject
%   is_correct - 0/1 indicating whether the subject got the right answer
%   rt - reaction time in milliseconds (time between stimulus onset time and button press time)
%
%   * In the case of catch trials for non-FIX tasks, all of the following are set to NaN
%     no matter what the subject does:
%       correct_response, change_mind, button_pressed, is_correct, rt
%
%   * If the subject did not press a button, is_correct is set to 0, and all of the 
%     following are set to NaN:
%       change_mind, button_pressed, rt
%
%   * If partial data, onset_abstime will be NaN for trials not conducted. In addition,
%     these trials not conducted will be scored as if subjects did not press a button.
%
% <results.summary> is a table with B rows, where B is the number of blocks in the run.
%   block_nr - only positive integer cases
%   crossing_nr - only positive integer cases
%   response_rate - out of non-catch trials, percentage with at least one button press
%   pct_correct - out of non-catch trials, percentage for which the correct answer was given
%   median_rt - median reaction time in milliseconds (ignoring catch trials and cases where
%               no button was pressed)
%
% Approach to handling all tasks other than fixation task:
% - We extract buttons in a response window that extends 0 to 4000 ms after
%   a given stimulus onset. We ignore any button that is not one of the choice buttons.
% - We process only the final button pressed (if there is more than one).
% - We identify the first instance of a button (if there is a series).
%
% Approach to handling fixation task:
% - We consider fixation circle frames starting from post-task-ITI and going
%   through the end of the block.
% - Each dot change is treated as if it is a trial.
% - We extract buttons in a response window that extends 0 to 1400 ms after
%   a given dot change event. We ignore any button that is not one of the choice buttons.
% - We process the first button pressed within the response window.
% - We do not allow the subject to change their mind.
% - The initial dot color is not treated as a change event.

%% Internal constants

deltatime = 40;  % buttons/triggers are automatically extended if less than this many milliseconds elapse while still the same
validkeys = {'1!' '2@' '3#' '4$' '5%' 'r' 'y' 'g' 'b' 't' 'absolutetimefor0' 'trigger' 'DONE'};  % things we expect to get
userkeys = {'1' '2' '3' '4' 'r' 'y' 'g' 'b'};  % user-driven buttons
choicebuttons = {'1' '2' '3' '4'};  % the official buttons we expect subjects to press
triggerkeys = {'5' 't'};  % this indicates buttons that are interpreted as a trigger
responsewindow_reg = [0 4000];  % we accept buttons in this range of milliseconds after stimulus onset (for all tasks other than fixation)
responsewindow_fix = [0 1400];  % for the fixation task, we accept buttons in this range of milliseconds after each dot change

%% Setup

% load the data
a1 = load(filename);

% initialize results
clear results;
results.trialinfo = table;
results.summary = table;

%% Expand simultaneous-keypress cases

% shorthand
timekeys = a1.data.timeKeys;

% expand multiple-keypress cases
timekeysB = {};
for p=1:size(timekeys,1)
  if iscell(timekeys{p,2})
    for pp=1:length(timekeys{p,2})
      timekeysB{end+1,1} = timekeys{p,1};
      timekeysB{end,2} = timekeys{p,2}{pp};
    end
  else
    timekeysB(end+1,:) = timekeys(p,:);
  end
end

%% Deal with basic timeframes and duration stuff

% shorthand
timeframes = a1.data.timing.timeframes;
presentationrate_hz = a1.params.stim.presentationrate_hz;  % experiment is intended to run at 60 frames per second
refresh_hz = a1.params.disp.refresh_hz;                    % the idealized display refresh rate, e.g. 60 Hz or 120 Hz
glitchcnt = a1.data.timing.glitchcnt;                      % number of timing glitches

% check for dropped frames
droppedcnt = sum(isnan(timeframes));                       % number of dropped frames

% check that the mean frame-to-frame difference is sane (it should differ from the idealized rate by less than 1 ms)
mdtf = nanmean(diff(timeframes));
assert(abs(mdtf*1000 - 1000/presentationrate_hz) < 1,'frame-to-frame differences are OFF!!');

% calc total empirical duration of the experiment (this includes the full completion of the last frame).
% in the case of partial data (run stopped early), we calculate duration up through the last valid recorded timeframe 
ix = find(~isnan(timeframes));
actualnum = ix(end);  % use the last valid recorded timeframe
totaldur = mdtf * actualnum;

% check that the "DONE" time (recorded after last frame) is close to the official total duration (within 50 ms)
ix = find(ismember(timekeysB(:,2),'DONE'));
if isempty(ix)
  warning('*** Did not find DONE in timekeys. This must be partial data??? BEWARE! ***');
  donetime = NaN;
else
  donetime = timekeysB{ix,1};
  assert(abs(donetime-totaldur) < 0.050,'the DONE time is mismatched!!');
end

% report to the command window
fprintf('==============================================================\n');
fprintf('Stimulus fps is %d frames per sec; display refresh rate is %d Hz.\n',presentationrate_hz,refresh_hz);
fprintf('Experiment duration (empirical / ideal) was %.3f / %.3f.\n',totaldur,length(timeframes)/presentationrate_hz);
fprintf('Frames per sec (empirical / ideal) was %.6f / %.6f.\n',1/mdtf,presentationrate_hz);
fprintf('Number of glitches: %d.\n',glitchcnt);
fprintf('Number of dropped frames: %d.\n',droppedcnt);
fprintf('==============================================================\n');

% record
results.totaldur = totaldur;
results.glitchcnt = glitchcnt;
results.droppedcnt = droppedcnt;
results.donetime = donetime;

%% Deal with dropped frames

% use linear interpolation to fill in the NaNs
ii0 = find(~isnan(timeframes));
ii1 = find(isnan(timeframes));
timeframes(ii1) = interp1(ii0,timeframes(ii0),ii1,'linear','extrap');

% record
results.timeframes = timeframes;  % user must use this version from behavioral analysis and not the version in the raw data!!

%% Deal with terrible button pre-processing stuff

% figure out absolute time for the first frame (this is determined via Matlab's now)
matlabnowtime = timekeysB{find(ismember(timekeysB(:,2),'absolutetimefor0')),1};

% figure out the time that the experiment was started. this time is relative to the
% time of the first frame. for example, -0.017 means we detected the trigger
% 17 ms before we were actually able to show the first stimulus frame.
mristarttime = timekeysB{find(ismember(timekeysB(:,2),'trigger')),1};

% clean up all the button stuff
oldkey = '';
oldkeytime = -Inf;
oldtriggertime = -Inf;
buttontimes = [];    % vector of times in seconds
buttonpressed = {};  % cell vector of characters
triggertimes = [];   % vector of times in seconds
for p=1:size(timekeysB,1)

  % warn if weird key found
  if ~ismember(timekeysB{p,2},validkeys)
    fprintf('*** Unknown key detected (%s); ignoring.\n',timekeysB{p,2});
    continue;
  end

  % figure out auxiliary and trigger events
  bad1a = ismember(timekeysB{p,2},{'absolutetimefor0' 'trigger' 'DONE'});  % auxiliary events
  bad1b = ~bad1a & ismember(timekeysB{p,2}(1),triggerkeys);                % trigger events (as specified by <triggerkeys>)
  bad = bad1a | bad1b;                                                     % either auxiliary or trigger events
  
  % is this a "held down" case? (is the current key a user-pressed key that is repeated and within <deltatime>?)
  bad2 = (isequal(timekeysB{p,2},oldkey) & timekeysB{p,1}-oldkeytime <= deltatime/1000);

  % if it appears to be a new key, we should do a special case check for simultaneity.
  % bad3 indicates if we should ignore the current key because it comes while some 
  % other key was originally held down first.
  bad3 = 0;
  if ~bad && ~isequal(timekeysB{p,2},oldkey)  % if not an auxiliary/trigger and appears to be a new user-pressed key
  
    % scan ahead...
    q = p+1;
    while 1
      if q > size(timekeysB,1)
        break;  % if we run out, just break
      end
      if timekeysB{q,1} > oldkeytime + deltatime/1000
        break;  % if we are past the window, just break
      end
      if isequal(timekeysB{q,2},oldkey)
        bad3 = 1;  % if we are within the deltatime AND it is the same as the old key, then mark the current one for ignoral
        break;
      end
      q = q + 1;
    end

  end
  
  % if this is a held-down button, just extend the time
  if bad2
    oldkeytime = timekeysB{p,1};
  end

  % if not bogus, then record the button time (for user-pressed keys)
  if ~(bad | bad2 | bad3)
    buttontimes = [buttontimes timekeysB{p,1}];
    buttonpressed = [buttonpressed {timekeysB{p,2}(1)}];
    oldkey = timekeysB{p,2};
    oldkeytime = timekeysB{p,1};
  end
  
  % deal with triggers
  if bad1b
    if timekeysB{p,1}-oldtriggertime <= deltatime/1000  % if we are within the delta, just extend the time
      oldtriggertime = timekeysB{p,1};
    else
      triggertimes = [triggertimes timekeysB{p,1}];     % otherwise, record it
      oldtriggertime = timekeysB{p,1};
    end
  end

end

% record
results.matlabnowtime = matlabnowtime;
results.mristarttime = mristarttime;
results.triggertimes = triggertimes;

% Note: we don't record the buttontimes/buttonpressed into results
%       because the plan is to fully process them.
% However, we do record simple counts of button presses (see below).

%% Do some basic counts of buttons and triggers

% count buttons
userkeycounts = [];
for pp=1:length(userkeys)
  userkeycounts(pp) = sum(ismember(buttonpressed,userkeys{pp}));
end

% record
results.userkeys = userkeys;
results.userkeycounts = userkeycounts;

%% Do some checks that we got sane triggers

% calc
totaln = length(triggertimes);  % total number of triggers detected

% if it seems that this is just a behavior experiment
if totaln < 5

  warning('*** Detected less than 5 triggers. This is probably a behavior-only experiment? ***');

% if it seems that this is a run with actual triggers from the scanner
else

  % calc
  expectedtriggers = ceil((length(timeframes)/presentationrate_hz)/a1.params.exp.TR);  % total number of TRs expected
  mdf = median(diff(triggertimes));                    % typical empirical TR diff
  
  % check that the mdf and the intended TR are within 50 ms
  assert(abs(mdf-a1.params.exp.TR) < 50/1000);  

  % find diffs that are more than 10% larger or 10% smaller than the mdf.
  % express these as a multiple of the mdf.
  temp = diff(triggertimes);
  mdferrs = temp(temp > mdf*1.1 | temp < mdf*.9) / mdf;
  
  % how many mdf errors are crazy? (i.e. more than .1 away from being a perfect multiple)
  numcrazymdferrs = sum(abs(mdferrs-round(mdferrs)) > .1);
  if numcrazymdferrs > 0
    warning('*** We encountered a crazy mdf error; all bets are off!!!***');
  end

  % if the total number of dropped triggers based on the mdf errors
  % is equal to the number of triggers that are missing, then all is fine.
  % otherwise, issue a warning!
  % NOTE: this check only really makes sense when numcrazymdferrs is 0.
  mdferrsok = sum(round(mdferrs)-1) == (expectedtriggers - totaln);
  if ~mdferrsok
    warning('*** Dropped triggers are not innocuous!! ***');
  end
  
  % record
  results.expectedtriggers = expectedtriggers;
  results.mdf = mdf;
  results.numcrazymdferrs = numcrazymdferrs;
  results.mdferrsok = mdferrsok;

end

%% Basic behavior checks

% General notes on the time table:
% - 1st frame in timeframes is presented (onset) at t=0.
% - duration is 60 means the end of the 60th frame is t=60.
% - the next event also starts at t=60.

% check that the total duration of all events is exactly the number of timeframes
assert(sum(a1.run_table.event_dur)==length(timeframes));

% check that the timing field in run_frames is exactly 0:X
assert(isequal(a1.run_frames.timing(:),(1:length(a1.run_frames.timing))'-1));

%% Now analyze the behavior!

% initialize
blockcnt = 0;    % which block number to process next
trialcnt = 1;    % which trial number to process next
rii = 1;         % the 1-based trial count within the current run

% check that the run_nr in the table is exactly the same as what is in params
assert(all(a1.run_table.run_nr==a1.params.run_nr));

% loop
ii = 1;  % this keeps track of the current row of interest
prevwarn = warning('query','MATLAB:table:RowsAddedExistingVars');
warning('off','MATLAB:table:RowsAddedExistingVars');
while 1

  % if we are at the end, get out
  if ii > length(a1.run_table.block_nr)
    break;
  end
  
  % if we encounter 0 in block_nr, it's time to go to the next block
  if a1.run_table.block_nr(ii) == 0
    blockcnt = blockcnt + 1;
    trialcnt = 1;  % reset to 1
    ii = ii + 1;
    continue;
  end
  
  % if the current row isn't what we want (the desired block and the desired trial), then just increment
  if ~(a1.run_table.block_nr(ii) == blockcnt & a1.run_table.trial_nr(ii) == trialcnt)
    ii = ii + 1;
    continue;
  end
  
  % OK: we have found the desired block and trial!
  
  % check that the crossing_nr is listed in a1.taskIDs
  assert(ismember(a1.run_table.crossing_nr(ii),a1.taskIDs));
  
  % handle tasks that aren't the FIX task
  if a1.run_table.task_class(ii) ~= 1

    % record
    results.trialinfo.session_nr(rii) =   a1.run_table.session_nr(ii);
    results.trialinfo.session_type(rii) = a1.run_table.session_type(ii);
    results.trialinfo.run_nr(rii) =       a1.run_table.run_nr(ii);
    results.trialinfo.block_nr(rii) =     a1.run_table.block_nr(ii);
    results.trialinfo.crossing_nr(rii) =  a1.run_table.crossing_nr(ii);
    results.trialinfo.trial_nr(rii) =     a1.run_table.trial_nr(ii);
    results.trialinfo.is_catch(rii) =     a1.run_table.is_catch(ii);
    
    % figure out if this is a single- or double-image trial
    if a1.run_table.trial_type(ii)==1
      idmatch = 94;  % for stim1
    else
      assert(a1.run_table.trial_type(ii)==2);
      idmatch = 95;  % for stim2
    end
    
    % go find the one row we are looking for and jump to that one
    ii2 = find( a1.run_table.block_nr == blockcnt & ...
                a1.run_table.trial_nr == trialcnt & ...
                a1.run_table.event_id == idmatch );
    assert(length(ii2)==1 & ii2 >= ii);
    ii = ii2;
    
    % calc some times
    stimonset = timeframes(a1.run_table.event_start(ii)+1);
    windowstart = stimonset + responsewindow_reg(1)/1000;
    windowend   = stimonset + responsewindow_reg(2)/1000;

    % record
    results.trialinfo.onset_start(rii)      = a1.run_table.event_start(ii);
    results.trialinfo.onset_abstime(rii)    = matlabnowtime + stimonset/60/60/24;  % convert from seconds to days and then add in
    results.trialinfo.correct_response(rii) = a1.run_table.correct_response(ii);  % NOTE: catch trials should be NaN

    % if this is a catch trial
    if a1.run_table.is_catch(ii) == 1
    
      results.trialinfo.change_mind(rii)       = NaN;
      results.trialinfo.button_pressed(rii)   = NaN;
      results.trialinfo.is_correct(rii)       = NaN;
      results.trialinfo.rt(rii)               = NaN;
      
    else
    
      % find buttons within response window. we only want buttons that are one of the choice buttons.
      okok = find( buttontimes > windowstart & ...
                   buttontimes <= windowend & ...
                   ismember(buttonpressed,choicebuttons) );

      % if no buttons pressed
      if isempty(okok)
    
        results.trialinfo.change_mind(rii)       = NaN;
        results.trialinfo.button_pressed(rii)   = NaN;
        results.trialinfo.is_correct(rii)       = 0;
        results.trialinfo.rt(rii)               = NaN;
      
      % if at least one button pressed
      else

        % if there is more than one unique button pressed, mark it as "change mind"
        if length(union(buttonpressed(okok),[])) > 1
          results.trialinfo.change_mind(rii) = 1;
        else
          results.trialinfo.change_mind(rii) = 0;
        end

        % handle the case of multiple buttons
        cnt0 = length(okok);
        while 1
          if cnt0-1 < 1
            break;  % if we are at the beginning, get out
          end
          if isequal(buttonpressed{okok(cnt0-1)},buttonpressed{okok(cnt0)})
            cnt0 = cnt0 - 1;  % if cnt0-1 is the same button, use that one!
          else
            break;  % if cnt0-1 is a different button, don't go any more into the past
          end
        end
        ok = okok(cnt0);  % score this button press (out of potentially many button presses)

        % calc
        timeofpress = buttontimes(ok);    % time that key was pressed
        thekey      = buttonpressed{ok};  % the key that was pressed, e.g. '1'
        
        % record
        results.trialinfo.button_pressed(rii) = str2double(thekey);
        results.trialinfo.is_correct(rii) = double((results.trialinfo.button_pressed(rii) == a1.run_table.correct_response(ii)));
        results.trialinfo.rt(rii) = (timeofpress - stimonset) * 1000;  % RT in ms

      end
      
    end

    % increment
    trialcnt = trialcnt + 1;
    rii = rii + 1;
    ii = ii + 1;

  % handle the FIX task
  else

    % NOTE: if this is a catch trial, the fixation behavior is still valid.
    % NOTE: the FIX task is always just single-image trials.
    
    % figure out when we should start looking at the fixation luminance
    assert(trialcnt==1);
    ii2 = find( a1.run_table.block_nr == blockcnt & ...
                a1.run_table.trial_nr == trialcnt & ...
                a1.run_table.event_id == 91 );  % look for the post-task-ITI
    assert(length(ii2)==1 & ii2 >= ii);
    fixstart = a1.run_table.event_start(ii2(1));  % we are in frame time
    
    % figure out when we should end
    ii2 = find( a1.run_table.block_nr ~= blockcnt & ...
                (1:length(a1.run_table.block_nr))' > ii );  % when block_nr switches to something different
    assert(length(ii2) > 0);
    fixend = a1.run_table.event_start(ii2(1)) - 1;  % notice the -1
    
    % calc
    df = diff(a1.run_frames.fix_abs_lum(fixstart+1 : fixend+1));  % vector of luminance diffs
    seq = find(df~=0) + 1;  % vector of indices indicating frames on which the dot changed. these are indices into the extracted range.

    % process each change
    for p=1:length(seq)
    
      % record
      results.trialinfo.session_nr(rii) =   a1.run_table.session_nr(ii);
      results.trialinfo.session_type(rii) = a1.run_table.session_type(ii);
      results.trialinfo.run_nr(rii) =       a1.run_table.run_nr(ii);
      results.trialinfo.block_nr(rii) =     a1.run_table.block_nr(ii);
      results.trialinfo.crossing_nr(rii) =  a1.run_table.crossing_nr(ii);
      results.trialinfo.trial_nr(rii) =     NaN;  % NOTE!
      results.trialinfo.is_catch(rii) =     0;    % NOTE!
  
      % calc some times
      stimonset = timeframes(fixstart-1+seq(p)+1);
      windowstart = stimonset + responsewindow_fix(1)/1000;
      windowend   = stimonset + responsewindow_fix(2)/1000;

      % record
      results.trialinfo.onset_start(rii)      = fixstart-1+seq(p);
      results.trialinfo.onset_abstime(rii)    = matlabnowtime + stimonset/60/60/24;  % convert from seconds to days and then add in
      results.trialinfo.correct_response(rii) = mod(sign(df(seq(p)-1)),3);  % 1 means positive luminance change, 2 means negative luminance change
      results.trialinfo.change_mind(rii)       = NaN;  % for FIX, no changing mind

      % find buttons within response window. we only want buttons that are one of the choice buttons.
      okok = find( buttontimes > windowstart & ...
                   buttontimes <= windowend & ...
                   ismember(buttonpressed,choicebuttons) );

      % if no buttons pressed
      if isempty(okok)

        results.trialinfo.button_pressed(rii)   = NaN;
        results.trialinfo.is_correct(rii)       = 0;
        results.trialinfo.rt(rii)               = NaN;

      % if at least one button pressed
      else

        % score the first button only!
        ok = okok(1);

        % calc
        timeofpress = buttontimes(ok);    % time that key was pressed
        thekey      = buttonpressed{ok};  % the key that was pressed, e.g. '1'

        % record
        results.trialinfo.button_pressed(rii) = str2double(thekey);
        results.trialinfo.is_correct(rii)     = double((results.trialinfo.button_pressed(rii) == results.trialinfo.correct_response(rii)));
        results.trialinfo.rt(rii)             = (timeofpress - stimonset) * 1000;  % RT in ms

      end
      
      % increment
      rii = rii + 1;

    end

    % increment
    ii = ii2(1);  % this skips ahead such that on the next iteration we will trigger going to the next block
    
  end

end
warning(prevwarn.state,'MATLAB:table:RowsAddedExistingVars');

%% Report counts of buttons and triggers

fprintf('Number of pressed buttons [expected number in brackets]:\n');
for pp=1:length(userkeys)
  fprintf('%s: % 5d   [% 5d]\n', ...
          userkeys{pp}, ...
          userkeycounts(pp), ...
          sum(results.trialinfo.correct_response == str2double(userkeys{pp})));  % NOTE: str2double is NaN for rygb
end
fprintf('\n');
fprintf('Number of triggers: %d\n',length(triggertimes));
fprintf('==============================================================\n');

%% Summarize behavioral performance

prevwarn = warning('query','MATLAB:table:RowsAddedExistingVars');
warning('off','MATLAB:table:RowsAddedExistingVars');
for p=1:max(results.trialinfo.block_nr)
  ii = find(results.trialinfo.block_nr==p);
  results.summary.block_nr(p)        = p;
  results.summary.crossing_nr(p)     = results.trialinfo.crossing_nr(ii(1));
  results.summary.response_rate(p)   = sum(~isnan(results.trialinfo.button_pressed(ii))) / sum(results.trialinfo.is_catch(ii)==0) * 100;  % how many buttons are not NaN / how many non-catch trials
  results.summary.pct_correct(p)     = nanmean(results.trialinfo.is_correct(ii))*100;
  results.summary.median_rt(p)       = nanmedian(results.trialinfo.rt(ii));
  fprintf('Block % 4d, Crossing % 4d, Response Rate % 4d%%, Percent Correct % 4d%%, Median RT % 7d ms\n', ...
          results.summary.block_nr(p), ...
          results.summary.crossing_nr(p), ...
          round(results.summary.response_rate(p)), ...
          round(results.summary.pct_correct(p)), ...
          round(results.summary.median_rt(p)));
end
fprintf('==============================================================\n');
warning(prevwarn.state,'MATLAB:table:RowsAddedExistingVars');
