function results = vcdeyetrackingpreprocessing(filename,behfilename,behresults,wantfigures,blinkpad,maxpolydeg);

% function results = vcdeyetrackingpreprocessing(filename,behfilename,behresults,wantfigures,blinkpad,maxpolydeg);
%
% <filename> is the .edf file for one run
% <behfilename> is the corresponding behavioral .mat file for that run
% <behresults> is the results of calling vcdbehavioralanalysis.m on <behfilename>
% <wantfigures> (optional) is
%   1 means to create figure windows
%   X means write out figures with X as the directory+filename prefix
%   Default: 1.
% <blinkpad> (optional) is [A B] with the number of milliseconds before and 
%   after each blink event to excise. Default: [200 200].
% <maxpolydeg> (optional) is the maximum polynomial degree to use for detrending.
%   Should be a non-negative integer. Default: round(L/2) where L is 
%   behresults.totaldur/60 (i.e. the total duration in minutes).
%
% Perform basic preprocessing of the eyetracking data. We perform the following
% steps: (1) convert the .edf file, (2) record a few important messages,
% (3) deal with synchronization and time-cropping, (4) convert to degrees
% of visual angle, (5) remove blinks, and (6) detrend and median-center both
% the x- and y-coordinates. We also include a number of sanity checks and asserts.
%
% Blinks are removed according to Eyelink's detected eblinks.
%
% We return <results> as a struct with:
%   <messages> - 1 x 4 cell vector with several specific messages from the .edf file.
%                These are: CAL CALIBRATION, CAL VALIDATION, GAZE_COORDS, MODE RECORD.
%                We extract the last instances found for the 1st and 2nd ones.
%                We check that there is exactly one instance of the 3rd and 4rd ones.
%   <samplingrate> - The sampling rate of the eyetracking data (in Hz) as reported in
%                    the .edf messages. Should be equal to 1000.
%   <synctimes> - The times corresponding to the SYNCTIME messages. These are whole
%                 numbers of milliseconds as reported by the Eyelink clock.
%                 In the case of partial data, only the first SYNCTIME will be found.
%   <eyedata> - 4 x samples with the eyetracking data. The rows are:
%               time - The time corresponding to each sample, where 0 corresponds
%                      to the onset of the first stimulus frame. The sampling rate
%                      is exactly 1000 Hz. The first sample is the one that is right 
%                      before the onset of the first stimulus frame. The last sample 
%                      is the one that is right after completion of the data 
%                      (according to behresults.totaldur) if it exists, but if not,
%                      the last sample is whatever the last recorded sample is.
%               x    - The x-coordinate of gaze position in degrees of visual angle
%               y    - The y-coordinate of gaze position in degrees of visual angle
%               pupilsize - Pupil size (in units provided by Eyelink)
%               Note that (0,0) corresponds to the center of the screen.
%               NaNs will exist in x, y, and pupilsize for time periods 
%               corresponding to removal of blinks.
%   <maxpolydeg> is the maximum polynomial degree used for detrending.
%   <targeterrs> is 1 x 5 (in the order of the event_id (991-995)) with the
%                   mean Euclidean distance from the target
%   <fixationradius> is the 75th percentile of Euclidean distances from the origin.
%                    A circle centered on the origin with radius <fixationradius>
%                    should encompass 75% of the data.
%
% We also generate several figures:
% 'timeseries' - This shows (1) raw data, (2) data after dva conversion, 
%                removal of blinks, and time synchronization, (3) data after
%                detrending and median-centering, (4) a version of 3 where
%                we zoom in to show the data from the beginning of the run
%                until the second task cue.
% 'targets' - This illustrates the eyetracking targets and associated 
%             eyetracking results. We plot samples obtained within 
%             300 to 1300 ms after target onset.
% 'total2dhist' - This shows a 2D histogram of gaze positions from the entire
%                 eyetracking run (after synchronization and time-cropping). 
%                 We use 0.25-deg bins. Various display elements of the VCD experiment
%                 are illustrated. We show the histogram with linear and log counts.
%                 The circle associated with <fixationradius> is plotted in yellow.
% 'locked' - We show pupil size and gaze position time-locked to stimuli and task cues.
%            Pupil size traces are extracted from 1000 ms before and 5000 ms after
%            onset of stimuli and task cues. Each trace is de-meaned, and the mean and
%            SE across traces are shown. Gaze positions are visualized using 2D 
%            histograms (0.25-deg bins).
%
% History:
% - 2025/07/21 - add results.fixationradius
% - 2025/07/21 - add results.targeterrs
% - 2025/07/21 - wantfigures==0 is no longer supported

%% Internal constants

targetwindow = [500 2000];    % number of milliseconds after target onset to extract from
onsetwindow = [-1000 5000];   % number of milliseconds for window around stim/taskcue
etloc = 4;                    % number of dva where eye tracking targets are placed (away from center)

%% Setup

% inputs
if ~exist('wantfigures','var') || isempty(wantfigures)
  wantfigures = 1;
end
if ~exist('blinkpad','var') || isempty(blinkpad)
  blinkpad = [200 200];
end
if ~exist('maxpolydeg','var') || isempty(maxpolydeg)
  maxpolydeg = [];
end
if isempty(maxpolydeg)
  maxpolydeg = round(behresults.totaldur/60/2);
end

% init
clear results;

% load the behavioral data file
a1 = load(behfilename);

%% Convert and load eyetracking data. Also, deal with timing and sanity checks.

% use edf2asc to convert the .edf file. note that this edf2asc call 
% includes metadata including events/messages.
% DO NOT REPORT TO COMMAND WINDOW.
t0 = maketempdir;
unix_wrapper(sprintf('edf2asc -y -p "%s" -miss NaN "%s"',t0,filename),0,0);  % -y = overwrite; -vel ?

% load in the .asc file and clean up
[~,tmpfile] = fileparts(stripext(filename));
tempascfile = fullfile(t0,[tmpfile '.asc']);
b1 = read_eyelink_asc_v3(tempascfile);
rmdirquiet(t0);

% NOTE: the edf2asc and read_eyelink_asc_v3 steps are slow. It could be cached/saved.

% extract some critical messages
%     {'MSG'}    {'3681821'}    {'!CAL CALIBRATION HV5 L LEFT GOOD  '                                                       }
%     {'MSG'}    {'3696191'}    {'!CAL VALIDATION HV5 L LEFT GOOD ERROR 0.54 avg. 1.07 max OFFSET 0.43 deg. -5.9,18.9 pix. '}
%     {'MSG'}    {'3701030'}    {'GAZE_COORDS 0.00 0.00 1919.00 1199.00 '                                                   }
%     {'MSG'}    {'3701031'}    {'!MODE RECORD CR 1000 2 0 L '                                                              }
%     {'MSG'}    {'3707917'}    {'SYNCTIME '                                                                                }
%     {'MSG'}    {'4145359'}    {'SYNCTIME '                                                                                }
results.messages = {};
  % get the last one of CAL CALIBRATION:
ok = b1.msg(cellfun(@(x) ~isempty(regexp(x,'!CAL CALIBRATION')),b1.msg(:,3)),3);
if isempty(ok)
  warning('*** No CAL CALIBRATION found!! ***');
  results.messages{1} = [];
else
  results.messages{1} = ok{end};
end
  % get the last one of CAL VALIDATION:
ok = b1.msg(cellfun(@(x) ~isempty(regexp(x,'!CAL VALIDATION')),b1.msg(:,3)),3);
if isempty(ok)
  warning('*** No CAL VALIDATION found!! ***');
  results.messages{2} = [];
else
  results.messages{2} = ok{end};
end
  % get the GAZE_COORDS
ok = b1.msg(cellfun(@(x) ~isempty(regexp(x,'GAZE_COORDS')),b1.msg(:,3)),3); assert(length(ok)==1);
results.messages{3} = ok{1};
  % get the MODE RECORD
ok = b1.msg(cellfun(@(x) ~isempty(regexp(x,'!MODE RECORD')),b1.msg(:,3)),3); assert(length(ok)==1);
results.messages{4} = ok{1};

% check that the screen coordinates are exactly as expected
coords = sscanf(results.messages{3}(length('GAZE_COORDS')+1:end),'%f');
assert(coords(1)==0);
assert(coords(2)==0);
assert(coords(3)==a1.params.disp.w_pix-1);
assert(coords(4)==a1.params.disp.h_pix-1);

% get the sampling rate and check that it is as expected
results.samplingrate = str2double(subscript(strsplit(results.messages{4},' '),4,1));
assert(results.samplingrate == 1000);

% extract the SYNCTIMEs.
% The first SYNCTIME occurs immediately after the detection of the trigger, and
%   corresponds to 'trigger' in timekeys. (behresults.mristarttime)
% The second SYNCTIME occurs immediately after the completion of the last frame, and
%   corresponds to 'DONE' in timekeys. (behresults.donetime)
% With partial data, we may not get the final 'DONE', so be careful!
ok = b1.msg(cellfun(@(x) ~isempty(regexp(x,'SYNCTIME')),b1.msg(:,3)),2); assert(length(ok)>=1);
if length(ok)==1
  warning('*** Did not get the second SYNCTIME. This may be partial data?? BEWARE. ***');
else
  assert(length(ok)==2);
end
results.synctimes = cellfun(@str2double,ok);

% cross-check the eyetracking data duration (according to synctimes) against the
% duration according to the stimulus computer's logs. assuming that this passes,
% we are going to assume that we can simply synchronize the first synctime with
% the 'trigger' in timekeys.
if length(results.synctimes) == 2
  assert(abs(diff(results.synctimes)/1000 - (behresults.donetime-behresults.mristarttime)) < 20/1000);
end

% get the eyetracking data (gaze and pupil data). at this point:
% x and y contain NaNs. pupilsize does not (it goes to 0).
% x goes from left to right (screen pixels).
% y goes from top to bottom (screen pixels).
results.eyedata = b1.dat;  % [4 (time/x/y/pupilsize) x samples]
origeyedata = results.eyedata;  % MARK: raw eyetracking data

% check that the time row consists of increasing integers (presumably reflecting milliseconds)
assert(isequal(results.eyedata(1,:),(1:size(results.eyedata,2)) + (results.eyedata(1,1)-1)));

%% Check for catastrophic failure

% this can happen if there is never an eyeball tracked (i think??)
if size(results.eyedata,1) ~= 4
  warning('*** Did not find 4 rows in eyedata!! Aborting processing and no figures will be made. ***');
  return;
end  

%% Do a bunch of preprocessing (convert to dva, remove blinks, synchronize and time-crop, detrend, median-center)

% convert to degrees of visual angle
cntr = [(0+(a1.params.disp.w_pix-1))/2 (0+(a1.params.disp.h_pix-1))/2];  % [width height]
wfun = @(x) (x - cntr(1)) / a1.params.disp.ppd;      % convert x in screen coordinates to dva
hfun = @(x) (cntr(2) - x) / a1.params.disp.ppd;      % convert y in screen coordinates to dva
results.eyedata(2,:) = wfun(results.eyedata(2,:));
results.eyedata(3,:) = hfun(results.eyedata(3,:));
xrng = wfun(coords([1 3])'+[-.5 .5]);  % a nice x-range (full screen)
yrng = hfun(coords([4 2])'+[.5 -.5]);  % a nice y-range (full screen)

% remove blinks by setting data to NaN
badix = [];
for p=1:length(b1.eblink.stime)
  badix = [badix b1.eblink.stime(p)-blinkpad(1):b1.eblink.etime(p)+blinkpad(2)];
end
ii = ismember(results.eyedata(1,:),badix);
results.eyedata(2:4,ii) = NaN;

% convert time in eyedata to the stimulus computer's time.
% after this step, 0 corresponds to the time of the first stimulus frame,
% and we are in units of seconds. note that we are making the strong
% assumption that the two clocks (eyelink and stimulus computer) keep time precisely.
results.eyedata(1,:) = (results.eyedata(1,:) - results.synctimes(1))/1000 + behresults.mristarttime;

% time-crop the data according to behresults.totaldur (note: partial data may result in shorter totaldur)
minix = firstel(find(results.eyedata(1,:)>0))-1;                  % first sample is right before first frame
temp = find(results.eyedata(1,:)>behresults.totaldur);
if ~isempty(temp)
  maxix = firstel(temp);                 % last sample is right after completion of last frame
else
  maxix = length(results.eyedata(1,:));  % last sample is just the last frame recorded
end
results.eyedata = results.eyedata(:,minix:maxix);

% check that the duration of the eyetracking is matched to the behavioral data
assert(abs(((results.eyedata(1,end) - results.eyedata(1,1)) + 1/1000) - behresults.totaldur) < 50/1000, ...
       'eyetracking data duration is mismatched to behresults.totaldur');

% figure out where we have baddata.
% check that the data are finite everywhere else.
baddata = isnan(results.eyedata(4,:));
assert(all(isfinite(flatten(results.eyedata(2:4,~baddata)))));

% detrend x and y by fitting low-order polynomials and subtracting
origeyedata2 = results.eyedata;  % MARK: prior to detrending
%polymatrix = constructpolynomialmatrix(size(results.eyedata,2),0:maxpolydeg);  % samples x polys
temp = linspace(0,1,size(results.eyedata,2));
polymatrix = [];
for p=0:maxpolydeg
  polymatrix(:,p+1) = temp .^ p;
end
X = polymatrix(~baddata,:);  % samples x polys
h = olsmatrix(X)*results.eyedata(2:3,~baddata)';  % weights x different-timeseries
% h = [];
% h(:,1) = fitl1line(X,results.eyedata(2,~baddata)')';
% h(:,2) = fitl1line(X,results.eyedata(3,~baddata)')';
thefit = polymatrix*h;  % samples x different-timeseries
results.eyedata(2:3,:) = results.eyedata(2:3,:) - thefit';

% median-center (after this, the medoid is the origin)
results.eyedata(2:3,:) = results.eyedata(2:3,:) - nanmedian(results.eyedata(2:3,:),2);

% record
results.maxpolydeg = maxpolydeg;

%% Prep for visualizations

% NO LONGER SUPPORTED:
% % if the user doesn't want figures, we are done
% if isequal(wantfigures,0)
%   return;
% end
if isequal(wantfigures,0)
  error('wantfigures==0 no longer supported');
end

% calc
wantfigwin = isequal(wantfigures,1);  % 1 means make figure windows; 0 means write figures to disk
if ~wantfigwin
  [figuredir,fileprefix] = stripfile(wantfigures);
end

%% Visualization of time-series

figureprep([100 100 1600 850],wantfigwin);
subplotresize(4,1);

% raw data
subplot(4,1,1); hold on;
plot(origeyedata(1,:),origeyedata(2,:),'r-');
plot(origeyedata(1,:),origeyedata(3,:),'b-');
plot(origeyedata(1,:),origeyedata(4,:),'k-');
xlim(origeyedata(1,[1 end]));
title('raw data (x=red, y=blue, pupil=black)');
xlabel('Time (Eyelink milliseconds)');
ylabel('Raw units');

% basic stuff
subplot(4,1,2); hold on;
yyaxis left;
plot(origeyedata2(1,:),origeyedata2(2,:),'r-');
plot(origeyedata2(1,:),origeyedata2(3,:),'b-');
ylabel('Degrees');
xlim([0-5 behresults.totaldur+5]);
%mx = 1.5*prctile(abs(flatten(origeyedata2(2:3,:))),99);
mx = 4;
ylim([-mx mx]);
uistack(straightline(0,'h','k-'),'bottom');
yyaxis right;
plot(origeyedata2(1,:),origeyedata2(4,:),'k-');
ylabel('Raw units');
xlabel('Time (s)');
title('dva, remove blinks, synchronization and time-crop');

% detrend
for spi=3:4

  subplot(4,1,spi); hold on;
  yyaxis left;
  plot(results.eyedata(1,:),results.eyedata(2,:),'r-');
  plot(results.eyedata(1,:),results.eyedata(3,:),'b-');
  ylabel('Degrees');
  xlim([0-5 behresults.totaldur+5]);
  ylim([-mx mx]);
  uistack(straightline(0,'h','k-'),'bottom');
  yyaxis right;
  plot(results.eyedata(1,:),results.eyedata(4,:),'k-');
  ylabel('Raw units');
  xlabel('Time (s)');
  title('detrend (fix=green, target=cyan, pupil=grays, taskcue=black, magneta=stim1/2 (stim2 is thick)');

  ii = find(ismember(a1.run_table.event_id,991:995));  % et_sac
  straightline(behresults.timeframes(a1.run_table.event_start(ii)+1),'v','c-');

  ii = find(ismember(a1.run_table.event_id,990));  % et_fix
  straightline(behresults.timeframes(a1.run_table.event_start(ii)+1),'v','g-');

  ii = find(ismember(a1.run_table.event_id,996));  % et_pupil_black
  set(straightline(behresults.timeframes(a1.run_table.event_start(ii)+1),'v','k-'),'Color',[.4 .4 .4]);

  ii = find(ismember(a1.run_table.event_id,997));  % et_pupil_white
  set(straightline(behresults.timeframes(a1.run_table.event_start(ii)+1),'v','k-'),'Color',[.8 .8 .8]);

  ii = find(ismember(a1.run_table.event_id,94));  % stim1
  straightline(behresults.timeframes(a1.run_table.event_start(ii)+1),'v','m-');

  ii = find(ismember(a1.run_table.event_id,95));  % stim2
  set(straightline(behresults.timeframes(a1.run_table.event_start(ii)+1),'v','m-'),'LineWidth',2);

  ii = find(ismember(a1.run_table.event_id,90));  % task-cue
  straightline(behresults.timeframes(a1.run_table.event_start(ii)+1),'v','k-');
  zoomend = ii(2);  % for the zoomed-in figure, let's stop at the second task cue
  
end

% zoom in
subplot(4,1,4);
xlim([0-1 behresults.timeframes(a1.run_table.event_start(zoomend)+1)+1]);
title('detrend (zoomed-in)');

% write to disk
if ~wantfigwin
  figurewrite([fileprefix '_timeseries'],[],[],figuredir);
end

%% Visualization of eyetracking targets

figureprep([100 100 500 500],wantfigwin); hold on;
ii = find(ismember(a1.run_table.event_id,991:995));  % et_sac (the targets)
assert(length(ii)==5);
targetxxs = [0 -etloc etloc 0 0];  % central, left, right, up, down.
targetyys = [0  0 0 etloc -etloc];
colors0 = cool(5);
results.targeterrs = [];
for p=1:length(ii)  % for each target, in the order they appear
  ix = a1.run_table.event_id(ii(p)) - 990;  % 1-5 indicating which target number
  
  targettime = behresults.timeframes(a1.run_table.event_start(ii(p))+1);
  windowstart = targettime + targetwindow(1)/1000;
  windowend   = targettime + targetwindow(2)/1000;
  
  tt = results.eyedata(1,:) >= windowstart & results.eyedata(1,:) < windowend;
  h = scatter(results.eyedata(2,tt),results.eyedata(3,tt),'r.');
  set(h,'CData',colors0(ix,:));

  scatter(targetxxs(ix),targetyys(ix),'ko','filled');
  
  % calculate error as mean Euclidean distance from target (mean over time samples)
  results.targeterrs(ix) = nanmean(sqrt(sum((results.eyedata(2:3,tt) - [targetxxs(ix); targetyys(ix)]).^2,1)));  % NOTE: nanmean
end
axis(2*[-etloc etloc -etloc etloc]);
axis square;
xlabel('X-position (deg)');
ylabel('Y-position (deg)');
title('Eyetracking targets');
colormap(colors0);
caxis([.5 5.5]);
h = colorbar;
set(h,'Ticks',1:5);
set(h,'TickLabels',{'center' 'left' 'right' 'up' 'down'});
if ~wantfigwin
  figurewrite([fileprefix '_targets'],[],[],figuredir);
end

%% Visualization of total 2D histogram

% quantify fixation extent
dists = sqrt(sum(results.eyedata(2:3,:).^2,1));    % 1 x samples with Euclidean distance from origin
results.fixationradius = prctile(dists,75);        % 75% of the time, the distances are this value or less

% proceed to figure
figureprep([100 100 900 900],wantfigwin);
for spi=1:2
  subplot(2,1,spi); hold on;
  xlim(xrng);
  ylim(yrng);
  xbins = [fliplr(-.25:-.25:xrng(1)) 0:.25:xrng(2)];
  ybins = [fliplr(-.25:-.25:yrng(1)) 0:.25:yrng(2)];
  n = hist2d(results.eyedata(2,:),results.eyedata(3,:),xbins,ybins);
  if spi==2
    n = log10(n);
  end
  imagesc(xbins,ybins,n);
  colormap(hot);
  caxis([0 max(n(:))+eps]);
  drawellipse(0,0,0,results.fixationradius,results.fixationradius,[],[],'y-');
  %%%% BUT, LET'S TRY QUICK REGRESS PUPIL? SEE HOW WELL.
  h = colorbar;
  if spi==2
    set(h,'Ticks',0:max(get(h,'Ticks')));
    set(h,'TickLabels',10.^get(h,'Ticks'));
  end
  plotrectangle([yrng xrng],'r-');
  axis equal tight;
  drawellipse(0, 0,0,4,4,[],[],'b-');
  drawellipse(4, 0,0,2,2,[],[],'g-');
  drawellipse(-4,0,0,2,2,[],[],'g-');
  plotrectangle([-4.2 4.2 -4.2 4.2],'w-');
  scatter(0,0,'kx');
  xlabel('X-position (deg)');
  ylabel('Y-position (deg)');
  for p=1:length(targetxxs)
    scatter(targetxxs(p),targetyys(p),'co');
  end
  title('display=red, 4degecc=blue, 4degstimaperture=green, targets=cyan, nssquare=white, fixextent=yellow');
end
subplotresize(2,1,.9,.85);
if ~wantfigwin
  figurewrite([fileprefix '_total2dhist'],[],[],figuredir);
end

%% Visualize pupil size and x,y as a function of stimuli and task-cue

% design:
%
% 1 is stim1/2 non-FIX; 2 is stim1/2 FIX; 3 is like 1 but catch; 4 is like 2 but catch; 5 is taskcue
%
% pupil:  1 2  3 4   5
% x,y:    6 7  8 9  10
figureprep([100 100 1500 700],wantfigwin);
taskclassfun = {@(x) x~=1 @(x) x==1 @(x) x~=1 @(x) x==1 [] ...
                @(x) x~=1 @(x) x==1 @(x) x~=1 @(x) x==1 []};
catchfun =     {@(x) x==0 @(x) x==0 @(x) x==1 @(x) x==1 [] ...
                @(x) x==0 @(x) x==0 @(x) x==1 @(x) x==1 []};
typ = [1 1 1 1 1 2 2 2 2 2];       % 1=pupil, 2=x,y
stimtask = [1 1 1 1 2 1 1 1 1 2];  % 1=stim, 2=task
titles = {'Stimulus-locked (regular tasks)' 'Stimulus-locked (fixation task)' ...
          'Catch trial (regular tasks)'     'Catch trial (fixation task)' ...
          'Task-cue-locked'};
for spi=1:10
  subplot(2,5,spi); hold on;
  switch stimtask(spi)
  case 1
    ii = find(ismember(a1.run_table.event_id,[94 95]) & taskclassfun{spi}(a1.run_table.task_class) & catchfun{spi}(a1.run_table.is_catch));
  case 2
    ii = find(ismember(a1.run_table.event_id,[90]));
  end
  tt = behresults.timeframes(a1.run_table.event_start(ii)+1);
  allx = [];
  ally = [];
  for p=1:length(tt)
    switch typ(spi)
    case 1
      finaltime0 = (onsetwindow(1):onsetwindow(2))/1000;
      tix = find(results.eyedata(1,:) >= tt(p)+(onsetwindow(1)-100)/1000 & results.eyedata(1,:) < tt(p)+(onsetwindow(2)+100)/1000);
      time0 = round(1000*(results.eyedata(1,tix)-tt(p)))/1000;  % round to nearest millisecond
      ok = ismember(time0,finaltime0);
      ally(end+1,1:sum(ok)) = zeromean(results.eyedata(4,tix(ok)));  % tricky business: there might be partial data
      ally(end,sum(ok)+1:end) = NaN;
    case 2
      tix = results.eyedata(1,:) >= tt(p)+onsetwindow(1)/1000 & results.eyedata(1,:) < tt(p)+onsetwindow(2)/1000;
      allx = [allx results.eyedata(2,tix)];
      ally = [ally results.eyedata(3,tix)];
    end
  end
  switch typ(spi)
  case 1
    if ~isempty(ally)
      plot(finaltime0,ally,'-');
      errorbar3(finaltime0,nanmean(ally,1),nanstd(ally,[],1)./sqrt(sum(~isnan(ally),1)),'v',[1 .7 .7]);
      plot(finaltime0,nanmean(ally,1),'r-','LineWidth',2);
    end
  case 2
    xbins = [fliplr(-.25:-.25:-3) 0:.25:3];
    ybins = [fliplr(-.25:-.25:-3) 0:.25:3];
    n = hist2d(allx,ally,xbins,ybins);
    imagesc(xbins,ybins,n);
    colormap(hot);
    caxis([0 max(n(:))+eps]);
    axis equal tight;
    scatter(0,0,'kx');
  end
  if spi==1
    ylim0 = ylim;
  end
  if spi>=1 && spi<=5
    xlim(onsetwindow/1000);
    ylim(ylim0);
    xlabel('Time relative to onset (s)');
    ylabel('Pupil size');
  end
  if spi>=6 && spi<=10
    xlim([-3 3]);
    ylim([-3 3]);
    xlabel('X-position (deg)');
    ylabel('Y-position (deg)');
  end
  if spi<=5
    title(titles{spi});
  end
  switch typ(spi)
  case 1
    straightline(0,'v','k-');
    straightline(0,'h','k-');
  case 2
  end
end
if ~wantfigwin
  figurewrite([fileprefix '_locked'],[],[],figuredir);
end

%%%%%%%%%%%%%%%%%

% JUNK:
%
% smooth and downsample?- perhaps smooth and downsample to 50ms?  or maybe the high frequency noise is low anyway
%   ff = constructbutterfilter1D(length(timenew),bpfilt * (length(timenew) * 1/resamplerate));
%   timenew = time(1):1/resamplerate:time(end);
%   xct = interp1(time,xc,timenew,'cubic');
% % band-pass filter [who cares?]
% ff = constructbutterfilter1D(length(timenew),bpfilt * (length(timenew) * 1/resamplerate));
% xct2 = tsfilter(xct,ff);
% yct2 = tsfilter(yct,ff);
