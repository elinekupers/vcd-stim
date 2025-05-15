function vcd_checkMonitorTiming(data)

dt = data.timing.mfi;
xtime             = [0:dt:((length(data.timing.timeframes)-1)*dt)];
xtime             = xtime(1:end-1);
desiredFlipTime   = data.timing.frameduration* data.timing.mfi *ones(1,length(xtime));
measuredFlipTime  = diff(data.timing.timeframes); 
medianFlipDur     = median(measuredFlipTime,'omitnan');
medianFlipTime    = medianFlipDur.*ones(1,length(measuredFlipTime));

% Set y range 
yrange = medianFlipDur + [-1,1].*(max(abs([min(measuredFlipTime) max(measuredFlipTime)]))-medianFlipDur);

% Plot!
figure(); set(gcf,'Position',[50 50 800 300],'color','w'); 
plot(xtime,measuredFlipTime, 'r-','LineWidth',2); hold on; 
plot(xtime,medianFlipTime, 'k-','LineWidth',2); hold on;  
plot(xtime,desiredFlipTime, 'g', 'LineWidth',3);


% Make plot pretty
xlabel('Time (s)')
ylabel('Measured flip duration (s)')
set(gca, 'TickDir', 'out', 'FontSize', 15)
ylim(yrange)
xlim([0 max(xtime)])
title('Inspection of timeframes');

% Get nr of frames between stimuli
allowable_delay = (medianFlipDur + ([-1,1].*0.0004));
glitches =  (measuredFlipTime < allowable_delay(1) | measuredFlipTime > allowable_delay(2));

% how many interstimulus frames differed from the median?
fprintf('Nr of glitches: \t %3.0f\n', sum(glitches))


% Expand in case we have multiple keypresses at once
timekeysB = {};
for kk = 1:size(data.timeKeys,1)
    
  if iscell(data.timeKeys{kk,2})
      
    for pp=1:length(data.timeKeys{kk,2})
        
        timekeysB{end+1,1} = data.timeKeys{kk,1};
        timekeysB{end,2}   = data.timeKeys{kk,2}{pp};
    end
  else
    timekeysB(end+1,:) = data.timeKeys(kk,:);
  end
end

% now process non-badkeys
oldkey = ''; oldkeytime = -Inf;
keytimes = [];
keybuttons = {};

deltatimeBAD = 0.25;
deltatime = 0.2;

for kk=1:size(timekeysB,1)
    
  if ~isequal(timekeysB{kk,2},'absolutetimefor0') && ...
     (timekeysB{kk,1}-oldkeytime > deltatime)
 
    keytimes = [keytimes timekeysB{kk,1}];  % record
    keybuttons = [keybuttons {timekeysB{kk,2}}];  % record

    oldkey     = timekeysB{kk,2};
    oldkeytime = timekeysB{kk,1};
  end
end

% add key presses to figure
hold on; plot(repmat(keytimes,2,1),repmat(yrange,length(keytimes),1)','k')
xlim([min(keytimes), max(keytimes)])
ax = gca;
legend(ax.Children(length(ax.Children):-1:(length(ax.Children)-4)),'Measured flip time','Median measured flip time','Expected flip time','button press'); legend box off
