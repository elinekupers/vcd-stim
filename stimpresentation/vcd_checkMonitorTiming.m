function fH = vcd_checkMonitorTiming(data)

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
fH = figure(); set(gcf,'Position',[50 50 800 300],'color','w'); 
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

ax = gca;
legend(ax.Children(length(ax.Children):-1:(length(ax.Children)-4)),'Measured flip time','Median measured flip time','Expected flip time'); legend boxoff
