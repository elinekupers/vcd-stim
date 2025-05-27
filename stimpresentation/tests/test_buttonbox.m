%% test_buttonbox.m

% Script for a quick and dirty check to see if button box/keyboard presses 
% come through and are registered by kbCheck.

vcd_startup

devices = PsychHID('Devices');
getoutearly = 0;
timekeys = {};
while 1
    [keyIsDown,secs,keyCode,~] = KbCheck(-3);  % previously -3 listen to all devices
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
    end
    
end