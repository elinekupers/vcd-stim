function el = vcd_setEyelinkParams(el)

% define Eyelink defaults (frequency, volume, duration);
% Optional: shrink the spread of the calibration/validation targets <x, y display proportion>
% if default outermost targets are not all visible in the bore.
% Default spread is 0.88, 0.83 (88% of the display horizontally and 83% vertically)
%   Eyelink('command', 'calibration_area_proportion 0.88 0.83');
%   Eyelink('command', 'validation_area_proportion 0.88 0.83');
el.cal_target_beep                  = [0 0 0]; % no beep
el.drift_correction_target_beep     = [0 0 0]; % no beep
el.targetbeep                       = false; % no beep
el.calibration_failed_beep          = [0 0 0]; % no beep
el.calibration_success_beep         = [0 0 0]; % no beep
el.drift_correction_failed_beep     = [0 0 0]; % no beep
el.drift_correction_success_beep    = [0 0 0]; % no beep
el.backgroundcolour                 = [127 127 127]; %same background as experiment
el.msgfontcolour                    = [0 0 0]; %same text as experimentel
el.calibrationtargetsize            = 1.5; % pixels ??
el.calibrationtargetwidth           = .6; % inner area pixels???
el.calibrationtargetcolor           = [0 0 0]; % black?

end