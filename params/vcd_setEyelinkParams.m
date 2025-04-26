function el = vcd_setEyelinkParams(el,backgroundcolour)

if ~exist('backgroundcolour','var') || isempty(backgroundcolour)
    backgroundcolour = [128 128 128];
else
    assert(numel(backgroundcolour)==3)
end

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
el.backgroundcolour                 = backgroundcolour; %same background as experiment
el.msgfontcolour                    = [0 0 0]; % black
el.calibrationtargetsize            = 2.0;     % size of calibration target as percentage of screen. 7TAS height is 1080, 2% = 22 pixels. (el default is 2.5).
el.calibrationtargetwidth           = 0.92;    % width of calibration target's border as percentage of screen. 7TAS height is 1080, 0.92% = 8 pixels. (el default is 1)
el.calibrationtargetcolor           = [0 0 0]; % black

% What about these parameters mentioned by Experiment Builder Manual???
% .outerSize -- Target Outer Size (Integer): The standard calibration and
% drift correction target is a filled circle (for peripheral detectability)
% with a central "hole" target (for accurate fixation). The disk is drawn in
% the calibration foreground color, and the hole is drawn in the calibration
% background color. The "Target Outer Size" property specifies the diameter
% of the outer disk of the default calibration target in pixels. This
% property is only available if "Use Custom Target" is not checked.
%
% .innerSize -- Target Inner Size (Integer): Diameter of the inner disk of
% the default calibration target in pixels). If hole size is 0, no central
% feature will be drawn. This property is only available if "Use Custom
% Target" is not checked.

end