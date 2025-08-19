function disp = vcd_getDisplayParams(dispname, varargin)
% VCD parameter function to get display parameters
%
%     disp = vcd_getDisplayParams(dispname, ['ppd_eccen_range_deg', <eccen>])
%
% INPUTS:
% * dispname  : (char) display name to load display params. Choose from: 
%                 '7TAS_BOLDSCREEN32'   : 27-inch iMac + 32-inch BOLDSCREEN LCD monitor @ CMRR's 7T MRI operator room 
%                                         - iMac (Retina 5K, 27-inch, 2017) 
%                                         - macOSX High Sierra, Version 10.13.6
%                                         - Processor 4.2 GHz Intel Core i7
%                                         - Memory 16 GB 2400 MHz DDR4
%                                         - Graphics Radeon Pro 575 4096 MB
%                                         - Eyelink 1000 eyetracking system
%                 'PPROOM_EIZOFLEXSCAN' : Mac Pro tower + 22-inch EIZO FlexScan SX2462W LCD monitor @ CMRR's psychophysics lab (room 1-149B)
%                                         - Mac Pro (late 2013)
%                                         - macOSX El Capitan, Version 10.11.6
%                                         - Processor: 3.50 GHz 6-core Intel Xeon E5
%                                         - Memory: 16 GB 1866 MHz DDR3 ECC
%                                         - Graphics: AMD FirePrio D500 3072 MB
%                                         - Cambridge Research Systems Bits# Stimulus Processor
%                                         - Eyelink 1000 eyetracking system
%                 'KKOFFICE_AOCQ3277'   : For debugging purposes only. 
%                                         Used as external monitor + EK's old MacbookPro laptop (macOS Mojave, build late 2013)            
%                 'EKHOME_ASUSVE247'    : For debugging purposes only. 
%                                         Used as external monitor + EK's old MacbookPro laptop (macOS Mojave, build late 2013) 
%                 'CCNYU_VIEWPIXX3D'    : A 23.5-inch ViewPixx display @ psychophysics room of Clay Curtis' lab @ NYU
%                                         - 24-inch iMac, M1, 2021. 
%
% * [ppd_eccen_range_deg]    : (double) eccentricity range in degrees visual
%                               angle to calculate the pixels per degree.
%                               Default: 4 degrees.
%
% OUTPUTS:
% * disp        : (struct) display parameters with the following fields:
%                               
%   disp.name                : (str) name of the display
%   disp.ppd_eccen_range_deg : (double) What extent of visual angle
%                               do we want to calculate the pixel per 
%                               degrees for? (in degrees).
%   disp.w_cm                : (double) width of the display in cm
%   disp.h_cm                : (double) height of the display in cm
%   disp.dist_cm             : (double) distance from eye to display in cm
%   disp.h_pix               : (double) total height of display in pixels
%   disp.w_pix               : (double) total width of display in pixels
%   disp.h_deg               : (double) total height of display in degrees
%   disp.w_deg               : (double) total width of display in degrees
%   disp.refresh_hz          : (double) desired (native) refresh rate of 
%                               the display in Hz
%   disp.ppd                 : (double) pixels per degree for a specified 
%                               eccentricity extent in deg relative to the 
%                               center of the screen. (See
%                               "ppd_eccen_range_deg" variable).
%   disp.xc                  : (double) x-center of the display relative to 
%                               upper left corner in pixels.
%   disp.yc                  : (double) y-center of the display relative to 
%                               upper left corner in pixels.
%   disp.clut                : (double) parameter used by knkutils functon 
%                               pton.m to set linear color look up table.
%   disp.fontsize            : (double) fontsize of text (in points???).
%   disp.el_monitor_size     : (double) eyelink display params used when
%                               initializing the eyelink. Numbers refer to 
%                               monitor size in millimeters (center to left, 
%                               top, right, and bottom). Field will be 
%                               empty for 'KKOFFICE_AOCQ3277', 
%                               'EKHOME_ASUSVE247', 'CCNYU_VIEWPIXX3D' 
%                               displays as we assume no eyetracking will
%                               be used for these display environments.
%   disp.el_screen_distance  : (double) eyelink display params used when
%                               initializing the eyelink. Number refer to 
%                               distance in millimeters from eye to top and
%                               bottom edge of the monitor. Field will be 
%                               empty for 'KKOFFICE_AOCQ3277', 
%                               'EKHOME_ASUSVE247', 'CCNYU_VIEWPIXX3D' 
%                               displays as we assume no eyetracking will
%                               be used for these display environments.
%
% We follow Psychtoolbox-3 (PTB) coordinate convention where [0,0] is the
% outer left and top edge of the first pixel in the upper left part of the
% screen. This means that [1,1] is the outer right and lower edge of the
% first pixel. Screen rect will be [0 0 w_pix h_pix], where xc and yc are
% in between two pixels. For BOLD screen, this is in between pixel 960 and
% 961.
%
% Therefore, we work with even sized number of pixels (this is enforced in
% vcd_getStimulusParams.m) and try to get the center on the stimulus
% support (which lies in between pixels) in the right spot using PTB
% CenterRectOnPoint.m. This function relies on PTB RectCenter.m, which
% rounds to the closest integer for a given [x,y] center coordinate.

% For pixels per degree calcuation, we do the following: 
% 1. Take the screen distance (in cm) which is the distance from the center
%   of the screen to 4 deg eccentricity:
%      screendistance = centerdistance_cm * tan(deg2rad(eccen_deg)); 
%   For example, for the pp room: screendistance = 99.0*tan(deg2rad(4)); 
%
% 2. We calculate the number of pixels on the screen occupying this
%   distance (for vCD we take the width of the screen): 
%      numpixels = (screendistance_cm) / screenwidth_cm) * screenwidth_pixels
%   For example, for the pp room: numpixels = (screendistance / 32.5) * 1200;
%
% 3. We calculate pixels per degree by taking the number of pixels occupying 
%   this distance, divided by the used eccentricity:
%      ppd = numpixels / eccen_deg;
%   For example, for the pp room: ppd = numpixels / 4;
%
% Written by E Kupers @ UMN
% 2024/11: initial version
% 2025/06: updated ppd calculation 
% 2025/06: added eyelink params
% 

%% Parse inputs

p = inputParser;
% MANDATORY INPUTS
p.addRequired('dispname', @(x) ismember(x,{'7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','PPROOM_EIZOFLEXSCAN','EKHOME_ASUSVE247','CCNYU_VIEWPIXX3D'}))
% OPTIONAL INPUTS
p.addParameter('ppd_eccen_range_deg' , 4, @isnumeric); % (degrees). 

% Parse inputs
p.parse(dispname, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    eval(sprintf('%s = p.Results.%s;', rename_me{ff},rename_me{ff}));
end
clear rename_me ff p

%% Create display parameter struct
disp      = struct(); 

% Add display name & eccentricity range to calculate pixels per degree
disp.name                = dispname;
disp.ppd_eccen_range_deg = ppd_eccen_range_deg;  % degrees

% Switch parameters according to input display name
switch dispname
    case '7TAS_BOLDSCREEN32'                    
        % 21-inch iMac MATLAB version 2016b and R2017b, psychtoolbox version 3.0.14 -- Flavor: Beta. SVN Revision 8301. December 30th 2016
        disp.w_cm        = 69.84;               % cm wide; note: beyond what subject can actaully see.
        disp.h_cm        = 39.29;               % cm high;
        disp.dist_cm     = 176+2+5.5;           % total of 183.5 cm: 176 cm from the mirror to glass of BOLDScreen +
                                                % 2 cm from glass to display + 5.5 from eye to mirror.
        disp.h_pix       = 1080;                % in pixels (BOLDscreen vertical, Nova1x32)
        disp.w_pix       = 1920;                % in pixels (BOLDscreen horizontal, Nova1x32)
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees (BOLDscreen vertical, Nova1x32) should be 12.2 deg
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees (BOLDscreen vertical, Nova1x32) should be 21.5 deg
        disp.refresh_hz  = 120;                 % monitor native refreshrate in Hz
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)
        disp.clut        = 0;                   % linear clut
        disp.fontsize    = 25;                  % fontsize of text
        % EYELINK DISPLAY PARAMS
        disp.el_monitor_size    = [-349.2, 196.45, 349.2, -196.45]; % monitor size in millimeters (center to left, top, right, and bottom). Numbers come from [39.29 cm height, 69.84 cm width] --> [392.9 mm height, 698.4 mm width] * 0.5.
        disp.el_screen_distance = [1845 1845];  % distance in millimeters from eye to top and bottom edge of the monitor. Given the 183.5 cm viewing distance, this is calculated as:  sqrt(183.5^2+(32.5/2)^2)*10 and then rounded to nearest integer.    
    case 'PPROOM_EIZOFLEXSCAN'                  
        % Mac tower has MATLAB version 2016b, psychtoolbox version 3.0.14 December 30th 2014
        disp.w_cm        = 52;                  % width in cm
        disp.h_cm        = 32.5;                % height in cm
        disp.dist_cm     = 99.0;                % eye to screen viewing distance in cm
        disp.h_pix       = 1200;                % height in pixels 
        disp.w_pix       = 1920;                % width in pixels 
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees 
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees
        disp.refresh_hz  = 60;                  % desired (native) refreshrate of monitor in Hz
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)
        disp.clut        = 0;                   % linear clut: amazingly, no lookup table needed!! when using user3 - gamma 2.2
        disp.fontsize    = 18;                  % fontsize of text
        % EYELINK DISPLAY PARAMS
        disp.el_monitor_size    = [-260.0, 162.5, 260.0, -162.5]; % monitor size in millimeters (center to left, top, right, and bottom). Numbers come from [32.5 cm height, 52 cm width] --> [325 cm height, 520 cm width] * 0.5
        disp.el_screen_distance = [1003 1003];  % distance in millimeters from eye to top and bottom edge of the monitor. Given the 99 cm viewing distance, this is calculated as:  sqrt(99^2+(32.5/2)^2)*10 and then rounded to nearest integer.
    case 'KKOFFICE_AOCQ3277'                    
        % Used as external monitor + EK stimulus laptop
        % MATLAB version 2018b, psychtoolbox version ??? Nov 17 2020 (git commit ef093cbf296115badddb995fa06452e34c8c7d02)
        disp.w_cm        = 71;                  % cm wide; 
        disp.h_cm        = 40;                  % cm high;
        disp.dist_cm     = 100;                 % from eye to display.
        disp.h_pix       = 1440;                % height in pixels 
        disp.w_pix       = 2560;                % width in pixels
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees 
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees
        disp.refresh_hz  = 60;                  % desired (native) refreshrate of monitor in Hz
        disp.ppd         = disp.h_pix/disp.h_deg; % pixels per degree
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)
        disp.clut        = 0;                   % linear clut
        disp.fontsize    = 12;                  % fontsize of text
        % EYELINK DISPLAY PARAMS
        disp.el_monitor_size    = [];           % assume no eyetracking
        disp.el_screen_distance = [];           % assume no eyetracking
     case 'EKHOME_ASUSVE247'                    
        % Used as external monitor + EK stimulus laptop
        % MATLAB version 2018b, Psychtoolbox-3 version ??? Nov 17 2020 (git commit ef093cbf296115badddb995fa06452e34c8c7d02)
        disp.w_cm        = 71;                  % cm wide; 
        disp.h_cm        = 40;                  % cm high;
        disp.dist_cm     = 100;                 % from eye to display.
        disp.h_pix       = 1080;                % height in pixels 
        disp.w_pix       = 1920;                % width in pixels
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees 
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees
        disp.refresh_hz  = 60;                  % desired (native) refreshrate of monitor in Hz
        disp.ppd         = disp.h_pix/disp.h_deg; % pixels per degree
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)       
        disp.clut        = 0;                   % linear clut
        disp.fontsize    = 12;                  % fontsize of text
        % EYELINK DISPLAY PARAMS
        disp.el_monitor_size    = [];           % assume no eyetracking
        disp.el_screen_distance = [];           % assume no eyetracking
    case 'CCNYU_VIEWPIXX3D'                     
        % A 23.5-inch ViewPixx display + iMac, 24in, M1, 2021. Psychtoolbox-3 version 3.0.19 on MATLAB 9.13
        disp.w_cm        = 30.1625;             % cm wide; (11.875in)
        disp.h_cm        = 53.34;               % cm high; (21 in)
        disp.dist_cm     = 100;                 % from eye to display.
        disp.h_pix       = 1080;                % height in pixels 
        disp.w_pix       = 1920;                % width in pixels
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees 
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees
        disp.refresh_hz  = 120;                 % desired (native) refreshrate of monitor in Hz
        disp.ppd         = disp.h_pix/disp.h_deg; % pixels per degree
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)       
        disp.clut        = 0;                   % assume linear clut (Check!!)
        disp.fontsize    = 12;                  % fontsize of text
        % EYELINK DISPLAY PARAMS
        disp.el_monitor_size    = [];           % assume no eyetracking
        disp.el_screen_distance = [];           % assume no eyetracking
end

% PIX PER DEGREE CALCULATION
% pixels per degree for 4 deg eccentricity distance from center
screendistance   = disp.dist_cm * tan(deg2rad(disp.ppd_eccen_range_deg));
numpixels        = (screendistance / disp.w_cm) * disp.w_pix;
disp.ppd         = numpixels / disp.ppd_eccen_range_deg; 

% assert xc and yc are integers
assert(isint(disp.xc)); assert(isint(disp.yc));

return