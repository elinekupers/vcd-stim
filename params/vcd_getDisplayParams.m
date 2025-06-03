function disp = vcd_getDisplayParams(dispname, varargin)
% VCD parameter function to get display parameters
%
%     disp = vcd_getDisplayParams(dispname, ['ppd_eccen_range_deg', <eccen>])
%
% INPUTS:
% * dispname
% * [ppd_eccen_range_deg]
%
% OUTPUTS:
% * disp          : (struct) display parameters with the following fields:
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
%
% We follow PTB coordinate convention where [0,0] is the outer left and top 
% edge of the first pixel in the upper left part of the screen. This means
% that [1,1] is the outer right and lower edge of the first pixel. 
% Screen rect will be [0 0 w_pix h_pix], where xc and yc are in between two
% pixels. For BOLD screen, this is in between pixel 960 and 961. 
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

%% Parse inputs

p = inputParser;
% MANDATORY INPUTS
p.addRequired('dispname', @(x) ismember(x,{'7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','PPROOM_EIZOFLEXSCAN','EKHOME_ASUSVE247'})) % display name to get the right display params. Choose from: 7TAS_BOLDSCREEN32, KKOFFICE_AOCQ3277, PPROOM_EIZOFLEXSCAN, 'EKHOME_ASUSVE247'
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

% Add display name
disp.name                = dispname;
disp.ppd_eccen_range_deg = 4;  % degrees

% Switch parameters according to input display name
switch dispname
    case '7TAS_BOLDSCREEN32'                    % MATLAB version 2016b, psychtoolbox version 3.0.14 December 30th 2016
        disp.w_cm        = 69.84;               % cm wide; note: beyond what subject can actaully see.
        disp.h_cm        = 39.29;               % cm high;
        disp.dist_cm     = 176+2+5.5;           % 176 cm from the mirror to glass of BOLDScreen,
                                                % 2 cm from glass to display, 5.5 from eye to mirror.
        disp.h_pix       = 1080;                % in pixels (BOLDscreen vertical, Nova1x32)
        disp.w_pix       = 1920;                % in pixels (BOLDscreen horizontal, Nova1x32)
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees (BOLDscreen vertical, Nova1x32) should be 12.2 deg
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees (BOLDscreen vertical, Nova1x32) should be 21.5 deg
        disp.refresh_hz  = 120;                 % monitor native refreshrate in Hz
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)
        disp.clut        = 0;                   % linear clut
        disp.fontsize    = 25;                  % fontsize of text
    
    case 'PPROOM_EIZOFLEXSCAN'                  % MATLAB version 2016b, psychtoolbox version 3.0.14 December 30th 2014
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
        
    case 'KKOFFICE_AOCQ3277'                    % EK stimlaptop uses MATLAB version 2018b, psychtoolbox version ??? Nov 17 2020 (git commit ef093cbf296115badddb995fa06452e34c8c7d02)
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
        
     case 'EKHOME_ASUSVE247'                    % EK stimlaptop uses MATLAB version 2018b, psychtoolbox version ??? Nov 17 2020 (git commit ef093cbf296115badddb995fa06452e34c8c7d02)
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
end

% PIX PER DEGREE CALCULATION
% pixels per degree for 4 deg eccentricity distance from center
screendistance   = disp.dist_cm * tan(deg2rad(disp.ppd_eccen_range_deg));
numpixels        = (screendistance / disp.w_cm) * disp.w_pix;
disp.ppd         = numpixels / disp.ppd_eccen_range_deg; 

% assert xc and yc are integers
assert(isint(disp.xc)); assert(isint(disp.yc));

return