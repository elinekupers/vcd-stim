function disp = vcd_getDisplayParams(dispname)

if ~exist('dispname','var') || isempty(dispname)
    dispname = '7TAS_BOLDSCREEN32';
end

disp = struct(); 
disp.name = dispname;

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

eccen_range = 4; % deg (what extent in dva do we want to calculate the pixel per degrees for?)

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
        
        tan(deg2rad(eccen_range));
        
        
        disp.ppd         = disp.h_pix/disp.h_deg; % pixels per degree
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
        disp.ppd         = disp.h_pix/disp.h_deg; % pixels per degree
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

% assert xc and yc are integers
assert(isint(disp.xc)); assert(isint(disp.yc));

return