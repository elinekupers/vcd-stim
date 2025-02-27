function disp = vcd_getDisplayParams(dispname)

if ~exist('dispname','var') || isempty(dispname)
    dispname = '7TASBOLDSCREEN32';
end

disp = struct(); 
disp.name = dispname;

switch dispname
    case '7TAS_BOLDSCREEN32'                     % MATLAB version 2016b, psychtoolbox version 3.0.14 December 30th 2016
        disp.w_cm        = 69.84;               % cm wide; beyond what subject can actaully see.
        disp.h_cm        = 39.29;               % cm high;
        disp.dist_cm     = 176+2+5.5;           % 176.5 cm from the mirror to glass of BOLDScreen,
%         disp.dist_cm     = 176.5+2+5.5;       % 176.5 cm from the mirror to glass of BOLDScreen, 
                                                % 2 cm from glass to display, 5.5 from eye to mirror (check!).
        disp.h_pix       = 1080;                % in pixels (BOLDscreen vertical, Nova1x32), 800 pixel x 800 pixel square is what people can see.
        disp.w_pix       = 1920;                % in pixels (BOLDscreen horizontal, Nova1x32), 1800 pixels is max extend
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees (BOLDscreen vertical, Nova1x32) should be 12.70 deg
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees (BOLDscreen vertical, Nova1x32)
        disp.refresh_hz  = 120;                 % refreshrate in Hz
        disp.ppd         = disp.h_pix/disp.h_deg;
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)
        disp.clut        = 0;                   % linear clut
        
    case 'KKOFFICE_AOCQ3277'                    % EK stimlaptop uses MATLAB version 2018b, psychtoolbox version ??? Nov 17 2020 (git commit ef093cbf296115badddb995fa06452e34c8c7d02)
        disp.w_cm        = 71;                  % cm wide; beyond what subject can actaully see.
        disp.h_cm        = 40;                  % cm high;
        disp.dist_cm     = 50;                  % from eye to display.
        disp.h_pix       = 1440;                % height in pixels 
        disp.w_pix       = 2560;                % width in pixels
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees 
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees
        disp.refresh_hz  = 60;                  % refreshrate in Hz
        disp.ppd         = disp.h_pix/disp.h_deg;
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)
        disp.clut        = 0;                   % linear clut
        
     case 'EKHOME_ASUSVE247'                    % EK stimlaptop uses MATLAB version 2018b, psychtoolbox version ??? Nov 17 2020 (git commit ef093cbf296115badddb995fa06452e34c8c7d02)
        disp.w_cm        = 71;                  % cm wide; beyond what subject can actaully see.
        disp.h_cm        = 40;                  % cm high;
        disp.dist_cm     = 50;                  % from eye to display.
        disp.h_pix       = 1080;                % height in pixels 
        disp.w_pix       = 1920;                % width in pixels
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees 
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees
        disp.refresh_hz  = 60;                  % refreshrate in Hz
        disp.ppd         = disp.h_pix/disp.h_deg;
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)       
        disp.clut        = 0;                   % linear clut
        
    case 'PPROOM_EIZOFLEXSCAN'                  % MATLAB version 2016b, psychtoolbox version 3.0.14 December 30th 2014
        disp.w_cm        = 52;                  % width in cm
        disp.h_cm        = 32.5;                % height in cm
        disp.dist_cm     = 99.0;                % eye to screen viewing distance in cm
        disp.h_pix       = 1200;                % height in pixels 
        disp.w_pix       = 1920;                % width in pixels 
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees 
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees
        disp.refresh_hz  = 60;                  % refreshrate in Hz
        disp.ppd         = disp.h_pix/disp.h_deg;
        disp.xc          = disp.w_pix/2;        % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;        % y-center relative to upper left corner (pixels)
        disp.clut        = 'gammaTable_EIZOFLEXSCAN.mat'; % lookup table
end

return