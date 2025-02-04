function disp = vcd_getDisplayParams(dispname)

if ~exist('dispname','var') || isempty(dispname)
    dispname = '7TASBOLDSCREEN32';
end

disp = struct(); 
disp.name = dispname;

switch dispname
    case '7TASBOLDSCREEN32'
        disp.w_cm        = 69.84;             % cm wide; beyond what subject can actaully see.
        disp.h_cm        = 39.29;             % cm high;
        disp.dist_cm     = 176.5+2+5.5;       % 176.5 cm from the mirror to glass of BOLDScreen, 
                                              % 2 cm from glass to display, 5.5 from eye to mirror (check!).
        disp.h_pix       = 1080;              % in pixels (BOLDscreen vertical, Nova1x32), 800 pixel x 800 pixel square is what people can see.
        disp.w_pix       = 1920;              % in pixels (BOLDscreen horizontal, Nova1x32), 1800 pixels is max extend
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees (BOLDscreen vertical, Nova1x32) should be 12.70 deg
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees (BOLDscreen vertical, Nova1x32)
        disp.refresh_hz  = 120;               % refreshrate in Hz
        disp.ppd         = disp.h_pix/disp.h_deg;
        disp.xc          = disp.w_pix/2;      % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;      % y-center relative to upper left corner (pixels)
    case 'KKOFFICEQ3277'
        disp.w_cm        = 71;                  % cm wide; beyond what subject can actaully see.
        disp.h_cm        = 40;                  % cm high;
        disp.dist_cm     = 50;                  % from eye to display.
        disp.h_pix       = 1080;              % in pixels 
        disp.w_pix       = 1920;              % in pixels 
        disp.h_deg       = pix2deg(disp.h_pix,disp.h_pix,disp.h_cm,disp.dist_cm); % in degrees 
        disp.w_deg       = pix2deg(disp.w_pix,disp.w_pix,disp.w_cm,disp.dist_cm); % in degrees
        disp.refresh_hz  = 60;               % refreshrate in Hz
        disp.ppd         = disp.h_pix/disp.h_deg;
        disp.xc          = disp.w_pix/2;      % x-center relative to upper left corner (pixels)
        disp.yc          = disp.h_pix/2;      % y-center relative to upper left corner (pixels)
        
end

return