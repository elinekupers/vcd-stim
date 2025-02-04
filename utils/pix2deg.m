function deg = pix2deg(stim_pix, res_pix, scr_cm, dist_cm) 
% Function to convert stimulus or screen nr of pixels to degree visual angle.
% 
% INPUTS:
% stim_pix  : size of stimulus in pixels
% res_pix   : resolution of screen in pixels (height or width) in pinxels
% scr_cm    : length of screen in centimeters (height or width, should
%               match res_pix)  in centimeters
% dist_cm   : distance from the screen to the eye in centimeters

deg = atan( (stim_pix/res_pix) * (scr_cm / 2) / dist_cm) / pi*180*2;

return