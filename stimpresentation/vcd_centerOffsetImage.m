function [offset] = vcd_centerOffsetImage(filename, displayname)
% VCD function to generate offset of central image on a display.
%
%       [offset] = vcd_centerOffsetImage(filename,displayname)
% 
% INPUTS: 
% * filename    : png image to load. If left empty, we will use VCD
%                 default: fullfile(vcd_rootPath, 'workspaces', 'stimuli'
%                 <dispname>,'calibrate','7TAS_BOLDSCREEN_pixelprecise_layout_cropped_white_square_matches_ring.png').
%                 On the 7TAS we rsynced three calibration images 
%                 to the vcd-stim code repo: 
%                 1) cropped_color.png - the union of the classic stimuli (in white) +
%                 central NS square (in green) and single dot ring (in red). 
%                 2) fullscreen_color.png - the union of the classic
%                 stimuli (in white) + central NS square (in green) and
%                 single dot ring (in red) on a full screen gray background
%                 (uint 128).
%                 3) cropped_white_square_matches_ring.png - the union of 
%                 the classic stimuli and central NS square that matches 
%                 the height of the single dot ring. Cutout is in white.
% * displayname : character string referring to the name of the monitor used. 
%                 If left empty, we will use '7TAS_BOLDSCREEN32'. 
%                 Options are: '7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277',
%                 'PPROOM_EIZOFLEXSCAN','EKHOME_ASUSVE247','CCNYU_VIEWPIXX3D'
% 
% This function relies on knkutils function ptviewimage.m:
% ptviewimage shows the images named by <files>, starting in the center
% of the screen. 
% 
% To move image position, you can use one of 3 options:
% 1. Use the arrow keys (increments of 1 pixel)
% 2. Use ijkl keys      (increments of 10 pixels)
% 3. Use wasd keys      (increments of 30 pixels)
% 
% If mode = 1, use spacebar or r to cycle through the different images.
%
% To exit: use escape
%
% ptviewimage inputs:
% * files       : full path to images you want to load. To load multiple
%                 images, you can use '/path/to/images/*.png'
% * imageflip   : (optional) is [J K] where J is whether to flip first
%                 dimension (flip up-down) and K is whether to flip second
%                 dimension (flip left-right).
%                 default: [0 0]. (no flipping)
% * mode        : (optional) 0 means normal operation,
%                 1 means automatically cycle through images (1 per
%                 second).  in this mode, the only possible key is escape,
%                 and you must hold it down until it is detected. 
%                 default: 0.
% ptviewimage outputs:
% * offset      : two signed integers indicating the final [x,y] position 
%                 of the image(s) relative to display center in pixels. 
%                 [0,0] means centered/no offset, 
%                 [10 20] means move 10-px right, 20-px down                 
% 
% Written by EK @ UMN 07/2025

% Check inputs
if ~exist('filename','var') || isempty(filename)
    filename = fullfile(vcd_rootPath, 'workspaces','stimuli',displayname,'calibrate','7TAS_BOLDSCREEN_pixelprecise_layout_cropped_white_square_matches_ring.png');
end

if ~exist('displayname','var') || isempty(displayname)
    displayname = '7TAS_BOLDSCREEN32';
end
assert(ismember(displayname,{'7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','PPROOM_EIZOFLEXSCAN','EKHOME_ASUSVE247','CCNYU_VIEWPIXX3D'}));

% Load display params
disp_params = vcd_getDisplayParams(displayname);

% Define ptviewimage params
imageflip = [0,0]; % flip up-down or left-right
mode      = 0;     % normal operation vs 1-s auto cycle


% Define pton params
skipsync = 0;

% pton input arguments are:
% 1: [width, height, framerate, bitdepth]
% 2: winsize (fraction: default is full extent)
% 3: clutfile -- 0 for linear CLUT (-2 for squaring CLUT for BOLDSCREEN to simulate normal monitors --> NOTE: we do this manually!)
% 4: skipsync (bool: 0 is false -- do not skip the text, 1 is true -- skip the test)
% 5: wantstereo (bool: default is false)
if strcmp(disp_params.name, 'PPROOM_EIZOFLEXSCAN')
    % apparently PP room monitor native refresh rate shows up as 0 (but is
    % actually 60 Hz), due to the Bits# stimulus processor.
    % PProom EIZOFLEXScan screen ptonparams are expected to be {[1920 1200 0 24],[], 0, 0}
    ptonparams = {[disp_params.w_pix disp_params.h_pix 0 24],[],disp_params.clut, skipsync};
else % '7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','EKHOME_ASUSVE247', 'CCNYU_VIEWPIXX3D'
    % Nova1x32 coil with BOLDscreen and big eye mirrors ptonparams are expected to be {[1920 1080 120 24],[], 0, 0}
    ptonparams = {[disp_params.w_pix disp_params.h_pix disp_params.refresh_hz 24],[],disp_params.clut, skipsync};
end

% Open PTB window
oldCLUT = pton(ptonparams{:});

% view image! log offset!
offset  = ptviewimage(filename,imageflip,mode); 

% Close textures
Screen('Close');

% Close PTB window
ptoff(oldCLUT);

return