function [floc_im, zebra_im, hmt_im] = vcd_create_afloc_stim(params, varargin)
% VCD function to create the localizer stimuli for VCD stimulus position.
%
%   afloc_stim = vcd_create_afloc_stim(params, varargin)
% 
% INPUTS
% 
%
% OUTPUTS
%   afloc_stim  
%
% 
% What we want:
% FoV: 8.4 x 8.4° for the majority, except for the Left, Right and single dot ring
% Presentation rate: 2 images per second (2 Hz)
% Presentation cycle: 4-s mini block (8 frames). 
% All (30?) conditions in one big condition table.
% 
% Stim type 1: Updated/adapted fLoc (w/ food items) + lame slow fixation task
% - Include phase scrambled images (to localize LOC)
% - Tools manipulative? ⇐ maybe?? Do people care??
% % 
% Stim type 2: Zebra stim pos localizer
% Thin ring to cover the dot locations - zebra patterns (or checkerboard?)
% Left location - zebra patterns (or checkerboard?)
% Right location -- zebra patterns (or checkerboard?)
% 'dotring'          = Thin (1-deg width) ring to cover the dot locations, iso eccentricity = 4 deg 
% 'single-dots'      = 16 images each containing one 1-deg diameter single dots, requires num to be 16.
% 'all-dots'         = 1 image with all 16 1-deg diameter single dot stimuli 
% 'circle-left'      = 4-deg diameter circle (for peripheral classic stimuli) positioned at 4 deg eccentricity, left from the fixation circle, on the horizontal meridian
% 'circle-right'     = 4-deg diameter circle (for peripheral classic stimuli) positioned at 4 deg eccentricity, right from the fixation circle, on the horizontal meridian
% 'circle-leftright' = two 4-deg diameter circles (for peripheral classic stimuli) positioned at 4 deg eccentricity, both left and right from the fixation circle, on the horizontal meridian
% 'square'           = Central 8.4-deg square location (for the NS stimuli) positioned in the center of the screen.
%
% Stim type 3: hMT+ localizer
% FULL SCREEN (1080x1920 pixels)
% Static dots (0% coherence?) vs. moving dots (100% coherence?) (expanding & contracting dots + Fixation task?). 
% Presentation rate: 60 Hz
%
% EXAMPLE:
% params.disp = vcd_getDisplayParams('7TAS_BOLDSCREEN32');
% params.stim = vcd_getStimParams('disp_name','7TAS_BOLDSCREEN32');
% [floc_im, zebra_im] = vcd_create_afloc_stim(params, 'num', 1)


%% Parse inputs
p = inputParser;
p.addRequired('params'          , @isstruct); % params struct
p.addParameter('pixoffset'      , [0 0]   , @isnumeric);                           % [x,y] offset of center in pixels  default = no offset: [0 0] 
p.addParameter('store_images'   , false   , @islogical);                           % store images as PNG or not
p.addParameter('save_dir'       , fullfile(vcd_rootPath,'workspaces', 'figs','afloc'), @ischar); % base folder to store PNG images
p.addParameter('stimfile'       , []      , @ischar);                              % mat file to store image

% Parse inputs
p.parse(params, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p

if ~isfield(params,'verbose') || isempty(params.verbose)
    params.verbose = false;
end    
    
% deal with center offset
if isequal(pixoffset,[0 0]) || ~isempty(pixoffset)
    xc = params.disp.xc + pixoffset(1);
    yc = params.disp.yc + pixoffset(2);
else
    xc = params.disp.xc;
    yc = params.disp.yc;
end

if store_images && ~exist('save_dir','dir')
    mkdir(save_dir);
end

%% ***** Stim type 1: Gray scale fLOC images ***** 
% Note: OG stimuli for food category are 1080x1080 pixels, all other categories are 1024x1024 pixels
% Stimuli are resized to NS square (741x741 pixels) and squared to account
% for linearized BOLD screen.
pth_to_floc = fullfile(vcd_rootPath,'workspaces','stimuli','RAW','vcd_afloc','floc_stimuli');

[floc_im,floc_cat] = vcd_floc_im(params,'square_pix_val',true,'store_imgs',true,'pth_to_floc',pth_to_floc);


%% ***** Stim type 2: VCD stim position zebra localizer *****
% zebra_im are made for 7TAS BOLD SCREEN: 1080x1080 pixels 
stim_loc_type = {'dotring','singledots','alldots','circleleft','circleright','circleleftright','square'};
num_images    = 144; % number of unique background images desired.  144 = same as fLoc
zebra_cpf     = 32;  % cycles per field of view for zebra pattern. One cycle = 1 black-white reversal. Higher numbers = finer pattern. Default: 32. 
zebra_im      = cell(1,length(stim_loc_type));
for ii = 1:length(stim_loc_type)
    zebra_im{ii} = vcd_zebra_stim_pos(params, 'stimtype',stim_loc_type{ii}, ...
                  'pixoffset', pixoffset, 'zebra_cpf', zebra_cpf, ...
                  'num', num_images, 'store_imgs', true);
end


%% ***** Stim type 3: hMT localizer 
% FULL SCREEN (1080x1920 pixels)
% Static dots (0% coherence?) vs. moving dots (100% coherence?) 
% (expanding & contracting dots + Fixation task?). 
% Presentation rate: 60 Hz for BOLDscreen native resolution 120 Hz
[hmt_im, dotlocs, rng_seed] = vcd_hMT_loc(params, verbose, store_imgs);


return

