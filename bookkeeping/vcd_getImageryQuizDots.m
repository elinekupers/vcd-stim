function [all_quiz_dots, params] = vcd_getImageryQuizDots(params, condition_master)
%
%  [all_quiz_dots, params] = vcd_getImageryQuizDots(params, condition_master)
%
% Purpose:
%   Create imagery quiz dot images for a specific display environment. Each
%   unique image has 20 specific quiz dot images that are called after the
%   delay period, in the second stimulus presentation window of the imagery
%   trial. The goal of the subject is to say how the quiz dots relate to
%   the mental image, according to the text prompt. For example, do the
%   quiz dots overlap the giraffe in the imagined giraffe scene? Yes/No
%    
%   Parameters for a specific display environment can be derived by calling
%   vcd_getStimParams.m.
%
% INPUTS:
%   params         : params struct, should contain a field for display params 
%                    (params.disp), and each stimulus class:
%                       - params.(stim_class).img_im : unique image nrs
%                       that will be used in imagery stim-task crossing 
%   condition_master :
%
% OUTPUTS:
%   all_quiz_dots.(stimClass) : struct with the following fields for each stimulus class: 
%       * images              : quiz dot images: height (pixels) by width (pixels) x 3 (rgb)
%       * masks               : alpha masks for each quiz dot image to crop out image edges
%       * prompts             : cell with individual text prompts related to each quiz image
%       * xy_coords_deg       : [x,y] coordinates of each dot location within quiz image
%       * xy_coords_pix       : [x,y] coordinates of each dot location within quiz image
%
% EXAMPLE:
% params = struct();
% params.verbose        = true; % visualize stimuli or not
% 
% % Get display params
% dispname = '7TAS_BOLDSCREEN32'; % Choose from: '7TAS_BOLDSCREEN32' or 'KKOFFICE_AOCQ3277' or 'EKHOME_ASUSVE247' or 'PPROOM_EIZOFLEXSCAN'
% params.disp   = vcd_getDisplayParams(dispname);
% 
% % Get stimulus parameters
% params.load_params                 = false; % get stim params
% params.store_params                = false; 
% params.overwrite_randomized_params = false; % DO NOT OVERWRITE PARAMS
% p.stim   = vcd_getStimParams('stim_class','all', 'disp_name',p.disp.name, ...
%     'load_params', params.load_params,...
%     'store_params', params.store_params, ...
%     'overwrite_randomized_params', params.overwrite_randomized_params);
%
% [all_quiz_dots, params] = vcd_getImageryQuizDots(params)


all_quiz_dots = struct();

for sc = 1:length(params.exp.stimClassLabels)
    
    stimClass = params.exp.stimClassLabels{sc};

switch stimClass
    
    case 'gabor'
        % Function to create uint8 quiz dot images related to mental image of gabor
        % Quiz image sz is a circular aperture with diameter: 5.658 degree diameter,
        % This diameter is 500 pixels for BOLDscreen or 365 pixels for EIZOFLEXSCAN.
        %
        % im (uint8) is a 6D array containing quiz images for 8 high contrast gabors (8 different orientations): 
        %   height (pixels) x width (pixels) x 3 (rgb) x 8 unique gabor images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % masks (uint8) is a 5D array with alpha transparency masks: 
        %   height (pixels) x width (pixels) x 8 unique gabor images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % text prompts (str) is a 3D cell array containing text prompt for each quiz image 
        %   8 unique gabor images x 10 unique quiz images x 2 quiz type (yes/no)
        %
        % xy_coords_deg (double) is a 6D array with [x,y] coordinates relative 
        % to center of stimulus [0,0] for the dots within a single quiz image in degrees visual angle: 
        %   x-pos (deg) x y-pos (deg) x 2 (dot nr) x 8 unique gabor images x unique 10 quiz images x 2 quiz types (yes/no)
        %
        % xy_coords_pix (double) is a 6D array with [x,y] coordinates relative 
        % to center of stimulus [0,0]) for the dots within a single quiz image in pixels: 
        %   x-pos (pixels) x y-pos (pixels) x 2 (dot nr) x 8 unique gabor images x 10 unique quiz images x 2 quiz types (yes/no)
        [im, mask, prompt,  xy_coords_deg, xy_coords_pix] = vcd_makeQuizDots_GBR(params);
         
    case 'rdk'
        % Function to create uint8 quiz dot images related to mental image of RDK
        % Quiz image sz is a circular aperture with diameter: 5.658 degree diameter,
        % This diameter is 500 pixels for BOLDscreen or 365 pixels for EIZOFLEXSCAN.
        %
        % im (uint8) is a 6D array containing quiz images for 8 high contrast gabors: 
        %   height (pixels) x width (pixels) x 3 (rgb) x 8 unique rdk images x 10 quiz unique images x 2 quiz types (yes/no)
        %
        % masks (uint8) is a 5D array with alpha transparency masks: 
        %   height (pixels) x width (pixels) x 8 unique rdk images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % text prompts (str) is a 3D cell array containing text prompt for each quiz image 
        %   8 unique rdk images x 10 unique quiz images x 2 quiz type (yes/no)
        %
        % xy_coords_deg (double) is a 6D array with [x,y] coordinates (relative to center of stimulus [0,0]) for the quiz dots in degrees visual angle: 
        %   x-pos (deg) x y-pos (deg) x 2 (dot nr) x 8 unique rdk images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % xy_coords_pix (double) is a 6D array with [x,y] coordinates (relative to center of stimulus [0,0]) for the quiz dots in pixels: 
        %   x-pos (pixels) x y-pos (pixels) x 2 (dot nr) x 8 unique rdk images x 10 unique quiz images x 2 quiz types (yes/no)
        [im, mask, prompt, xy_coords_deg, xy_coords_pix] = vcd_makeQuizDots_RDK(params);
        
    case 'dot'
        % Function to create uint8 quiz dot images related to mental image of single dot
        % Quiz image sz is a circular aperture with diameter: 5.658 degree diameter,
        % This diameter is 500 pixels for BOLDscreen or 365 pixels for EIZOFLEXSCAN.

        % im (uint8) is a 6D array containing quiz images for 8 single dots: 
        %   height (pixels) x width (pixels) x 3 (rgb) x 8 unique dot images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % masks (uint8) is a 5D array with alpha transparency masks: 
        %   height (pixels) x width (pixels) x 8 unique dot images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % text prompts (str) is a 3D cell array containing text prompt for each quiz image 
        %   8 unique dot images x 10 unique quiz images x 2 quiz type (yes/no)
        %
        % xy_coords_deg (double) is a 6D array with [x,y] coordinates (relative to center of stimulus [0,0]) for the quiz dots in degrees visual angle: 
        %   x-pos (deg) x y-pos (deg) x 2 (dot nr) x 8 unique dot images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % xy_coords_pix (double) is a 6D array with [x,y] coordinates (relative to center of stimulus [0,0]) for the quiz dots in pixels: 
        %   x-pos (pixels) x y-pos (pixels) x 2 (dot nr) x 8 unique dot images x 10 unique quiz images x 2 quiz types (yes/no)
        [im, mask, prompt, xy_coords_deg, xy_coords_pix] = vcd_makeQuizDots_DOT(params);
        
    case 'obj'
        % Function to create uint8 quiz dot images related to mental image of single dot
        % Quiz image sz is a circular aperture with diameter: 5.658 degree diameter,
        % This diameter is 500 pixels for BOLDscreen or 365 pixels for EIZOFLEXSCAN.

        % im (uint8) is a 6D array containing quiz images for 8 high contrast gabors: 
        %   height (pixels) x width (pixels) x 3 (rgb) x 8 unique obj images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % masks (uint8) is a 5D array with alpha transparency masks: 
        %   height (pixels) x width (pixels) x 8 unique obj images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % text prompts (str) is a 3D cell array containing text prompt for each quiz image 
        %   8 unique obj images x 10 quiz images x 2 quiz type (yes/no)
        %
        % xy_coords_deg (double) is a 6D array with [x,y] coordinates (relative to center of stimulus [0,0]) for the quiz dots in degrees visual angle: 
        %   x-pos (deg) x y-pos (deg) x 2 (dot nr) x 8 unique obj images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % xy_coords_pix (double) is a 6D array with [x,y] coordinates (relative to center of stimulus [0,0]) for the quiz dots in pixels: 
        %   x-pos (pixels) x y-pos (pixels) x 2 (dot nr) x 8 unique obj images x 10 unique quiz images x 2 quiz types (yes/no)        
        [im, mask, prompt, xy_coords_deg, xy_coords_pix] = vcd_gmakeQuizDots_OBJ(params);
        
    case 'ns'
        % Function to create uint8 quiz dot images and prompts to probe and quiz subject's mental image
        % of naturalistic scene. Quiz image sz is a square aperture with width/height: 8.4 degrees,
        % This diameter is 742 pixels for BOLDscreen or 541 pixels for EIZOFLEXSCAN.

        % im (uint8) is a 6D array containing quiz images for 8 high contrast gabors: 
        %   height (pixels) x width (pixels) x 3 (rgb) x 15 unique ns images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % NO alpha transparency mask needed 
        %  (or can be 5D array with NaNs: height (pixels) x width (pixels) x 15 unique ns images x 10 unique quiz images x 2 quiz types)
        %
        % text prompts (str) is a 3D cell array containing text prompt for each quiz image 
        %   15 unique ns images x 10 unique quiz images x 2 quiz type (yes/no)
        %
        % xy_coords_deg (double) is a 6D array with [x,y] coordinates (relative to center of stimulus [0,0]) for the quiz dots in degrees visual angle: 
        %   x-pos (deg) x y-pos (deg) x 2 (dot nr) x 15 unique ns images x 10 unique quiz images x 2 quiz types (yes/no)
        %
        % xy_coords_pix (double) is a 6D array with [x,y] coordinates (relative to center of stimulus [0,0]) for the quiz dots in pixels: 
        %   x-pos (pixels) x y-pos (pixels) x 2 (dot nr) x 15 unique ns images x 10 unique quiz images x 2 quiz types (yes/no)
        [im, ~, prompt, xy_coords_deg, xy_coords_pix] = vcd_makeQuizDots_NS(params);
end
    
all_quiz_dots.(stimClass).images        = im;
all_quiz_dots.(stimClass).masks         = mask;
all_quiz_dots.(stimClass).prompts       = prompt;
all_quiz_dots.(stimClass).xy_coords_deg = xy_coords_deg;
all_quiz_dots.(stimClass).xy_coords_pix = xy_coords_pix;

end
