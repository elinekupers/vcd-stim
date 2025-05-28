function [output_im, output_mask] = vcd_applyContrastDecrement(params, c_onset, stim_class,  input_im, varargin)
% VCD function to apply a temporal contrast modulation to a static image
% (or multiple frames of a movie). In the case of VCD, we use an inverted
% Gaussian window, where the mean image contrast dips slightly over time
% starting at the time point (c_onset).
%
%    output_im, output_mask] = vcd_applyContrastDecrement(params, c_onset, ...
%                                           stim_class, input_im, input_mask)
%
%  NOTE: We expect different luminance ranges depending on the stimulus class: 
%   * Gabors    : 1-255 (created for linearized display)
%   * RDKs      : 1-255 (created for linearized display)
%   * Dots      : 1-255 (created for linearized display)
%   * Objects   : 1-255 (initially made with a gamma=2.0 display, but manually converted to work for a linearized display)
%   * NS        : 0-255 (initially made with a gamma=2.0 display, but manually converted to work for a linearized display)
%  
%
%
% INPUTS:
% * params                          : parameter struct. Requires the fields:
%    stim.(<stim_class>).duration     duration of the stimulus in nr of
%                                     presentation frames. E.g., a duration
%                                     of 60 frames at 60 Hz is 1 second.
%    stim.cd.t_gauss                : Temporal contrast modulation
%                                     function, that describes the fraction
%                                     by which the mean contrast of the
%                                     input image will be altered.
%                                     decrement in frames.
%                                     The values of the temporal function
%                                     envelope describe the relative
%                                     decrement of the given input image
%                                     relative to contrast level of the
%                                     input image. 1 means no change (input
%                                     image contrast is unaltered) 0 means
%                                     mean luminance (all contrast is removed).
% * c_onset                         : The first time point of the temporal
%                                     contrast modulation (in frames,
%                                     relative to the first frame of the
%                                     stimulus). Onset time needs to be
%                                     early enough such that duration of
%                                     temporal contrast modulation function
%                                     ends before or at the same time as
%                                     the stimulus duration.
% * stim_class                      : stimulus class of the given input image.
%                                     Needs to be one of the following:
%                                     'gabor', 'rdk', 'dot', 'obj', ns'.
% * input_im                        : input image in pixels (uint8). Can be
%                                      a 2D achromatic image
%                                      (height,width), a 3D RGB image
%                                      without alpha transparency mask
%                                      (height,width,3), or a 3D RGB image
%                                      with alpha transparency mask
%                                      (height,width,4). If there is an
%                                      alpha transparency mask, it will be
%                                      ignored when applying the contast
%                                      decrement.
% * input_mask                      : (optional): separated alpha
%                                      transparency mask(s) for the
%                                      input_im(s). input_mask has the same
%                                      dimensions as input_im in pixels
%                                      (uint8), and is expected to be a 2D
%                                      image (height,width) or a 3D image
%                                      (height,width,3). If input_im has an
%                                      alpha transparency mask, input_mask
%                                      will be ignored. input_mask values
%                                      can be binary [0 means mask pixels,
%                                      1 means show pixels] or continuous
%                                      values between [0-255], where 0 =
%                                      mask, 255 = show.
%
% OUTPUTS:
% * output_im                       : (uint8) images with temporal contrast
%                                      modulation. Each
%                                      cell is a time frame with a single
%                                      image with the same image
%                                      dimensions as the input_im (h x w or
%                                      h x w x 3 or h x w x 4).
%                                      The number of cells are determined
%                                      by the total stimulus duration, not
%                                      by the temporal contrast modulation
%                                      function, unless the stimulus
%                                      duration and modulation have the
%                                      same duration. This allows us to use
%                                      variable onset times while keeping
%                                      the number of frames in the output
%                                      image the same.
% * output_mask                     : (uint8) alpha transparency masks
%                                      corresponding to output_im. Each
%                                      cell is a time frame with a single
%                                      mask image with the same image
%                                      dimensions and range as the 
%                                      input_mask (h x w or h x w x 3),
%                                      range between [1-255]. We duplicate 
%                                      the masks such that output_mask has 
%                                      the same size as output_im.
%
% Written by Eline Kupers @ UMN 2025/05
%%

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'       , @isstruct);
p0.addRequired('c_onset'      , @(x) isnumeric(x) & (x>0));
p0.addRequired('stim_class'   , @(x) ischar(x) & ismember(x, {'gabor','rdk','dot','obj','ns'}));
p0.addRequired('input_im'     , @iscell);
p0.addParameter('input_mask'  , [], @iscell); 

% Parse inputs
p0.parse(params,c_onset,stim_class,input_im,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

% Check if user defined the temporal modulation function and stim duration
assert(isfield(params.stim.cd, 't_gauss'))
assert(~isempty(params.stim.cd.t_gauss))
assert(isfield(params.stim.(stim_class), 'duration'))
assert(~isempty(params.stim.(stim_class).duration))

% Check if onset of contrast decrement falls within stimulus duration
assert(c_onset < params.stim.(stim_class).duration)

% Check if temporal modulation function fits within stimulus duration
assert(length(params.stim.cd.t_gauss) <= params.stim.(stim_class).duration)

% Check if onset of contrast decrement allows for temporal modulation function
% to fall entirely within stimulus duration
assert((c_onset+length(params.stim.cd.t_gauss)-1) <= params.stim.(stim_class).duration)

% Check if we deal with a static image (one frame or a movie (multiple
% frames).
if size(input_im,1) == 1
    input_im0 = repmat(input_im, params.stim.(stim_class).duration, 1);
    input_im = input_im0; clear input_im0;
    
    if ~isempty(input_mask) %#ok<NODEF>
        input_mask0 = repmat(input_mask, params.stim.(stim_class).duration, 1);
        input_mask = input_mask0; clear input_mask0;
    end
end


% Copy the input image
output_im = input_im;
output_mask = input_mask;

% loop over time points of the temporal contrast modulation function
for tt = 1:length(params.stim.cd.t_gauss)
    
    % Get image
    tmp_im = output_im{c_onset+tt-1}; % subtract one because we want to start at c_onset
    
    if ~isempty(input_mask) 
        tmp_mask = output_mask{c_onset+tt-1}; % subtract one because we want to start at c_onset
    end
    
    % Get image size
    sz0 = size(tmp_im);
    
    % Check if there is an alpha mask in fourth dim of input_im or
    % separately defined by input_mask)
    
    if ndims(tmp_im)==3 && size(tmp_im,3)==4
        tmp_im = tmp_im(:,:,1:3);
        if isempty(tmp_mask)
           tmp_mask = double(tmp_im(:,:,4));
        end
    end
    
    % apply alpha transparency mask
    if exist('tmp_mask','var') && ~isempty(tmp_mask) 
        
        if strcmp(stim_class,'obj')
            mask_idx = (double(tmp_mask)./255)>0.5; % note simple hack
        else
            mask_idx = (double(tmp_mask)./255)>0;
        end
        if ndims(mask_idx)==2 && ndims(tmp_im)==3
            mask_idx = repmat(mask_idx, 1,1,3);
        end
        tmp_im0   = double(tmp_im);
        tmp_im0  = tmp_im0(mask_idx);
    else
        tmp_im0  = double(tmp_im);
    end
    
    if strcmp(stim_class, 'ns')
        % SCENES are in color, where RGB channels range between values of 0-255.
        tmp_im_g        = rgb2gray(tmp_im0);                                % convert rgb to gray
        tmp_im_g_norm   = (tmp_im_g./255).^2;                       % convert grayscale image range from [0-255] to [0-1] (normalized image)
        tmp_im_norm     = (tmp_im0./255).^2;                        % convert color image range from [0-255] to [0-1] (normalized image)
        mn_im           = mean(tmp_im_g_norm(:));                           % compute the mean of grayscale image
        % subtract the mean luminance of this scene from the normalized
        % color image, then scale contrast for the given time frame in the
        % contrast modulation function, and add back the mean.
        tmp_im_c        = ((tmp_im_norm-mn_im).*params.stim.cd.t_gauss(tt)) + mn_im;  
        tmp_im_c        = (255.*sqrt(tmp_im_c));                            % bring back to 0-255
        
    elseif strcmp(stim_class,'obj')
        % Objects are grayscale, where RGB channels range between luminance values [1-255]   
        tmp_im_norm     = ((tmp_im0-1)./254).^2;                            % convert image luminance range [0-1]
        mn_im           = mean(tmp_im_norm(:));                             % compute the mean of image
        % center around mean, scale contrast for the given time frame in 
        % the contrast modulation function, bring min/max range back to [0-1]  
        tmp_im1         = ((tmp_im_norm - mn_im)*params.stim.cd.t_gauss(tt)) + mn_im;                  
        tmp_im_c        = ( (254.*sqrt(tmp_im1)) +1 );                      % bring min/max range back to [1-255]
        
    else % Gabors, RDKs, Dots are gray scale, where RGB channels range between luminance values [1-255], and mean luminance is 128
        tmp_im_norm = (tmp_im0./255)-0.5;                                      % center around 0, range [-0.5 0.5]
        tmp_im1     = tmp_im_norm.*params.stim.cd.t_gauss(tt);                    % scale contrast for the given time frame in the contrast modulation function
        tmp_im_c    = ( (255.*(tmp_im1+0.5)) );                               % bring back to 1-255
    end
    
    % Create full uint8 image
    if exist('tmp_mask','var') && ~isempty(tmp_mask) 
        tmp_im_full = tmp_im;                           % No need to resize
        tmp_im_full(mask_idx) = uint8(tmp_im_c);        % Convert from double to uint8
    else
        tmp_im_c    = reshape(tmp_im_c,sz0);            % Resize to 2D (or 3D)
        tmp_im_full = uint8(tmp_im_c);                  % Convert from double to uint8
    end

    % accumulate contrast modulated image
    output_im{c_onset+tt-1} = tmp_im_full;
    
    
    
end % tt

return