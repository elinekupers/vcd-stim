function output_im = vcd_applyContrastDecrement(params, c_onset, stim_class,  input_im)
% VCD function to apply a temporal contrast modulation to a static image
% (or multiple frames of a movie). In the case of VCD, we use an inverted
% Gaussian window, where the mean image contrast dips slightly over time
% starting at the time point (c_onset).
%
%    output_im = vcd_applyContrastDecrement(params, c_onset, stim_class, input_im)
%
% INPUTS:
% * params                          : parameter struct. Requires the fields:
%    stim.(<stim_class>).duration     duration of the stimulus in nr of
%                                     presentation frames. E.g., a duration
%                                     of 60 frames at 60 Hz is 1 second.
%    params.stim.cd.t_gauss         : Temporal contrast modulation
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
%
% OUTPUTS:
% * output_im                       : (uint8) output images with temporal contrast
%                                      modulation. Each cell is a frame.
%                                      The number of cells are determined
%                                      by the total stimulus duration, not
%                                      by the temporal contrast modulation
%                                      function, unless the stimulus
%                                      duration and modulation have the
%                                      same duration. This allows us to use
%                                      variable onset times while keeping
%                                      the number of frames in the output
%                                      image the same.

% Check if we deal with a static image (one frame or a movie (multiple
% frames).
if size(input_im,1) == 1
    input_im0 = repmat(input_im, params.stim.(stim_class).duration, 1);
    input_im = input_im0; clear input_im0;
end

% Copy the input image
output_im = input_im;

% loop over time points of the temporal contrast modulation function
for tt = 1:length(params.stim.cd.t_gauss)
    
    % Get image
    tmp_im = output_im{c_onset+tt-1}; % subtract one because we want to start at c_onset
    
    % Get image size
    sz0 = size(tmp_im);
    
    % Check if there is an alpha mask in fourth dim, we will ignore it
    if ndims(tmp_im)==3 && size(tmp_im,3)==4
        tmp_im = tmp_im(:,:,1:3);
    end
    
    if strcmp(stim_class, 'ns')
        % SCENES are in color, where RGB channels range between values of 0-255.
        tmp_im_g = rgb2gray(tmp_im);  % rgb to gray
        tmp_im_g_norm  = (double(tmp_im_g)./255).^2; % range [0-1];
        tmp_im_norm = (double(tmp_im)./255).^2;
        mn_g = mean(tmp_im_g_norm(:));
        
        % subtract the mean luminance of this scene
        tmp_im_c = ((tmp_im_norm-mn_g).*params.stim.cd.t_gauss(tt)) + mn_g;
        tmp_im_c = uint8(255.*sqrt(tmp_im_c)); % bring back to 0-255
        
    elseif strcmp(stim_class,'obj')
        % Objects (range between 1-255)
        tmp_im_c = ((tmp_im-1)./254).^2;                                % convert luminance range [0 1]
        tmp_im_c = tmp_im_c-0.5;                                        % center around 0, min/max range [-0.5 0.5]
        tmp_im_c = tmp_im_c.*params.stim.cd.t_gauss(tt);                % scale contrast
        tmp_im_c = uint8(254.*sqrt(tmp_im_c))+1;                        % bring back to 1-255
        % resize to 3D
        tmp_im_c = reshape(tmp_im_c,sz0);
        
    else % Gabors, RDKs, Dots
        tmp_im_c = (tmp_im./255)-0.5;                                   % center around 0, range [-0.5 0.5]
        tmp_im_c = tmp_im_c.*params.stim.cd.t_gauss(tt);                % scale contrast
        tmp_im_c = uint8( bsxfun(@plus, (255.*tmp_im_c), double(params.stim.bckgrnd_grayval))); % bring back to 1-255
        % resize to 3D
        tmp_im_c = reshape(tmp_im_c,sz0);
    end
    
    % accumulate contrast modulated image
    output_im{c_onset+tt-1} = tmp_im_c;
    
end % tt