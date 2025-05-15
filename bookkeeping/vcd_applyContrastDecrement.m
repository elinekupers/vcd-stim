function output_im = vcd_applyContrastDecrement(params, c_onset, stim_class,  input_im)


if size(input_im,1) == 1
    input_im0 = repmat(input_im, params.stim.(stim_class).duration, 1);
    input_im = input_im0; clear input_im0;
end

output_im = input_im;


if strcmp(stim_class,'ns')
    
    % only one image, in first column
    assert(all(cellfun(@isempty, input_im(:,2))));
    
    for tt = 1:length(params.stim.cd.t_gauss)
        
        % get image
        tmp_im = output_im{(c_onset-1)+tt,1};
        sz0 = size(tmp_im);
        
        if ndims(tmp_im)==3 && size(tmp_im,3)==4
            tmp_im = tmp_im(:,:,1:3);
        end
        
        tmp_im_g = rgb2gray(tmp_im);  % rgb to gray
        tmp_im_g_norm  = (double(tmp_im_g)./255).^2; % range [0-1];
        tmp_im_norm = (double(tmp_im)./255).^2;
        mn_g = mean(tmp_im_g_norm(:));

        % subtract the mean luminance of this scene
        tmp_im_c = ((tmp_im_norm-mn_g).*params.stim.cd.t_gauss(tt)) + mn_g;
        tmp_im_c = uint8(255.*sqrt(tmp_im_c)); % bring back to 0-255
        
        output_im{(c_onset-1)+tt,1} = tmp_im_c;
    end
    
else % Gabors, RDKs, Single dot, Objects
    
    % check right/left stimulus locations
    nsides = find(~cellfun(@isempty, input_im(1,:)));
    
    for side = nsides
        
        % Get onset of contrast decrement within the
        % stimulus period
%         c_onset(side) = feval(cdsoafun);

        for tt = 1:length(params.stim.cd.t_gauss)
            
            tmp_im = double(output_im{c_onset+tt-1,side});
            sz0 = size(tmp_im);
            
            % Check if there is an alpha mask in fourth dim
            if ndims(tmp_im)==3 && size(tmp_im,3)==4
                tmp_im = tmp_im(:,:,1:3);
            end
            
            % Objects
            if strcmp(stim_class,'obj')
                
                tmp_im_c = ((tmp_im-1)./254).^2; % range [0 1]
                tmp_im_c = tmp_im_c-0.5; % center around 0, min/max range [-0.5 0.5]
                tmp_im_c = tmp_im_c.*params.stim.cd.t_gauss(tt); % scale
                tmp_im_c = uint8(254.*sqrt(tmp_im_c))+1; % bring back to 1-255
                
                tmp_im_c_rz = reshape(tmp_im_c,sz0);
                
            else % Gabors, RDKs, Dots
                tmp_im_c = (tmp_im./255)-0.5; % center around 0, range [-0.5 0.5]
                tmp_im_c = tmp_im_c.*params.stim.cd.t_gauss(tt); % scale
                tmp_im_c = uint8( bsxfun(@plus, (255.*tmp_im_c), double(params.stim.bckgrnd_grayval))); % bring back to 1-255
                
                tmp_im_c_rz = reshape(tmp_im_c,sz0);
            end
            
            output_im{c_onset(side)+tt-1,side} = tmp_im_c_rz;
        end
        
    end % sides
end  % if ns