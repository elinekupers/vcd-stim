function [floc_im,floc_cat] = vcd_floc_im(params,varargin)
% VCD function to rescale OG fLOC images to be the same size as VCD scenes
% for 7TAS BOLD SCREEN

% Example:
% [floc_im,floc_cat] = vcd_floc_im(params,'square_pix_val',true,'store_imgs',true)

%% Parse inputs
p = inputParser;
p.addRequired('params'          , @isstruct); % params struct
p.addParameter('num_floc_im'    , 144, @isnumeric);
p.addParameter('square_pix_val' , true, @islogical); % for BOLD display
p.addParameter('pth_to_floc'    , fullfile(vcd_rootPath,'workspaces','stimuli','RAW','vcd_afloc','floc_stimuli'), @ischar);
p.addParameter('store_imgs'     , false, @islogical);
p.addParameter('save_dir'       , fullfile(vcd_rootPath,'figs',params.disp.name,'localizer','afloc'), @ischar);

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

%%
floc_cat    = {'adult_face','body','car','child_face','corridor','food','house','instrument','limb','number','scrambled','word'};

% Preallocate space, each high-level category gets dim
num_floc_im   = 144; % we anticipate 144 original images
floc_im       = uint8(ones(params.stim.ns.img_sz_pix,params.stim.ns.img_sz_pix, 3, length(floc_cat), num_floc_im));


% Give an update to experimenter
fprintf('[%s]: Resampling the stimuli; this may take a while..\n',mfilename);
tic;

for ff = 1:length(floc_cat)
    
   
    for ii = 1:num_floc_im
        statusdots(ii,num_floc_im);
        
         % get OG images
        d = dir(fullfile(pth_to_floc,floc_cat{ff}, sprintf('*-%d.jpg',ii)));
        
        % Read in image --> NOTE: initial image lum vals ranges from [0-255]
        tmp = imread(fullfile(d.folder,d.name));
        
        % calculate scale factor (we want the final size to be the same
        % as the VCD scenes: 741x741
        dres = params.stim.ns.img_sz_pix/ size(tmp,1);  % scale factor to apply
        
        % before resizing convert to double
        sc0  = double(tmp); % image lum vals are now double
        
        % Resize
        sc1  = imresize(sc0,dres); % resize
        
        % If needed, square pixels for linearized display
        if params.stim.ns.square_pix_val
            sc2 = 255.*((sc1./255).^2);
        else % skip step
            sc2 = sc1;
        end
        
        % Convert to uint8
        sc3 = uint8(sc2);
        
        % check for clipping
        assert( min(sc3(:))>=0 && max(sc3(:))<=255)
        
        % replicate for RGB
        sc4 = repmat(sc3, [1 1 3]);
        
        % Add image to the array
        floc_im(:,:,:,ff,ii) = sc4;

        % save image if requested
        if store_imgs
            if square_pix_val
                save_dir2 = fullfile(save_dir,'resized_squared',floc_cat{ff});
            else
                save_dir2 = fullfile(save_dir,'resized',floc_cat{ff});
            end
            
            if ~exist(save_dir2,'dir')
                mkdir(save_dir2);
            end

            imwrite(sc4, fullfile(save_dir2, sprintf('%03d_vcd_afloc_%s.png', ii, floc_cat{ff})));
        end
        
        % clean up
        clear sc4 sc3 sc2 sc1 sc0 tmp
    end
end

fprintf('[%]: done!',mfilename); toc

