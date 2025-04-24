function [objects, masks, im_order, info] = vcd_objects(params)
% VCD function:
%   [objects, masks, info] = vcd_objects(params)
%
% Purpose:
%   Luminance-correct, load, and store object images for VCD experimental display.
%   Requires two files:
%    1. a set of individual object image files, where the folder is defined
%       by p.stim.obj.indivobjfile, e.g.:
%       fullfile(vcd_rootPath,'workspaces/stimuli/vcd_complex_objects/*_rot*.png)
%    2. A csv file defined in p.stim.obj.infofile, e.g.: 
%       fullfile(vcd_rootPath,'workspaces/objects_info.csv')
%   Given that 0, 90, and 180 degrees are right, forwards, and left facing, 
%   respectively. Facing directions of objects in the following files are:
%       * 1-23 (0-44 deg): face sideways
%       * 24-46 (46-90 deg): face forward
%       * 47-68 (92-134 deg): face forward
%       * 69-91 (136-180 deg): face sideways
%   
%   NOTE: If params.stim.obj.dres~=1, function rescales images by
%   converting to double, [0-1] lum range (assume range is 1-255), and
%   square as images were made with gamma = 2.
% 
%   When params.store_params = true, this function will store generated
%   object images as a single mat file in params.stim.obj.stimfile (e.g.:
%   fullfile(vcd_rootPath,'workspaces',....
%   'stimuli',<disp_name>,'objects.mat').
%
%
% INPUTS:
%  params            :   struct with stimulus params
%
% OUTPUTS:
%  objects      :   uint8 images of complex images, dimensions are
%                   width x height x 3 (rgb) x object's subordinate
%                   category x num of rotations
%  masks        :   uint8 alpha transparency masks for ptb to remove the 
%                   background. Dimensions are width (pixels) x height
%                   (pixels) x object's subordinate category x num of
%                   rotations
%  im_order     :   (struct) with fields containing a subset of info about 
%                   objects in the order they were loaded:
%                       object_name: {1×80 cell} - same as filename in info
%                       base_rot: [1×80 double]  - base rotation of core object
%                       abs_rot: [1×80 double]   - absolute rotation of object (after adding delta)
%                       rel_rot: [1×80 double]   - rel rotation of object (i.e. the delta)
%                       unique_im_nr: [1×80 double] - same as unique_im_nr in info
%                       facing_dir: {1×80 cell}  - same as facing_dir in info
%                       rel_rot_name: {1×80 cell} - same as rel_rot_name in info
%  info         :   Loaded csv from workspaces/stimuli/ with png file names
%                   and object category information. 
%              
%      filename                         unique_im_nr   stim_pos  stim_pos_i  superordinate    superordinate_i     basic      basic_i  subordinate   subordinate_i  rot_abs  rot_rel   facing_dir   is_in_img  rel_rot_name
%  {'faces_damonwayans_rot31.png'  }        65         {'left' }    1          {'human'}           1           {'facemale'}     1      {'damon'}        1            60       0       {'forward'}    1       {'rightwards'}              
%      ...
%  {'places_watertower_rot75.png'  }        430       {'right'}     2         {'place' }           4           {'tower'}        2      {'watertower'}   3            148      8       {'sideways'}   0       {'rightwards'}    
%
% Written by Eline Kupers 2024/12, updated 2025/04

%% If you haven't done so yet, first preprocess object images:
% s_preprocessObjects.m


%% Prepare images

% Get preprocessed image folder
d = dir(sprintf('%s*',fullfile(params.stim.obj.infofile)));

% Read in table with filenames
info = readtable(fullfile(d(end).folder,d(end).name));

% Define superordinate and basic categories, and number of exemplars per basic category
idx0            = info.rot_rel == 0;
subordinate     = info.subordinate;
base_rot        = info.rot_abs(idx0);  % alternating 2 view ± 4 steps spaced 2 deg apart

% Define all rotations
rotations       = [0, params.stim.obj.delta_from_ref];

% reshape unique image nr for WM test images
unique_ref_im   = reshape(params.stim.obj.unique_im_nrs_WM, [],sum(idx0));

% Get (rescaled) image extent
extent   = params.stim.obj.img_sz_pix.*params.stim.obj.dres;

% Preallocate space
objects = uint8(ones(extent, extent,3, length(subordinate(idx0)),length(rotations)));
masks = uint8(ones(extent, extent, length(subordinate(idx0)),length(rotations)));

counter = 1; % we use the lazy counter way to keep track of image

% Loop over each image
for sub = 1:length(subordinate(idx0))
    
    % loop over each rotation
    for rr = 1:length(rotations)
        
        % Get file name
        rot = 1+((base_rot(sub) + rotations(rr))/2);
        fname = sprintf('*%s*_rot%02d.png', subordinate{sub},rot);
        d = dir(fullfile(params.stim.obj.indivobjfile,fname));
        if isempty(d)
            error('[%s]: Can''t find image file!')
        end
        
        if rr==1
            assert(strcmp(info.filename{sub},d.name))
        else
            assert(strcmp(cell2mat(info.filename(length(subordinate(idx0)) + (rr-1) + (length(params.stim.obj.delta_from_ref)*(sub-1)))),d.name));
        end
        
        % Load image
        [im0, ~, alpha0] = imread(fullfile(d.folder,d.name),'BackgroundColor','none');
        
        % If we need to rescale images:
        % convert to double, [0-1] lum range (assume range is 1-255),
        % and square as images were made with gamma = 2.
        if ~isequal(params.stim.obj.dres,1) && ~isempty(params.stim.obj.dres)
            fprintf('[%s]: Resampling the stimuli..',mfilename);
            
            im0     = (double(im0-1)./254).^2; % convert to [0-1] and square
            alpha0  = (double(alpha0./255)).^2; % convert to [0-1] and square
            
            im1     = imresize(im0, params.stim.obj.dres); % scale if needed
            alpha1  = imresize(alpha0,params.stim.obj.dres); % scale if needed
            
            im = uint8(1+(sqrt(im1)*254));
            alpha_im = uint8((sqrt(alpha1)*255));
        else
            im = im0;
            alpha_im = alpha0;
        end
        clear im0 alpha0
        
        % in case we want to debug:
        % figure(99); clf; imagesc(im,'AlphaData',alpha_im); colormap gray; colorbar; axis image
        
        objects(:,:,:,sub,rr) = repmat(im, [1 1 3]);
        masks(:,:,sub,rr) = alpha_im;
        
        % Keep track of image order
        im_order.object_name(counter)    = subordinate(sub);
        im_order.base_rot(counter)       = base_rot(sub);
        im_order.abs_rot(counter)        = base_rot(sub)+rotations(rr);
        im_order.rel_rot(counter)        = rotations(rr);
        if rr == 1
            im_order.unique_im_nr(counter)   = params.stim.obj.unique_im_nrs(sub);
        else
            im_order.unique_im_nr(counter)   = unique_ref_im(rr-1,sub);
        end
        
        % check facing direction        
        if (im_order.abs_rot(counter) < 45) || (im_order.abs_rot(counter) > 135)
            facing_dir = {'sideways'};
        elseif (im_order.abs_rot(counter) > 45) || (im_order.abs_rot(counter) < 135)
            facing_dir = {'forward'};
        end
        im_order.facing_dir(counter)     = facing_dir;
        
        % we define rotations relative to the reference object rotation
        if (im_order.abs_rot(counter)-90) < 0 && im_order.rel_rot(counter) < 0
            im_order.rel_rot_name(counter) = {'rightward'};
        elseif (im_order.abs_rot(counter)-90) < 0 && im_order.rel_rot(counter) > 0
            im_order.rel_rot_name(counter) = {'leftward'};
        elseif (im_order.abs_rot(counter)-90) > 0 && im_order.rel_rot(counter) > 0
            im_order.rel_rot_name(counter) = {'rightward'};
        elseif (im_order.abs_rot(counter)-90) > 0 && im_order.rel_rot(counter) < 0
            im_order.rel_rot_name(counter) = {'leftward'};
        elseif (im_order.abs_rot(counter)-90) == 0 && im_order.base_rot(counter) < 0 && im_order.rel_rot(counter) < 0
            im_order.rel_rot_name(counter) = {'rightward'};
        elseif (im_order.abs_rot(counter)-90) == 0 && im_order.base_rot(counter) < 0 && im_order.rel_rot(counter) > 0
            im_order.rel_rot_name(counter) = {'leftward'};
        elseif (im_order.abs_rot(counter)-90) == 0 && im_order.base_rot(counter) > 0 && im_order.rel_rot(counter) > 0
            im_order.rel_rot_name(counter) = {'rightward'};
        elseif (im_order.abs_rot(counter)-90) == 0 && im_order.base_rot(counter) > 0 && im_order.rel_rot(counter) < 0
            im_order.rel_rot_name(counter) = {'leftward'};
        elseif im_order.rel_rot(counter) == 0
            im_order.rel_rot_name(counter) = {'none'};
        end
        
        % Update counter
        counter = counter+1;
    end
end

if params.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(params.stim.obj.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',params.stim.obj.stimfile,datestr(now,30))),'objects','masks','info','im_order','-v7.3');
end

%
%     figure(1); clf;
%     figure(2); clf;
%     figure(3); clf;
%     figure(4); clf;
%     for ii = 1:size(objects,5)
%         figure(1); subplot(3,3,ii); hold all;
%         I = imshow(objects(:,:,:,1,ii),[1 255]);
%         I.AlphaData = masks(:,:,1,ii)>0;
%         axis image
%         
%         figure(2); subplot(3,3,ii); hold all;
%         I = imshow(objects(:,:,:,2,ii),[1 255]);
%         I.AlphaData = masks(:,:,2,ii)>0;
%         axis image
%         
%         figure(3); subplot(3,3,ii); hold all;
%         I = imshow(objects(:,:,:,3,ii),[1 255]);
%         I.AlphaData = masks(:,:,3,ii)>0;
%         axis image
%         
%         figure(4); subplot(3,3,ii); hold all;
%         I = imshow(objects(:,:,:,4,ii),[1 255]);
%         I.AlphaData = masks(:,:,4,ii)>0;
%         axis image
%     end
return

