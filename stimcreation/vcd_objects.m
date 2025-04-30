function [objects, masks, info] = vcd_objects(params)
% VCD function:
%   [objects, masks, info] = vcd_objects(params)
%
% Purpose:
%   Luminance-correct, load, and store object images for VCD experimental
%   display. By default, we do not preprocess images and assume this has
%   already been done.
%
%   After preprocessing, you have two files required to execute the rest of
%   the function:
%    1. a set of individual object image files, where the folder is defined
%       by params.stim.obj.indivobjfile, e.g.:
%       fullfile(vcd_rootPath,'workspaces/stimuli/vcd_complex_objects/*_rot*.png)
%    2. A csv file defined in params.stim.obj.infofile, e.g.: 
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
%   fullfile(vcd_rootPath,'workspaces','stimuli',<disp_name>,'objects.mat').
%
%   We also create a struct called "im_order", which contains subset of 
%   the info table in the order they were loaded, which makes it easier to
%   keep track and double check indexing.
%    - object_name: {1×80 cell} - same as filename in info
%    - base_rot: [1×80 double]  - base rotation of core object
%    - abs_rot: [1×80 double]   - absolute rotation of object (after adding delta)
%    - rel_rot: [1×80 double]   - rel rotation of object (i.e. the delta)
%    - unique_im_nr: [1×80 double] - same as unique_im_nr in info
%    - facing_dir: {1×80 cell}  - same as facing_dir in info
%    - rel_rot_name: {1×80 cell} - same as rel_rot_name in info
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
%  info         :   Loaded csv from workspaces/stimuli/ with png file names
%                   and object category information. 
%      filename         : (cell) filename of original png, e.g.: {'faces_damonwayans_rot31.png'} or {'places_watertower_rot75.png'}                
%      unique_im        : (double) unique image nr for each object: 
%                         range 65-80, 367-430. 
%      stim_pos_i       : (double) stimulus position index. 1=left, 2=right
%      stim_pos         : (cell) stimulus position, same as stim_loc_i but
%                           human readable ({'left'} or {'right'})
%      superordinate    : (cell) superordinate semantic category label:
%                           {'human','animal','object','place'}
%      superordinate_i  : (double) same as superordinate, but indexed by nr
%                           1:'human' through 4: 'place'    
%      basic            : (cell) basic semantic category label for each
%                          superordinate semantic categorycategory:
%                           {'facefemale','facefemale', ...
%                           'small','large',...
%                           'tool',food','vehicle',...
%                           'building'};
%      basic_i          : (double) same as basic, but indexed by nr
%                           1: facefemale, 2: facefemale
%                           1: small, 2: large
%                           1: tool, 2: food, 3: vehicle
%                           1: building.
%      subordinate      : (cell) subordinate semantic category label, name 
%                           for each individual core object: 
%                           {'damon','lisa','sophia','parrot','cat','bear',...
%                            'giraffe','drill','brush','pizza','banana',...
%                            'bus','suv''church','house','watertower'};
%      subordinate_i    : (double) same as subordinate, but indexed by nr
%                           1:damon through 16: watertower. 
%      rot_abs          : (double) absolute rotation (deg) of object.
%                          ranges from 10-170, in evenly steps of 10 deg
%                          (excluding 90 degrees). 
%                           0 deg = profile view facing to the right. 
%                           90 deg = forward facing view. 
%                           180 deg = profile view facing to the left. 
%      rot_rel          : (double) relative rotation (deg) of object from
%                           core object rotation, ranges from -8, -4, +4,
%                           +8 deg.
%      facing_dir       : (cell) facing direction, 'forward' directions are
%                           between 46-134 degree. 'sideways' directions
%                           are between 0-44 and 136-180 degrees. Rotations
%                           of 45 and 135 are ill-defined and does not occur.
%      rel_rot_name     : (cell) relative rotation name for WM test images
%                           whether the relative rotation results in the
%                           object facing more left or rightwards from the
%                           perspective of the observer. 
%                           'rightwards' is when absolute rotation is
%                           between 0-90 and relative rotation is > 0, OR
%                           when absolute rotation is between 91-180 and
%                           relative rotation < 0 deg. 
%                           'leftwards' is when absolute rotation is
%                           between 0-90 and relative rotation is < 0, OR
%                           when absolute rotation is between 91-180
%                           relative rotation is > 0 deg.
%      is_in_img_ltm    : (logical) whether the object is part of the
%                           subselected stimuli used in imagery and
%                           long-term memory task.
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
unique_ref_im   = reshape(params.stim.obj.unique_im_nrs_wm_test, [],sum(idx0));

% Get (rescaled) image extent
extent   = params.stim.obj.img_sz_pix.*params.stim.obj.dres;

% Preallocate space
objects = uint8(ones(extent, extent,3, length(subordinate(idx0)),length(rotations)));
masks   = uint8(ones(extent, extent, length(subordinate(idx0)),length(rotations)));

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
            im_order.unique_im(counter)   = params.stim.obj.unique_im_nrs_core(sub);
        else 
            im_order.unique_im(counter)  = unique_ref_im(rr-1,sub);
        end
        
        % check facing direction, when the offset may technically result in
        % a flipping of the facing direction, go by original core facing
        % direction, as subjects will never perform the PC task on working
        % memory test images.
        if im_order.rel_rot(counter) == 0
            if (im_order.abs_rot(counter) < 45) || (im_order.abs_rot(counter) > 135)
                im_order.facing_dir(counter) = {'sideways'};
            elseif (im_order.abs_rot(counter) > 45) || (im_order.abs_rot(counter) < 135)
                im_order.facing_dir(counter) = {'forward'};
            end
        elseif im_order.rel_rot(counter) ~= 0
            og_rot = im_order.base_rot(counter);
            if (og_rot < 45) || (og_rot > 135)
                im_order.facing_dir(counter) = {'sideways'};
            elseif (og_rot > 45) || (og_rot < 135)
                im_order.facing_dir(counter) = {'forward'};
            end
        end

        % we define rotations relative to the reference object rotation
            % if abs rotation is between 0-90 and rel rotation is negative
        if (im_order.abs_rot(counter)-90) < 0 && im_order.rel_rot(counter) < 0 
            im_order.rel_rot_name(counter) = {'rightward'};
            % if abs rotation is between 0-90 and rel rotation is positive
        elseif (im_order.abs_rot(counter)-90) < 0 && im_order.rel_rot(counter) > 0 
            im_order.rel_rot_name(counter) = {'leftward'};
            % if abs rotation is between 90-180 and rel rotation is negative
        elseif (im_order.abs_rot(counter)-90) > 0 && im_order.rel_rot(counter) < 0 
            im_order.rel_rot_name(counter) = {'rightward'};
            % if abs rotation is between 90-180 and rel rotation is positive
        elseif (im_order.abs_rot(counter)-90) > 0 && im_order.rel_rot(counter) > 0 
            im_order.rel_rot_name(counter) = {'leftward'};
            % CORNER CASE 1: if abs rotation is exactly 90, but core rotation was between 0-90, and rel rotation is negative
        elseif (im_order.abs_rot(counter)-90) == 0 && im_order.base_rot(counter) < 0 && im_order.rel_rot(counter) < 0 
            im_order.rel_rot_name(counter) = {'rightward'};
            % CORNER CASE 2:  if abs rotation is exactly 90, but core rotation was between 0-90, and rel rotation is positive
        elseif (im_order.abs_rot(counter)-90) == 0 && im_order.base_rot(counter) < 0 && im_order.rel_rot(counter) > 0
            im_order.rel_rot_name(counter) = {'leftward'};
            % CORNER CASE 3:  if abs rotation is exactly 90, but core rotation was between 90-180, and rel rotation is negative
        elseif (im_order.abs_rot(counter)-90) == 0 && im_order.base_rot(counter) > 0 && im_order.rel_rot(counter) < 0
            im_order.rel_rot_name(counter) = {'rightward'};
            % CORNER CASE 4:  if abs rotation is exactly 90, but core rotation was between 90-180, and rel rotation is positive
        elseif (im_order.abs_rot(counter)-90) == 0 && im_order.base_rot(counter) > 0 && im_order.rel_rot(counter) > 0
            im_order.rel_rot_name(counter) = {'leftward'};
        elseif im_order.rel_rot(counter) == 0
            im_order.rel_rot_name(counter) = {'none'};
        end
        
        % Update counter
        counter = counter+1;
    end
end

%% Store stimuli if requested
if params.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(params.stim.obj.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',params.stim.obj.stimfile,datestr(now,30))),'objects','masks','info','im_order','-v7.3');
end

%% Visualize stimuli if requested
if params.verbose
   
    counter = 1;
    if params.stim.store_imgs
        saveFigDir = fullfile(vcd_rootPath,'figs',params.disp.name, 'obj','visual_checks');
        if ~exist(saveDir,'dir'); mkdir(saveDir); end
    end
    
    figure(99); set(gcf, 'Position', [300   584   868   753]);
    for objectNr = 1:size(objects,4)
        for rot = 1:size(objects,5)
            
            if rot == 1
                im_nr = info.unique_im_nr(objectNr);
            else 
                im_nr = info.unique_im_nr(size(objects,4) + ((objectNr-1)*4) + rot-1);
            end
            cla;
            I = imshow(objects(:,:,:,objectNr,rot),[1 255]);
            I.AlphaData = masks(:,:,objectNr,rot);
            axis image
            
            title(sprintf('Object:%02d Rot:%02d Delta:%02d',objectNr,im_order.abs_rot(counter),im_order.rel_rot(counter)), 'FontSize',20);
            
            if rot == 1, dd = 0;
            else, dd = rot; end
            
            if params.stim.store_imgs
                % Print figure with axes
                print(fullfile(saveFigDir,sprintf('%02d_object%02d_delta%02d', im_nr,objectNr, dd)),'-dpng','-r150');
                % Export PNG
                imwrite(objects(:,:,:,objectNr,rot), fullfile(vcd_rootPath,'figs',params.disp.name,...
                    'obj',sprintf('%02d_object%02d_delta%02d.png',im_nr, objectNr, dd)),'Alpha',masks(:,:,objectNr,rot));
            end
            
            % update counter
            counter = counter+1;
        end
    end
end

return

