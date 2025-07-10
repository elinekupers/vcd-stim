function [objects, masks, objects_catch, masks_catch, info] = vcd_objects(params, verbose, store_imgs)
% VCD function:
%   [objects, masks, objects_catch, mask_catch, info] = vcd_objects(params, verbose, store_imgs)
%
% Purpose:
%   Luminance-controlled, load, and store rotated object images for VCD
%   experimental display. By default, we do not preprocess images and
%   assume this has already been done.
%
%   If you haven't done so already, you can control object's luminance by
%   running the script: s_preprocessObjects.m. This line is commented out
%   by default, assuming the user will only do this once at the beginning
%   and does not want to redo this step automatically everytime this
%   function is called. Preprocessing step requires "raw" png files living
%   in the folder:
%    ./vcd-stim/workspaces/stimuli/RAW/vcd_objects/all_to_process/*
%   Once preprocess, png files are expected to live in the folder defined
%   by "params.stim.obj.indivobjfile", which is by default:
%    ./vcd-stim/workspaces/stimuli/<dispname>/vcd_objects_2degstep_lumcorrected/*'
%
%   This function will resize a subselection of preprocessed
%   object pngs if requested and reorder the files according to their
%   unique image number. It will also create a .csv info file, which
%   contains information about the object stimuli used in VCD. See OUTPUTS
%   below for more details. The csv file will be stored in the folder
%   defined by params.stim.obj.infofile, e.g.:
%       fullfile(vcd_rootPath,'workspaces/objects_info.csv')
%
%   Objects can rotate between 0 and 180 degrees, where 0, 90, and 180
%   degrees correspond to right, forwards, and left facing, respectively.
%   Facing directions of objects in the following files are:
%       * 1-23  (0-44 deg)   : face sideways
%       * 24-46 (46-90 deg)  : face forward
%       * 47-68 (92-134 deg) : face forward
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
% INPUTS:
%  params       :   struct with stimulus params
%    *** this function requires the following struct fields ***
%    bckgrnd_grayval                : (double, integral) background gray value (128)
%    stim.obj.img_sz_pix            : (double, integral) size of image
%    support (pixels) f
%                                       . Must be an even number.
%    stim.obj.crop_img_sz_pix       : (double, integral) size of cropped image support (pixels) of the
%                                       dot. Must be an even number.
%    stim.obj.super_cat             : (double, integral) superordinate object category (ranges between 1-5)
%    stim.obj.basic_cat             : (double, integral) basic object category (ranges between 1-2)
%    stim.obj.sub_cat               : (double, integral) subordinate object category (ranges between 1-3)
%    stim.obj.facing_dir_deg        : (double, integral) list of facing directions for each object (deg)
%    stim.obj.affordance            : (double, integral) number associated with each affordance label for each object,
%                                        where 1: greet, 2: grasp, 3: enter, and 4: observe.
%    stim.obj.x0_pix                : (double, integral) horizontal center positions for left/right object apertures (pixels)
%    stim.obj.y0_pix                : (double, integral) vertical center positions for left/right object apertures (pixels)
%    stim.obj.unique_im_nrs_core    : (double, integral) unique stimulus numbers for core images
%    stim.obj.unique_im_nrs_wm_test : (double, integral) unique stimulus numbers for working memory test images
%  verbose           : (logical) show debug figures
%  store_imgs        : (logical) store stimuli and debug figures as pngs
%
% OUTPUTS:
%  objects      :   uint8 5D array with images of individual object stimuli
%                   dimensions are:
%                   width (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x height (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x 3 (rgb)
%                   x 16 object's subordinate categories (e.g., damon, lisa, etc.)
%                   x 5 rotations (0, -24, -12, +12, +24 degrees).
%  masks        :   uint8 4D array with alpha transparency masks for
%                   Psychtoolbox to remove the object from the background.
%                   Dimensions are:
%                   width (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x height (BOLDscreen: 512 pixels; PProom: 512 pixels)
%                   x 16 object's subordinate categories (e.g., damon, lisa, etc.)
%                   x 5 rotations (0, -24, -12, +12, +24 degrees).
% objects_catch :   uint8 5D array with object catch images for each
%                   individual object, dimension are:
%                   height (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x width (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x 3 (rgb)
%                   x 16 object subordinate categories (e.g., damon, lisa, etc.)
%                   x 18 catch rotations sampled from all rotation objects
%                   (0-180 degrees), except the core rotation of that
%                   particular object.
% masks_catch :   uint8 4D array with alpha transparency masks for
%                   Psychtoolbox to remove the object from the background.
%                   Dimensions are:
%                   width (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x height (BOLDscreen: 512 pixels; PProom: 512 pixels) 
%                   x 16 object's subordinate categories (e.g., damon, lisa, etc.)
%                   x 18 catch rotations sampled from all rotation objects
%                   (0-180 degrees), except the core rotation of that
%                   particular object.
%  info         :   386x12 table with information about object stimuli.
%                   Columns include:
%      filename             : (cell) filename of original png, e.g.: 
%                             {'faces_damonwayans_rot30.png'} or {'places_watertower_rot66.png'}
%      unique_im            : (double) unique image nr for each object:
%                              range 65-80 (core images), 239-302 (wm test images).
%      stim_pos_name        : (cell) stimulus position ({'left'} or {'right'}),
%                             same as stim_loc_i but human readable.
%      stim_pos             : (double) stimulus position index. 1=left, 2=right
%      super_cat_name       : (cell) superordinate semantic category label:
%                             {'human','animal','object','food','place'}
%      super_cat            : (double) same as superordinate, but indexed by nr
%                             1:'human' through 5: 'place'
%      basic_cat_name       : (cell) basic semantic category label for each
%                             superordinate semantic categorycategory:
%                             {'facefemale','facefemale', ...
%                             'small','large',...
%                             'tool','vehicle',...
%                             'man-made','produce',...
%                             'building'};
%      basic_cat           : (double) same as basic, but indexed by nr
%                             1: facefemale, 2: facefemale
%                             1: small, 2: large
%                             1: tool, 2: vehicle
%                             1: man-made, 2: produce
%                             1: building.
%      sub_cat_name        : (cell) subordinate semantic category label, name
%                            for each individual core object:
%                            {'damon','lisa','sophia','parrot','cat','bear',...
%                            'giraffe','drill','brush','pizza','banana',...
%                            'bus','suv''church','house','watertower'};
%      sub_cat             : (double) same as subordinate, but indexed by nr
%                            1:damon, 1:lisa, 2:sophia,
%                            1:parrot, 2:cat, 1:bear, 2: giraffe,
%                            1:drill, 2: brush,
%                            1:bus, 2: suv,
%                            1:pizza, 1: banana,
%                            1:church, 2: house, 3: watertower
%      affordance_name     : (cell) affordance label for each object, names
%                            are one of 4: {'greet','grasp','enter','observe'};
%      affordance_cat      : (double) same as affordance_name, but
%                            computer-readable, where 1: greet, 2: grasp,
%                            3: enter, and 4: observe.
%      is_specialcore     : (logical) whether the object is part of the
%                            subselected stimuli used in imagery and
%                            long-term memory task.
%      base_rot           : (double) rotation (in deg) of object shown in
%                            the first stimulus array (prior to applying
%                            relative rotation ("rel_rot). base_rot is the
%                            same for core images, and only relevant for wm
%                            test images only. Values range from 26-154 deg,
%                            in evenly steps of 8 deg (excluding 90 deg).
%      rot_abs            : (double) absolute rotation (deg) of object.
%                            ranges from 26-154, in evenly steps of 8 deg
%                            (excluding 90 degrees).
%                            0 deg = profile view facing to the right.
%                            90 deg = forward facing view.
%                            180 deg = profile view facing to the left.
%      rot_rel          : (double) relative rotation (deg) of object from
%                           core object rotation, ranges from -24, -12, +12,
%                           +24 deg. Only relevant for wm test images.
%      facing_dir       : (cell) facing direction, 'forward' directions are
%                           between 46-134 degree. 'sideways' directions
%                           are between 0-44 and 136-180 degrees. Rotations
%                           of 45, 90, and 135 are ill-defined and don't occur.
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
%     is_objectcatch    : (logical) whether the object is part of the
%                            catch object rotations used in PC-OBJ task.
%
% Written by Eline Kupers 2024/12, updated 2025/04 and 2025/06

%% !!! IMPORTANT !!! 
% If you haven't done so yet, first preprocess object images (once):
% s_preprocessObjects; <-- Check variable names inside the script before you run this!!!


%% Prepare info table

% Read in table with filenames
info = table();
[info.filename, info.unique_im] = vcd_getOBJfilenames;

% Define superordinate and basic categories, and number of exemplars per basic category
[~,idx0]     = intersect(info.unique_im, params.stim.obj.unique_im_nrs_core);
core_im_name = info.filename(idx0);

% Define superordinate and basic categories, and number of exemplars per basic category
superordinate = cat(2, repmat(params.stim.obj.super_cat(1),1,3), ...
                       repmat(params.stim.obj.super_cat(2),1,4),...
                       repmat(params.stim.obj.super_cat(3),1,4),...
                       repmat(params.stim.obj.super_cat(4),1,2),...
                       repmat(params.stim.obj.super_cat(5),1,3))';                    % 5 superordinate categories
basic               = catcell(2,params.stim.obj.basic_cat)';         % 1,2,or 3 basic categories
subordinate         = catcell(2,params.stim.obj.sub_cat)';           % 16 subordiante categories: each object
affordance          = catcell(2,params.stim.obj.affordance)';
[~,superordinate_i] = ismember(superordinate, params.stim.obj.super_cat);

basic_i = []; subordinate_i = [];
for bb = 1:length(params.stim.obj.super_cat)
    [~,bi] = ismember(basic,unique(params.stim.obj.basic_cat{bb},'stable'));
    basic_i = cat(1,basic_i,bi(bi>0));
    
    sb = 1:length(params.stim.obj.sub_cat{bb});
    subordinate_i = cat(2,subordinate_i,sb);
end

% Define affordance index nr
[~,affordance_i]    = ismember(affordance,{'greet','grasp','enter','observe'});

% Get info about wm test images
n_obj        = length(core_im_name);
n_wm_changes = length(params.stim.obj.delta_from_ref);

% Get info about object catch images
n_catch_images = size(params.stim.obj.catch_rotation,2);

% Predefine info table
info.stim_pos_name   = cat(1, repmat({'left';'right'}, length(core_im_name)/2,1) , ... core im
    repmat(repelem({'left';'right'},n_wm_changes), length(core_im_name)/2,1), ... wm im
    repmat(repelem({'left';'right'},n_catch_images), length(core_im_name)/2,1)); % catch im
info.stim_pos        = cat(1, repmat([1;2],length(core_im_name)/2,1) , ....core im
    repmat(repelem([1;2],n_wm_changes), length(core_im_name)/2,1), ...   wm im
    repmat(repelem([1;2],n_catch_images), length(core_im_name)/2,1)); % catch im
info.super_cat_name  = cat(1,superordinate, ... core im
    repelem(superordinate,n_wm_changes), ... wm im
    repelem(superordinate,n_catch_images)); % catch im
info.super_cat       = cat(1,superordinate_i, ... core im
    repelem(superordinate_i,n_wm_changes), ... wm im
    repelem(superordinate_i,n_catch_images)); % catch im
info.basic_cat_name  = cat(1,basic, ... core im
    repelem(basic,n_wm_changes),... wm im
    repelem(basic,n_catch_images)); % catch im
info.basic_cat       = cat(1,basic_i, ... core im
    repelem(basic_i,n_wm_changes),... wm im
    repelem(basic_i,n_catch_images)); % catch im
info.sub_cat_name    = cat(1,subordinate, ... core im
    repelem(subordinate,n_wm_changes), ...wm im
    repelem(subordinate,n_catch_images)); % catch im
info.sub_cat         = cat(1,subordinate_i', ... core im
    repelem(subordinate_i', n_wm_changes),... % wm im
    repelem(subordinate_i', n_catch_images)); % wm im
info.affordance_name = cat(1,affordance,repelem(affordance,n_wm_changes),repelem(affordance,n_catch_images));
info.affordance_cat  = cat(1,affordance_i,repelem(affordance_i,n_wm_changes),repelem(affordance_i,n_catch_images));
info.is_specialcore  = ismember(info.unique_im, params.stim.obj.unique_im_nrs_specialcore);  % is_specialcore (logical)
info.is_objectcatch  = ismember(info.unique_im, params.stim.obj.unique_im_nrs_objcatch);  % is_specialcore (logical)

% Define all rotations
catch_rotations = params.stim.obj.catch_rotation';
delta_rotations = [0, params.stim.obj.delta_from_ref];
info.base_rot   = cat(1, params.stim.obj.facing_dir_deg', repelem(params.stim.obj.facing_dir_deg,n_wm_changes)', catch_rotations(:));
rot_abs_wm      = params.stim.obj.facing_dir_deg + params.stim.obj.delta_from_ref';
info.abs_rot    = cat(1, params.stim.obj.facing_dir_deg', rot_abs_wm(:), catch_rotations(:));  % wm: alternating 2 view ± 4 steps spaced 2 deg apart
info.rel_rot    = cat(1, repmat(delta_rotations(1),n_obj,1),repmat(delta_rotations(2:end)',n_obj,1), repmat(delta_rotations(1),n_obj*n_catch_images,1));

% Get cropped rescaled image to reduce glitches when drawing stimuli on the screen
crop_im = zeros(params.stim.obj.crop_img_sz_pix, params.stim.obj.crop_img_sz_pix); % 512 pixels x 512 pixels

% Preallocate space
extent        = params.stim.obj.crop_img_sz_pix.*params.stim.obj.dres;
objects       = uint8(ones(extent, extent, 3, n_obj, length(delta_rotations)));
masks         = uint8(ones(extent, extent, n_obj, length(delta_rotations)));
objects_catch = uint8(ones(extent, extent, 3, n_obj, n_catch_images));
masks_catch   = uint8(ones(extent, extent, n_obj, n_catch_images));

% keep track of individual objects
obj_sub      = cat(2,[1:16],repelem([1:16],n_wm_changes), repelem([1:16],n_catch_images))';
objcatch_sub = cat(2,zeros(1,n_obj),zeros(1,n_obj*n_wm_changes), repmat([1:n_catch_images],1,n_obj))';

% loop over each object rotation
fprintf('[%s]: Load and store object stimuli',mfilename); tic
for rr = 1:length(info.rel_rot)
    fprintf('.')

    sub = obj_sub(rr);
    
    % Get file name
    rot = 1+((info.base_rot(rr) + info.rel_rot(rr))/2);
    fname = sprintf('*%s*_rot%02d.png', subordinate{sub},rot);
    d = dir(fullfile(params.stim.obj.indivobjfile,fname));
    if isempty(d)
        error('[%s]: Can''t find image file!')
    end
    
    if info.rel_rot(rr) == 0
        assert(strcmp(info.filename{rr},d.name))
    else
        assert(strcmp(cell2mat(info.filename(rr)),d.name));
    end
    
    % Load image
    [im0, ~, alpha0] = imread(fullfile(d.folder,d.name),'BackgroundColor','none');
    
    
    % If needed, crop image and alpha transparency mask
    if size(extent,1) < size(im0,1) || size(extent,2) < size(im0,2)
        im1    = double(im0);
        alpha1 = double(alpha0);
        im0    = uint8(placematrix(crop_im, im1, []));    % empty third input argument means enter original stimuli with respect to cropping image.
        alpha0 = uint8(placematrix(crop_im, alpha1, [])); % empty third input argument means enter original stimuli with respect to cropping image.
        clear im1 alpha1
    end
    
    % If we need to rescale images:
    % convert to double, [0-1] lum range (assume range is 1-255),
    % and square as images were made with gamma = 2.
    if ~isequal(params.stim.obj.dres,1) && ~isempty(params.stim.obj.dres)
        error('[%s]: wtf. we don''t resize object images anymore. Please check params.stim.obj.dres value',mfilename)
    end
   
    % in case we want to debug:
    % figure(99); clf; imagesc(im,'AlphaData',alpha_im); colormap gray; colorbar; axis image
    if info.is_objectcatch(rr)
        objects_catch(:,:,:,sub,objcatch_sub(rr)) = repmat(im0, [1 1 3]);
        masks_catch(:,:,sub,objcatch_sub(rr))     = alpha0;
    else
        objects(:,:,:,sub,info.rel_rot(rr) == delta_rotations) = repmat(im0, [1 1 3]);
        masks(:,:,sub,info.rel_rot(rr) == delta_rotations)     = alpha0;
    end
    
    clear im0 alpha0 
    
    % check facing direction, when the offset may technically result in
    % a flipping of the facing direction, go by original core facing
    % direction, as subjects will never perform the PC task on working
    % memory test images.
    if info.rel_rot(rr) == 0
        if (info.abs_rot(rr) < 45) || (info.abs_rot(rr) > 135)
            info.facing_dir(rr) = {'sideways'};
        elseif (info.abs_rot(rr) > 45) || (info.abs_rot(rr) < 135)
            info.facing_dir(rr) = {'forward'};
        end
    elseif info.rel_rot(rr) ~= 0
        og_rot = info.base_rot(rr);
        if (og_rot < 45) || (og_rot > 135)
            info.facing_dir(rr) = {'sideways'};
        elseif (og_rot > 45) || (og_rot < 135)
            info.facing_dir(rr) = {'forward'};
        end
    end
    
    % we define rotations relative to the reference object rotation
    % if abs rotation is between 0-90 and rel rotation is negative
    if (info.abs_rot(rr)-90) < 0 && info.rel_rot(rr) < 0
        info.rel_rot_name(rr) = {'rightward'};
        % if abs rotation is between 0-90 and rel rotation is positive
    elseif (info.abs_rot(rr)-90) < 0 && info.rel_rot(rr) > 0
        info.rel_rot_name(rr) = {'leftward'};
        % if abs rotation is between 90-180 and rel rotation is negative
    elseif (info.abs_rot(rr)-90) > 0 && info.rel_rot(rr) < 0
        info.rel_rot_name(rr) = {'rightward'};
        % if abs rotation is between 90-180 and rel rotation is positive
    elseif (info.abs_rot(rr)-90) > 0 && info.rel_rot(rr) > 0
        info.rel_rot_name(rr) = {'leftward'};
        % CORNER CASE 1: if abs rotation is exactly 90, but core rotation was between 0-90, and rel rotation is negative
    elseif (info.abs_rot(rr)-90) == 0 && info.base_rot(rr) < 0 && info.rel_rot(rr) < 0
        info.rel_rot_name(rr) = {'rightward'};
        % CORNER CASE 2:  if abs rotation is exactly 90, but core rotation was between 0-90, and rel rotation is positive
    elseif (info.abs_rot(rr)-90) == 0 && info.base_rot(rr) < 0 && info.rel_rot(rr) > 0
        info.rel_rot_name(rr) = {'leftward'};
        % CORNER CASE 3:  if abs rotation is exactly 90, but core rotation was between 90-180, and rel rotation is negative
    elseif (info.abs_rot(rr)-90) == 0 && info.base_rot(rr) > 0 && info.rel_rot(rr) < 0
        info.rel_rot_name(rr) = {'rightward'};
        % CORNER CASE 4:  if abs rotation is exactly 90, but core rotation was between 90-180, and rel rotation is positive
    elseif (info.abs_rot(rr)-90) == 0 && info.base_rot(rr) > 0 && info.rel_rot(rr) > 0
        info.rel_rot_name(rr) = {'leftward'};
    elseif info.rel_rot(rr) == 0
        info.rel_rot_name(rr) = {'none'};
    end
    
end

fprintf('\n[%s]: Done! ', mfilename)
toc
fprintf('\n');

%% Store stimuli if requested
if store_imgs
    fprintf('[%s]:Storing images..\n',mfilename)
    saveDir = fileparts(fullfile(params.stim.obj.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',params.stim.obj.stimfile,datestr(now,30))),'objects','masks','objects_catch','masks_catch','info','-v7.3');
    
    saveDir = fileparts(fullfile(params.stim.obj.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',params.stim.obj.infofile,datestr(now,30))))
end

%% Visualize stimuli if requested

if store_imgs
    saveFigDir = fullfile(vcd_rootPath,'figs',params.disp.name, 'obj');
    if ~exist(saveFigDir,'dir'); mkdir(saveFigDir); end
    
    % Export PNG -- APPLY WEIGHTED ALPHA MASK BEFORE STORING.
    % WE DO NOT STORE ALPHA MASK SEPARATELY IN PNG
    for objectNr = 1:size(objects,4)
        for rot = 1:size(objects,5)
            
            if rot == 1
                im_nr = info.unique_im(objectNr);
            else
                im_nr = info.unique_im(size(objects,4) + ((objectNr-1)*4) + rot-1);
            end
            
            im_w_alpha = uint8( (double(objects(:,:,:,objectNr,rot)) .* double(masks(:,:,objectNr,rot))/255)  +  (128 .* (1 - double(masks(:,:,objectNr,rot))/255)) );
            imwrite(im_w_alpha, fullfile(vcd_rootPath,'figs',params.disp.name,...
                'obj',sprintf('%04d_object%02d_delta%02d.png',im_nr, objectNr, rot-1)));
        end
    end
    
    catch_images = find(info.is_objectcatch);
    for objectNr = 1:size(objects_catch,4)
        for rot = 1:size(objects_catch,5)
            
            idx   = catch_images(((objectNr-1)*n_catch_images) + rot);
            im_nr = info.unique_im(idx);

            im_w_alpha = uint8( (double(objects_catch(:,:,:,objectNr,rot)) .* double(masks_catch(:,:,objectNr,rot))/255)  +  (128 .* (1 - double(masks_catch(:,:,objectNr,rot))/255)) );
            imwrite(im_w_alpha, fullfile(vcd_rootPath,'figs',params.disp.name,...
                'obj',sprintf('%04d_object%02d_catch%02d.png',im_nr, objectNr, rot)));
        end
    end
end

% if verbose
%
%     debug figure for wm test images
%     counter = 1;
%     figure(99); set(gcf, 'Position', [300   584   868   753]);
%     for objectNr = 1:size(objects,4)
%         for rot = 1:size(objects,5)
%
%             if rot == 1
%                 im_nr = info.unique_im(objectNr);
%             else
%                 im_nr = info.unique_im(size(objects,4) + ((objectNr-1)*4) + rot-1);
%             end
%             cla;
%             I = imshow(objects(:,:,:,objectNr,rot),[1 255]);
%             I.AlphaData = masks(:,:,objectNr,rot);
%             axis image
%
%             title(sprintf('Object:%04d Rot:%02d Delta:%02d',objectNr,info.abs_rot(counter),info.rel_rot(counter)), 'FontSize',20);
%
%             if rot == 1, dd = 0;
%             else, dd = rot; end
%
%             if store_imgs
%               saveFigDir = fullfile(vcd_rootPath,'figs',params.disp.name, 'obj','visual_checks');
%               if ~exist(saveFigDir,'dir'); mkdir(saveFigDir); end
%                 % Print figure with axes
%                 print(fullfile(saveFigDir,sprintf('%02d_object%02d_delta%02d', im_nr,objectNr, dd)),'-dpng','-r150');
%             end
%
%             update counter
%             counter = counter+1;
%         end
%     end
% end

return

