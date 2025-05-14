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
%   We also create a struct called "info", which contains subset of 
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
%                           {'human','animal','object','food','place'}
%      superordinate_i  : (double) same as superordinate, but indexed by nr
%                           1:'human' through 5: 'place'    
%      basic            : (cell) basic semantic category label for each
%                          superordinate semantic categorycategory:
%                           {'facefemale','facefemale', ...
%                           'small','large',...
%                           'tool','vehicle',...
%                           'man-made','produce',...
%                           'building'};
%      basic_i          : (double) same as basic, but indexed by nr
%                           1: facefemale, 2: facefemale
%                           1: small, 2: large
%                           1: tool, 2: vehicle
%                           1: man-made, 2: produce
%                           1: building.
%      subordinate      : (cell) subordinate semantic category label, name 
%                           for each individual core object: 
%                           {'damon','lisa','sophia','parrot','cat','bear',...
%                            'giraffe','drill','brush','pizza','banana',...
%                            'bus','suv''church','house','watertower'};
%      subordinate_i    : (double) same as subordinate, but indexed by nr
%                           1:damon, 1:lisa, 2:sophia, 
%                           1:parrot, 2:cat, 1:bear, 2: giraffe,
%                           1:drill, 2: brush,
%                           1:bus, 2: suv,
%                           1:pizza, 1: banana,
%                           1:church, 2: house, 3: watertower
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
%      is_specialcore    : (logical) whether the object is part of the
%                           subselected stimuli used in imagery and
%                           long-term memory task.
%
% Written by Eline Kupers 2024/12, updated 2025/04

%% If you haven't done so yet, first preprocess object images:
% s_preprocessObjects.m


%% Prepare images

% Read in table with filenames
info = table();
[info.filename, info.unique_im] = vcd_getOBJfilenames;

% Define superordinate and basic categories, and number of exemplars per basic category
[~,idx0] = intersect(info.unique_im, params.stim.obj.unique_im_nrs_core);
[~,idx1] = intersect(info.unique_im,params.stim.obj.unique_im_nrs_wm_test);
% [~,idx2] = intersect(info.unique_im,params.stim.obj.unique_im_nrs_img_test);

core_im_name     = info.filename(idx0);
wm_test_im_name  = info.filename(idx1);
% img_test_im_name = info.filename(idx2);

% Define superordinate and basic categories, and number of exemplars per basic category
superordinate = cat(2, repmat(params.stim.obj.super_cat(1),1,3), ...
                          repmat(params.stim.obj.super_cat(2),1,4),...
                          repmat(params.stim.obj.super_cat(3),1,4),...
                          repmat(params.stim.obj.super_cat(4),1,2),...
                          repmat(params.stim.obj.super_cat(5),1,3))';                    % 5 superordinate categories
basic         = catcell(2,params.stim.obj.basic_cat)';         % 1,2,or 3 basic categories
subordinate   = catcell(2,params.stim.obj.sub_cat)';           % 16 subordiante categories: each object
affordance    = catcell(2,params.stim.obj.affordance)';
[~,superordinate_i] = ismember(superordinate, params.stim.obj.super_cat);

basic_i = []; subordinate_i = [];
for bb = 1:length(params.stim.obj.super_cat)
    [~,bi] = ismember(basic,unique(params.stim.obj.basic_cat{bb},'stable'));
    basic_i = cat(1,basic_i,bi(bi>0));
    
    sb = 1:length(params.stim.obj.sub_cat{bb});
    subordinate_i = cat(2,subordinate_i,sb);
end

[~,affordance_i]    = ismember(affordance,{'greet','grasp','enter','observe'});

% Get info about images
n_obj        = length(core_im_name);
n_wm_im      = length(wm_test_im_name);
n_wm_changes   = length(params.stim.obj.delta_from_ref);

% Predefine info table
info.stim_pos_name   = cat(1, repmat({'left';'right'},length(core_im_name)/2,1) , ....
                              repmat(repelem({'left';'right'},n_wm_changes), length(core_im_name)/2,1));
info.stim_pos        = cat(1, repmat([1;2],length(core_im_name)/2,1) , ....
                              repmat(repelem([1;2],n_wm_changes), length(core_im_name)/2,1));
                          
info.super_cat_name  = cat(1,superordinate, ... core im
                                repelem(superordinate,n_wm_changes));... wm im
info.super_cat       = cat(1,superordinate_i, ... core im
                                repelem(superordinate_i,n_wm_changes));... wm im
info.basic_cat_name  = cat(1,basic, ... core im
                                repelem(basic,n_wm_changes));... wm im
info.basic_cat       = cat(1,basic_i, ... core im
                                repelem(basic_i,n_wm_changes));... wm im
info.sub_cat_name    = cat(1,subordinate, ... core im
                                repelem(subordinate,n_wm_changes));... wm im
info.sub_cat         = cat(1,subordinate_i', ... core im
                                repelem(subordinate_i', n_wm_changes)); % wm im
info.affordance_name = cat(1,affordance,repelem(affordance,n_wm_changes)); 
info.affordance_cat  = cat(1,affordance_i,repelem(affordance_i,n_wm_changes)); 
info.is_specialcore  = ismember(info.unique_im,params.stim.obj.unique_im_nrs_specialcore);  % is_specialcore (logical)

% Define all rotations
rotations     = [0, params.stim.obj.delta_from_ref];
info.base_rot = cat(1,params.stim.obj.facing_dir_deg',repelem(params.stim.obj.facing_dir_deg,n_wm_changes)'); 
rot_abs_wm    = params.stim.obj.facing_dir_deg + params.stim.obj.delta_from_ref';
info.abs_rot  = cat(1,params.stim.obj.facing_dir_deg', rot_abs_wm(:));  % alternating 2 view ± 4 steps spaced 2 deg apart
info.rel_rot  = cat(1, repmat(rotations(1),n_obj,1),repmat(rotations(2:end)',n_obj,1));

% reshape unique image nr for WM test images
unique_wm_im   = reshape(params.stim.obj.unique_im_nrs_wm_test, [],n_obj);

% Get (rescaled) image extent
extent   = params.stim.obj.img_sz_pix.*params.stim.obj.dres;

% Preallocate space
objects = uint8(ones(extent, extent, 3, n_obj, length(rotations)));
masks   = uint8(ones(extent, extent, n_obj, length(rotations)));

counter = 1; % we use the lazy counter way to keep track of image

% Loop over each image
for sub = 1:n_obj
    
    % loop over each rotation
    for rr = 1:length(rotations)
        
        % Get file name
        rot = 1+((info.base_rot(sub) + rotations(rr))/2);
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

        % check facing direction, when the offset may technically result in
        % a flipping of the facing direction, go by original core facing
        % direction, as subjects will never perform the PC task on working
        % memory test images.
        if info.rel_rot(counter) == 0
            if (info.abs_rot(counter) < 45) || (info.abs_rot(counter) > 135)
                info.facing_dir(counter) = {'sideways'};
            elseif (info.abs_rot(counter) > 45) || (info.abs_rot(counter) < 135)
                info.facing_dir(counter) = {'forward'};
            end
        elseif info.rel_rot(counter) ~= 0
            og_rot = info.base_rot(counter);
            if (og_rot < 45) || (og_rot > 135)
                info.facing_dir(counter) = {'sideways'};
            elseif (og_rot > 45) || (og_rot < 135)
                info.facing_dir(counter) = {'forward'};
            end
        end

        % we define rotations relative to the reference object rotation
            % if abs rotation is between 0-90 and rel rotation is negative
        if (info.abs_rot(counter)-90) < 0 && info.rel_rot(counter) < 0 
            info.rel_rot_name(counter) = {'rightward'};
            % if abs rotation is between 0-90 and rel rotation is positive
        elseif (info.abs_rot(counter)-90) < 0 && info.rel_rot(counter) > 0 
            info.rel_rot_name(counter) = {'leftward'};
            % if abs rotation is between 90-180 and rel rotation is negative
        elseif (info.abs_rot(counter)-90) > 0 && info.rel_rot(counter) < 0 
            info.rel_rot_name(counter) = {'rightward'};
            % if abs rotation is between 90-180 and rel rotation is positive
        elseif (info.abs_rot(counter)-90) > 0 && info.rel_rot(counter) > 0 
            info.rel_rot_name(counter) = {'leftward'};
            % CORNER CASE 1: if abs rotation is exactly 90, but core rotation was between 0-90, and rel rotation is negative
        elseif (info.abs_rot(counter)-90) == 0 && info.base_rot(counter) < 0 && info.rel_rot(counter) < 0 
            info.rel_rot_name(counter) = {'rightward'};
            % CORNER CASE 2:  if abs rotation is exactly 90, but core rotation was between 0-90, and rel rotation is positive
        elseif (info.abs_rot(counter)-90) == 0 && info.base_rot(counter) < 0 && info.rel_rot(counter) > 0
            info.rel_rot_name(counter) = {'leftward'};
            % CORNER CASE 3:  if abs rotation is exactly 90, but core rotation was between 90-180, and rel rotation is negative
        elseif (info.abs_rot(counter)-90) == 0 && info.base_rot(counter) > 0 && info.rel_rot(counter) < 0
            info.rel_rot_name(counter) = {'rightward'};
            % CORNER CASE 4:  if abs rotation is exactly 90, but core rotation was between 90-180, and rel rotation is positive
        elseif (info.abs_rot(counter)-90) == 0 && info.base_rot(counter) > 0 && info.rel_rot(counter) > 0
            info.rel_rot_name(counter) = {'leftward'};
        elseif info.rel_rot(counter) == 0
            info.rel_rot_name(counter) = {'none'};
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
    save(fullfile(sprintf('%s_%s.mat',params.stim.obj.stimfile,datestr(now,30))),'objects','masks','info','-v7.3');
    
    saveDir = fileparts(fullfile(params.stim.obj.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',params.stim.obj.infofile,datestr(now,30))))
end

%% Visualize stimuli if requested

if params.stim.store_imgs
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
end

if params.verbose
    
    % debug figure
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
%             if params.stim.store_imgs
%               saveFigDir = fullfile(vcd_rootPath,'figs',params.disp.name, 'obj','visual_checks');
%               if ~exist(saveFigDir,'dir'); mkdir(saveFigDir); end
%                 % Print figure with axes
%                 print(fullfile(saveFigDir,sprintf('%02d_object%02d_delta%02d', im_nr,objectNr, dd)),'-dpng','-r150');
%             end
%             
            % update counter
%             counter = counter+1;
%         end
%     end
end

return

