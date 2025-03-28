function [objects, masks, im_order, info, p] = vcd_complexobjects(p)
% VCD function to generate complex object images:
% 
%   [objects, masks, im_order, info, p] = vcd_complexobjects(p)
%
% Purpose:
%   Load complex object images for VCD experimental display.
%   Requires two files:
%    1. a set of individual object image files, where the folder is defined
%       by p.stim.cobj.indivobjfile, e.g.:
%       fullfile(vcd_rootPath,'workspaces/stimuli/vcd_complex_objects/*_rot*.png)
%    2. A csv file defined in p.stim.cobj.infofile, e.g.: 
%       fullfile(vcd_rootPath,'workspaces/objects_info.csv')
%
%   This function will store generated object images as a single mat file 
%   in p.stim.cobj.stimfile (e.g.: fullfile(vcd_rootPath,'workspaces','stimuli','objects.mat')
%   when p.store_params = true;
%
%   We reorder and select unique object images from all_unique_im.cobj:
%   unique_im loc  super  basic  sub  facing_dir
%      1     1     1     1     1    56
%      2     2     1     2     2   124
%      3     1     1     2     3    56
%      4     2     2     1     1   124
%      5     1     2     1     2    56
%      6     2     2     2     3   124
%      7     1     2     2     4    56
%      8     2     3     1     1   124
%      9     1     3     1     2    56
%     10     2     3     2     3   124
%     11     1     3     2     4    56
%     12     2     3     3     5   124
%     13     1     3     3     6    56
%     14     2     4     1     1   124
%     15     1     4     1     2    56
%     16     2     4     1     3   124
%
% INPUTS:
%  p            :   struct with stimulus params
%
% OUTPUTS:
%  objects      :   uint8 images of complex images, dimensions are
%                   width x height x 3 (rgb) x object's subordinate
%                   category x num of rotations
%  masks        :   uint8 alpha transparency masks for ptb to remove the 
%                   background. Dimensions are width (pixels) x height
%                   (pixels) x object's subordinate category x num of
%                   rotations
%  im_order     :   table with information about object rotation and unique
%                   image number
%  info         :   Loaded csv from workspaces/stimuli/ with png file names
%                   and object category information. 
%  p            :   updated params struct
%
% Written by Eline Kupers 2024/12

%% load images
d = dir(sprintf('%s*',fullfile(p.stim.cobj.infofile)));
info = readtable(fullfile(d(end).folder,d(end).name));

% Define superordinate and basic categories, and number of exemplars per basic category
mask            = info.unique_im>0;
subordinate     = unique(info.subordinate,'stable');
canonical_view  = info.rot_abs(mask);                % alternating 2 canonical view ± 5 2-deg steps

rotation_names  = [info.Properties.VariableNames(~cellfun(@isempty, regexp(info.Properties.VariableNames, 'rot_minus*'))), ...
    info.Properties.VariableNames(~cellfun(@isempty, regexp(info.Properties.VariableNames, 'rot_plus*')))];

rotations = 0;
for ii = 1:length(rotation_names)
    rotations = [rotations, choose(regexp(rotation_names{ii}, 'rot_minus'),-1,1).*str2num(cell2mat(regexp(rotation_names{ii}, '\d+','match')))];
end

%% EK HACK --> process more images to delete this line
rotations = rotations([1, 3:end-1]);
rotation_names2 = catcell(2,{{'none'}, rotation_names(2:end-1)});
%% EK HACK END

% Preallocate info table
n_cols   = length(canonical_view)*length(rotations);
im_order = table(cell(n_cols,1),NaN(n_cols,1),NaN(n_cols,1),cell(n_cols,1),NaN(n_cols,1));
im_order.Properties.VariableNames = {'object_name','abs_rot','rel_rot','rot_name','unique_im'};

% Get (rescaled) image extent
extent   = p.stim.cobj.og_res_stim.*p.stim.cobj.dres;

% Preallocate space
objects = uint8(ones(extent, extent,3,...
    length(subordinate),length(rotations)));
masks = uint8(ones(extent, extent, ...
    length(subordinate),length(rotations)));

counter = 1; % we use the lazy counter way to keep track of image

for sub = 1:length(subordinate)
    
    for rr = 1:length(rotations)
        
        % Get file name
        rot = (canonical_view(sub)/2) + mod(rotations(rr),2);
        d = dir(fullfile(p.stim.cobj.indivobjfile,sprintf('*%s*_rot%02d.png', subordinate{sub},rot)));
        
        % Load image
        [im0, ~, alpha0] = imread(fullfile(d.folder,d.name),'BackgroundColor','none');
        
        % convert to double, [0-1] lum range (assume range is 1-255),
        % and square as images were made with gamma = 2.
        if ~isequal(p.stim.cobj.dres,1) && ~isempty(p.stim.cobj.dres)
            fprintf('[%s]: Resampling the stimuli..',mfilename);
            
            im0     = (double(im0-1)./254).^2; % convert to [0-1] and square
            alpha0  = (double(alpha0./255)).^2; % convert to [0-1] and square
            
            im1     = imresize(im0, p.stim.cobj.dres); % scale if needed
            alpha1  = imresize(alpha0,p.stim.cobj.dres); % scale if needed
            
            im = uint8(1+(sqrt(im1)*254));
            alpha_im = uint8((sqrt(alpha1)*255));
        else
            im = im0;
            alpha_im = alpha0;
        end
        clear im0 alpha0
        
        %             figure(99); clf; imagesc(im,'AlphaData',alpha_im); colormap gray; colorbar; axis image
        
        objects(:,:,:,sub,rr) = repmat(im, [1 1 3]);
        masks(:,:,sub,rr) = alpha_im;
        
        % Keep track of image order
        im_order.object_name(counter)    = subordinate(sub);
        im_order.canonical_view(counter) = canonical_view(sub);
        im_order.abs_rot(counter)        = canonical_view(sub)+rotations(rr);
        im_order.rel_rot(counter)        = rotations(rr);
        im_order.rot_name(counter)       = rotation_names2(rr);
        im_order.unique_im(counter)      = sub;
        
        % Update counter
        counter = counter+1;
    end
end

if p.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(p.stim.cobj.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',p.stim.cobj.stimfile,datestr(now,30))),'objects','masks','info','im_order','-v7.3');
end

%
%     figure(1); clf;
%     figure(2); clf;
%     figure(3); clf;
%     figure(4); clf;
%     for ii = 1:size(objects,5)
%         figure(1); subplot(3,3,ii); hold all;
%         imshow(objects(:,:,1,1,ii),[1 255]);
%
%         figure(2); subplot(3,3,ii); hold all;
%         imshow(objects(:,:,1,2,ii),[1 255]);
%
%         figure(3); subplot(3,3,ii); hold all;
%         imshow(objects(:,:,1,3,ii),[1 255]);
%
%         figure(4); subplot(3,3,ii); hold all;
%         imshow(objects(:,:,1,4,ii),[1 255]);
%     end
return

