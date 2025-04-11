function [scenes,lures,cblind, im_order,info,p] = vcd_naturalscenes(p)
% VCD function to load, and if needed resize/square pix values for natural scenes:
%
%  [scenes,lures,cblind,im_order,info,p] = vcd_naturalscenes(p)
%
% Purpose:
%   Load subset of natural scene images used in the Natural Scenes dataset
%   for specified VCD experimental display. If p.store_imgs is set to
%   true then this function will store the loaded images as a single mat
%   file as defined by p.stim.ns.stimfile, e.g.: 
%   fullfile(vcd_rootPath,'workspaces','stimuli','scenes.mat')
%
% Requirements:
%   1. csv file with png file names defined by p.stim.ns.infofile, e.g.: 
%       fullfile(vcd_rootPath,'workspaces','scenes_info.csv')
%   2. Folder with png images, defined by p.stim.ns.indivscenefiles, e.g.: 
%       fullfile(vcd_rootPath,'workspaces','stimuli','vcd_natural_scenes')
%   
%
% INPUTS:
%  p                  : struct with stimulus params, with the following fields
%   store_imgs        : store finalized scene images and info or not.
%   disp.name         : name of used monitor (to check for squaring pixel values) (see vcd_getDisplayParams.m)
%   stim.ns.infofile  : csv file name (see vcd_getStimulusParams.m)
%   stim.ns.lure_im   : header names for lure images in csv file
%   stim.ns.change_im : header names for change blindness images in csv file
%   stim.ns.og_res_stim : original size of scene images (in pixels), can be height or width as we assume image is square
%   stim.ns.dres      : scale factor (fraction) applied to OG scene resolution
%   stim.ns.stimfile  : name of matfile final images will be stored
%
% OUTPUTS:
%  scenes       :   uint8 images of natural scenes, dimensions are
%                   width (pixels) x height (pixels) x 3 (rgb) x 5
%                   superordinate categories x 2 scene locations
%                   (indoor/outdoor) x 3 obj locations (left/middle/right)
%  lures        :   uint8 lure images of 30 VCD natural scenes, dimensions are
%                   width (pixels) x height (pixels) x 3 (rgb) x 5
%                   superordinate categories x 2 scene locations
%                   (indoor/outdoor) x 3 obj locations (left/middle/right)
%                   x 4 lure types (1: most similar, 4 least similar)
%  cblind       :   uint8 change blindness images of VCD natural scenes, 
%                   dimensions are width (pixels) x height (pixels) x 3
%                   (rgb) x 5 superordinate categories x 2 scene locations
%                   (indoor/outdoor) x 3 obj locations (left/middle/right)
%                   x 4 change blindness types (1: easy_add, 2: hard_add,
%                   3: easy_remove, 4: hard_remove).
%  im_order     :   270x8 table with information about scene super and basic
%                   level semantic category, affordance type, scene loc,
%                   object location in scene, lure name, change blindness
%                   name, unique image number.
%  info         :   Loaded csv from workspaces/stimuli/ with png file names
%                   and scene/object category information. Has some overlap 
%                   of information with im_order, such as category info.
%  p            :   updated params struct
%
% Written by Eline Kupers 2024/12
%

%% Load images from stim info file and resize if requested
if ~isfield(p.stim.ns, 'infofile')
    p.stim.ns.infofile = fullfile(vcd_rootPath,'workspaces','info','scene*');
end

% find info file
d = dir(sprintf('%s*.csv',p.stim.ns.infofile));

if isempty(d)
    error('[%s]: Stimfile or infofile cannot be found or loaded',mfilename)
end

% Read info table
info = readtable(fullfile(d(end).folder, d(end).name));

% Define superordinate and basic categories, and number of exemplars per basic category
superordinate = unique(info.superordinate,'stable'); % 4 superordinate categories
basic         = unique(info.basic,'stable');         % basic
ns_loc        = unique(info.ns_loc,'stable');        % indoor/outdoor
obj_loc       = unique(info.obj_loc,'stable');       % dominant object location

% Get info about images
n_images          = length(superordinate)*length(ns_loc)*length(obj_loc);
n_lures           = length(p.stim.ns.lure_im);
n_changeblindness = length(p.stim.ns.change_im);

% Predefine im_order table
n_cols   = n_images + (n_images*n_lures) + (n_images*n_changeblindness);
im_order = table(cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),NaN(n_cols,1));
im_order.Properties.VariableNames = {'super_cat','basic_cat','affordance','scene_loc','obj_loc','lure_im','change_im','unique_im'};

% Preallocate space for images
scenes0 = uint8(ones(p.stim.ns.og_res_stim,p.stim.ns.og_res_stim,3,...
    length(superordinate), length(ns_loc),length(obj_loc)));
cblind0 = uint8(ones(p.stim.ns.og_res_stim,p.stim.ns.og_res_stim,3,...
    length(superordinate), length(ns_loc),length(obj_loc),n_lures));
lures0 = uint8(ones(p.stim.ns.og_res_stim,p.stim.ns.og_res_stim,3,...
    length(superordinate), length(ns_loc),length(obj_loc),n_changeblindness));


%% Loop over superordinate, inside/outside, and varying object location images
counter = 1;

for ss = 1:length(superordinate)
    for ex = 1:length(ns_loc)
        for bb = 1:length(obj_loc)
            
            % Get scene file name from info table
            im_png = info.filename{ ...
                strcmp(info.superordinate,superordinate(ss)) & ...
                strcmp(info.ns_loc,ns_loc(ex)) & ...
                strcmp(info.obj_loc,obj_loc(bb))};
            
            d = dir(fullfile(vcd_rootPath,'workspaces','stimuli','vcd_natural_scenes', im_png));
            
            % Read in image
            scenes0(:,:,:,ss,ex,bb) = imread(fullfile(d.folder,d.name));

            % Keep track of category etc
            im_order.super_cat(counter) = superordinate(ss);
            im_order.basic_cat(counter) = info.basic(bb);
            im_order.affordance(counter)= info.obj_act(bb);
            im_order.scene_loc(counter) = ns_loc(ex);
            im_order.obj_loc(counter)   = obj_loc(bb);
            im_order.unique_im(counter) = (ex-1)*length(obj_loc) + bb + (ss-1)*(length(ns_loc)*length(obj_loc));
            im_order.lure_im(counter)   = {NaN};
            im_order.change_im(counter) = {NaN};
            
            counter = counter+1;
            
            % Find corresponding change blindness images
            for cb_idx = 1:n_changeblindness
                
                % Find image filename in info table
                fn = sprintf('change_img%02d',cb_idx);
                og_scene = strsplit(im_png,'.png'); og_scene = og_scene{1};
                
                cb_im_name = info.(fn){strcmp(info.superordinate,superordinate(ss)) & ...
                                        strcmp(info.ns_loc,ns_loc(ex)) & ...
                                        strcmp(info.obj_loc,obj_loc(bb))};
                
                d = dir(fullfile(vcd_rootPath,'workspaces','stimuli','vcd_natural_scenes',...
                    'changes', og_scene, cb_im_name));
                if isempty(d)
                    error('[%s]: can''t find change blindness image file',mfilename)
                end                
                % Load image into array
                cblind0(:,:,:,ss,ex,bb,cb_idx) = imread(fullfile(d.folder,d.name));

                % Keep track of category, etc.
                im_order.super_cat(counter) = superordinate(ss);
                im_order.basic_cat(counter) = info.basic(bb);
                im_order.affordance(counter)= info.obj_act(bb);
                im_order.scene_loc(counter) = ns_loc(ex);
                im_order.obj_loc(counter)   = obj_loc(bb);
                im_order.unique_im(counter) = (ex-1)*length(obj_loc) + bb + (ss-1)*(length(ns_loc)*length(obj_loc));
                im_order.lure_im(counter)   = {NaN};
                im_order.change_im(counter) = {fn};
                
                counter = counter+1;
                
            end
            
            
            % Now do the same for LTM lures (can't use the same look
            % because nr of lures and change blindness varies per image
            for lure_idx = 1:n_lures

                % Find image filename in info table
                fn = sprintf('lure_img%02d',lure_idx);
                lure_im_name = info.(fn){ ...
                    strcmp(info.superordinate,superordinate(ss)) & ...
                    strcmp(info.ns_loc,ns_loc(ex)) & ...
                    strcmp(info.obj_loc,obj_loc(bb))};
                
                d = dir(fullfile(vcd_rootPath,'workspaces','stimuli','vcd_natural_scenes',...
                    'lures', og_scene, lure_im_name));
                if isempty(d)
                    error('[%s]: can''t find lure image file',mfilename)
                end
                    
                % Load image into array
                lures0(:,:,:,ss,ex,bb,lure_idx) = imread(fullfile(d.folder,d.name));
                
                % Keep track of category, etc.
                im_order.super_cat(counter) = superordinate(ss);
                im_order.basic_cat(counter) = info.basic(bb);
                im_order.affordance(counter)= info.obj_act(bb);
                im_order.scene_loc(counter) = ns_loc(ex);
                im_order.obj_loc(counter)   = obj_loc(bb);
                im_order.unique_im(counter) = (ex-1)*length(obj_loc) + bb +  (ss-1)*(length(ns_loc)*length(obj_loc));
                im_order.lure_im(counter)   = {fn};
                im_order.change_im(counter) = {NaN};
                
                counter = counter+1;
            end
            
        end
    end
end

%% Resize image if needed
if ~isempty(p.stim.ns.dres) && ~isequal(p.stim.ns.dres,1)
    fprintf('[%s]: Resampling the stimuli; this may take a while',mfilename);
    tic;

    % collapse og scenes
    scenes_rz = reshape(scenes0,size(scenes0,1),size(scenes0,2),size(scenes0,3), ...
        size(scenes0,4)*size(scenes0,5)*size(scenes0,6));
    
    % collapse lure conditions
    lures_rz = reshape(lures0,size(lures0,1),size(lures0,2),size(lures0,3), ...
        size(lures0,4)*size(lures0,5)*size(lures0,6),size(lures0,7));
    
    % collapse change blindness conditions    
    cblind_rz = reshape(cblind0,size(cblind0,1),size(cblind0,2),size(cblind0,3), ...
        size(cblind0,4)*size(cblind0,5)*size(cblind0,6),size(cblind0,7));
    
    % create temp var
    temp_s = cast([],class(scenes_rz));
    temp_l = cast([],class(lures_rz));
    temp_c = cast([],class(cblind_rz));
    
    % loop over scenes, lures, clindness images
    for pp = 1:size(scenes_rz,4)
        statusdots(pp,size(scenes_rz,4));
        
        % before resizing convert to double
        sc0  = double(scenes_rz(:,:,:,pp)); % image lum vals are now double
        sc1  = uint8(imresize(sc0,p.stim.ns.dres)); % resize, convert to uint8     
        
        % check for clipping
        assert( min(sc1(:))>=0 && max(sc1(:))<=255)
        
        temp_s(:,:,:,pp) = sc1;
        
        % Loop over lures
        for ll = 1:n_lures
            l0   = double(lures_rz(:,:,:,pp,ll)); % image lum vals ranges from [0-255]
            l1   = uint8(imresize(l0,p.stim.ns.dres));

            % check for clipping
            assert( min(l1(:))>=0 && max(l1(:))<=255)
            
            temp_l(:,:,:,pp,ll) = l1;
        end
        
        % Loop over cblindness images
        for kk = 1:n_changeblindness
            c0 = double(cblind_rz(:,:,:,pp,kk)); % image lum vals ranges from [0-255]
            c1 = uint8(imresize(c0,p.stim.ns.dres));
            
            % check for clipping
            assert( min(c1(:))>=0 && max(c1(:))<=255)
            
            temp_c(:,:,:,pp,kk) = c1;
        end
    end
    
    % accumulate images
    scenes_rz = temp_s;
    lures_rz  = temp_l;
    cblind_rz = temp_c;
    
    % Reshape back
    scenes = reshape(scenes_rz, size(scenes_rz,1),size(scenes_rz,2),size(scenes_rz,3),...
        size(scenes0,4),size(scenes0,5),size(scenes0,6));
    
    lures = reshape(lures_rz, size(lures_rz,1),size(lures_rz,2),size(lures_rz,3),...
        size(lures0,4),size(lures0,5),size(lures0,6),size(lures0,7));
    
    cblind = reshape(cblind_rz, size(cblind_rz,1),size(cblind_rz,2),size(cblind_rz,3),...
        size(cblind0,4),size(cblind0,5),size(cblind0,6),size(cblind0,7));
    
    fprintf('[%]: done!',mfilename); toc
    
    % Clean up
    clear scenes0 lures0 cblind0 scenes_rz lures_rz cblind_rz temp_s tmp_l tmp_c
    
else % do nothing
    scenes = scenes0;
    lures  = lures0;
    cblind = cblind0; 
end


%% Visualize resized images
makeprettyfigures;
figure; set(gcf,'Position', [156    91   881   706],'color','w');
scenes0 = reshape(scenes,size(scenes,1),size(scenes,2),size(scenes,3),size(scenes,4)*size(scenes,5)*size(scenes,6));
lures0  = reshape(lures,size(lures,1),size(lures,2),size(lures,3),size(lures,4)*size(lures,5)*size(lures,6),size(lures,7));
cblind0 = reshape(cblind,size(cblind,1),size(cblind,2),size(cblind,3),size(cblind,4)*size(cblind,5)*size(cblind,6),size(cblind,7));

for ss = 1:size(scenes0,4)
    
    clf;
    imagesc(scenes0(:,:,:,ss));
    title(sprintf('Im %02d resized',ss), 'FontSize',20);
    axis image; box off
    set(gca,'CLim',[1 255]);
    if p.stim.store_imgs
        saveDir = fullfile(vcd_rootPath,'figs',p.disp.name,'ns','resized');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        print(fullfile(saveDir, sprintf('ns_%02d', ss)),'-dpng','-r150');
    end
    
    for ll = 1:size(lures0,5)
        clf;
        imagesc(lures0(:,:,:,ss,ll));
        title(sprintf('Im %02d, lure %02d resized ',ss,ll), 'FontSize',20);
        axis image; box off
        set(gca,'CLim',[1 255]);
        if p.stim.store_imgs
            print(fullfile(saveDir, sprintf('ns_%02d_lure%02d', ss,ll)),'-dpng','-r150');
        end
        
        clf;
        imagesc(cblind0(:,:,:,ss,ll)); 
        title(sprintf('Im %02d, cblindness %02d resized',ss,ll), 'FontSize',20);
        axis image; box off
        set(gca,'CLim',[1 255]);
        if p.stim.store_imgs
            print(fullfile(saveDir, sprintf('ns_%02d_cblind%02d', ss,ll)),'-dpng','-r150');
        end
    end
end

%% For monitors acting linearly, square pixel values 
if any(strcmp(p.disp.name, {'7TAS_BOLDSCREEN32', 'PPROOM_EIZOFLEXSCAN'}))
    fprintf('[%]: Square pixel values given using %s with linear gamma CLUT',mfilename, p.disp.name);
    
    % collapse og scenes category dim (w x h x 3 x 30)
    scenes_sq = reshape(scenes,size(scenes,1),size(scenes,2),size(scenes,3), ...
        size(scenes,4)*size(scenes,5)*size(scenes,6));
    
    % collapse lures category dim (w x h x 3 x 30 x 4)
    lures_sq = reshape(lures,size(lures,1),size(lures,2),size(lures,3), ...
        size(lures,4)*size(lures,5)*size(lures,6),size(lures,7));
    
    % collapse change blindness category dim (w x h x 3 x 30 x 4)
    cblind_sq = reshape(cblind,size(cblind,1),size(cblind,2),size(cblind,3), ...
        size(cblind,4)*size(cblind,5)*size(cblind,6),size(cblind,7));
    
    % create temp var
    temp_s = cast([],class(scenes_sq));
    temp_l = cast([],class(lures_sq));
    temp_c = cast([],class(cblind_sq));
    
    for ss = 1:size(scenes_sq,4)
        
        temp_s(:,:,:,ss) = uint8(255.*((double(scenes_sq(:,:,:,ss))./255).^2));
        
        for ll = 1:size(lures_sq,5)
            temp_l(:,:,:,ss,ll) = uint8(255.*((double(lures_sq(:,:,:,ss,ll))./255).^2));
            temp_c(:,:,:,ss,ll) = uint8(255.*((double(cblind_sq(:,:,:,ss,ll))./255).^2));
        end
    end
    
    figure; set(gcf,'Position', [156    91   881   706],'color','w');
    for ss = 1:size(scenes,4)
        clf;
        imagesc(temp_s(:,:,:,ss));
        title(sprintf('Im %02d resized & squared',ss), 'FontSize',20);
        axis image; box off
        set(gca,'CLim',[1 255]);
        if p.stim.store_imgs
            saveDir = fullfile(vcd_rootPath,'figs',p.disp.name,'ns','resized_and_squared');
            if ~exist(saveDir,'dir'), mkdir(saveDir); end
            print(fullfile(saveDir, sprintf('ns_%02d', ss)),'-dpng','-r150');
        end
        
        for ll = 1:size(lures,5)
            clf;
            imagesc(temp_l(:,:,:,ss,ll));
            title(sprintf('Im %02d, lure %02d resized & squared',ss,ll), 'FontSize',20);
            axis image; box off
            set(gca,'CLim',[1 255]);
            if p.stim.store_imgs
                print(fullfile(saveDir, sprintf('ns_%02d_lure%02d', ss,ll)),'-dpng','-r150');
            end
            
            clf;
            imagesc(temp_c(:,:,:,ss,ll));
            title(sprintf('Im %02d, cblindness %02d resized & squared',ss,ll), 'FontSize',20);
            axis image; box off
            set(gca,'CLim',[1 255]);
            
            if p.stim.store_imgs
                print(fullfile(saveDir, sprintf('ns_%02d_cblind%02d', ss,ll)),'-dpng','-r150');
            end
        end
    end

    %% Reshape back 
    % Scene dims: (w x h x 3 x 5 super cat x 2 scene loc x 3 obj loc)
    scenes = reshape(temp_s, size(scenes,1),size(scenes,2),size(scenes,3),...
                             size(scenes,4),size(scenes,5),size(scenes,6));
    % Lure dims: (w x h x 3 x 5 super cat x 2 scene loc x 3 obj loc x 4 lure types)
    lures  = reshape(temp_l, size(lures,1),size(lures,2),size(lures,3),...
                            size(lures,4),size(lures,5),size(lures,6),size(lures,7));
    % Change blindness dims: (w x h x 3 x 5 super cat x 2 scene loc x 3 obj loc x 4 blindness types)
    cblind = reshape(temp_c, size(cblind,1),size(cblind,2),size(cblind,3),...
                             size(cblind,4),size(cblind,5),size(cblind,6),size(cblind,7));

    % Clean up
    clear temp_s temp_l temp_c
end

%% Store images if requested
if p.stim.store_imgs
    fprintf('[%]: Storing images',mfilename);
    saveDir = fileparts(fullfile(p.stim.ns.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s_%s.mat',p.stim.ns.stimfile,p.disp.name,datestr(now,30))),'scenes','lures','cblind','info','im_order','-v7.3');
end


return

