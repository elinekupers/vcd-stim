function [scenes,ltm_lures,wm_im, im_order,info] = vcd_naturalscenes(params)
% VCD function to load and resize/square pix values for natural scenes:
%
%  [scenes,lures,cblind,im_order,info,params] = vcd_naturalscenes(params)
%
% Purpose:
%   Load subset of natural scene images used in the Natural Scenes Dataset
%   for specified VCD experimental display. If params.store_imgs is set to
%   true then this function will store the loaded images as a single mat
%   file as defined by p.stim.ns.stimfile, e.g.: 
%   fullfile(vcd_rootPath,'workspaces','stimuli','scenes.mat')
%
% Requirements:
%   1. csv file with png file names defined by p.stim.ns.infofile, e.g.: 
%       fullfile(vcd_rootPath,'workspaces','scenes_info.csv')
%   2. Folder with png images, defined by p.stim.ns.indivscenefiles, e.g.: 
%       fullfile(vcd_rootPath,'workspaces','stimuli','RAW','vcd_natural_scenes')
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
%  ltm_lures    :   uint8 4 LTM lure images for the 15 chosen VCD natural 
%                   scenes, dimensions are:
%                       height (pixels) x width (pixels) x 3 (rgb) x 5
%                       superordinate categories x 2 scene locations
%                       (indoor/outdoor) x 3 obj locations (left/center/right)
%                       x 4 lure types (1: most similar, 4 least similar)
%  wm_im        :   uint8 wm test images for each of the 30 VCD natural 
%                   scenes (think of a change blindness experiment)
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
%
% Written by Eline Kupers 2024/12, updated 2025/04
%

%% Load images from stim info file and resize if requested
if ~isfield(params.stim.ns, 'infofile')
    params.stim.ns.infofile = fullfile(vcd_rootPath,'workspaces','info','scene*');
end

% find info file
d = dir(sprintf('%s*.csv',params.stim.ns.infofile));
if isempty(d)
    error('[%s]: Infofile cannot be found!',mfilename)
end

% Read info table, we import as char, other the ns_loc and obj_loc columns
% are converted to NaNs. One annoying part is that we have to convert image
% numbers back to integers (and using a loop because matlab won't zero pad)
opts = detectImportOptions(fullfile(d(end).folder, d(end).name));
opts = setvartype(opts,'char');  % or 'string'
info = readtable(fullfile(d(end).folder, d(end).name),opts);

for ii = 1:length(info.unique_im_nr)
    unique_im_nrs(ii) = str2num(info.unique_im_nr{ii});
end
[~,idx0] = intersect(params.stim.ns.unique_im_nrs,unique_im_nrs);

% Define superordinate and basic categories, and number of exemplars per basic category
superordinate = unique(info.superordinate,'stable'); % 4 superordinate categories
basic         = unique(info.basic,'stable');         % basic
ns_loc        = unique(info.ns_loc,'stable');        % indoor/outdoor
obj_loc       = unique(info.obj_loc,'stable');       % dominant object location
ns_loc(find(~cell2mat(cellfun(@isempty, regexp(ns_loc,'NaN'), 'UniformOutput', false))))=[];
obj_loc(find(~cell2mat(cellfun(@isempty, regexp(obj_loc,'NaN'), 'UniformOutput', false))))=[];

% WM image filenames
wm_test_im_name = repmat(params.stim.ns.change_im,1,length(idx0));

% Get info about images
n_scenes          = length(idx0)*length(ns_loc)*length(obj_loc);
n_ltm_lures       = length(params.stim.ns.lure_im);
n_wm_im = length(params.stim.ns.change_im);

% Predefine im_order table
n_cols   = n_scenes + (n_scenes*n_ltm_lures) + (n_scenes*n_wm_im);
im_order = table(NaN(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),NaN(n_cols,1),NaN(n_cols,1));
im_order.Properties.VariableNames = {'unique_im','super_cat','basic_cat','affordance','scene_loc','obj_loc','change_im','lure_im'};

% Preallocate space for images
scenes0 = uint8(ones(params.stim.ns.og_res_stim,params.stim.ns.og_res_stim,3,...
    length(superordinate), length(ns_loc),length(obj_loc)));
wmtest0 = uint8(ones(params.stim.ns.og_res_stim,params.stim.ns.og_res_stim,3,...
    length(superordinate), length(ns_loc),length(obj_loc),n_ltm_lures));
ltmlures0 = uint8(ones(params.stim.ns.og_res_stim,params.stim.ns.og_res_stim,3,...
    length(superordinate), length(ns_loc),length(obj_loc),n_wm_im));


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
            
            d = dir(fullfile(params.stim.ns.indivscenefile, im_png));
            
            % Read in image
            scenes0(:,:,:,ss,ex,bb) = imread(fullfile(d.folder,d.name));

            % Keep track of category etc
            im_order.super_cat(counter) = superordinate(ss);
            im_order.basic_cat(counter) = info.basic(bb);
            im_order.affordance(counter)= info.obj_act(bb);
            im_order.scene_loc(counter) = ns_loc(ex);
            im_order.obj_loc(counter)   = obj_loc(bb);
            im_order.unique_im(counter) = unique_im_nrs((ex-1)*length(obj_loc) + bb + (ss-1)*(length(ns_loc)*length(obj_loc)));
            im_order.lure_im(counter)   = NaN;
            im_order.change_im(counter) = NaN;
            
            row_nr = find(strcmp(info.filename,im_png));
            
            counter = counter+1;
            
            % Find corresponding change blindness images
            for cb_idx = 1:n_wm_im
                
                % Find image filename in info table
                og_scene = strsplit(im_png,'.png'); og_scene = og_scene{1};
                fn = sprintf('%s/%s.png',og_scene,wm_test_im_name{cb_idx});

                cb_im_name = info.filename{find(strcmp(info.filename,fn))};
                
                d = dir(fullfile(params.stim.ns.indivscenefile,...
                    'changes', cb_im_name));
                if isempty(d)
                    error('[%s]: can''t find change blindness image file',mfilename)
                end                
                % Load image into array
                wmtest0(:,:,:,ss,ex,bb,cb_idx) = imread(fullfile(d.folder,d.name));

                % Keep track of category, etc.
                im_order.super_cat(counter) = superordinate(ss);
                im_order.basic_cat(counter) = info.basic(bb);
                im_order.affordance(counter)= info.obj_act(bb);
                im_order.scene_loc(counter) = ns_loc(ex);
                im_order.obj_loc(counter)   = obj_loc(bb);
                im_order.unique_im(counter) = params.stim.ns.unique_im_nrs_WM(cb_idx+((ex-1)*length(obj_loc) + bb + (ss-1)*(length(ns_loc)*length(obj_loc))));
                im_order.lure_im(counter)   = NaN;
                im_order.change_im(counter) = cb_idx;
                
                counter = counter+1;
                
            end
            
            
            % Now do the same for LTM lures (can't use the same look up,
            % because nr of lures and change blindness varies per image)
            if str2num(info.is_in_img_ltm{row_nr})
            
                for lure_idx = 1:n_ltm_lures
                    
                    % Find image filename in info table
                    fn = find(~cell2mat(cellfun(@isempty, regexp(info.filename,sprintf('%s',og_scene),'ONCE'), 'UniformOutput', false)) & ...
                            ~cell2mat(cellfun(@isempty, regexp(info.filename,sprintf('lure%02d.png',lure_idx),'ONCE'), 'UniformOutput', false)));
                    lure_im_name = info.filename{fn};
                    
                    d = dir(fullfile(params.stim.ns.indivscenefile,...
                        'lures', lure_im_name));
                    if isempty(d)
                        error('[%s]: can''t find lure image file',mfilename)
                    end
                    
                    % Load image into array
                    ltmlures0(:,:,:,ss,ex,bb,lure_idx) = imread(fullfile(d.folder,d.name));
                    
                    % Keep track of category, etc.
                    im_order.super_cat(counter) = superordinate(ss);
                    im_order.basic_cat(counter) = info.basic(bb);
                    im_order.affordance(counter)= info.obj_act(bb);
                    im_order.scene_loc(counter) = ns_loc(ex);
                    im_order.obj_loc(counter)   = obj_loc(bb);
                    im_order.unique_im(counter) = params.stim.ns.unique_im_nrs_LTM_lures((ex-1)*length(obj_loc) + bb +  (ss-1)*(length(ns_loc)*length(obj_loc)));
                    im_order.lure_im(counter)   = lure_idx;
                    im_order.change_im(counter) = NaN;
                    
                    counter = counter+1;
                end
            end 
        end
    end
end

%% Resize image if needed
if ~isempty(params.stim.ns.dres) && ~isequal(params.stim.ns.dres,1)
    fprintf('[%s]: Resampling the stimuli; this may take a while',mfilename);
    tic;

    % collapse og scenes
    scenes_rz = reshape(scenes0,size(scenes0,1),size(scenes0,2),size(scenes0,3), ...
        size(scenes0,4)*size(scenes0,5)*size(scenes0,6));
    
    % collapse lure conditions
    lures_rz = reshape(ltmlures0,size(ltmlures0,1),size(ltmlures0,2),size(ltmlures0,3), ...
        size(ltmlures0,4)*size(ltmlures0,5)*size(ltmlures0,6),size(ltmlures0,7));
    
    % collapse change blindness conditions    
    wm_rz = reshape(wmtest0,size(wmtest0,1),size(wmtest0,2),size(wmtest0,3), ...
        size(wmtest0,4)*size(wmtest0,5)*size(wmtest0,6),size(wmtest0,7));
    
    % create temp var
    temp_sc = cast([],class(scenes_rz));
    temp_ltm = cast([],class(lures_rz));
    temp_wm = cast([],class(wm_rz));
    
    % loop over scenes, lures, clindness images
    for pp = 1:size(scenes_rz,4)
        statusdots(pp,size(scenes_rz,4));
        
        % before resizing convert to double
        sc0  = double(scenes_rz(:,:,:,pp)); % image lum vals are now double
        sc1  = uint8(imresize(sc0,params.stim.ns.dres)); % resize, convert to uint8     
        
        % check for clipping
        assert( min(sc1(:))>=0 && max(sc1(:))<=255)
        
        temp_sc(:,:,:,pp) = sc1;
        
        % Loop over lures
        for ll = 1:n_ltm_lures
            l0   = double(lures_rz(:,:,:,pp,ll)); % image lum vals ranges from [0-255]
            l1   = uint8(imresize(l0,params.stim.ns.dres));

            % check for clipping
            assert( min(l1(:))>=0 && max(l1(:))<=255)
            
            temp_ltm(:,:,:,pp,ll) = l1;
        end
        
        % Loop over cblindness images
        for kk = 1:n_wm_im
            c0 = double(wm_rz(:,:,:,pp,kk)); % image lum vals ranges from [0-255]
            c1 = uint8(imresize(c0,params.stim.ns.dres));
            
            % check for clipping
            assert( min(c1(:))>=0 && max(c1(:))<=255)
            
            temp_wm(:,:,:,pp,kk) = c1;
        end
    end
    
    % accumulate images
    scenes_rz = temp_sc;
    lures_rz  = temp_ltm;
    wm_rz = temp_wm;
    
    % Reshape back
    scenes = reshape(scenes_rz, size(scenes_rz,1),size(scenes_rz,2),size(scenes_rz,3),...
        size(scenes0,4),size(scenes0,5),size(scenes0,6));
    
    ltm_lures = reshape(lures_rz, size(lures_rz,1),size(lures_rz,2),size(lures_rz,3),...
        size(ltmlures0,4),size(ltmlures0,5),size(ltmlures0,6),size(ltmlures0,7));
    
    wm_im = reshape(wm_rz, size(wm_rz,1),size(wm_rz,2),size(wm_rz,3),...
        size(wmtest0,4),size(wmtest0,5),size(wmtest0,6),size(wmtest0,7));
    
    fprintf('[%]: done!',mfilename); toc
    
    % Clean up
    clear scenes0 ltmlures0 wmtest0 scenes_rz lures_rz wm_rz temp_sc tmp_ltm tmp_wm
    
else % do nothing
    scenes     = scenes0;
    ltm_lures  = ltmlures0;
    wm_im      = wmtest0; 
end

if params.verbose
    %% Visualize resized images
    makeprettyfigures;
    figure; set(gcf,'Position', [156    91   881   706],'color','w');
    scenes0 = reshape(scenes,size(scenes,1),size(scenes,2),size(scenes,3),size(scenes,4)*size(scenes,5)*size(scenes,6));
    ltmlures0  = reshape(ltm_lures,size(ltm_lures,1),size(ltm_lures,2),size(ltm_lures,3),size(ltm_lures,4)*size(ltm_lures,5)*size(ltm_lures,6),size(ltm_lures,7));
    wmtest0 = reshape(wm_im,size(wm_im,1),size(wm_im,2),size(wm_im,3),size(wm_im,4)*size(wm_im,5)*size(wm_im,6),size(wm_im,7));
    
    for ss = 1:size(scenes0,4)
        
        clf;
        imagesc(scenes0(:,:,:,ss));
        title(sprintf('Im %02d resized',ss), 'FontSize',20);
        axis image; box off
        set(gca,'CLim',[1 255]);
        if params.stim.store_imgs
            saveDir = fullfile(vcd_rootPath,'figs',params.disp.name,'ns','resized');
            if ~exist(saveDir,'dir'), mkdir(saveDir); end
            print(fullfile(saveDir, sprintf('ns_%02d', ss)),'-dpng','-r150');
        end
        
        for ll = 1:size(ltmlures0,5)
            clf;
            imagesc(ltmlures0(:,:,:,ss,ll));
            title(sprintf('Im %02d, LTM lure im %02d resized ',ss,ll), 'FontSize',20);
            axis image; box off
            set(gca,'CLim',[1 255]);
            if params.stim.store_imgs
                print(fullfile(saveDir, sprintf('ns_%02d_ltm_lure%02d', ss,ll)),'-dpng','-r150');
            end
            
            clf;
            imagesc(wmtest0(:,:,:,ss,ll));
            title(sprintf('Im %02d, WM test im %02d resized',ss,ll), 'FontSize',20);
            axis image; box off
            set(gca,'CLim',[1 255]);
            if params.stim.store_imgs
                print(fullfile(saveDir, sprintf('ns_%02d_wm_test%02d', ss,ll)),'-dpng','-r150');
            end
        end
    end
end

%% For monitors acting linearly, square pixel values 
if any(strcmp(params.disp.name, {'7TAS_BOLDSCREEN32', 'PPROOM_EIZOFLEXSCAN'}))
    fprintf('[%]: Square pixel values given using %s with linear gamma CLUT',mfilename, params.disp.name);
    
    % collapse og scenes category dim (w x h x 3 x 30)
    scenes_sq = reshape(scenes,size(scenes,1),size(scenes,2),size(scenes,3), ...
        size(scenes,4)*size(scenes,5)*size(scenes,6));
    
    % collapse lures category dim (w x h x 3 x 30 x 4)
    lures_sq = reshape(ltm_lures,size(ltm_lures,1),size(ltm_lures,2),size(ltm_lures,3), ...
        size(ltm_lures,4)*size(ltm_lures,5)*size(ltm_lures,6),size(ltm_lures,7));
    
    % collapse change blindness category dim (w x h x 3 x 30 x 4)
    wm_sq = reshape(wm_im,size(wm_im,1),size(wm_im,2),size(wm_im,3), ...
        size(wm_im,4)*size(wm_im,5)*size(wm_im,6),size(wm_im,7));
    
    % create temp var
    temp_sc = cast([],class(scenes_sq));
    temp_ltm = cast([],class(lures_sq));
    temp_wm = cast([],class(wm_sq));
    
    for ss = 1:size(scenes_sq,4)
        
        temp_sc(:,:,:,ss) = uint8(255.*((double(scenes_sq(:,:,:,ss))./255).^2));
        
        for ll = 1:size(lures_sq,5)
            temp_ltm(:,:,:,ss,ll) = uint8(255.*((double(lures_sq(:,:,:,ss,ll))./255).^2));
            temp_wm(:,:,:,ss,ll) = uint8(255.*((double(wm_sq(:,:,:,ss,ll))./255).^2));
        end
    end
    
    if params.verbose
        figure; set(gcf,'Position', [156    91   881   706],'color','w');
        for ss = 1:size(scenes_sq,4)
            clf;
            imagesc(temp_sc(:,:,:,ss));
            title(sprintf('Im %02d resized & squared',ss), 'FontSize',20);
            axis image; box off
            set(gca,'CLim',[1 255]);
            if params.stim.store_imgs
                saveDir = fullfile(vcd_rootPath,'figs',params.disp.name,'ns','resized_and_squared');
                if ~exist(saveDir,'dir'), mkdir(saveDir); end
                print(fullfile(saveDir, sprintf('ns_%02d', ss)),'-dpng','-r150');
            end
            
            for ll = 1:size(lures_sq,5)
                clf;
                imagesc(temp_ltm(:,:,:,ss,ll));
                title(sprintf('Im %02d, ltm lure %02d resized & squared',ss,ll), 'FontSize',20);
                axis image; box off
                set(gca,'CLim',[1 255]);
                if params.stim.store_imgs
                    print(fullfile(saveDir, sprintf('ns_%02d_ltm_lure%02d', ss,ll)),'-dpng','-r150');
                end
                
                clf;
                imagesc(temp_wm(:,:,:,ss,ll));
                title(sprintf('Im %02d, wm test im %02d resized & squared',ss,ll), 'FontSize',20);
                axis image; box off
                set(gca,'CLim',[1 255]);
                
                if params.stim.store_imgs
                    print(fullfile(saveDir, sprintf('ns_%02d_wm_im%02d', ss,ll)),'-dpng','-r150');
                end
            end
        end
    end
    
    %% Reshape back 
    % Scene dims: (w x h x 3 x 5 super cat x 2 scene loc x 3 obj loc)
    scenes = reshape(temp_sc, size(scenes,1),size(scenes,2),size(scenes,3),...
                             size(scenes,4),size(scenes,5),size(scenes,6));
    % Lure dims: (w x h x 3 x 5 super cat x 2 scene loc x 3 obj loc x 4 lure types)
    ltm_lures  = reshape(temp_ltm, size(ltm_lures,1),size(ltm_lures,2),size(ltm_lures,3),...
                            size(ltm_lures,4),size(ltm_lures,5),size(ltm_lures,6),size(ltm_lures,7));
    % Change blindness dims: (w x h x 3 x 5 super cat x 2 scene loc x 3 obj loc x 4 blindness types)
    wm_im = reshape(temp_wm, size(wm_im,1),size(wm_im,2),size(wm_im,3),...
                             size(wm_im,4),size(wm_im,5),size(wm_im,6),size(wm_im,7));

    % Clean up
    clear temp_s temp_l temp_c
end

%% Store images if requested
if params.stim.store_imgs
    fprintf('[%]: Storing images',mfilename);
    saveDir = fileparts(fullfile(params.stim.ns.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s_%s.mat',params.stim.ns.stimfile,params.disp.name,datestr(now,30))),'scenes','ltm_lures','wm_im', 'im_order','info,'-v7.3');
end


return

