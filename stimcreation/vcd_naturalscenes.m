function [scenes,ltm_lures,wm_im,info] = vcd_naturalscenes(params)
% VCD function to load and resize/square pix values for natural scenes:
%
%  [scenes,ltm_lures,wm_im,info] = vcd_naturalscenes(params)
%
% Purpose:
%   Load subset of natural scene images used in the Natural Scenes Dataset
%   for specified VCD experimental display. If params.store_imgs is set to
%   true then this function will store the loaded images as a single mat
%   file as defined by p.stim.ns.stimfile, e.g.: 
%   fullfile(vcd_rootPath,'workspaces','stimuli',<disp_name>,'scenes_<disp_name>_YYYYMMDDTHHMMSS.mat')
%
% Requirements:
%   1. csv file with png file names defined by p.stim.ns.infofile, e.g.: 
%       fullfile(vcd_rootPath,'workspaces','scenes_info.csv')
%   2. Folder with png images, defined by p.stim.ns.indivscenefiles, e.g.: 
%       fullfile(vcd_rootPath,'workspaces','stimuli','RAW','vcd_natural_scenes')
%
% While loading and resizing scenes, we also create a struct "im_order" with
% information in the order that images are loaded, including scene super 
% and basic level semantic category, affordance type, scene loc, object 
% location in scene, lure name, wm test im name, unique image number.
%
% INPUTS:
%  params             : struct with stimulus params, with the following fields
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
%                       x 4 lure image types (1: most similar, 4 least similar)
%  wm_im        :   uint8 wm test images for each of the 30 VCD natural 
%                   scenes (think of a change blindness experiment)
%                   dimensions are width (pixels) x height (pixels) x 3
%                   (rgb) x 5 superordinate categories x 2 scene locations
%                   (indoor/outdoor) x 3 obj locations (left/middle/right)
%                   x 4 wm image types (1: easy_add, 2: hard_add,
%                   3: easy_remove, 4: hard_remove).
%  info         :   Loaded csv from workspaces/stimuli/ with png file names
%                   and scene/object category information. Has some overlap 
%                   of information with im_order, such as category info.
%      filename         : (cell) filename of original png.                
%      unique_im        : (double) unique image nr for each object: range 65-430   
%      superordinate    : (cell) superordinate semantic category label:
%                           {'human','animal','object','food','place'}
%      superordinate_i  : (double) same as superordinate, but indexed by nr
%                           1:'human' through 5: 'place'    
%      exemplar         : (cell) name for each image:
%                           {'face','cat','giraffe','donut',...
%                           'banana',vase','bus',...
%                           'bathroom','building'};
%      exemplar_i       : (double) index for each individual image within a
%                           superordinate class: 1-6 faces, 1-6 animals,
%                           1-6 food items, 1-6 objects, 1-6 places
%      basic           : (double) basic semantic category label:
%                           {'indoor', 'outdoor'}
%      basic_i          : (double) same as basic, but indexed by nr
%                           1: indoor, 2: outdoor
%      subordinate      : (cell) subordinate semantic category label, one
%                           of three object spatial locations relative to
%                           central fixation circle: {'left','center','right'};
%      subordinate_i    : (double) same as subordinate, but indexed by nr
%                           1:left, 2: center, 3: right. 
%      obj_act          : (cell) object affordance inferred by the experimenter 
%                          {'greet','grasp','walk','observe'}
%      obj_act_1        : (double) same as obj_act, but indexed by nr
%                           1:greet, 2: grasp, 3: walk, 4:observe
%      is_in_img_ltm    : (logical) whether the scene is part of the
%                           subselected stimuli used in imagery and
%                           long-term memory task.
%      is_lure          : (double) whether the scenes are lures for long-
%                           term memory task: ordered 1-6. Where 1 is most
%                           similar to core scene and 6 is least similar.
%                           Note that some core scenes only have 5 lures.
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

% Read info table
info = readtable(fullfile(d(end).folder, d(end).name));
[~,idx0] = intersect(params.stim.ns.unique_im_nrs_core,info.unique_im);

% Define superordinate and basic categories, and number of exemplars per basic category
superordinate = unique(info.superordinate,'stable'); % 5 superordinate categories
basic         = unique(info.basic,'stable');         % indoor/outdoor
subordinate   = unique(info.subordinate,'stable');   % dominant object location

% Delete empties
basic(cellfun(@isempty, basic)) = [];
subordinate(cellfun(@isempty, subordinate)) = [];

% WM image filenames
wm_test_im_name = repmat(params.stim.ns.change_im_name,1,length(idx0));

% unique image nrs
unique_im       = params.stim.ns.unique_im_nrs_core;

% Get info about images
n_scenes     = length(idx0)*length(basic)*length(subordinate);
n_ltm_lures  = length(params.stim.ns.lure_im);
n_wm_im      = length(params.stim.ns.change_im);

% Predefine im_order table
n_cols   = n_scenes + (n_scenes*n_ltm_lures) + (n_scenes*n_wm_im);
im_order = table(NaN(n_cols,1), ...unique_im
                cell(n_cols,1), ...super_cat
                cell(n_cols,1), ...exemplar
                cell(n_cols,1), ...affordance
                cell(n_cols,1), ...basic_cat
                cell(n_cols,1), ...sub_cat
                NaN(n_cols,1), ... wm change_im
                NaN(n_cols,1)); ...lure_im
im_order.Properties.VariableNames = {'unique_im','super_cat','exemplar','affordance','basic_cat','sub_cat','change_im','lure_im'};

% Preallocate space for images
scenes0 = uint8(ones(params.stim.ns.og_res_stim,params.stim.ns.og_res_stim,3,...
    length(superordinate), length(basic),length(subordinate)));
wmtest0 = uint8(ones(params.stim.ns.og_res_stim,params.stim.ns.og_res_stim,3,...
    length(superordinate), length(basic),length(subordinate),n_ltm_lures));
ltmlures0 = uint8(ones(params.stim.ns.og_res_stim,params.stim.ns.og_res_stim,3,...
    length(superordinate), length(basic),length(subordinate),n_wm_im));


%% Loop over superordinate, basic (inside/outside), and subordiante (object location) images
counter = 1;

for ss = 1:length(superordinate)
    for ex = 1:length(basic)
        for bb = 1:length(subordinate)
            
            % Get scene file name from info table
            im_png = info.filename{ ...
                strcmp(info.superordinate,superordinate(ss)) & ...
                strcmp(info.basic,basic(ex)) & ...
                strcmp(info.subordinate,subordinate(bb))};
            
            d = dir(fullfile(params.stim.ns.indivscenefile, im_png));
            
            % Read in image
            scenes0(:,:,:,ss,ex,bb) = imread(fullfile(d.folder,d.name));

            % Keep track of category etc
            im_order.super_cat(counter) = superordinate(ss);
            im_order.basic_cat(counter) = info.basic(bb);
            im_order.affordance(counter)= info.obj_act(bb);
            im_order.unique_im(counter) = unique_im((ex-1)*length(subordinate) + bb + (ss-1)*(length(basic)*length(subordinate)));
            im_order.lure_im(counter)   = NaN;
            im_order.change_im(counter) = NaN;
            
            row_nr = find(strcmp(info.filename,im_png));
            
            counter = counter+1;
            
            % Find corresponding wm test image images
            for cb_idx = 1:n_wm_im
                
                % Find image filename in info table
                og_scene = strsplit(im_png,'.png'); og_scene = og_scene{1};
                fn = sprintf('%s/%s.png',og_scene,wm_test_im_name{cb_idx});

                cb_im_name = info.filename{find(strcmp(info.filename,fn))};
                
                d = dir(fullfile(params.stim.ns.indivscenefile,...
                    'changes', cb_im_name));
                if isempty(d)
                    error('[%s]: can''t find wm test image file',mfilename)
                end                
                % Load image into array
                wmtest0(:,:,:,ss,ex,bb,cb_idx) = imread(fullfile(d.folder,d.name));

                % Keep track of category, etc.
                im_order.super_cat(counter) = superordinate(ss);
                im_order.exemplar(counter)  = info.exemplar(bb);
                im_order.affordance(counter)= info.obj_act(bb);
                im_order.basic_cat(counter) = basic(ex);
                im_order.sub_cat(counter)   = subordinate(bb);
                im_order.unique_im(counter) = params.stim.ns.unique_im_nrs_wm_test(cb_idx+((ex-1)*length(subordinate) + bb + (ss-1)*(length(basic)*length(subordinate))));
                im_order.lure_im(counter)   = false;
                im_order.change_im(counter) = cb_idx;
                
                counter = counter+1;
                
            end
            
            
            % Now do the same for LTM lures (can't use the same look up,
            % because nr of lures and wm test images varies)
            if info.is_in_img_ltm(row_nr)
            
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
                    im_order.exemplar(counter)  = info.exemplar(bb);
                    im_order.affordance(counter)= info.obj_act(bb);
                    im_order.basic_cat(counter) = basic(ex);
                    im_order.sub_cat(counter)   = subordinate(bb);
                    im_order.unique_im(counter) = params.stim.ns.unique_im_nrs_ltm_lures((ex-1)*length(subordinate) + bb +  (ss-1)*(length(basic)*length(subordinate)));
                    im_order.lure_im(counter)   = lure_idx;
                    im_order.change_im(counter) = false;
                    
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
    
    % collapse ltm lure conditions
    lures_rz = reshape(ltmlures0,size(ltmlures0,1),size(ltmlures0,2),size(ltmlures0,3), ...
        size(ltmlures0,4)*size(ltmlures0,5)*size(ltmlures0,6),size(ltmlures0,7));
    
    % collapse wm im conditions    
    wm_rz = reshape(wmtest0,size(wmtest0,1),size(wmtest0,2),size(wmtest0,3), ...
        size(wmtest0,4)*size(wmtest0,5)*size(wmtest0,6),size(wmtest0,7));
    
    % create temp var
    temp_sc = cast([],class(scenes_rz));
    temp_ltm = cast([],class(lures_rz));
    temp_wm = cast([],class(wm_rz));
    
    % loop over scenes, lures, wm images
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
        
        % Loop over wm images
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
    scenes0    = reshape(scenes,size(scenes,1),size(scenes,2),size(scenes,3),size(scenes,4)*size(scenes,5)*size(scenes,6));
    ltmlures0  = reshape(ltm_lures,size(ltm_lures,1),size(ltm_lures,2),size(ltm_lures,3),size(ltm_lures,4)*size(ltm_lures,5)*size(ltm_lures,6),size(ltm_lures,7));
    wmtest0    = reshape(wm_im,size(wm_im,1),size(wm_im,2),size(wm_im,3),size(wm_im,4)*size(wm_im,5)*size(wm_im,6),size(wm_im,7));
    
    for ss = 1:size(scenes0,4)
        
        % Plot scene
        clf;
        imagesc(scenes0(:,:,:,ss));
        
        title(sprintf('Im %02d resized',ss), 'FontSize',20);
        axis image; box off
        set(gca,'CLim',[1 255]);
        if params.stim.store_imgs
            saveFigDir1 = fullfile(vcd_rootPath,'figs',params.disp.name,'ns','visual_checks','resized_with_axes');
            if ~exist(saveFigDir1,'dir'), mkdir(saveFigDir1); end
            print(fullfile(saveFigDir1, sprintf('%03d_vcd_ns%02d', params.stim.ns.unique_im_nrs_core(ss), ss)),'-dpng','-r150');
            
            saveFigDir2 = fullfile(vcd_rootPath,'figs',params.disp.name,'ns','visual_checks','resized');
            if ~exist(saveFigDir2,'dir'), mkdir(saveFigDir2); end
            imwrite(scenes0(:,:,:,ss), fullfile(saveFigDir2, sprintf('%03d_vcd_ns%02d.png', params.stim.ns.unique_im_nrs_core(ss), ss)));
        end
        
        for ll = 1:size(ltmlures0,5)
            if ~isequal(ltmlures0(:,:,:,ss,ll),ones(size(ltmlures0,1),size(ltmlures0,2),size(ltmlures0,3)))
                % Plot LTM lure image
                clf; imagesc(ltmlures0(:,:,:,ss,ll));
                title(sprintf('Im %02d, LTM lure im %02d resized ',ss,ll), 'FontSize',20);
                axis image; box off
                set(gca,'CLim',[1 255]);
                
                if params.stim.store_imgs
                    print(fullfile(saveFigDir1, sprintf('%03d_vcd_ns%02d_ltm_lure%02d', params.stim.ns.unique_im_nrs_ltm_lures(ss), ss,ll)),'-dpng','-r150');
                    imwrite(ltmlures0(:,:,:,ss,ll), fullfile(saveFigDir2, sprintf('%03d_vcd_ns%02d_ltm_lure%02d.png', params.stim.ns.unique_im_nrs_ltm_lures(ss), ss,ll)));
                end
            end
            % Plot WM image
            clf;
            imagesc(wmtest0(:,:,:,ss,ll));
            title(sprintf('Im %02d, WM test im %02d resized',ss,ll), 'FontSize',20);
            axis image; box off
            set(gca,'CLim',[1 255]);
            
            if params.stim.store_imgs
                print(fullfile(saveFigDir1, sprintf('%03d_vcd_ns%02d_wm_im%02d',params.stim.ns.unique_im_nrs_wm_test(ss), ss,ll)),'-dpng','-r150');
                imwrite(wmtest0(:,:,:,ss,ll), fullfile(saveFigDir2, sprintf('%03d_vcd_ns%02d_wm_im%02d.png', params.stim.ns.unique_im_nrs_wm_test(ss), ss,ll)));
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
    
    % collapse wm category dim (w x h x 3 x 30 x 4)
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
                saveFigDir1 = fullfile(vcd_rootPath,'figs',params.disp.name,'ns','visual_checks','resized_and_squared_with_axes');
                if ~exist(saveFigDir1,'dir'), mkdir(saveFigDir1); end
                print(fullfile(saveFigDir1, sprintf('%03d_vcd_ns%02d', params.stim.ns.unique_im_nrs_core(ss), ss)),'-dpng','-r150');
                
                saveFigDir2 = fullfile(vcd_rootPath,'figs',params.disp.name,'ns','visual_checks','resized_and_squared');
                if ~exist(saveFigDir2,'dir'), mkdir(saveFigDir2); end
                imwrite(temp_sc(:,:,:,ss), fullfile(saveFigDir2, sprintf('%03d_vcd_ns%02d.png', params.stim.ns.unique_im_nrs_core(ss), ss)));
            end
            
            for ll = 1:size(lures_sq,5)
                if ~isequal(temp_ltm(:,:,:,ss,ll),ones(size(ltmlures0,1),size(ltmlures0,2),size(ltmlures0,3)))
                    clf;
                    imagesc(temp_ltm(:,:,:,ss,ll));
                    title(sprintf('Im %02d, ltm lure %02d resized & squared',ss,ll), 'FontSize',20);
                    axis image; box off
                    set(gca,'CLim',[1 255]);
                    if params.stim.store_imgs
                        print(fullfile(saveFigDir1, sprintf('%03d_vcd_ns%02d_ltm_lure%02d', params.stim.ns.unique_im_nrs_ltm_lures(ss),ss,ll)),'-dpng','-r150');
                        imwrite(temp_ltm(:,:,:,ss,ll), fullfile(saveFigDir2, sprintf('%03d_vcd_ns%02d_ltm_lure%02d.png', params.stim.ns.unique_im_nrs_ltm_lures(ss), ss,ll)));
                    end
                end
                clf;
                imagesc(temp_wm(:,:,:,ss,ll));
                title(sprintf('Im %02d, wm test im %02d resized & squared',ss,ll), 'FontSize',20);
                axis image; box off
                set(gca,'CLim',[1 255]);
                
                if params.stim.store_imgs
                    print(fullfile(saveFigDir1, sprintf('%03d_vcd_ns%02d_wm_im%02d', params.stim.ns.unique_im_nrs_wm_test(ss),ss,ll)),'-dpng','-r150');
                    imwrite(temp_wm(:,:,:,ss,ll), fullfile(saveFigDir2, sprintf('%03d_vcd_ns%02d_wm_im%02d.png', params.stim.ns.unique_im_nrs_wm_test(ss), ss,ll)));
                end
            end
        end
    end
    
    %% Reshape back 
    % Scenes 6 dims: (w x h x 3 x 5 super cat x 2 scene loc x 3 obj loc)
    scenes = reshape(temp_sc, size(scenes,1),size(scenes,2),size(scenes,3),...
                             size(scenes,4),size(scenes,5),size(scenes,6));
    % LTM lures 7 dims: (w x h x 3 x 5 super cat x 2 basic cat x 3 sub cat x 4 ltm lure images)
    ltm_lures  = reshape(temp_ltm, size(ltm_lures,1),size(ltm_lures,2),size(ltm_lures,3),...
                            size(ltm_lures,4),size(ltm_lures,5),size(ltm_lures,6),size(ltm_lures,7));
    % WM ims 7 dims: (w x h x 3 x 5 super cat x 2 basic cat x 3 sub cat x 4 wm test images)
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
    save(fullfile(sprintf('%s_%s_%s.mat',params.stim.ns.stimfile,params.disp.name,datestr(now,30))),'scenes','ltm_lures','wm_im', 'im_order','info','-v7.3');
end


return

