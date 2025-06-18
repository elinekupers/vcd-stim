function [scenes,ltm_lures,wm_im,info] = vcd_naturalscenes(params, verbose, store_imgs)
% VCD function to load and resize/square pix values for natural scenes:
%
%  [scenes,ltm_lures,wm_im,info] = vcd_naturalscenes(params, verbose, store_imgs)
%
% Purpose:
%   Load subset of natural scene images used in the Natural Scenes Dataset
%   for specified VCD experimental display. If store_imgs is set to
%   true then this function will store the loaded images as a single mat
%   file as defined by p.stim.ns.stimfile, e.g.:
%   fullfile(vcd_rootPath,'workspaces','stimuli',<disp_name>,'scenes_<disp_name>_YYYYMMDDTHHMMSS.mat')
%
% Requirements:
%   Folder with png images, defined by p.stim.ns.indivscenefiles, e.g.:
%       fullfile(vcd_rootPath,'workspaces','stimuli','RAW','vcd_natural_scenes')
%
% While loading and resizing scenes, we also create an "info" table with
% information in the order that images are loaded, including scene super
% and basic level semantic category, affordance type, scene loc, object
% location in scene, lure name, wm test im name, unique image number.
% this info table will be stored as a csv file defined by
% params.stim.ns.infofile (default is
% fullfile(vcd_rootPath,'workspaces','scenes_info.csv')).
%
% INPUTS:
%  params             : (struct) stimulus params with the following fields:
%    disp.name                  : name of used monitor (to check for squaring pixel values) (see vcd_getDisplayParams.m)
%    stim.ns.infofile           : csv file name (see vcd_getStimulusParams.m)
%    stim.ns.lure_im            : header names for lure images in csv file
%    stim.ns.change_im          : header names for change blindness images in csv file
%    stim.ns.og_res_stim        : original size of scene images (in pixels), can be height or width as we assume image is square
%    stim.ns.dres               : scale factor (fraction) applied to OG scene resolution
%    stim.ns.stimfile           : name of matfile final images will be stored
%    stim.ns.core_png_folder    : folder, where to find core scene pngs?
%    stim.ns.wmtest_png_folder  : folder, where to find individual wm test pngs?
%    stim.ns.ltmlure_png_folder : folder, where to find individual novel ltm lure pngs?
%    stim.ns.imgtest_png_folder : folder, where to find individual img test pngs?
%  verbose           : (logical) show debug figures
%  store_imgs        : (logical) store stimuli and debug figures as pngs 
%
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
%      stim_pos_name    : (cell) name of stimulus position on the screen.
%                           All scenes are {'center'}. (human readable)
%      stim_pos         : (double) index of "stim_pos"
%                           All scenes are 3. (machine readable)
%      super_cat_name   : (cell) superordinate semantic category label:
%                           {'human','animal','object','food','place'}
%      super_cat        : (double) same as superordinate, but indexed by nr
%                           1:'human' through 5: 'place'
%      basic_cat_name   : (double) basic semantic category label:
%                           {'indoor', 'outdoor'}
%      basic_cat        : (double) same as basic, but indexed by nr
%                           1: indoor, 2: outdoor
%      sub_cat_name     : (cell) subordinate semantic category label, one
%                           of three object spatial locations relative to
%                           central fixation circle: {'left','center','right'};
%      sub_cat          : (double) same as subordinate, but indexed by nr
%                           1:left, 2: center, 3: right.
%      affordance_name  : (cell) object affordance inferred by the experimenter
%                          {'greet','grasp','walk','observe'}
%      affordance_cat   : (double) same as obj_act, but indexed by nr
%                           1:greet, 2: grasp, 3: walk, 4:observe
%      is_specialcore   : (logical) whether the scene is part of the
%                           subselected stimuli used in imagery and
%                           long-term memory task.
%      is_lure          : (logical) whether the scenes are lures for long-
%                           term memory task: ordered 1-4. Where 1 is most
%                           similar to core scene and 4 is least similar.
%                           Note that we have additional lure images that
%                           we do not use.
%
% Written by Eline Kupers 2024/12, updated 2025/04
%

%% Load images from stim info file and resize if requested
if ~isfield(params.stim.ns, 'infofile')
    params.stim.ns.infofile = fullfile(vcd_rootPath,'workspaces','info','scene*');
end

% Read info table
% info = readtable(fullfile(d(end).folder, d(end).name));
info = table();
[info.filename, info.unique_im] = vcd_getNSfilenames;

% get image filenames
[~,idx0] = intersect(info.unique_im,params.stim.ns.unique_im_nrs_core);
[~,idx1] = intersect(info.unique_im,params.stim.ns.unique_im_nrs_wm_test);
[~,idx2] = intersect(info.unique_im,params.stim.ns.unique_im_nrs_ltm_lures);
% [~,idx3] = intersect(info.unique_im,params.stim.ns.unique_im_nrs_img_test);
core_im_name     = info.filename(idx0);
core_im_nr       = info.unique_im(idx0);
wm_test_im_name  = info.filename(idx1);
wm_test_im_nr    = info.unique_im(idx1);
ltm_lure_im_name = info.filename(idx2);
ltm_lure_im_nr    = info.unique_im(idx2);
% img_test_im_name = info.filename(idx3);

% Define superordinate and basic categories, and number of exemplars per basic category
superordinate = unique(params.stim.ns.super_cat, 'stable');                 % 5 superordinate categories
basic         = unique(catcell(2,params.stim.ns.basic_cat),'stable');       % 2 basic categories: indoor/outdoor
sub_tmp       = catcell(2,cellfun(@(x) strsplit(x,'_'), params.stim.ns.sub_cat{1}, 'UniformOutput', false));
sub_tmp       = reshape(sub_tmp',2,[]);
subordinate   = sub_tmp(2,:);                                               % 3 subordiante categories: dominant object location
clear sub_tmp

affordance = {};
for ii = 1:length(superordinate)
    affordance = cat(1,affordance, reshape(catcell(1,params.stim.ns.affordance(ii,:)),1,[])');
end
[~,affordance_i] = ismember(affordance,{'greet','grasp','walk','observe'});


% Get info about images
n_scenes     = length(core_im_name);
n_wm_im      = length(wm_test_im_name);
n_ltm_lures  = length(ltm_lure_im_name);
n_wm_changes = length(params.stim.ns.change_im_name);
n_ltm_lure_types = length(params.stim.ns.lure_im);

% Predefine im_order table
info.stim_pos_name   = repmat('center',size(info,1),1);
info.stim_pos        = repmat(3,size(info,1),1);
info.super_cat_name  = cat(1,repelem(superordinate,length(basic)*length(subordinate))', ... core im
    repelem(superordinate,n_wm_changes*length(basic)*length(subordinate))',... wm test
    repelem(superordinate,n_ltm_lures/length(superordinate))'); % ltm lures
info.super_cat       = cat(1,repelem(1:length(superordinate),length(basic)*length(subordinate))', ... core im
    repelem(1:length(superordinate),n_wm_changes*length(basic)*length(subordinate))',... wm test
    repelem(1:length(superordinate),n_ltm_lures/length(superordinate))'); % ltm lures
info.basic_cat_name  = cat(1,repmat(basic(:),length(superordinate)*length(subordinate),1), ... core im
    repelem(repmat(basic(:),length(superordinate)*length(subordinate),1),4),... wm test
    repmat({NaN},n_ltm_lures,1)); % ltm lures
info.basic_cat       = cat(1,repmat([1:length(basic)]',length(superordinate)*length(subordinate),1), ... core im
    repelem(repmat([1:length(basic)]',length(superordinate)*length(subordinate),1),4),... wm test
    NaN(n_ltm_lures,1)); % ltm lures
info.sub_cat_name    = cat(1,repmat(repelem(subordinate(:),length(basic)),length(superordinate),1), ... core im
    repmat(repelem(subordinate(:),length(basic)*n_wm_changes),length(superordinate),1),... wm test
    repmat({NaN},n_ltm_lures,1)); % ltm lures
info.sub_cat         = cat(1,repmat(repelem([1:length(subordinate)]',length(basic)),length(superordinate),1), ... core im
    repmat(repelem([1:length(subordinate)]',length(basic)*n_wm_changes),length(superordinate),1),... wm test
    NaN(n_ltm_lures,1)); % ltm lures
info.affordance_name = cat(1,affordance,repelem(affordance,n_wm_changes), repmat({NaN},n_ltm_lures,1));
info.affordance_cat  = cat(1,affordance_i,repelem(affordance_i,n_wm_changes), NaN(n_ltm_lures,1));
info.is_specialcore  = ismember(info.unique_im,params.stim.ns.unique_im_nrs_specialcore);  % is_specialcore (logical)
info.is_lure         = ismember(info.unique_im,params.stim.ns.unique_im_nrs_ltm_lures);  % is lure (logical)

% Preallocate space for images
scenes0 = uint8(ones(params.stim.ns.og_res_stim,params.stim.ns.og_res_stim,3,...
    length(superordinate), length(basic),length(subordinate)));
wmtest0 = uint8(ones(params.stim.ns.og_res_stim,params.stim.ns.og_res_stim,3,...
    length(superordinate), length(basic),length(subordinate),n_wm_changes));
ltmlures0 = uint8(ones(params.stim.ns.og_res_stim,params.stim.ns.og_res_stim,3,...
    length(superordinate), length(basic),length(subordinate),n_ltm_lure_types));


%% Loop over superordinate, basic (inside/outside), and subordiante (object location) images
specialcore_counter = 0;
counter = 0;

for ii = 1:n_scenes
    
    [~,ss] = ismember(info.super_cat_name(ii),superordinate);
    [~,tt] = ismember(info.basic_cat_name(ii), basic);
    [~,uu] = ismember(info.sub_cat_name(ii),subordinate);
    im_nr(ii) = info.unique_im(ii);
    
    % Get scene file name from info table
    fn_scn = info.filename{ii};
    
    d = dir(fullfile(params.stim.ns.core_png_folder, fn_scn));
    
    % Read in image
    scenes0(:,:,:,ss,tt,uu) = imread(fullfile(d.folder,d.name));
    
    % Find corresponding wm test image images
    for wm_idx = 1:n_wm_changes
        
        % Find image filename in info table
        fn_wm = wm_test_im_name{(n_wm_changes*(ii-1))+wm_idx};
        
        d = dir(fullfile(params.stim.ns.core_png_folder,'wm_test', fn_wm));
        if isempty(d)
            error('[%s]: can''t find wm test image file',mfilename)
        end
        % Load image into array
        wmtest0(:,:,:,ss,tt,uu,wm_idx) = imread(fullfile(d.folder,d.name));
    end
    
    % Now do the same for LTM lures (can't use the same counter,
    % because nr of lures and wm test images varies)
    if info.is_specialcore(ii)
        
        specialcore_counter = specialcore_counter +1;
        
        for lure_idx = 1:n_ltm_lure_types
            
            % Find image filename in info table
            fn_lure = ltm_lure_im_name{(n_ltm_lure_types*(specialcore_counter-1))+lure_idx};
            
            d = dir(fullfile(params.stim.ns.core_png_folder,'ltm_novel_lures', fn_lure));
            if isempty(d)
                error('[%s]: can''t find lure image file',mfilename)
            end
            
            % Load image into array
            ltmlures0(:,:,:,ss,tt,uu,lure_idx) = imread(fullfile(d.folder,d.name));
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
        for ll = 1:n_ltm_lure_types
            l0   = double(lures_rz(:,:,:,pp,ll)); % image lum vals ranges from [0-255]
            l1   = uint8(imresize(l0,params.stim.ns.dres));
            
            % check for clipping
            assert( min(l1(:))>=0 && max(l1(:))<=255)
            
            temp_ltm(:,:,:,pp,ll) = l1;
        end
        
        % Loop over wm images
        for kk = 1:n_wm_changes
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

if verbose
    %% Visualize resized images
    makeprettyfigures;
    figure; set(gcf,'Position', [156    91   881   706],'color','w');
    
    specialcore_counter = 0;
    
    for ii = 1:n_scenes
        
        [~,ss] = ismember(info.super_cat_name(ii),superordinate);
        [~,tt] = ismember(info.basic_cat_name(ii), basic);
        [~,uu] = ismember(info.sub_cat_name(ii),subordinate);
        im_nr(ii) = info.unique_im(ii);
        
        % Plot scene
        clf;
        imagesc(scenes(:,:,:,ss,tt,uu)); drawnow;
        
        title(sprintf('Im %02d resized',ii), 'FontSize',20);
        axis image; box off
        set(gca,'CLim',[1 255]);
        if store_imgs
            saveFigDir2 = fullfile(vcd_rootPath,'figs',params.disp.name,'ns','resized');
            if ~exist(saveFigDir2,'dir'), mkdir(saveFigDir2); end
            imwrite(scenes(:,:,:,ss,tt,uu), fullfile(saveFigDir2, sprintf('%04d_vcd_ns%02d.png', im_nr(ii), ii)));
        end

        for wm_idx = 1:size(wm_im,7)
            % Plot WM image
            clf;
            imagesc(wm_im(:,:,:,ss,tt,uu,wm_idx));
            title(sprintf('Im %02d, WM test im %02d resized',ii,wm_idx), 'FontSize',20);
            axis image; box off
            set(gca,'CLim',[1 255]);
            
            if store_imgs
                imwrite(wm_im(:,:,:,ss,tt,uu,wm_idx), fullfile(saveFigDir2, sprintf('%04d_vcd_ns%02d_wm_im%02d.png', wm_test_im_nr(((ii-1)*n_wm_changes)+ll), ii,wm_idx)));
            end
        end
        
        
        if info.is_specialcore(ii)
            
            specialcore_counter = specialcore_counter +1;
            
            for lure_idx = 1:n_ltm_lure_types

                % Plot LTM lure image
                clf; imagesc(ltm_lures(:,:,:,ss,tt,uu,lure_idx));
                title(sprintf('NS %02d, LTM lure im %02d resized ',ii,lure_idx), 'FontSize',20);
                axis image; box off
                set(gca,'CLim',[1 255]);
                
                if store_imgs
                    imwrite(ltm_lures(:,:,:,ss,tt,uu,lure_idx), ...
                        fullfile(saveFigDir2, sprintf('%04d_vcd_ns%02d_ltm_lure%02d.png', ltm_lure_im_nr(((specialcore_counter-1)*n_ltm_lure_types)+lure_idx), ii,lure_idx)));
                end

            end
        end
    end
end

%% For monitors acting linearly, square pixel values
if any(strcmp(params.disp.name, {'7TAS_BOLDSCREEN32', 'PPROOM_EIZOFLEXSCAN','CCNYU_VIEWPIXX3D'}))
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
    
    
    % Reshape back
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
    clear temp_sc temp_ltm temp_wm
end

% Visualize squared and resized images
if verbose && any(strcmp(params.disp.name, {'7TAS_BOLDSCREEN32', 'PPROOM_EIZOFLEXSCAN','CCNYU_VIEWPIXX3D'}))
    
    specialcore_counter = 0;
    
    figure; set(gcf,'Position', [156    91   881   706],'color','w');
    
    for ii = 1:n_scenes
        
        [~,ss] = ismember(info.super_cat_name(ii),superordinate);
        [~,tt] = ismember(info.basic_cat_name(ii), basic);
        [~,uu] = ismember(info.sub_cat_name(ii),subordinate);
        im_nr(ii) = info.unique_im(ii);
        
        clf;
        imagesc(scenes(:,:,:,ss,tt,uu)); drawnow
        title(sprintf('Im %02d resized & squared',ss), 'FontSize',20);
        axis image; box off
        set(gca,'CLim',[1 255]);
        if store_imgs
            saveFigDir2 = fullfile(vcd_rootPath,'figs',params.disp.name,'ns','resized_and_squared');
            if ~exist(saveFigDir2,'dir'), mkdir(saveFigDir2); end
            imwrite(scenes(:,:,:,ss,tt,uu), fullfile(saveFigDir2, sprintf('%04d_vcd_ns%02d.png',im_nr(ii), ii)));
        end
        
        for wm_idx = 1:size(wm_im,7)
            
            clf;
            imagesc(wm_im(:,:,:,ss,tt,uu,wm_idx)); drawnow;
            title(sprintf('Im %02d, wm test im %02d resized & squared',ii,wm_idx), 'FontSize',20);
            axis image; box off
            set(gca,'CLim',[1 255]);
            
            if store_imgs
                imwrite(wm_im(:,:,:,ss,tt,uu,wm_idx), fullfile(saveFigDir2, sprintf('%04d_vcd_ns%02d_wm_im%02d.png', wm_test_im_nr(((ii-1)*n_wm_changes)+wm_idx), ii,wm_idx)));
            end
        end
        
        if info.is_specialcore(ii)
            
            specialcore_counter = specialcore_counter +1;
            
            for lure_idx = 1:n_ltm_lure_types
                title(sprintf('Im %02d, ltm lure %02d resized & squared',ii,lure_idx), 'FontSize',20);
                axis image; box off
                set(gca,'CLim',[1 255]);
                if store_imgs
                    imwrite(ltm_lures(:,:,:,ss,tt,uu,lure_idx), ...
                        fullfile(saveFigDir2, sprintf('%04d_vcd_ns%02d_ltm_lure%02d.png', ltm_lure_im_nr(((specialcore_counter-1)*n_ltm_lure_types)+lure_idx), ii,lure_idx)));
                end
                
            end
        end
    end
end



%% Store images if requested
if store_imgs
    fprintf('[%]: Storing images',mfilename);
    saveDir = fileparts(fullfile(params.stim.ns.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',params.stim.ns.stimfile,datestr(now,30))),'scenes','ltm_lures','wm_im','info','-v7.3');
    
    saveDir = fileparts(fullfile(params.stim.ns.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',params.stim.ns.infofile,datestr(now,30))))
end


return

