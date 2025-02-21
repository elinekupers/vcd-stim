function scenes = vcd_naturalscenes(p)
%
%  scenes = vcd_naturalscenes(p)
%
% Purpose:
%   Load natural scene images for VCD experimental display.
%   Requires a valid stim.ns.stimfile, e.g.: fullfile(vcd_rootPath,'workspaces','scenes.mat')
%   or valid csv file in stim.ns.infofile, e.g.: fullfile(vcd_rootPath,'workspaces','scenes_info.csv')
%
% INPUTS:
%  p            : struct with stimulus params
%
% OUTPUTS:
% scenes        :   uint8 images of natural scenes, dimensions are
%                   width x height x 3 (rgb) x num exemplars, num of basic
%                   categories, num of superordinate categories
%
% Written by Eline Kupers 2024/12
%

%% Load existing resized images
% if isfield(p.stim.ns, 'stimfile') && exist(p.stim.ns.stimfile,'file')
%     load(p.stim.ns.stimfile);

%% Load images from stim info file and resize if requested
if ~isfield(p.stim.ns, 'infofile')
    p.stim.ns.infofile = fullfile(vcd_rootPath,'workspaces','info','scene*');
end

d = dir(sprintf('%s*.csv',p.stim.ns.infofile));

if isempty(d)
    error('[%s]: Stimfile or infofile cannot be found or loaded',mfilename)
end


info = readtable(fullfile(d(end).folder, d(end).name));

% Define superordinate and basic categories, and number of exemplars per basic category
superordinate = unique(info.superordinate,'stable'); % 4 superordinate categories
basic         = unique(info.basic,'stable');  % basic
ns_loc        = unique(info.ns_loc,'stable'); % indoor/outdoor
obj_loc       = unique(info.obj_loc,'stable'); % dominant object location

n_images    = length(superordinate)*length(ns_loc)*length(obj_loc);
n_lures     = 4;
n_changeblindness = 4;

% Preallocate table about image order
n_cols = n_images + (n_images*n_lures) + (n_images*n_changeblindness);
im_order = table(cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),cell(n_cols,1),NaN(n_cols,1));
im_order.Properties.VariableNames = {'super_cat','basic_cat','affordance','scene_loc','obj_loc','lure_im','change_im','unique_im'};


if p.stim.ns.iscolor
    scenes0 = uint8(ones(p.stim.ns.og_res_stim,p.stim.ns.og_res_stim,3,...
        length(superordinate), length(ns_loc),length(obj_loc)));
    cblind0 = uint8(ones(p.stim.ns.og_res_stim,p.stim.ns.og_res_stim,3,...
        length(superordinate), length(ns_loc),length(obj_loc),n_lures));
    lures0 = uint8(ones(p.stim.ns.og_res_stim,p.stim.ns.og_res_stim,3,...
        length(superordinate), length(ns_loc),length(obj_loc),n_changeblindness));
else
    scenes0 = uint8(ones(p.stim.ns.og_res_stim,p.stim.ns.og_res_stim, ...
        length(superordinate), length(ns_loc),length(obj_loc)));
end

counter = 1;

for ss = 1:length(superordinate)
    for ex = 1:length(ns_loc)
        for bb = 1:length(obj_loc)
            
            im_png = info.filename{ ...
                strcmp(info.superordinate,superordinate(ss)) & ...
                strcmp(info.ns_loc,ns_loc(ex)) & ...
                strcmp(info.obj_loc,obj_loc(bb))};
            d = dir(fullfile(vcd_rootPath,'workspaces','stimuli','vcd_natural_scenes', im_png));
            
            if p.stim.ns.iscolor
                scenes0(:,:,:,ss,ex,bb) = imread(fullfile(d.folder,d.name));
            else
                scenes0(:,:,ss,ex,bb) = imread(fullfile(d.folder,d.name));
            end
            
            im_order.super_cat(counter) = superordinate(ss);
            im_order.basic_cat(counter) = info.basic(bb);
            im_order.affordance(counter)= info.obj_act(bb);
            im_order.scene_loc(counter) = ns_loc(ex);
            im_order.obj_loc(counter)   = obj_loc(bb);
            im_order.unique_im(counter) = (ex-1)*length(obj_loc) + bb + (ss-1)*(length(ns_loc)*length(obj_loc));
            im_order.lure_im(counter)   = {NaN};
            im_order.change_im(counter) = {NaN};
            
            counter = counter+1;
            
            for cb_idx = 1:n_changeblindness
                
                % change_blindness
                fn = sprintf('change_img%d',cb_idx);
                og_scene = strsplit(im_png,'.png'); og_scene = og_scene{1};
                
                cb_im_name = info.(fn){strcmp(info.superordinate,superordinate(ss)) & ...
                                        strcmp(info.ns_loc,ns_loc(ex)) & ...
                                        strcmp(info.obj_loc,obj_loc(bb))};
                
                d = dir(fullfile(vcd_rootPath,'workspaces','stimuli','vcd_natural_scenes',...
                    'changes', og_scene, cb_im_name));
                
                if ~isempty(d)
                    if p.stim.ns.iscolor
                        cblind0(:,:,:,ss,ex,bb,cb_idx) = imread(fullfile(d.folder,d.name));
                    else
                        cblind0(:,:,ss,ex,bb,cb_idx) = imread(fullfile(d.folder,d.name));
                    end
                else
                    if p.stim.ns.iscolor
                        cblind0(:,:,:,ss,ex,bb,cb_idx) = NaN(size(cblind0,1),size(cblind0,2),3);
                    else
                        cblind0(:,:,ss,ex,bb,cb_idx) = NaN(size(cblind0,1),size(cblind0,2));
                    end
                end
                
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
            
            for lure_idx = 1:n_lures
                % LTM lures
                fn = sprintf('lure_img%d',lure_idx);
                lure_im_name = info.(fn){ ...
                    strcmp(info.superordinate,superordinate(ss)) & ...
                    strcmp(info.ns_loc,ns_loc(ex)) & ...
                    strcmp(info.obj_loc,obj_loc(bb))};
                d = dir(fullfile(vcd_rootPath,'workspaces','stimuli','vcd_natural_scenes',...
                    'lures', og_scene, lure_im_name));
                if ~isempty(d)
                    if p.stim.ns.iscolor
                        lures0(:,:,:,ss,ex,bb,lure_idx) = imread(fullfile(d.folder,d.name));
                    else
                        lures0(:,:,ss,ex,bb,lure_idx) = imread(fullfile(d.folder,d.name));
                    end
                else
                    if p.stim.ns.iscolor
                        lures0(:,:,:,ss,ex,bb,lure_idx) = NaN(size(lures0,1),size(lures0,2),3);
                    else
                        lures0(:,:,ss,ex,bb,lure_idx) = NaN(size(lures0,1),size(lures0,2));
                    end
                end
                
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

% resize if desired
if ~isempty(p.stim.ns.dres) && ~isequal(p.stim.ns.dres,1)
    disp('Resampling the stimuli; this may take a while');
    tic;
    if p.stim.ns.iscolor
        scenes_rz = reshape(scenes0,size(scenes0,1),size(scenes0,2),size(scenes0,3), ...
            size(scenes0,4)*size(scenes0,5)*size(scenes0,6));
        
        lures_rz = reshape(lures0,size(lures0,1),size(lures0,2),size(lures0,3), ...
            size(lures0,4)*size(lures0,5)*size(lures0,6),size(lures0,7));
        
        cblind_rz = reshape(cblind0,size(cblind0,1),size(cblind0,2),size(cblind0,3), ...
            size(cblind0,4)*size(cblind0,5)*size(cblind0,6),size(cblind0,7));
        
        temp_s = cast([],class(scenes_rz));
        temp_l = cast([],class(lures_rz));
        temp_c = cast([],class(cblind_rz));
        
        
        for pp = 1:size(scenes_rz,4)
            statusdots(pp,size(scenes_rz,4));
            temp_s(:,:,:,pp) = imresize(scenes_rz(:,:,:,pp),p.stim.ns.dres);
            
            for ll = 1:n_lures
                temp_l(:,:,:,pp,ll) = imresize(lures_rz(:,:,:,pp,ll),p.stim.ns.dres);
            end
            for kk = 1:n_changeblindness
                temp_c(:,:,:,pp,kk) = imresize(cblind_rz(:,:,:,pp,kk),p.stim.ns.dres);
            end
        end
        
    else
        scenes_rz = reshape(scenes0,size(scenes0,1),size(scenes0,2), ...
            size(scenes0,3)*size(scenes0,4)*size(scenes0,5));
        lures_rz = reshape(lures0,size(lures0,1),size(lures0,2), ...
            size(lures0,3)*size(lures0,4)*size(lures0,5),size(lures0,7));
        
        cblind_rz = reshape(cblind0,size(cblind0,1),size(cblind0,2), ...
            size(cblind0,3)*size(cblind0,4)*size(cblind0,5),size(cblind0,7));
        
        temp_s = cast([],class(scenes_rz));
        temp_l = cast([],class(lures_rz));
        temp_c = cast([],class(cblind_rz));
        
        for pp = 1:size(scenes_rz,3)
            statusdots(pp,size(scenes_rz,3));
            temp_s(:,:,pp) = imresize(scenes_rz(:,:,pp),p.stim.ns.dres);
            
            for ll = 1:n_lures
                temp_l(:,:,pp,ll) = imresize(lures_rz(:,:,pp,ll),p.stim.ns.dres);
            end
            for kk = 1:n_changeblindness
                temp_c(:,:,pp,kk) = imresize(cblind_rz(:,:,pp,kk),p.stim.ns.dres);
            end
        end
    end
    
    
    toc
    fprintf('[%]: done!\n',mfilename);
    scenes_rz = temp_s;
    lures_rz  = temp_l;
    cblind_rz = temp_c;
else
    scenes_rz = scenes0;
    lures_rz  = lures0;
    cblind_rz = cblind0;
    
end

scenes = reshape(scenes_rz, size(scenes_rz,1),size(scenes_rz,2),size(scenes_rz,3),...
    size(scenes0,4),size(scenes0,5),size(scenes0,6));

lures = reshape(lures_rz, size(lures_rz,1),size(lures_rz,2),size(lures_rz,3),...
    size(lures0,4),size(lures0,5),size(lures0,6),size(lures0,7));

cblind = reshape(cblind_rz, size(cblind_rz,1),size(cblind_rz,2),size(cblind_rz,3),...
    size(cblind0,4),size(cblind0,5),size(cblind0,6),size(cblind0,7));

clear scenes0 lures0 cblind0 scenes_rz lures_rz cblind_rz

if p.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(p.stim.ns.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',p.stim.ns.stimfile,datestr(now,30))),'scenes','lures','cblind','info','im_order','-v7.3');
end


return

