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
if isfield(p.stim.ns, 'stimfile') && exist(p.stim.ns.stimfile,'file')
    load(p.stim.ns.stimfile,'scenes');

%% Load images from stim info file and resize if requested
elseif isfield(p.stim.ns, 'infofile') && exist(p.stim.ns.infofile,'file')
    
    p.stim.ns.stimfile = fullfile(vcd_rootPath,'workspaces','stimuli','scenes.mat');
    
    t = readtable(p.stim.ns.infofile);
    
    % Define superordinate and basic categories, and number of exemplars per basic category
    superordinate = unique(t.superordinate,'stable'); % 4 superordinate categories
    basic         = unique(t.basic,'stable');  % basic 
    ns_loc        = unique(t.ns_loc,'stable'); % indoor/outdoor
    obj_loc       = unique(t.obj_loc,'stable'); % dominant object location

    n_images = length(superordinate)*length(ns_loc)*length(obj_loc);
    
    if p.stim.ns.iscolor
        scenes0 = uint8(ones(p.stim.ns.og_res_stim,p.stim.ns.og_res_stim,3,...
            length(superordinate), length(ns_loc),length(obj_loc)));
    else
        scenes0 = uint8(ones(p.stim.ns.og_res_stim,p.stim.ns.og_res_stim, ...
             length(superordinate), length(ns_loc),length(obj_loc)));
    end
    
    for ss = 1:length(superordinate)
        for ex = 1:length(ns_loc)
            for bb = 1:length(obj_loc)
                d = dir(fullfile(vcd_rootPath,'workspaces','stimuli','vcd_natural_scenes', ...
                    t.filename{ ...
                    strcmp(t.superordinate,superordinate(ss)) & ...
                    strcmp(t.ns_loc,ns_loc(ex)) & ...
                    strcmp(t.obj_loc,obj_loc(bb))}));
            
                if p.stim.ns.iscolor
                    scenes0(:,:,:,ss,ex,bb) = imread(fullfile(d.folder,d.name));
                else
                    scenes0(:,:,ss,ex,bb) = imread(fullfile(d.folder,d.name));
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
                
            temp = cast([],class(scenes_rz));
            for pp = 1:size(scenes_rz,4)
                statusdots(pp,size(scenes_rz,4));
                temp(:,:,:,pp) = imresize(scenes_rz(:,:,:,pp),p.stim.ns.dres);
            end
            
        else
            scenes_rz = reshape(scenes0,size(scenes0,1),size(scenes0,2), ...
                size(scenes0,3)*size(scenes0,4)*size(scenes0,5));

            temp = cast([],class(scenes_rz(:,:,pp)));
            for pp = 1:size(scenes_rz,4)
                statusdots(pp,size(scenes_rz,3));
                temp(:,:,pp) = imresize(scenes_rz(:,:,pp),p.stim.ns.dres);
            end
        end
        
        
        toc
        fprintf('[%]: done!\n',mfilename);
        scenes_rz = temp;
    else
        scenes_rz = scenes0;
    end
    
    scenes = reshape(scenes_rz, size(scenes_rz,1),size(scenes_rz,2),size(scenes_rz,3),...
        size(scenes0,4),size(scenes0,5),size(scenes0,6));
    
    clear scenes0 scenes_rz
    
    if p.stim.store_imgs
        save(p.stim.ns.stimfile, 'scenes','-v7.3')
    end
else
    error('[%s]: Stimfile or infofile cannot be found or loaded',mfilename)
end

return