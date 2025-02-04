function scenes = vcd_naturalscenes(stim)
%
%  scenes = vcd_naturalscenes(stim)
%
% Purpose:
%   Load natural scene images for VCD experimental display.
%   Requires a valid stim.ns.stimfile, e.g.: fullfile(vcd_rootPath,'workspaces','scenes.mat')
%   or valid csv file in stim.ns.infofile, e.g.: fullfile(vcd_rootPath,'workspaces','scenes_info.csv')
%
% INPUTS:
%  stim         : struct with stimulus params 
%
% OUTPUTS:
% scenes        :   uint8 images of natural scenes, dimensions are
%                   width x height x 3 (rgb) x num exemplars, num of basic
%                   categories, num of superordinate categories
%
% Written by Eline Kupers 2024/12
%

%% Load existing resized images
if isfield(stim.ns, 'stimfile') && exist(stim.ns.stimfile,'file')
    load(stim.ns.stimfile,'scenes');

%% Load images from stim info file and resize if requested
elseif isfield(stim.ns, 'infofile') && exist(stim.ns.infofile,'file')
    
    stim.ns.stimfile = fullfile(vcd_rootPath,'workspaces','scenes.mat');
    
    t = readtable(stim.ns.infofile);
    
    % Define superordinate and basic categories, and number of exemplars per basic category
    superordinate = unique(t.superordinate,'stable'); % 4 superordinate categories
    basic = unique(t.basic,'stable');  % basic ,'kid_face', 'dog','bird','cow','tool','transport','sign','streets_trails','bathrooms','kitchens'
    n_exemplars = 7; %length(unique(t.exemplar)); % for now we have 6 exemplars (subordinate category level)

    
    if stim.iscolor
        scenes0 = uint8(ones(stim.ns.og_res_stim,stim.ns.og_res_stim,3,...
            n_exemplars,length(basic),length(superordinate)));
    else
        scenes0 = uint8(ones(stim.ns.og_res_stim,stim.ns.og_res_stim, ...
            n_exemplars*n_views,length(superordinate)));
    end
    
    for ss = 1:length(superordinate)
        bb = ss; % for now we only have one basic category per superordinate 
        for ex = 1:n_exemplars
            d = dir(fullfile(vcd_rootPath,'workspaces','ns', superordinate{ss},...
                basic{ss}, t.filename{ ...
                strcmp(t.superordinate,superordinate(ss)) & ...
                strcmp(t.basic,basic(bb)) & ...
                (t.exemplar==ex)}));
            if stim.iscolor
                scenes0(:,:,:,ex,bb,ss) = imread(fullfile(d.folder,d.name));
            else
                scenes0(:,:,ex,bb,ss) = imread(fullfile(d.folder,d.name));
            end
        end
    end
    
    % resize if desired
    if ~isempty(stim.ns.dres) && ~isequal(stim.ns.dres,1)
        tic;
        fprintf('[%]: resampling the stimuli; this may take a while',mfilename);
        if stim.iscolor
            scenes_rz = reshape(scenes0,size(scenes0,1),size(scenes0,2),size(scenes0,3), ...
                    size(scenes0,4)*size(scenes0,5)*size(scenes0,6));
                
            temp = cast([],class(scenes_rz));
            for pp = 1:size(scenes_rz,4)
                statusdots(pp,size(scenes_rz,4));
                temp(:,:,:,pp) = imresize(scenes_rz(:,:,:,pp),stim.ns.dres);
            end
            
        else
            scenes_rz = reshape(scenes0,size(scenes0,1),size(scenes0,2), ...
                size(scenes0,3)*size(scenes0,4)*size(scenes0,5));

            temp = cast([],class(scenes_rz(:,:,pp)));
            for pp = 1:size(scenes_rz,4)
                statusdots(pp,size(scenes_rz,3));
                temp(:,:,pp) = imresize(scenes_rz(:,:,pp),stim.ns.dres);
            end
        end
        
        fprintf('[%]: done!\n',mfilename);
        toc
        scenes_rz = temp;
    else
        scenes_rz = scenes0;
    end
    
    scenes = reshape(scenes_rz, size(scenes_rz,1),size(scenes_rz,2),size(scenes_rz,3),...
        size(scenes0,4),size(scenes0,5),size(scenes0,6));
    
    clear scenes0 scenes_rz
    
    if stim.store_imgs
        save(stim.ns.stimfile, 'scenes','-v7.3')
    end
else
    error('[%s]: Stimfile or infofile cannot be found or loaded',mfilename)
end

return