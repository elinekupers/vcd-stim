function images = vcd_loadAllUniqueImages(params)
% VCD function to load in all unique, pregenerated RGB images for each
% class
% Note image files are stored locally (workspaces > stimuli), and need to
% be created prior to running this function (see s_createStim.m) for every
% monitor display setup, such that it has the <params.disp.name> in the
% stimulus file name.
       
%% Predefine cell arrays for stimuli
run_images = cell(size(subj_run,1),2); % second dim represent left and right side stimuli
run_alpha_masks = run_images;

% Check if we need to load images
if isfield(params, 'images')
    images = params.images;  params = rmfield(params, 'images');
end

if ~exist('images','var') || isempty(fieldnames(images))
    images = struct('gabor',[],'rdk',[],'dot',[],'obj',[],'ns',[],...
        'fix',[], 'info',[], 'image_order',[],'alpha',[]);
end

if isempty(images.gabor)
    fprintf('[%s]: Loading gabor images..',mfilename);
    % GABORS: 6D array: [x,y,8 orient, 4 phase,3 contrast, og + 4 delta]
    d = dir(sprintf('%s*.mat', params.stim.gabor.stimfile));
    load(fullfile(d(end).folder,d(end).name), 'gabors','masks','info');
    
    images.gabor = gabors; clear gabors;
    images.info.gabor  = info; clear info;
    images.alpha.gabor = masks; clear masks;
end

if isempty(images.rdk{1})
    fprintf('[%s]: Loading rdk images..',mfilename);
    % RDKs: 130 mat files: 8 mot dir x 3 coherence levels x 5 deltas (0 + 4 deltas)
    stimDir = dir(fullfile(sprintf('%s*',params.stim.rdk.stimfile)));
    
    rdk_files = dir(fullfile(stimDir,'*.mat'));
    
    for dd = 1:length(rdk_files)
        stimfile = fullfile(rdk_files(dd).folder,rdk_files(dd).name);
        load(stimfile, 'frames','mask','info');
        
        images.rdk{dd}       = frames; clear frames
        images.alpha.rdk{dd} = mask; clear mask;
        images.info.rdk{dd}  = info; clear info;
    end
end

if isempty(images.dot)
    % Simple dot: 2D array: [x,y]
    d = dir(sprintf('%s*.mat', params.stim.dot.stimfile));
    load(fullfile(d(end).folder,d(end).name), 'simple_dot','mask','info');
    images.dot = simple_dot; clear simple_dot;
    images.alpha.dot = mask; clear mask;
    images.info.dot = info; clear info;
end


if isempty(images.obj)
    fprintf('[%s]: Loading object images..',mfilename);
    % Complex objects: 4D array: [x,y,16 object, og + 10 rotation]
    d = dir(sprintf('%s*.mat', params.stim.obj.stimfile));
    load(fullfile(d(end).folder,d(end).name), 'objects','masks','info');
    images.obj = objects; clear objects;
    images.alpha.obj = masks; clear masks;
    images.info.obj = info; clear info;
end


if isempty(images.ns)
    fprintf('[%s]: Loading scene images..',mfilename);
    % NS: 6D array: [x,y,3, 5 superordinate cat, 2 ns_loc, 3 obj_loc]
    % CBlind: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 change images];
    % Lures: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 lure images];
    d = dir(sprintf('%s*.mat', params.stim.ns.stimfile));
    load(fullfile(d(end).folder,d(end).name), 'scenes','lures','cblind','info');
    images.ns = scenes; clear scenes;
    images.lures = lures; clear lures;
    images.cblind = cblind; clear cblind;
    images.info.ns = info; clear info;
end


fprintf('Done! ');  toc;
fprintf('\n');

return

