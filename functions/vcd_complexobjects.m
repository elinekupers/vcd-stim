function objects = vcd_complexobjects(stim)
%
%  objects = vcd_complexobjects(stim)
%
% Purpose:
%   Load complex object images for VCD experimental display.
%   Requires a valid stim.cobj.stimfile, e.g.: fullfile(vcd_rootPath,'workspaces','objects.mat')
%   or valid csv file in stim.cobj.infofile, e.g.: fullfile(vcd_rootPath,'workspaces','objects_info.csv')
%
% INPUTS:
%  stim         : struct with stimulus params
%
% OUTPUTS:
% objects       :   uint8 images of complex images, dimensions are
%                   width x height x 3 (rgb) x num views (or poses), 
%                   num exemplars, num of basic categories,
%                   num of superordinate categories
%
% Written by Eline Kupers 2024/12
%

% load images
if isfield(stim.cobj, 'stimfile') && exist(stim.cobj.stimfile,'file')
    load(stim.cobj.stimfile,'objects');
else
    
    stim.cobj.stimfile = fullfile(vcd_rootPath,'workspaces','objects.mat');
    stim.cobj.infofile = fullfile(vcd_rootPath,'workspaces','object_info.csv');
    
    t = readtable(fullfile(stim.cobj.infofile));
    
    % Define superordinate and basic categories, and number of exemplars per basic category
    superordinate = unique(t.superordinate,'stable'); % 4 superordinate categories
    basic = unique(t.basic,'stable');  % 4 basic categories
    n_exemplars = 6; %length(unique(t.view)); % for now we have 6 exemplars (subordinate category level)
    n_views = 7;
    
    if stim.iscolor
        objects = uint8(ones(stim.cobj.img_sz_pix, stim.cobj.img_sz_pix,3,...
            n_views,n_exemplars,length(basic),length(superordinate)));
    else
        objects = uint8(ones(stim.cobj.img_sz_pix, stim.cobj.img_sz_pix, ...
            n_views,n_exemplars,length(basic),length(superordinate)));
    end
    
    for ss = 1:length(superordinate)
        bb = ss; %         for bb = 1:length(basic); % for now we only have one basic category per superordinate
        for ex = 1:n_exemplars
            for vv = 1:n_views
                % Get file name
                d = dir(fullfile(vcd_rootPath,'workspaces','cobj', superordinate{ss},...
                    basic{bb}, t.filename{ ...
                    strcmp(t.superordinate,superordinate(ss)) & ...
                    strcmp(t.basic,basic(bb)) & ...
                    (t.basic_i==vv) & ...
                    (t.view==ex)}));
                
                % Load image
                [im0, ~, transparency] = imread(fullfile(d.folder,d.name),'BackgroundColor','none');
                
                if stim.iscolor
                    % Add RGB dim when grayscale image
                    if ismatrix(im0)
                        im0 = uint8(repmat(im0,[1 1 3]));
                    end
                    
                    % If color, no background, then add gray background
                    if ~isempty(transparency)
                        transparency_3d = repmat(transparency,[1 1 3]);
                        im0(~transparency_3d) = 127;
                    end
                    
                else
                    % convert to gray if no color
                    if ~stim.iscolor && ndim(im0)==3
                        im0 = rgb2gray(im0);
                    end
                    % Add background to 2D grayscale image
                    if ~isempty(transparency)
                        im0(~transparency) = 127;
                    end
                end
                
                
                im = imresize(im0,[stim.cobj.img_sz_pix,stim.cobj.img_sz_pix]);
                clear im0
                
                if stim.iscolor
                    objects(:,:,:,vv,ex,bb,ss) = im;
                else
                    objects(:,:,vv,ex,bb,ss) = im;
                end
            end
        end
        
    end
    
    if stim.store_imgs    
        save(stim.cobj.stimfile, 'objects','-v7.3')
    end
end

%     figure(1); clf;
%     figure(2); clf;
%     figure(3); clf;
%     figure(4); clf;
%     for ii = 1:24
%         figure(1); subplot(4,6,ii); hold all;
%         imshow(objects(:,:,:,ii,1),[1 255]);
%
%         figure(2); subplot(4,6,ii); hold all;
%         imshow(objects(:,:,:,ii,2),[1 255]);
%
%         figure(3); subplot(4,6,ii); hold all;
%         imshow(objects(:,:,:,ii,3),[1 255]);
%
%         figure(4); subplot(4,6,ii); hold all;
%         imshow(objects(:,:,:,ii,4),[1 255]);
%     end
return

