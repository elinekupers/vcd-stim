function objects = vcd_complexobjects(p)
%
%  objects = vcd_complexobjects(p)
%
% Purpose:
%   Load complex object images for VCD experimental display.
%   Requires a valid stim.cobj.stimfile, e.g.: fullfile(vcd_rootPath,'workspaces','objects.mat')
%   or valid csv file in stim.cobj.infofile, e.g.: fullfile(vcd_rootPath,'workspaces','objects_info.csv')
%
% INPUTS:
%  p            : struct with stimulus params
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
if isfield(p.stim.cobj, 'stimfile') && exist(p.stim.cobj.stimfile,'file')
    load(p.stim.cobj.stimfile,'objects');
else

    info = readtable(fullfile(p.stim.cobj.infofile));
    
    % Define superordinate and basic categories, and number of exemplars per basic category
    superordinate = unique(info.superordinate,'stable'); % 4 superordinate categories
    basic         = unique(info.basic,'stable');                 % 8 basic categories (4x2)
    subordinate   = unique(info.subordinate,'stable');             % 16 sub categories (8x2)
    rotation      = unique(info.rot);
    n_views       = length(rotation);
    
    % Preallocate space
    if p.stim.cobj.iscolor
        objects = uint8(ones(p.stim.cobj.img_sz_pix, p.stim.cobj.img_sz_pix,3,...
            length(subordinate),n_views));
    else
        objects = uint8(ones(p.stim.cobj.img_sz_pix, p.stim.cobj.img_sz_pix, ...
            length(subordinate),n_views));
    end
    
    for sub = 1:length(subordinate)
        for rr = 1:length(rotation)
            
            % Get file name
            d = dir(fullfile(vcd_rootPath,'workspaces','stimuli','vcd_complex_objects_2degstep',...
                sprintf('*%s*_rot%d.png', subordinate{sub},(rotation(rr)/2)+1)));
            
            % Load image
            [im0, ~, transparency] = imread(fullfile(d.folder,d.name),'BackgroundColor','none');
            
            if p.stim.cobj.iscolor
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
                if ~p.stim.cobj.iscolor && ndims(im0)==3
                    im0 = rgb2gray(im0);
                end
                % Add background to 2D grayscale image
                if ~isempty(transparency)
                    im0(~transparency) = 127;
                end
            end
            
            
            im = imresize(im0,[p.stim.cobj.img_sz_pix,p.stim.cobj.img_sz_pix]);
            clear im0
            
            if p.stim.cobj.iscolor
                objects(:,:,:,sub,rr) = im;
            else
                objects(:,:,sub,rr) = im;
            end
        end
        
        
    end
    
    if p.stim.store_imgs    
        fprintf('\nStoring images..')
        saveDir = fileparts(fullfile(p.stim.cobj.stimfile));
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(sprintf('%s_%s.mat',p.stim.cobj.stimfile,datestr(now,30))),'objects','info','-v7.3');
    end
end
% 
%     figure(1); clf;
%     figure(2); clf;
%     figure(3); clf;
%     figure(4); clf;
%     for ii = 1:size(objects,4)
%         figure(1); subplot(2,11,ii); hold all;
%         imshow(objects(:,:,1,ii),[1 255]);
% 
%         figure(2); subplot(2,11,ii); hold all;
%         imshow(objects(:,:,2,ii),[1 255]);
% 
%         figure(3); subplot(2,11,ii); hold all;
%         imshow(objects(:,:,3,ii),[1 255]);
% 
%         figure(4); subplot(2,11,ii); hold all;
%         imshow(objects(:,:,4,ii),[1 255]);
%     end
return

