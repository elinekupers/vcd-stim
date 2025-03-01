function [objects, masks] = vcd_complexobjects(p)
%
%  objects = vcd_complexobjects(p)
%
% Purpose:
%   Load complex object images for VCD experimental display.
%   Requires a valid stim.cobj.stimfile, e.g.: fullfile(vcd_rootPath,'workspaces','objects.mat')
%   or valid csv file in stim.cobj.infofile, e.g.: fullfile(vcd_rootPath,'workspaces','objects_info.csv')
%
% INPUTS:
%  p            :   struct with stimulus params
%
% OUTPUTS:
%  objects      :   uint8 images of complex images, dimensions are
%                   width x height x 3 (rgb) x num views (or poses), 
%                   num exemplars, num of basic categories,
%                   num of superordinate categories
%  masks        :   uint8 alpha masks of complex objects, removing the
%                   background: width x height x 3 (rgb) x num views (or poses), 
%                   num exemplars, num of basic categories,
%                   num of superordinate categories
%
% Written by Eline Kupers 2024/12
%

% load images
if isfield(p.stim.cobj, 'stimfile') && exist(p.stim.cobj.stimfile,'file')
    load(p.stim.cobj.stimfile,'objects');
else
    d = dir(sprintf('%s*',fullfile(p.stim.cobj.infofile)));
    info = readtable(fullfile(d(end).folder,d(end).name));
    
    % Define superordinate and basic categories, and number of exemplars per basic category
    mask            = info.unique_im>0;
    subordinate     = unique(info.subordinate,'stable');
    filename        = info.filename(mask);
    canonical_view  = info.rot_abs(mask);                        % alternating 2 canonical view ± 5 2-deg steps
    rotation_names  = [info.Properties.VariableNames(~cellfun(@isempty, regexp(info.Properties.VariableNames, 'rot_minus*'))), ...
                        info.Properties.VariableNames(~cellfun(@isempty, regexp(info.Properties.VariableNames, 'rot_plus*')))];
    rotations = 0;
    for ii = 1:length(rotation_names)
        rotations = [rotations, choose(regexp(rotation_names{ii}, 'rot_minus'),-1,1).*str2num(cell2mat(regexp(rotation_names{ii}, '\d+','match')))];
    end             
    
    rotations = rotations([1, 3:end-1]);
    rotation_names2 = catcell(2,{{'none'}, rotation_names(2:end-1)});
    
    % Preallocate info table
    n_cols = length(canonical_view)*length(rotations);
    im_order = table(cell(n_cols,1),NaN(n_cols,1),NaN(n_cols,1),cell(n_cols,1),NaN(n_cols,1));
    im_order.Properties.VariableNames = {'object_name','abs_rot','rel_rot','rot_name','unique_im'};
    
    % 
    og_im_sz = 1024;
    extent  = og_im_sz/2 + ([-1 1].*((p.stim.cobj.img_sz_pix+1)/1.5));
    new_im_z = diff(extent)+1;
    
    % Preallocate space
    if p.stim.cobj.iscolor
        objects = uint8(ones(new_im_z, new_im_z,3,...
            length(subordinate),length(rotations)));
        masks = objects;
    else
        objects = uint8(ones(new_im_z, new_im_z,2,...
            length(canonical_view),length(rotations)));
        masks = objects;
    end
    

    counter = 1;
    for sub = 1:length(subordinate)
        
        for rr = 1:length(rotations)
            
            % Get file name
            rot = (canonical_view(sub)/2) + mod(rotations(rr),2);
            d = dir(fullfile(vcd_rootPath,'workspaces','stimuli','vcd_complex_objects_2degstep_lumcorrected',...
                sprintf('*%s*_rot%d.png', subordinate{sub},rot)));
           
            % Load image
            [im0, ~, alpha0] = imread(fullfile(d.folder,d.name),'BackgroundColor','none');
            
            
            % we have done preprocessing already
%             if p.stim.cobj.iscolor
%                 % Gray --> Color: Add RGB dim when grayscale image
%                 if ismatrix(im0)
%                     im0 = uint8(repmat(im0,[1 1 3]));
%                 end
%                 
%                 % If color, no background, then add gray background
%                 if ~isempty(transparency)
%                     transparency_3d = repmat(transparency,[1 1 3]);
%                     im0(~transparency_3d) = p.stim.bckgrnd_grayval;
%                 end
%                 
%             else
%                 % convert to gray if no color
%                 if ~p.stim.cobj.iscolor && ndims(im0)==3
%                     im0 = rgb2gray(im0);
%                 end
%                 % Add background to 2D grayscale image
%                 if ~isempty(transparency)
%                     im0(~transparency) = p.stim.bckgrnd_grayval;
%                 end
%             end
            
            % assume range is 1-255
            im0     = (double(im0-1)./254).^2;
            alpha0  = (double(alpha0./255)).^2;
            
            % crop image
            cropme = extent(1):extent(2);
            
            alpha0_cropped = alpha0(cropme,cropme);
            im0_cropped = im0(cropme,cropme);

            
            im = uint8(1+(sqrt(im0_cropped)*254));
            alpha_im = uint8((sqrt(alpha0_cropped)*255));

%             im = imresize(im0_cropped,[p.stim.cobj.img_sz_pix,p.stim.cobj.img_sz_pix]);
%             aplha = imresize(alpha0_cropped,[p.stim.cobj.img_sz_pix,p.stim.cobj.img_sz_pix]);
            
            clear im0 alpha0 im0_cropped alpha0_cropped
            
%             figure(99); clf; imagesc(im0_cropped,'AlphaData',alpha0_cropped); colormap gray; colorbar; axis image
%             figure(99); clf; imagesc(im,'AlphaData',alpha_im); colormap gray; colorbar; axis image

            if p.stim.cobj.iscolor
                objects(:,:,:,sub,rr) = im;
                masks(:,:,:,sub,rr) = alpha_im;
            else
                objects(:,:,sub,rr) = im;
                masks(:,:,sub,rr) = alpha_im;
            end
            
            im_order.object_name(counter) = subordinate(sub);
            im_order.canonical_view(counter) = canonical_view(sub);
            im_order.abs_rot(counter) = canonical_view(sub)+rotations(rr);
            im_order.rel_rot(counter) = rotations(rr);
            im_order.rot_name(counter) = rotation_names2(rr);
            im_order.unique_im(counter) = sub;
            
            counter = counter+1;
        end
    end
    
    % reorder and select unique object images
    %     all_unique_im.cobj:
    % unique_im loc  super  basic  sub  facing_dir
    %      1     1     1     1     1    56
    %      2     2     1     2     2   124
    %      3     1     1     2     3    56
    %      4     2     2     1     1   124
    %      5     1     2     1     2    56
    %      6     2     2     2     3   124
    %      7     1     2     2     4    56
    %      8     2     3     1     1   124
    %      9     1     3     1     2    56
    %     10     2     3     2     3   124
    %     11     1     3     2     4    56
    %     12     2     3     3     5   124
    %     13     1     3     3     6    56
    %     14     2     4     1     1   124
    %     15     1     4     1     2    56
    %     16     2     4     1     3   124
    
    
    
    if p.stim.store_imgs    
        fprintf('\nStoring images..')
        saveDir = fileparts(fullfile(p.stim.cobj.stimfile));
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(sprintf('%s_%s.mat',p.stim.cobj.stimfile,datestr(now,30))),'objects','masks','info','im_order','-v7.3');
    end
end
% 
%     figure(1); clf;
%     figure(2); clf;
%     figure(3); clf;
%     figure(4); clf;
%     for ii = 1:size(objects,5)
%         figure(1); subplot(3,3,ii); hold all;
%         imshow(objects(:,:,1,1,ii),[1 255]);
% 
%         figure(2); subplot(3,3,ii); hold all;
%         imshow(objects(:,:,1,2,ii),[1 255]);
% 
%         figure(3); subplot(3,3,ii); hold all;
%         imshow(objects(:,:,1,3,ii),[1 255]);
% 
%         figure(4); subplot(3,3,ii); hold all;
%         imshow(objects(:,:,1,4,ii),[1 255]);
%     end
return

