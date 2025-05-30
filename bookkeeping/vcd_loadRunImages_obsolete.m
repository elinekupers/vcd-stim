function [run_images, run_alpha_masks, images] = vcd_loadRunImages_obsolete(run_image_order, block, params)

%% %%%%%%%%%%%%% PRE-LOAD STIMULI %%%%%%%%%%%%%
if isfield(params, 'images')
    images = params.images;  params = rmfield(params, 'images');
end

if ~exist('images','var') || isempty(fieldnames(images))
    images = struct('gabor',[],'rdk',[],'dot',[],'obj',[],'ns',[],...
        'fix',[], 'info',[], 'image_order',[],'alpha',[]);
end

fprintf('[%s]: Loading images..',mfilename); tic;
for ii = 1:length(run_image_order)
    statusdots(ii,length(block));
    
    st_crossing = block(ii).name;
    tmp = strsplit(st_crossing,'-');
    taskClass = tmp{1};
    stimClass = tmp{2};
    
    if any(strcmp(taskClass,{'pre','post'}))
        % do nothing
    else
        % check we assigned the string to the right label
        assert(any(strcmp(taskClass,params.exp.taskClassLabels)));
        assert(any(strcmp(stimClass,params.exp.stimClassLabels)));
        
        fprintf('\n\t%s: ',stimClass)
        switch stimClass
            case 'gabor'
                sz = size(block(ii).trial);
                if sz(1)>sz(2) && any(sz~=1)
                    numTrials = sz(1);
                    numSides = sz(2);
                elseif sz(1)<sz(2) && any(sz~=1)
                    numTrials = sz(2);
                    numSides = sz(1);
                end
                
                if isempty(images.gabor)
                    % GABORS: 6D array: [x,y,8 orient, 4 phase,3 contrast, og + 4 delta]
                    d = dir(sprintf('%s*.mat', params.stim.gabor.stimfile));
                    load(fullfile(d(end).folder,d(end).name), 'gabors','masks','info');
                    
                    images.gabor = gabors; clear gabors;
                    images.info.gabor  = info; clear info;
                    images.alpha.gabor = masks; clear masks;
                end
                
                for jj = 1:numTrials
                    for nn = 1:numSides
                        i3 = run_image_order{ii}{jj,1}(nn);
                        run_images{ii,jj,nn,1} = images.gabor(:,:,:,i3,1);
                        run_alpha_masks{ii,jj,nn,1} = images.alpha.gabor(:,:,i3,1);
                        
                        if strcmp(taskClass,'wm')
                            i4 = run_image_order{ii}{jj,2}(nn);
                            run_images{ii,jj,nn,2} = images.gabor(:,:,:,i3,i4+1);
                            run_alpha_masks{ii,jj,nn,1} = images.alpha.gabor(:,:,i3,i4+1);
                        end
                        %                 if strcmp(taskClass,'ltm')
                        %
                        %                     run_images{ii,jj,1,1} = images.ns(:,:,:,i4,i5,i6);
                        %                 end
                    end
                end
                
            case 'rdk'
                sz = size(block(ii).trial);
                if sz(1)>sz(2) && any(sz~=1)
                    numTrials = sz(1);
                    numSides = sz(2);
                elseif sz(1)<sz(2) && any(sz~=1)
                    numTrials = sz(2);
                    numSides = sz(1);
                end
                
                for jj = 1:numTrials
                    for nn = 1:numSides
                        % RDKs: 130 mat files: 8 directions x 3 coherence levels x 5 deltas (0 + 4 deltas)
                        stimDir = dir(fullfile(sprintf('%s*',params.stim.rdk.stimfile)));
                        filename = sprintf('%d_rdk_ori%d_coh%d_delta%d',...
                            block(ii).trial(jj,nn).unique_im_nr, ...
                            find(block(ii).trial(jj,nn).motdir==params.stim.rdk.dots_direction), ...
                            find(block(ii).trial(jj,nn).coh==params.stim.rdk.dots_coherence), ...
                            0);
                        
                        stimfile = fullfile(stimDir(1).folder,stimDir(1).name,sprintf('%s.mat', filename));
                        if exist(stimfile,'file')
                            load(stimfile, 'frames','mask');
                        else
                            error('[%s]: Can''t find RDK stim file!')
                        end
                        
                        % each file contains a 4D array: [x,y,3,frames]
                        run_images{ii,jj,nn,1} = frames;
                        run_alpha_masks{ii,jj,nn,1} = mask;
                        
                        if strcmp(taskClass,'wm')
                            delta_idx = find((block(ii).trial(jj,nn).motdir - block(ii).trial(jj,nn).ref_delta) == params.stim.rdk.delta_from_ref);
                            
                            filename_delta = sprintf('%d_rdk_ori%d_coh%d_delta%d',...
                                block(ii).trial(jj,nn).unique_im_nr, ...
                                find(block(ii).trial(jj,nn).motdir==params.stim.rdk.dots_direction), ...
                                find(block(ii).trial(jj,nn).coh==params.stim.rdk.dots_coherence), ...
                                delta_idx);
                            
                            stimfile_delta = fullfile(stimDir(1).folder,stimDir(1).name,sprintf('%s.mat', filename_delta));
                            
                            load(stimfile_delta, 'frames');
                            run_images{ii,jj,nn,2} = frames;
                            run_alpha_masks{ii,jj,nn,2} = run_alpha_masks{ii,jj,nn,1};
                        end
                        %                 if strcmp(taskClass,'ltm')
                        %
                        %                     run_images{ii,jj,1,1} = images.ns(:,:,:,i4,i5,i6);
                        %                 end
                    end
                    
                end
                
                
            case 'dot'
                sz = size(block(ii).trial);
                if sz(1)>sz(2) && any(sz~=1)
                    numTrials = sz(1);
                    numSides = sz(2);
                elseif sz(1)<sz(2) && any(sz~=1)
                    numTrials = sz(2);
                    numSides = sz(1);
                end
                
                if isempty(images.dot)
                    % Simple dot: 2D array: [x,y]
                    d = dir(sprintf('%s*.mat', params.stim.dot.stimfile));
                    load(fullfile(d(end).folder,d(end).name), 'simple_dot','mask','info');
                    images.dot = simple_dot; clear simple_dot;
                    images.alpha.dot = mask; clear mask;
                    images.info.dot = info; clear info;
                end
                
                if ndims(images.alpha.dot)==3
                    images.alpha.dot = images.alpha.dot;
                end
                for jj = 1:numTrials
                    for nn = 1:numSides
                        
                        % regular stim interval
                        run_images{ii,jj,nn,1} = images.dot;
                        run_alpha_masks{ii,jj,nn,1} = images.alpha.dot;
                        
                        if strcmp(taskClass,'wm')
                            run_images{ii,jj,nn,2} = images.dot;
                            run_alpha_masks{ii,jj,nn,2} = images.alpha.dot;
                        end
                        %                 if strcmp(taskClass,'ltm')
                        %
                        %                     run_images{ii,jj,1,1} = images.ns(:,:,:,i4,i5,i6);
                        %                 end
                    end
                end
                
            case 'obj'
                sz = size(block(ii).trial);
                if sz(1)>sz(2) && any(sz~=1)
                    numTrials = sz(1);
                    numSides = sz(2);
                elseif sz(1)<sz(2) && any(sz~=1)
                    numTrials = sz(2);
                    numSides = sz(1);
                end
                
                if isempty(images.obj)
                    % Complex objects: 4D array: [x,y,16 object, og + 10 rotation]
                    d = dir(sprintf('%s*.mat', params.stim.obj.stimfile));
                    load(fullfile(d(end).folder,d(end).name), 'objects','masks','info','im_order');
                    images.obj = objects; clear objects;
                    images.alpha.obj = masks; clear masks;
                    images.info.obj = info; clear info;
                    images.image_order.obj = im_order; clear im_order;
                end
                
                for jj = 1:numTrials
                    for nn = 1:numSides
                        [i4,i5] = ind2sub([size(images.obj,4),size(images.obj,5)],run_image_order{ii}{jj,1}(nn));
                        
                        run_images{ii,jj,nn,1} = images.obj(:,:,:,i4,i5);
                        run_alpha_masks{ii,jj,nn,1} = images.alpha.obj(:,:,i4,i5);
                        
                        if strcmp(taskClass,'wm')
                            [i4,i5] = ind2sub([size(images.obj,4),size(images.obj,5)],run_image_order{ii}{jj,2}(nn));
                            run_images{ii,jj,nn,2} = images.obj(:,:,:,i4,i5);
                            run_alpha_masks{ii,jj,nn,2} = images.alpha.obj(:,:,i4,i5);
                        end
                        %                 if strcmp(taskClass,'ltm')
                        %
                        %                     run_images{ii,jj,1,1} = images.ns(:,:,:,i4,i5,i6);
                        %                 end
                        
                    end
                end
                
            case 'ns'
                sz = size(block(ii).trial);
                if sz(1)>sz(2) && any(sz~=1)
                    numTrials = sz(1);
                elseif sz(1)<sz(2) && any(sz~=1)
                    numTrials = sz(2);
                end
                
                if isempty(images.ns)
                    % NS: 6D array: [x,y,3, 5 superordinate cat, 2 ns_loc, 3 obj_loc]
                    % CBlind: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 change images];
                    % Lures: 7D array: [x,y,5 superordinate cat, 2 ns_loc, 3 obj_loc, 4 lure images];
                    d = dir(sprintf('%s*.mat', params.stim.ns.stimfile));
                    load(fullfile(d(end).folder,d(end).name), 'scenes','lures','cblind','info','im_order');
                    images.ns = scenes; clear scenes;
                    images.lures = lures; clear lures;
                    images.cblind = cblind; clear cblind;
                    images.info.ns = info; clear info;
                    images.image_order.ns = im_order; clear im_order;
                end
                
                for jj = 1:numTrials
                    curr_im = images.info.ns(run_image_order{ii}{jj,1},:);
                    i4 = curr_im.superordinate_i;
                    i5 = curr_im.ns_loc_i;
                    i6 = curr_im.obj_loc_i;
                    
                    % third dim has image and alpha mask
                    run_images{ii,jj,1,1} = images.ns(:,:,:,i4,i5,i6);
                    
                    if strcmp(taskClass,'wm')
                        i7 = run_image_order{ii}{jj,2};
                        run_images{ii,jj,1,2} = images.cblind(:,:,:,i4,i5,i6,i7);
                    end
                    %                 if strcmp(taskClass,'ltm')
                    %                     i7 = run_image_order{ii}{jj,2};
                    %                     run_images{ii,jj,1,2} = images.cblind(:,:,:,i4,i5,i6,i7);
                    %                 end
                end
                
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% CHECK CLUT & COLOR & SIZE   %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         if ~isempty(params.stim.(stimClass).dres) && ...
        %                 sum(params.stim.(stimClass).dres)~=0
        fprintf('check stim size..');
        for jj = 1:numTrials
            
            statusdots(jj,length(numTrials));
            
            tmp0 = squeeze(run_images(ii,jj,:,:));
            tmp0_mask = squeeze(run_alpha_masks(ii,jj,:,:));

            sz0 = size(tmp0);
            
            tmp0 = reshape(tmp0,1,[]);
            tmp0_mask = reshape(tmp0_mask,1,[]);
            
            im_rz = []; im_mask_rz = [];
            for kk = 1:length(tmp0)
                
                if ~isempty(tmp0{kk})
                    sz = cell2mat(reshape(cellfun(@size, tmp0(kk),'UniformOutput',false),[],1));
                    
                    if params.stim.(stimClass).iscolor && (length(sz) == 3); sz = sz([1,2]); end
                    if strcmp(stimClass,'rdk'); sz = sz([1,2]); end
                    if isequal(sz(1),sz(2)), sz = sz(1); end
                    
                     
                    if any(params.stim.(stimClass).img_sz_pix ~= sz)
                        
                        % recompute scale factor in case we didn't adjust
                        % stim size to the current monitor
                        scale_factor = params.stim.(stimClass).img_sz_pix/sz;
                        params.stim.(stimClass).dres = scale_factor; % overwrite scale factor
                        
                        % reshape to get one vector of pixels by 1 (or time)
                        tmp0_im = double(tmp0{kk});
                        
                        if ~isempty(tmp0_mask{kk})
                            tmp0_im_mask = tmp0_mask{kk};
                        end
                        %                         % GRAY SCALE has 2 or 3 dim
                        %                         if ~params.stim.(stimClass).iscolor
                        %                             tmp_im = reshape(tmp0_im, ...
                        %                                 size(tmp0_im,1),... x
                        %                                 size(tmp0_im,2),... y
                        %                                 []);
                        %
                        %                             % COLOR has 3 or 4 dim
                        %                         elseif params.stim.(stimClass).iscolor
                        %                             tmp_im = reshape(tmp0_im, ...
                        %                                 size(tmp0_im,1),... x
                        %                                 size(tmp0_im,2),... y
                        %                                 3, ... % rgb
                        %                                 []);
                        %                         end
                        if ndims(tmp0_im)==4
                            numFrames = size(tmp0_im,4);
                        
                            for ll = 1:numFrames
                                im_rz(:,:,:,ll) = imresize(tmp0_im(:,:,:,ll),scale_factor); %% DEFAULT IS BICUBIC
                            end
                            
                            if exist('tmp0_im_mask','var') % assume there is only one mask, the same for every time point 
                                im_mask_rz = imresize(tmp0_im_mask,scale_factor); %% DEFAULT IS BICUBIC
                            end
                        
                        else
                            im_rz = imresize(tmp0_im,scale_factor);
                            if exist('tmp0_im_mask','var')
                                im_mask_rz = imresize(tmp0_im_mask,scale_factor); %% DEFAULT IS BICUBIC
                            end
                            
                        end
                        
                        run_im_tmp{kk} = im_rz;
                        run_im_tmp_mask{kk} = im_mask_rz;
                    else
                        run_im_tmp{kk} = tmp0{kk};
                        run_im_tmp_mask{kk} = tmp0_mask{kk};
                    end
                end
            end
            
            % Reshape trial images back to 2xN
            run_im_tmp = reshape(run_im_tmp,sz0(1),sz0(2));
            run_images(ii,jj,:,:) = run_im_tmp;
            
            if exist('run_im_tmp_mask','var')
                run_im_tmp_mask = reshape(run_im_tmp_mask,sz0(1),sz0(2));
                run_alpha_masks(ii,jj,:,:) = run_im_tmp_mask;
            end
        end
    end
    
    
    % If using BOLDSCREEN, we want to square pixel values for objects and
    % scenes
    if strcmp(params.disp.name,'7TAS_BOLDSCREEN32') && ...
            ~any(strcmp(taskClass,{'pre','post'}))  && ...
        any(strcmp(stimClass,{'obj','ns'}))
        
        assert(params.stim.(stimClass).square_pix_val==true)
        
        fprintf('square pix values for CLUT..');
        % reshape to get one vector of images
        run_im_tmp0 = squeeze(run_images(ii,jj,:,:));
        sz0 = size(run_im_tmp0);
        run_im_tmp = reshape(run_im_tmp0,[],1);
        
        for mm = 1:length(run_im_tmp)
            statusdots(mm,length(run_im_tmp));
            img = run_im_tmp{mm};
            
            if ~isempty(img)
                if ~params.stim.(stimClass).iscolor
                    img = double(img);
                    pix_range = [min(img(:)),max(img(:))];
                    
                    % if pix lum range is 0-1
                    if pix_range(2)<=1
                        img  = floor((img*255)+params.stim.bckgrnd_grayval);
                        img = img.^2;
                        
                        % if pix lum range is 1-255
                    elseif pix_range(2)<=255
                        if ismatrix(img)
                            img = img.^2;
                        elseif ndims(img)==3 % rdks
                            img_vec = reshape(img,size(img,1)*size(img,2),[]);
                            img_vec = img_vec.^2;
                            img = reshape(img_vec,size(img,1),size(img,2),size(img,3));
                        end
                    end
                    
                    % (Q: if pix lum range is 0-255, do we shift to 1-255??)
                    % elseif pix_range(1)<1 && pix_range(2)<=255
                    %   img = img.^2;
                    % end
                    img = uint8(img);
                    
                    % if image contains color
                elseif params.stim.(stimClass).iscolor
                    img = double(img); % convert to double
                    img = ((img./255).^2)*255; % rescale range to [0-1], square img, return 255 range
                    img = uint8(img); % convert back to uint8
                end
            end
            
            run_im_tmp{mm} = img;
        end
        
        run_images(ii,jj,:,:) = reshape(run_im_tmp, sz0(1),sz0(2));
        
    end
end

fprintf('done! '); toc
fprintf('\n');

return





