function [run_images, images] = vcd_loadRunImages(run_image_order, block, params)

%% %%%%%%%%%%%%% PRE-LOAD STIMULI %%%%%%%%%%%%%
if isfield(params, 'images')
    images = params.images;  params = rmfield(params, 'images');
end

if ~exist('images','var') || isempty(fieldnames(images))
    images = struct('gabor',[],'rdk',[],'dot',[],'cobj',[],'ns',[],...
        'fix',[], 'info',[], 'image_order',[]);
end

fprintf('[%s]: Loading images..',mfilename); tic;
for ii = 1:length(run_image_order)
    statusdots(ii,length(block));
    
    st_crossing = block(ii).name;
    tmp = strsplit(st_crossing,'-');
    taskClass = tmp{1};
    stimClass = tmp{2};
    
    if ~any(strcmp(taskClass,{'pre','post'}))
        % check we assigned the string to the right label
        assert(any(strcmp(taskClass,params.exp.taskClassLabels)));
        assert(any(strcmp(stimClass,params.exp.stimClassLabels)));
    end
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
                load(fullfile(d(end).folder,d(end).name), 'gabors','info');
                
                images.gabor = gabors; clear gabors;
                images.info.gabor = info; clear info;
            end
            
            for jj = 1:numTrials
                for nn = 1:numSides
                    i3 = run_image_order{ii}{jj,1}(nn);
                    run_images{ii,jj,nn,1} = images.gabor(:,:,i3,1);
                    
                    if strcmp(taskClass,'wm')
                        i4 = run_image_order{ii}{jj,2}(nn);
                        run_images{ii,jj,nn,2} = images.gabor(:,:,i3,i4+1);
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
                    
                    stimfile = fullfile(stimDir.folder,stimDir.name,sprintf('%s.mat', filename));
                    if exist(stimfile,'file')
                        load(stimfile, 'frames');
                    else
                        error('[%s]: Can''t find RDK stim file!')
                    end
                    % each file contains a 4D array: [x,y,3,frames]
                    run_images{ii,jj,nn,1} = frames;
                   
                    if strcmp(taskClass,'wm')
                        filename_delta = sprintf('%s_rdk_ori%d_coh%d_delta%d',...
                            block(ii).trial(jj,nn).unique_im_nr, ...
                            block(ii).trial(jj,nn).motdir, ...
                            block(ii).trial(jj,nn).coh, ...
                            run_image_order{ii}{jj,2}(nn));
                        
                        stimfile_delta = fullfile(stimDir,sprintf('%s.mat', filename_delta));
                        
                        load(stimfile_delta, 'frames');
                        run_images{ii,jj,nn,2} = frames;
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
                load(fullfile(d(end).folder,d(end).name), 'simple_dot','info');
                images.dot = simple_dot; clear simple_dot;
                images.info.dot = info; clear info;
            end
            for jj = 1:numTrials
                for nn = 1:numSides
                    
                    run_images{ii,jj,nn,1} = images.dot;
                    if strcmp(taskClass,'wm')
                        
                        run_images{ii,jj,nn,1} = images.dot;
                    end
                    %                 if strcmp(taskClass,'ltm')
                    %
                    %                     run_images{ii,jj,1,1} = images.ns(:,:,:,i4,i5,i6);
                    %                 end
                end
            end
            
        case 'cobj'
            sz = size(block(ii).trial);
            if sz(1)>sz(2) && any(sz~=1)
                numTrials = sz(1);
                numSides = sz(2);
            elseif sz(1)<sz(2) && any(sz~=1)
                numTrials = sz(2);
                numSides = sz(1);
            end
            
            if isempty(images.cobj)
                % Complex objects: 4D array: [x,y,16 object, og + 10 rotation]
                d = dir(sprintf('%s*.mat', params.stim.cobj.stimfile));
                load(fullfile(d(end).folder,d(end).name), 'objects','info','im_order');
                images.cobj = objects; clear objects;
                images.info.cobj = info; clear info;
                images.image_order.cobj = im_order; clear im_order;
            end
            
            for jj = 1:numTrials
                for nn = 1:numSides
                    [i3,i4] = ind2sub([size(images.cobj,3),size(images.cobj,4)],run_image_order{ii}{jj,1}(nn));
                    
                    run_images{ii,jj,nn,1} = images.cobj(:,:,i3,i4);
                    
                    if strcmp(taskClass,'wm')
                        [i3,i4] = ind2sub([size(images.cobj,3),size(images.cobj,4)],run_image_order{ii}{jj,2}(nn));
                        run_images{ii,jj,nn,2} = images.cobj(:,:,i3,i4);
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
                numTrials = sz(2);
            elseif sz(1)<sz(2) && any(sz~=1)
                numTrials = sz(1);
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
                i4 = images.info.ns.superordinate_i(run_image_order{ii}{jj,1},:);
                i5 = images.info.ns.ns_loc_i(run_image_order{ii}{jj,1},:);
                i6 = images.info.ns.obj_loc_i(run_image_order{ii}{jj,1},:);
                
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
    
    % If using BOLDSCREEN, we want to square pixel values
    if strcmp(params.disp.name,'7TAS_BOLDSCREEN32') && ...
            ~any(strcmp(taskClass,{'pre','post'}))
        
        fprintf('square pix values for CLUT..',mfilename);
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
                        img = img.^2;
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


           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% CHECK CLUT & COLOR %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            % WE ALREADY RESIZED IMAGES!!!
            
%             if ~any(strcmp(taskClass,{'pre','post'})) && ...
%                     ~isempty(params.stim.(stimClass).dres) && ...
%                     sum(params.stim.(stimClass).dres)~=0
%                 
%                 fprintf('[%]: Resampling..\n',mfilename);
%                 
%                 % reshape to get one vector of images
%                 run_im_tmp0 = run_images{ii,jj,:,:};
%                 sz0 = size(run_im_tmp0);
%                 run_im_tmp = reshape(run_im_tmp0,1, []);
%                 
%                 
%                 for mm = 1:length(run_im_tmp)
%                     statusdots(ii,length(run_im_tmp));
%                     if ~isempty(run_im_tmp{mm})
%                         
%                         
%                         % GRAY SCALE has 4 dim
%                         if ~params.stim.(stimClass).iscolor
%                             tmp_im = reshape(run_im_tmp{mm}, ...
%                                 size(run_im_tmp{mm},1),... x
%                                 size(run_im_tmp{mm},2),... y
%                                 []);
%                             
%                             % COLOR has 3 dim
%                         elseif params.stim.(stimClass).iscolor
%                             tmp_im = reshape(run_im_tmp{mm}, ...
%                                 size(run_im_tmp{mm},1),... x
%                                 size(run_im_tmp{mm},2),... y
%                                 3, ... % rgb
%                                 []);
%                         end
%                         
%                         numIn = size(tmp_im);
%                         numIn = numIn(end);
%                         
%                         for kk = 1:numIn
%                             
%                             if params.stim.(stimClass).iscolor
%                                 im_rz(:,:,:,kk) = imresize(tmp_im(:,:,:,kk),params.stim.(stimClass).dres); %% DEFAULT IS BICUBIC
%                             else
%                                 im_rz(:,:,kk) = imresize(tmp_im(:,:,kk),params.stim.(stimClass).dres);
%                             end
%                         end
%                         run_im_tmp{mm} = im_rz;
%                     end
%                     
%                     
%                 end
%                 run_im_tmp1 = reshape(run_im_tmp,sz0(1),sz0(2));
%                 run_images{ii,jj,:,:} = run_im_tmp1;
%             end
%             

