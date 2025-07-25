function [fix_im, mask, info] = vcd_fixationDot(params, verbose, store_imgs)
% VCD function:
%  [fix_im, info] = vcd_fixationDot(params, verbose, store_imgs)
%
% Purpose:
%   Create a set of fixation dot images for experimental display. We want a
%   dot that will change between luminance between 5 levels up 
%   or down + 6th level: mean luminance. See vcd_setStimParams.m for 
%   fixation parameters.
%
% INPUTS:
% params  : params stuct (see vcd_setStimParams.m)
%               * stim.store_imgs
%               * stim.fix.dotcenterdiam_pix
%               * stim.fix.dotthickcenterdiam_pix
%               * stim.fix.dotthincenterdiam_pix
%               * stim.fix.dotlum
% verbose           : (logical) show debug figures
% store_imgs        : (logical) store stimuli and debug figures as pngs 
%
% OUTPUTS:
% fix_im  : fixation dot images, 5-dim array:
%               w (pixels) x h (pixels) x 3 x 6 luminance levels x 6 rim widths
% mask    : alpha mask for ptb: w (pixels) x h (pixels) x 2 (one gray
%                   layer, one transparency layer)
% info    : table with information about fixation dot conditions 
%
% Written by Eline Kupers 2025/02


%% Create support width dims: 2*fixation diam x 2*fixation diam x 3 x 5 luminance levels x 5 dot rims
support_x = 2*params.stim.fix.dotcenterdiam_pix;
support_y = support_x; 

fix_im  = uint8(zeros([support_x, support_y, 3, length(params.stim.fix.dotlum), 6])); 

x = [1:(2*params.stim.fix.dotcenterdiam_pix)]-params.stim.fix.dotcenterdiam_pix;
[XX,~] = meshgrid(x,x);

% Where to insert luminance val? (divide diam by 2 to get radius, which is
% expected by makecircleimage)
fixationmask_inner    = find(makecircleimage(support_x, params.stim.fix.dotcenterdiam_pix/2));
fixationmask_rimthin  = find(makecircleimage(support_x, params.stim.fix.dotthinborderdiam_pix/2) - ...
                           makecircleimage(support_x, params.stim.fix.dotcenterdiam_pix/2));  
fixationmask_rimthick = find(makecircleimage(support_x, params.stim.fix.dotthickborderdiam_pix/2) - ...
                           makecircleimage(support_x, params.stim.fix.dotcenterdiam_pix/2));  

% Find left & right half moons (for spatial cue)
left_idx = find((XX <= 0));
right_idx = find((XX > 0));
[~,li] = intersect(fixationmask_rimthick,left_idx);
[~,ri] = intersect(fixationmask_rimthick,right_idx);

fixationmask_rimthick_left = fixationmask_rimthick(li);
fixationmask_rimthick_right = fixationmask_rimthick(ri);
              
% Loop over center dot luminance values
for lum = 1:length(params.stim.fix.dotlum)
    
    % everything is initially gray
    temp0 = params.stim.bckgrnd_grayval*ones(support_x*support_y, 3);        
    fixation_rimthin  = temp0; 
    fixation_rimthick_white = temp0;
    fixation_rimthick_black = temp0;
    
    % full thin and thick rim base
    fixation_rimthin(fixationmask_rimthin,:)   = params.stim.fix.color(1,:) .* ones(size(fixationmask_rimthin,1), 3);  % add white rim
    fixation_rimthick_white(fixationmask_rimthick,:) = params.stim.fix.color(1,:) .* ones(size(fixationmask_rimthick,1), 3);  % add white rim
    fixation_rimthick_black(fixationmask_rimthick,:) = params.stim.fix.color(3,:) .* ones(size(fixationmask_rimthick,1), 3);  % add white rim

    % left half red rim base
    fixation_rimthick_left = fixation_rimthick_white;
    fixation_rimthick_left(fixationmask_rimthick_left,:) = params.stim.fix.color(2,:) .* ones(size(fixationmask_rimthick_left,1), 3); 
    
    % right half red rim base
    fixation_rimthick_right = fixation_rimthick_white;
    fixation_rimthick_right(fixationmask_rimthick_right,:) = params.stim.fix.color(2,:) .* ones(size(fixationmask_rimthick_right,1),3);  
    
    % both left&right sides black rim base (for ns)
    fixation_rimthick_both = fixation_rimthick_white;
    fixation_rimthick_both([fixationmask_rimthick_right;fixationmask_rimthick_left],:) = ...
                        params.stim.fix.color(2,:) .* ones(size(fixationmask_rimthick_right,1)+size(fixationmask_rimthick_left,1), 3);  

    % Fill inner circle with the specified luminance 
    fixation_rimthin(fixationmask_inner,:)          = repmat(params.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    fixation_rimthick_white(fixationmask_inner,:)   = repmat(params.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    fixation_rimthick_black(fixationmask_inner,:)   = repmat(params.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    fixation_rimthick_right(fixationmask_inner,:)   = repmat(params.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    fixation_rimthick_left(fixationmask_inner,:)    = repmat(params.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    fixation_rimthick_both(fixationmask_inner,:)    = repmat(params.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    
    % Repmat to get RGB, reshape back to square image
    fixation_rimthin1        = reshape(fixation_rimthin, [support_x, support_y, 3]);
    fixation_rimthick_white1 = reshape(fixation_rimthick_white,[support_x, support_y, 3]);
    fixation_rimthick_black1 = reshape(fixation_rimthick_black,[support_x, support_y, 3]);
    fixation_rimthick_left1  = reshape(fixation_rimthick_left, [support_x, support_y, 3]);
    fixation_rimthick_right1 = reshape(fixation_rimthick_right,[support_x, support_y, 3]);
    fixation_rimthick_both1  = reshape(fixation_rimthick_both, [support_x, support_y, 3]);
    
    % Convert to uint8 images
    fix_im(:,:,:,lum,1) = uint8(fixation_rimthin1);
    fix_im(:,:,:,lum,2) = uint8(fixation_rimthick_white1); % white
    fix_im(:,:,:,lum,3) = uint8(fixation_rimthick_left1);
    fix_im(:,:,:,lum,4) = uint8(fixation_rimthick_right1);
    fix_im(:,:,:,lum,5) = uint8(fixation_rimthick_both1);
    fix_im(:,:,:,lum,6) = uint8(fixation_rimthick_black1); % black
    
    % clean up
    clear temp0 
end

%% Create alpha mask

% Create support to be the same size as dot image
alpha_mask0   = zeros(support_x, support_y); % initially everything is invisible

% flatten image
alpha_mask0_thin   = alpha_mask0(:);
alpha_mask0_thick   = alpha_mask0(:);

% create a circle mask (divide by 4 because makecircleimage expects radius,
% number is used to pad on both left and right side)
alpha_idx_thick     = find(makecircleimage(support_x, (params.stim.fix.dotthickborderdiam_pix/2)));  % add 2 pixels on either size to avoid the rim being cut
alpha_idx_thin      = find(makecircleimage(support_x, (params.stim.fix.dotthinborderdiam_pix/2))); % add 2 pixels on either size to avoid the rim being cut

% everything inside the alpha mask is 50% opacity, outside mask is invisible
alpha_mask0_thin(alpha_idx_thin)   = ones(length(alpha_idx_thin),1)*255*params.stim.fix.dotopacity;
alpha_mask0_thick(alpha_idx_thick) = ones(length(alpha_idx_thick),1)*255*params.stim.fix.dotopacity;

% resize alpha mask 
alpha_mask_thin = reshape(alpha_mask0_thin,support_x,support_y);
alpha_mask_thick = reshape(alpha_mask0_thick,support_x,support_y);

% convert to uint8 and concat
alpha_mask_thin  = uint8(alpha_mask_thin);
alpha_mask_thick = uint8(alpha_mask_thick);
mask = cat(3, alpha_mask_thin, repmat(alpha_mask_thick,[1 1 5]));

%% create info table
lum_info        = repmat(params.stim.fix.dotlum',6,1);
width_info      = repelem({'thin','thick','thick','thick','thick','thick'},length(params.stim.fix.dotlum));
color_side_info = repelem({'both','both','left','right','both','both'},length(params.stim.fix.dotlum));
rim_color       = cat(1, repmat(params.stim.fix.color(1,:),2*length(params.stim.fix.dotlum),1),...
                    repmat(params.stim.fix.color(2,:),3*length(params.stim.fix.dotlum),1),...
                    repmat(params.stim.fix.color(3,:),length(params.stim.fix.dotlum),1));
rim_color  = mat2cell(rim_color,ones(size(rim_color,1),1));
info       = table(lum_info,width_info',color_side_info',rim_color,'VariableNames',{'luminance','rim_width','color_side','rim_color'});

%% Store images and info table if requested
if store_imgs
    fprintf('[%s]: Storing images..',mfilename)
    tic
    saveMatStimFileDir = fileparts(fullfile(params.stim.fix.stimfile));
    if ~exist(saveMatStimFileDir,'dir'), mkdir(saveMatStimFileDir); end
    save(fullfile(sprintf('%s_%s.mat',params.stim.fix.stimfile,datestr(now,30))),'fix_im','mask','info','-v7.3');

    saveMatStimFileDir = fileparts(fullfile(params.stim.fix.infofile));
    if ~exist(saveMatStimFileDir,'dir'), mkdir(saveMatStimFileDir); end
    writetable(info,fullfile(sprintf('%s_%s.csv',params.stim.fix.infofile,datestr(now,30))));
    fprintf('done! '); toc
end

%% Visualize fixation circle if requested
if verbose
    makeprettyfigures;
    
    saveFigDir = fullfile(params.saveFigsFolder,'fix');
    if ~exist(saveFigDir,'dir'); mkdir(saveFigDir); end
    
    fH = figure(101); clf;
    set(fH, 'Position', [1   400   750   578], 'color','w')
    counter = 1;
    for ii = 1:size(fix_im,4)
        for jj = 1:size(fix_im,5)
            subplot(6,6,counter)
            imagesc(squeeze(fix_im(:,:,:,ii,jj)), 'AlphaData',params.stim.fix.dotopacity)
            axis square
            axis off
            set(gca,'CLim',[1 255])
            title(sprintf('%01d, %01d',ii, jj))
            counter = counter+1;
            if store_imgs
                imwrite(fix_im(:,:,:,ii,jj), fullfile(saveFigDir,sprintf('vcd_fixcircle_w_alphamask_%01d_%01d.png', ii, jj)), 'Alpha' , double(mask(:,:,jj))/255);
            end
        end
    end
    if store_imgs
        filename = sprintf('vcd_fixcircle_all_w_alphamask.png');
        print(fH,'-dpng','-r300',fullfile(saveFigDir,filename));
    end
    
    fH = figure(101); clf;
    set(fH, 'Position', [1   400   750   578], 'color','w')
    counter = 1;
    for ii = 1:size(fix_im,4)
        for jj = 1:size(fix_im,5)
            subplot(6,6,counter)
            imagesc(squeeze(fix_im(:,:,:,ii,jj)))
            axis square
            axis off
            set(gca,'CLim',[1 255])
            title(sprintf('%01d, %01d',ii, jj))
            counter = counter+1;
            if store_imgs
                imwrite(fix_im(:,:,:,ii,jj), fullfile(saveFigDir,sprintf('vcd_fixcircle_%01d_%01d.png', ii, jj)));
            end
        end
    end
    if store_imgs
        filename = sprintf('vcd_fixcircle_all.png');
        print(fH,'-dpng','-r300',fullfile(saveFigDir,filename));
    end
end

return


