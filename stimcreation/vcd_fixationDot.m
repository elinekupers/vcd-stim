function [fix_im, mask, info, p] = vcd_fixationDot(p)
% VCD function:
%  [fix_im, info, p] = vcd_fixationDot(p)
%
% Purpose:
%   Create a set of fixation dot images for experimental display. We want a
%   dot that will change between luminance between 5 (or more?) levels up 
%   or down. See vcd_setStimParams.m for fixation parameters.
%
% INPUTS:
%   p       : params stuct (see vcd_setStimParams.m)
%               * stim.store_imgs
%               * stim.fix.dotcenterdiam_pix
%               * stim.fix.dotthickcenterdiam_pix
%               * stim.fix.dotthincenterdiam_pix
%               * stim.fix.dotlum
%
% OUTPUTS:
%   fix_im  : fixation dot images, 5-dim array:
%               w (pixels) x h (pixels) x 3 x 5 luminance levels x 2 rim widths
%   mask    : alpha mask for ptb: w (pixels) x h (pixels) x 2 (one gray
%                   layer, one transparency layer)
%   info    : table with information about fixation dot conditions 
%   p       : updated params struct
%
% Written by Eline Kupers 2025/02


%% Create support width dims: 2*fixation diam x 2*fixation diam x 3 x 5 luminance levels x 5 dot rims
support_x = 2*p.stim.fix.dotcenterdiam_pix;
support_y = support_x; 

fix_im  = uint8(zeros([support_x, support_y, 3, length(p.stim.fix.dotlum), 5])); 

x = [1:(2*p.stim.fix.dotcenterdiam_pix)]-p.stim.fix.dotcenterdiam_pix;
[XX,~] = meshgrid(x,x);

% Where to insert luminance val? (divide diam by 2 to get radius, which is
% expected by makecircleimage)
fixationmask_inner    = find(makecircleimage(support_x, p.stim.fix.dotcenterdiam_pix/2));
fixationmask_rimthin  = find(makecircleimage(support_x, p.stim.fix.dotthinborderdiam_pix/2) - ...
                           makecircleimage(support_x, p.stim.fix.dotcenterdiam_pix/2));  
fixationmask_rimthick = find(makecircleimage(support_x, p.stim.fix.dotthickborderdiam_pix/2) - ...
                           makecircleimage(support_x, p.stim.fix.dotcenterdiam_pix/2));  

% Find left & right half moons (for spatial cue)
left_idx = find((XX <= 0));
right_idx = find((XX > 0));
[~,li] = intersect(fixationmask_rimthick,left_idx);
[~,ri] = intersect(fixationmask_rimthick,right_idx);

fixationmask_rimthick_left = fixationmask_rimthick(li);
fixationmask_rimthick_right = fixationmask_rimthick(ri);
              
% Loop over center dot luminance values
for lum = 1:length(p.stim.fix.dotlum)
    
    % everything is initially gray
    temp0 = p.stim.bckgrnd_grayval*ones(support_x*support_y, 3);        
    fixation_rimthin  = temp0; 
    fixation_rimthick = temp0;

    % full thin and thick rim base
    fixation_rimthin(fixationmask_rimthin,:)   = p.stim.fix.color(1,:) .* ones(size(fixationmask_rimthin,1), 3);  % add white rim
    fixation_rimthick(fixationmask_rimthick,:) = p.stim.fix.color(1,:) .* ones(size(fixationmask_rimthick,1), 3);  % add white rim
    
    % left half red rim base
    fixation_rimthick_left = fixation_rimthick;
    fixation_rimthick_left(fixationmask_rimthick_left,:) = p.stim.fix.color(2,:) .* ones(size(fixationmask_rimthick_left,1), 3); 
    
    % right half red rim base
    fixation_rimthick_right = fixation_rimthick;
    fixation_rimthick_right(fixationmask_rimthick_right,:) = p.stim.fix.color(2,:) .* ones(size(fixationmask_rimthick_right,1),3);  
    
    % both left&right sides black rim base (for ns)
    fixation_rimthick_both = fixation_rimthick;
    fixation_rimthick_both([fixationmask_rimthick_right;fixationmask_rimthick_left],:) = ...
                        p.stim.fix.color(2,:) .* ones(size(fixationmask_rimthick_right,1)+size(fixationmask_rimthick_left,1), 3);  

    % Fill inner circle with the specified luminance 
    fixation_rimthin(fixationmask_inner,:)          = repmat(p.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    fixation_rimthick(fixationmask_inner,:)         = repmat(p.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    fixation_rimthick_right(fixationmask_inner,:)   = repmat(p.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    fixation_rimthick_left(fixationmask_inner,:)    = repmat(p.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    fixation_rimthick_both(fixationmask_inner,:)    = repmat(p.stim.fix.dotlum(lum),[length(fixationmask_inner) 3]);
    
    % Repmat to get RGB, reshape back to square image
    fixation_rimthin1        = reshape(fixation_rimthin, [support_x, support_y, 3]);
    fixation_rimthick1       = reshape(fixation_rimthick,[support_x, support_y, 3]);
    fixation_rimthick_left1  = reshape(fixation_rimthick_left, [support_x, support_y, 3]);
    fixation_rimthick_right1 = reshape(fixation_rimthick_right,[support_x, support_y, 3]);
    fixation_rimthick_both1  = reshape(fixation_rimthick_both, [support_x, support_y, 3]);
    
    % Convert to uint8 images
    fix_im(:,:,:,lum,1) = uint8(fixation_rimthin1);
    fix_im(:,:,:,lum,2) = uint8(fixation_rimthick1);
    fix_im(:,:,:,lum,3) = uint8(fixation_rimthick_left1);
    fix_im(:,:,:,lum,4) = uint8(fixation_rimthick_right1);
    fix_im(:,:,:,lum,5) = uint8(fixation_rimthick_both1);
    
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
alpha_idx_thick     = find(makecircleimage(support_x, p.stim.fix.dotthickborderdiam_pix/4)); 
alpha_idx_thin      = find(makecircleimage(support_x, p.stim.fix.dotthinborderdiam_pix/4)); 

% everything inside the alpha mask is 50% opacity, outside mask is invisible
alpha_mask0_thin(alpha_idx_thin)   = ones(length(alpha_idx_thin),1)*255*p.stim.fix.dotopacity;
alpha_mask0_thick(alpha_idx_thick) = ones(length(alpha_idx_thick),1)*255*p.stim.fix.dotopacity;

% resize alpha mask 
alpha_mask_thin = reshape(alpha_mask0_thin,support_x,support_y);
alpha_mask_thick = reshape(alpha_mask0_thick,support_x,support_y);

% convert to uint8 and concat
alpha_mask_thin  = uint8(alpha_mask_thin);
alpha_mask_thick = uint8(alpha_mask_thick);
mask = cat(3, alpha_mask_thin, repmat(alpha_mask_thick,[1 1 4]));

%% create info table
lum_info   = repmat(p.stim.fix.dotlum',5,1);
width_info = repelem({'thin','thick','thick_left','thick_right','thick_both'},length(p.stim.fix.dotlum));
info       = table(lum_info,width_info','VariableNames',{'luminance','rim_width'});

%% Store images and info table if requested
if p.stim.store_imgs
    fprintf('[%s]: Storing images..',mfilename)
    tic
    saveDir = fileparts(fullfile(p.stim.fix.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s_%s.mat',p.stim.fix.stimfile,p.disp.name,datestr(now,30))),'fix_im','mask','info','-v7.3');

    saveDir = fileparts(fullfile(p.stim.fix.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info,fullfile(sprintf('%s_%s.csv',p.stim.fix.infofile,datestr(now,30))));
    fprintf('done! '); toc
end

%% Visualize fixation circle if requested
if params.verbose
    makeprettyfigures;
    
    saveFigDir = fullfile(vcd_rootPath,'figs',params.disp.name,'fix','visual_check');
    if ~exist(saveFigDir,'dir'); mkdir(saveFigDir); end
    
    fH = figure(101); clf;
    set(fH, 'Position', [1   400   750   578], 'color','w')
    counter = 1;
    for ii = 1:size(fix_im,4)
        for jj = 1:size(fix_im,5)
            subplot(5,5,counter)
            imagesc(squeeze(fix_im(:,:,:,ii,jj)), 'AlphaData',params.stim.fix.dotopacity)
            axis square
            axis off
            set(gca,'CLim',[1 255])
            title(sprintf('%01d, %01d',ii, jj))
            counter = counter+1;
            if params.store_imgs
                imwrite(fix_im(:,:,:,ii,jj), fullfile(vcd_rootPath,'figs',params.disp.name,'fix',sprintf('vcd_fixdots_w_alphamask_%01d_%01d.png', ii, jj)));
            end
        end
    end
    if params.store_imgs
        filename = sprintf('vcd_fixdots_w_alphamask.png');
        print(fH,'-dpng','-r300',fullfile(saveFigDir,filename));
    end
    
    fH = figure(101); clf;
    set(fH, 'Position', [1   400   750   578], 'color','w')
    counter = 1;
    for ii = 1:size(fix_im,4)
        for jj = 1:size(fix_im,5)
            subplot(5,5,counter)
            imagesc(squeeze(fix_im(:,:,:,ii,jj)))
            axis square
            axis off
            set(gca,'CLim',[1 255])
            title(sprintf('%01d, %01d',ii, jj))
            counter = counter+1;
            if params.store_img
                imwrite(fix_im(:,:,:,ii,jj), fullfile(vcd_rootPath,'figs',params.disp.name,'fix',sprintf('vcd_fixdots_%01d_%01d.png', ii, jj)));
            end
        end
    end
    if params.store_imgs
        filename = sprintf('vcd_fixdots.png');
        print(fH,'-dpng','-r300',fullfile(saveFigDir,filename));
    end
end

return


