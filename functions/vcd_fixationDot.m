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
%   info    : table with fix image information
%   p       : updated params struct
%
% Written by Eline Kupers 2025/02


% Create support: 2*fixation diam x 2*fixation diam x 3 x 5 luminance levels x 2 dot rims
fix_im   = zeros([2*p.stim.fix.dotcenterdiam_pix, 2*p.stim.fix.dotcenterdiam_pix, 3, length(p.stim.fix.dotlum), 2]); 

x = [1:(2*p.stim.fix.dotcenterdiam_pix)]-p.stim.fix.dotcenterdiam_pix;
[XX,~] = meshgrid(x,x);

% Where to insert luminance val
fixationmask_inner    = find(makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotcenterdiam_pix/4));
fixationmask_rimthin  = find( makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotthinborderdiam_pix/4) - ...
                           makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotcenterdiam_pix/4));  
fixationmask_rimthick = find( makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotthickborderdiam_pix/4) - ...
                           makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotcenterdiam_pix/4));  

% Find left & right half moons
left_idx = find((XX <= 0));
right_idx = find((XX > 0));
[~,li] = intersect(fixationmask_rimthick,left_idx);
[~,ri] = intersect(fixationmask_rimthick,right_idx);

fixationmask_rimthick_left = fixationmask_rimthick(li);
fixationmask_rimthick_right = fixationmask_rimthick(ri);

                       
                       
for lum = 1:length(p.stim.fix.dotlum)
    temp0 = double(p.stim.bckgrnd_grayval)*ones(2*p.stim.fix.dotcenterdiam_pix*2*p.stim.fix.dotcenterdiam_pix, 1);        % everything is initially black
    fixation_rimthin = temp0; 
    fixation_rimthick = temp0;

    % full rim
    fixation_rimthin(fixationmask_rimthin,:)   = repmat(double(255),[size(fixationmask_rimthin,1) 1]);  % add rim
    fixation_rimthick(fixationmask_rimthick,:) = repmat(double(255),[size(fixationmask_rimthick,1) 1]);  % add rim
    
    % half black rim
    fixation_rimthick_left = fixation_rimthick;
    fixation_rimthick_left(fixationmask_rimthick_left,:) = repmat(double(1),[size(fixationmask_rimthick_left,1) 1]);  % add left black rim
    
    fixation_rimthick_right = fixation_rimthick;
    fixation_rimthick_right(fixationmask_rimthick_right,:) = repmat(double(1),[size(fixationmask_rimthick_right,1) 1]);  % add left black rim
    
    % but we fill it with the specified color
    fixation_rimthin(fixationmask_inner,:)     = repmat(double(p.stim.fix.dotlum(lum)),[length(fixationmask_inner) 1]);
    fixation_rimthick(fixationmask_inner,:)    = repmat(double(p.stim.fix.dotlum(lum)),[length(fixationmask_inner) 1]);
    
    fixation_rimthick_right(fixationmask_inner,:)     = repmat(double(p.stim.fix.dotlum(lum)),[length(fixationmask_inner) 1]);
    fixation_rimthick_left(fixationmask_inner,:)    = repmat(double(p.stim.fix.dotlum(lum)),[length(fixationmask_inner) 1]);
    
    % repmat to get RGB, reshape back to square image
    fix_im(:,:,:,lum,1) = repmat(reshape(fixation_rimthin,[2*p.stim.fix.dotcenterdiam_pix, 2*p.stim.fix.dotcenterdiam_pix, 1]),[1 1 3]);
    fix_im(:,:,:,lum,2) = repmat(reshape(fixation_rimthick,      [2*p.stim.fix.dotcenterdiam_pix, 2*p.stim.fix.dotcenterdiam_pix, 1]),[1 1 3]);
    fix_im(:,:,:,lum,3) = repmat(reshape(fixation_rimthick_left, [2*p.stim.fix.dotcenterdiam_pix, 2*p.stim.fix.dotcenterdiam_pix, 1]),[1 1 3]);
    fix_im(:,:,:,lum,4) = repmat(reshape(fixation_rimthick_right,[2*p.stim.fix.dotcenterdiam_pix, 2*p.stim.fix.dotcenterdiam_pix, 1]),[1 1 3]);
    
    clear temp0 
end

% create alpha mask
alpha_mask0   = zeros(2*p.stim.fix.dotcenterdiam_pix, 2*p.stim.fix.dotcenterdiam_pix);
alpha_mask0   = alpha_mask0(:);
alpha_idx     = find(makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotthickborderdiam_pix/2));

alpha_mask0(alpha_idx) = ones(length(alpha_idx),1).*255;

alpha_mask_rz = reshape(alpha_mask0,2*p.stim.fix.dotcenterdiam_pix,2*p.stim.fix.dotcenterdiam_pix);
mask = uint8(alpha_mask_rz);

% figure; 
% for ii = 1:4
%     subplot(2,2,ii);
%     imagesc(uint8(squeeze(fix_im(:,:,:,5,ii))))
%     set(gca,'CLim',[0 255])
%     colormap gray
% end

% create info table
lum_info = repmat(p.stim.fix.dotlum',4,1);
width_info = repelem({'thick','thin','thin_left','thin_right',},length(p.stim.fix.dotlum));
info = table(lum_info,width_info','VariableNames',{'luminance','rim_width'});

if p.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(p.stim.fix.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',p.stim.fix.stimfile,datestr(now,30))),'fix_im','mask','info','-v7.3');

    saveDir = fileparts(fullfile(p.stim.fix.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info,fullfile(sprintf('%s_%s.csv',p.stim.fix.infofile,datestr(now,30))));


end


return