function [fix_im, info, p] = vcd_fixationDot(p)
% VCD function:
%  fixIm = vcd_fixationDot(p)
%
% Purpose:
%   Create a set of fixation dot images for experimental display.
%   See vcd_setStimParams.m for fixation parameters.
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
%   info    : table with fix image information
%   p       : updated params struct
%
% Written by Eline Kupers 2025/02

% we want a dot that will change between luminance between 5 (or more?) levels
% up or down:

% 2*fixation diam x 2*fixation diam x 3 x 5 luminance levels x 2 dot rims
fix_im   = zeros([2*p.stim.fix.dotcenterdiam_pix, 2*p.stim.fix.dotcenterdiam_pix, 3, length(p.stim.fix.dotlum), 2]); 

% Where to insert luminance val
fixationmask_inner    = find(makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotcenterdiam_pix/2));

fixationmask_rimthin    = find( makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotthinborderdiam_pix/2) - ...
                           makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotcenterdiam_pix/2));  

fixationmask_rimthick   = find( makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotthickborderdiam_pix/2) - ...
                           makecircleimage(2*p.stim.fix.dotcenterdiam_pix, p.stim.fix.dotcenterdiam_pix/2));  

for lum = 1:length(p.stim.fix.dotlum)
    temp0 = zeros(2*p.stim.fix.dotcenterdiam_pix*2*p.stim.fix.dotcenterdiam_pix, 3);        % everything is initially black
    fixation_rimthin = temp0; fixation_rimthick = temp0;

    fixation_rimthin(fixationmask_rimthin,:) = repmat(double(255),[size(fixationmask_rimthin,1) 3]);  % add rim
    fixation_rimthick(fixationmask_rimthick,:) = repmat(double(255),[size(fixationmask_rimthick,1) 3]);  % add rim
%     temp0 = repmat(double(3),[size(temp0,1) 1]);                 % but we fill it with the specified color

    fixation_rimthin(fixationmask_inner,:)  = repmat(double(p.stim.fix.dotlum(lum)),[length(fixationmask_inner) 3]);               % but we fill it with the specified color
    fixation_rimthick(fixationmask_inner,:)  = repmat(double(p.stim.fix.dotlum(lum)),[length(fixationmask_inner) 3]);               % but we fill it with the specified color



    fix_im(:,:,:,lum,1) = reshape(fixation_rimthin,[2*p.stim.fix.dotcenterdiam_pix, 2*p.stim.fix.dotcenterdiam_pix, 3]);
    fix_im(:,:,:,lum,2) = reshape(fixation_rimthick,[2*p.stim.fix.dotcenterdiam_pix, 2*p.stim.fix.dotcenterdiam_pix, 3]);

    clear temp0 
end

% fixationalphamask_thin = makecircleimage(2*p.stim.fix.dotcenterdiam_pix,p.stim.fix.dotthinborderdiam_pix/2);  % 2*fixationsize x 2*fixationsize; double [0,255] alpha values (50% of 255 in circle, 0 outside)
% fixationalphamask_thick = makecircleimage(2*p.stim.fix.dotcenterdiam_pix,p.stim.fix.dotthickborderdiam_pix/2);  % 2*fixationsize x 2*fixationsize; double [0,255] alpha values (50% of 255 in circle, 0 outside)
% fixationimage(:,:,:,:,1) = bsxfun(@times, fixationimage(:,:,:,:,1), fixationalphamask_thin);
% fixationimage(:,:,:,:,2) = bsxfun(@times, fixationimage(:,:,:,:,2), fixationalphamask_thick);


% create info table
lum_info = repmat(p.stim.fix.dotlum',2,1);
width_info = repelem({'thin','thick'},length(p.stim.fix.dotlum));
info = table(lum_info,width_info','VariableNames',{'luminance','rim_width'});

if p.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(p.stim.fix.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',p.stim.fix.stimfile,datestr(now,30))),'fix_im','info','-v7.3');

    saveDir = fileparts(fullfile(p.stim.fix.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info,fullfile(sprintf('%s_%s.csv',p.stim.fix.infofile,datestr(now,30))));


end


return