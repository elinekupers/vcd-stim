function bckgrnd_im = vcd_pinknoisebackground(params, varargin)
% VCD function:
%   bckgrnd_im = vcd_pinknoisebackground(params, [type,] [borderwidth,] [num,] [pixoffset])
% 
% Purpose:
%   Create a full-field, pink noise background with mean luminance blank 
%   gap for stimuli in the center. This function is modular in terms of the
%   monitor size.
%   
% INPUTS:
%   params       : struct with params, should contain the following fields:
%                   - params.disp : monitor display params (struct) with 
%                                  pixels per degree in field of view (pix) 
%                   - params.disp.w_pix  : width of screen resolution (pix)
%                   - params.disp.h_pix  : height of screen resolution (pix)
%                   - params.disp.xc_pix : x-center of screen (pix)
%                   - params.disp.yc_pix : y-center of screen (pix)
%                   - params.disp.ppd    : pixels per degree (num)
%                   - params.stim.bckground : (struct) background params 
%                   - params.stim.bckground.stimfile : (str) where to store images
%                   - params.stim.store_imgs : (bool) store images or not?
%   gaptype      : (str) what type of gap do you want. Choose from:
%                   - puzzle: center square overlayed with 2 parafoveal circular patches.
%                   - dotring: simple dot iso-eccentricity ring.
%                   - comb: combination of puzzle piece and dotring. 
%                   - circle: a 4.2 or 5.2 deg radius circle (default)
%   borderwidth  : (str) how wide should the gap between stimuli and gap be? 
%                     Choose from:
%                   - skinny: aligned to stimulus extend
%                   - fat: stimulus extend + 2 deg
%   num          : (int) number of generated noise images
%   pixoffset 	 : (int) offset of [x,y]-center in pixels. Default is [0,0] 
%   
% OUTPUTS:
%   bckgrnd_im   : background images, array:  w (pixels) x h (pixels) x 3 x num
%
%
% Written by Eline Kupers 2024/12, updated 2025/03

%% Parse inputs
p = inputParser;
p.addRequired('params'          , @isstruct); % params struct
p.addParameter('gaptype'        , 'comb' , @(x) any(strcmp(x, {'puzzle','dotring','comb', 'circle'})))
p.addParameter('borderwidth'    , 'fat'  , @(x) any(strcmp(x, {'skinny','fat'}))) % skinny=no gap or fat=+2 deg extra
p.addParameter('num'            , 1      , @isnumeric);                           % number of generated noise images
p.addParameter('pixoffset'      , [0 0]  , @isnumeric);                           % [x,y] offset of center in pixels

% Parse inputs
p.parse(params, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p


%% Generate pink noise image
pixelrange = [1 255]; % pixel range
scale_images_fixedrng = @(x,b) uint8((x - min(b)) / (max(b) - min(b)) * diff(pixelrange) + min(pixelrange));

im1 = NaN(params.disp.h_pix,params.disp.w_pix,num);
for ii = 1:num
    
    % generate pink noise image
    im0 = generatepinknoise(params.disp.w_pix,1,1,0); % KNK function: mode 0 means fixed amplitude spectrum + random phase
    
    % trim edges (initial generated pinknoise image is a square)
    screen_edge_pix = (params.disp.w_pix - params.disp.h_pix)/2;
    trimMe          = [1:screen_edge_pix; (params.disp.w_pix-screen_edge_pix+1):params.disp.w_pix];
    
    % Clip the range based on the predefined std
    std_norm        = params.stim.bckground.std_clip_range*std(im0(:));

    im0(trimMe,:) = [];
    im1(:,:,ii) = scale_images_fixedrng(im0,std_norm.*[-1 1]);
end

% deal with center offset
if isequal(pixoffset,[0 0]) || ~isempty(pixoffset)
    xc = params.disp.xc + pixoffset(1);
    yc = params.disp.yc + pixoffset(2);
else
    xc = params.disp.xc;
    yc = params.disp.yc;
end

% create support for addition stim masks
bckground_mask = true(size(im1,1),size(im1,2));
x = (1:params.disp.w_pix)-xc;
y = (1:params.disp.h_pix)-yc;
[XX,YY] = meshgrid(x,y);

% Do we want a skinny or fat border around the mask?
if strcmp(borderwidth,'fat')
    rim = 2*params.disp.ppd;
elseif strcmp(borderwidth,'skinny')
    rim = 0;
end

% Create the stimulus elements, take the union of the elements to create the cutout
switch gaptype
    
    case 'circle'
        
        % Create the circle mask.
        radius  = ((params.stim.ns.img_sz_pix/2) + rim) *params.disp.ppd;

        circMask = (YY - 0).^2 ...
            + (XX - 0).^2 <= radius.^2;

        bckground_mask = logical(circMask);
    
    case 'puzzle'
        
        % Create center square
        x_square = [1:(params.stim.ns.img_sz_pix+rim)] - ((params.stim.ns.img_sz_pix+rim)/2);
        bckground_mask(XX < min(x_square))=false;
        bckground_mask(XX > max(x_square))=false;
        bckground_mask(YY < min(x_square))=false;
        bckground_mask(YY > max(x_square))=false;
        
        % Next create the circles in the image.
        radius  = (params.stim.gabor.img_sz_pix + rim)/2;
        
        circle_left = (YY - params.stim.gabor.y0_pix(1)).^2 ...
            + (XX - params.stim.gabor.x0_pix(1)).^2 <= radius.^2;

        circle_right = (YY - params.stim.gabor.y0_pix(2)).^2 ...
            + (XX - params.stim.gabor.x0_pix(3)).^2 <= radius.^2;
        
        bckground_mask = bckground_mask + circle_left;
        bckground_mask = bckground_mask + circle_right;
        bckground_mask = logical(bckground_mask);

    case 'dotring'
        
        % Create the circles in the image.
        radius_inner  = (params.stim.dot.iso_eccen - (params.stim.dot.radius_deg) + rim) *params.disp.ppd;
        radius_outer  = (params.stim.dot.iso_eccen + (params.stim.dot.radius_deg) + rim) *params.disp.ppd;

        circle_inner = (YY - 0).^2 ...
            + (XX - 0).^2 <= radius_inner.^2;

        circle_outer = (YY - 0).^2 ...
            + (XX - 0).^2 <= radius_outer.^2;

        circMask = circle_outer-circle_inner;

        bckground_mask = logical(circMask);
        
    case 'comb'
        
        % Create center square
        x_square = [1:(params.stim.ns.img_sz_pix+rim)] - ((params.stim.ns.img_sz_pix+rim)/2);
        bckground_mask(XX < min(x_square))=false;
        bckground_mask(XX > max(x_square))=false;
        bckground_mask(YY < min(x_square))=false;
        bckground_mask(YY > max(x_square))=false;
        
        % Peripheral circles
        radius  = (params.stim.gabor.img_sz_pix + rim)/2;
        
        circle_left = (YY - params.stim.gabor.y0_pix(1)).^2 ...
            + (XX - params.stim.gabor.x0_pix(1)).^2 <= radius.^2;

        circle_right = (YY - params.stim.gabor.y0_pix(2)).^2 ...
            + (XX - params.stim.gabor.x0_pix(2)).^2 <= radius.^2;

        % Dot circle image.
        radius_inner  = ((params.stim.dot.iso_eccen - params.stim.dot.radius_deg)*params.disp.ppd) - (rim/2);
        radius_outer  = ((params.stim.dot.iso_eccen + params.stim.dot.radius_deg)*params.disp.ppd) + (rim/2);

        circle_inner = (YY - 0).^2 ...
            + (XX - 0).^2 <= radius_inner.^2;

        circle_outer = (YY - 0).^2 ...
            + (XX - 0).^2 <= radius_outer.^2;

        circMask = circle_outer-circle_inner;
        
        bckground_mask = bckground_mask + circle_left + circle_right + circMask;
        bckground_mask = logical(bckground_mask);
        

end

% Reshape back into an image with a gap
tmp = reshape(im1,size(im1,1)*size(im1,2),[]);
tmp(bckground_mask(:),:) = params.stim.bckgrnd_grayval;
tmp = uint8(tmp);

bckgrnd_im = reshape(tmp, size(bckground_mask,1),size(bckground_mask,2),num);
bckgrnd_im = repmat(bckgrnd_im, [1 1 1 3]);
bckgrnd_im = permute(bckgrnd_im, [1 2 4 3]);

% Store images if requested
if params.stim.store_imgs
    fprintf('[%s]:Storing background image(s)..',mfilename);
    tic
    saveDir = fileparts(fullfile(params.stim.bckground.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s_%s_%s.mat',params.stim.bckground.stimfile,gaptype,borderwidth,datestr(now,30))),'bckgrnd_im','-v7.3');
    toc
    fprintf('\n')
    
    saveFigDir = fullfile(vcd_rootPath,'figs',params.disp.name,'background');
    if ~exist(saveFigDir,'var') || isempty(saveFigDir)
        mkdir(saveFigDir); end
    
    for jj = 1:size(bckgrnd_im,4)
        filename = sprintf('%04d_vcd_background_%s%s.png',jj,gaptype, borderwidth);
        imwrite(bckgrnd_im(:,:,:,jj), fullfile(saveFigDir,filename));
    end
end

if params.verbose
    vcd_visualizeBackground(params, bckgrnd_im, [], gaptype, borderwidth)
end

return
