function bckgrnd_im = vcd_pinknoisebackground(p, type, borderwidth, num)

if ~exist('type','var') || isempty(type)
    type = 'puzzle'; % can dotring or comb (puzzle+dotring)
end

if ~exist('borderwidth','var') || isempty(borderwidth)
    borderwidth = 'skinny'; % can also be fat (+2 deg extra)
end

if ~exist('num','var') || isempty(num)
    num = 1; % number of generated noise images
end

% generate pink noise image
pixelrange = [1 255]; % pixel range
scale_images_fixedrng = @(x,b) uint8((x - min(b)) / (max(b) - min(b)) * diff(pixelrange) + min(pixelrange));
im0 = generatepinknoise(p.disp.w_pix,1,num,0);

% trim edges (pinknoise image is square)
screen_edge_pix = (p.disp.w_pix - p.disp.h_pix)/2;
trimMe          = [1:screen_edge_pix; (p.disp.w_pix-screen_edge_pix+1):p.disp.w_pix];
std_norm        = p.stim.bckground.std_clip_range*std(im0(:));

im1 = NaN(p.disp.h_pix,p.disp.w_pix,num);
for ii = 1:num
    tmp = im0(:,:,ii);
    tmp(trimMe,:) = [];
    im1(:,:,ii) = scale_images_fixedrng(tmp,std_norm.*[-1 1]);
end

% create support for addition stim masks
bckground_mask = true(size(im1,1),size(im1,2));
x = (1:p.disp.w_pix)-p.disp.xc;
y = (1:p.disp.h_pix)-p.disp.yc;
[XX,YY] = meshgrid(x,y);

% Do we want a skinny or fat border around the mask?
if strcmp(borderwidth,'fat')
    rim = 2*p.disp.ppd;
elseif strcmp(borderwidth,'skinny')
    rim = 0;
end

% 
switch type
    
    case 'puzzle'
        
        % Create center square
        x_square = [1:(p.stim.ns.img_sz_pix+rim)] - ((p.stim.ns.img_sz_pix+rim)/2);
        bckground_mask(XX < min(x_square))=false;
        bckground_mask(XX > max(x_square))=false;
        bckground_mask(YY < min(x_square))=false;
        bckground_mask(YY > max(x_square))=false;
        
        % Next create the circles in the image.
        radius  = (p.stim.gabor.img_sz_pix + rim)/2;
        
        circle_left = (YY - p.stim.gabor.y0_pix(1)).^2 ...
            + (XX - p.stim.gabor.x0_pix(1)).^2 <= radius.^2;

        circle_right = (YY - p.stim.gabor.y0_pix(3)).^2 ...
            + (XX - p.stim.gabor.x0_pix(3)).^2 <= radius.^2;
        
        bckground_mask = bckground_mask + circle_left;
        bckground_mask = bckground_mask + circle_right;
        bckground_mask = logical(bckground_mask);

    case 'dotring'
        
        % Create the circles in the image.
        radius_inner  = (p.stim.dot.iso_eccen - (p.stim.dot.radius_deg) + rim) *p.disp.ppd;
        radius_outer  = (p.stim.dot.iso_eccen + (p.stim.dot.radius_deg) + rim) *p.disp.ppd;

        circle_inner = (YY - 0).^2 ...
            + (XX - 0).^2 <= radius_inner.^2;

        circle_outer = (YY - 0).^2 ...
            + (XX - 0).^2 <= radius_outer.^2;

        circMask = circle_outer-circle_inner;

        bckground_mask = logical(circMask);
        
    case 'comb'
        
        % Create center square
        x_square = [1:(p.stim.ns.img_sz_pix+rim)] - ((p.stim.ns.img_sz_pix+rim)/2);
        bckground_mask(XX < min(x_square))=false;
        bckground_mask(XX > max(x_square))=false;
        bckground_mask(YY < min(x_square))=false;
        bckground_mask(YY > max(x_square))=false;
        
        % Peripheral circles
        radius  = (p.stim.gabor.img_sz_pix + rim)/2;
        
        circle_left = (YY - p.stim.gabor.y0_pix(1)).^2 ...
            + (XX - p.stim.gabor.x0_pix(1)).^2 <= radius.^2;

        circle_right = (YY - p.stim.gabor.y0_pix(3)).^2 ...
            + (XX - p.stim.gabor.x0_pix(3)).^2 <= radius.^2;

        % Dot circle image.
        radius_inner  = ((p.stim.dot.iso_eccen - p.stim.dot.radius_deg)*p.disp.ppd) - rim;
        radius_outer  = ((p.stim.dot.iso_eccen + p.stim.dot.radius_deg)*p.disp.ppd) + rim;

        circle_inner = (YY - 0).^2 ...
            + (XX - 0).^2 <= radius_inner.^2;

        circle_outer = (YY - 0).^2 ...
            + (XX - 0).^2 <= radius_outer.^2;

        circMask = circle_outer-circle_inner;
        
        bckground_mask = bckground_mask + circle_left + circle_right + circMask;
        bckground_mask = logical(bckground_mask);
        

end

% bckgrnd_im = bsxfun(@(x,y) x.*y, im1, bckground_mask);
tmp = reshape(im1,size(im1,1)*size(im1,2),[]);
tmp(bckground_mask(:),:) = p.stim.bckgrnd_grayval;
bckgrnd_im = reshape(tmp, size(bckground_mask,1),size(bckground_mask,2),num);

% if p.verbose
%     figure; clf;
%     for jj = 1:num
%         subplot(1,num,jj); imagesc(bckgrnd_im(:,:,jj)); colormap gray
%         axis image
%     end
% end


if p.stim.store_imgs
    fprintf('\nStoring images..')
    saveDir = fileparts(fullfile(p.stim.gabor.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(saveDir,sprintf('background_%s_%s_%s_%s.mat',type,borderwidth,p.disp.name,datestr(now,30))),'bckgrnd_im','-v7.3');
end


return