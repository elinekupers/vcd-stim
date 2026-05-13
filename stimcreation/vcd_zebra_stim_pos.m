function zebra_im = vcd_zebra_stim_pos(params, varargin)
% VCD function to create the zebra pattern localizer stimuli for VCD stimulus position.
%
%   zebra_im = vcd_zebra_stim_pos(params, varargin)
% 
% OUTPUTS
%   zebra_im        : 4D array (uint8) pixels x pixels x 3 (RGB) x stimtypes x num of images
%
% ** Stim loc localizer **
% 'dotring'         = 1-deg width ring at 4 deg eccentricity 
% 'singledots'      = 1-deg diameter single dot (1 per image), requires num to be 16.
% 'alldots'         = all 16 1-deg diameter single dot stimuli in one image (16 per image)
% 'circleleft'      = 4-deg diameter circle (for peripheral classic stimuli) positioned at 4 deg eccentricity, left from the fixation circle, on the horizontal meridian
% 'circleright'     = 4-deg diameter circle (for peripheral classic stimuli) positioned at 4 deg eccentricity, right from the fixation circle, on the horizontal meridian
% 'circleleftright' = two 4-deg diameter circles (for peripheral classic stimuli) positioned at 4 deg eccentricity, both left and right from the fixation circle, on the horizontal meridian
% 'square'          = 8.4-deg square (for the NS stimuli) positioned in the center of the screen.

% params.disp = vcd_getDisplayParams('7TAS_BOLDSCREEN32');
% params.stim = vcd_getStimParams('disp_name','7TAS_BOLDSCREEN32');
% stimfile = fullfile(vcd_rootPath,'workspaces','stimuli','7TAS_BOLDSCREEN32',sprintf('afloc_zebra_pos_7TAS_BOLDSCREEN32_%s.mat',datestr(now,30)));
% zebra_im = vcd_create_zebra_stim_loc(params, 'stimtype', 'square', 'num', 3, 'zebra_cpf', 32, 'stimfile', stimfile)

%% Parse inputs
p = inputParser;
p.addRequired('params'          , @isstruct); % params struct
p.addParameter('stimtype'       , 'square', @(x) any(strcmp(x, {'none','dotring','singledots','alldots','circleleft','circleright','circleleftright','square'})))
p.addParameter('num'            , 1       , @isnumeric);                           % number of unique background images desired.  default: 1.
p.addParameter('pixoffset'      , [0 0]   , @isnumeric);                           % [x,y] offset of center in pixels  default = no offset: [0 0] 
p.addParameter('zebra_cpf'      , 32      , @isnumeric);                           % cycles per field of view for zebra pattern. One cycle = 1 black-white reversal. Higher numbers = finer pattern. Default: 32. 
p.addParameter('store_imgs'     , false, @islogical);
p.addParameter('save_dir'       , fullfile(vcd_rootPath,'figs',params.disp.name,'localizer','afloc'), @ischar);

% Parse inputs
p.parse(params, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p

if ~isfield(params,'verbose') || isempty(params.verbose)
    params.verbose = false;
end

% deal with center offset
if isequal(pixoffset,[0 0]) || ~isempty(pixoffset)
    xc = params.disp.xc + pixoffset(1);
    yc = params.disp.yc + pixoffset(2);
else
    xc = params.disp.xc;
    yc = params.disp.yc;
end


% For STIMULUS POSITION
if strcmp(stimtype,'single-dots')
    if num < 16
        num = 16;
        warning('[%s]: Number of requested images is less than number of single dots! Will change num to 16', mfilename)
    else
        % do nothing
    end
end

%%%%%%%%%%%%%% CREATE zebra patterns full-field %%%%%%%%%%%%%%

% preallocate space:
% start with square with 1.5 * width screen resolution (zebra pattern can
% only be square), then we shave off top/bottom.
zebra_im0  = NaN(params.disp.h_pix,params.disp.h_pix);  % square 1080x1080
zebra_im   = NaN(params.disp.h_pix,params.disp.h_pix, 3, num); % dims height pixels x height pixels x 3 (RBG), num requested images

for pp = 1:num
    
    % Create zebra pattern background
    tmp = (imagefilter(rand(1.5*params.disp.h_pix,1.5*params.disp.h_pix)-0.5,constructbutterfilter(1.5*params.disp.h_pix,1.5*zebra_cpf,5)) > 0) - 0.5;  % values either -.5 or .5
    tmp_zebra = placematrix(zebra_im0,tmp,[]); % crop zebra pattern to 1080x1080 screen height resolution
    
    
    % Create support for addition stim masks
    bckground_mask = true(params.disp.h_pix,params.disp.w_pix);
    x = (1:params.disp.w_pix)-xc;
    y = (1:params.disp.h_pix)-yc;
    [XX,YY] = meshgrid(x,y);
   
   % Create the stimulus elements (take the union of the elements if multiple elements)
   switch stimtype
       
       case 'none'
           % do nothing, keep it gray
           bckground_mask = [];
           
       case 'circleleft'
           % Next create the circles in the image.
           radius  = (params.stim.gabor.img_sz_pix)/2;
           
           circle_left = (YY - params.stim.gabor.y0_pix(1)).^2 ...
               + (XX - params.stim.gabor.x0_pix(1)).^2 <= radius.^2;
           
           bckground_mask = circle_left;
           bckground_mask = logical(bckground_mask);
           
       case 'circleright'
           % Next create the circles in the image.
           radius  = (params.stim.gabor.img_sz_pix)/2;
           
           circle_right = (YY - params.stim.gabor.y0_pix(2)).^2 ...
               + (XX - params.stim.gabor.x0_pix(2)).^2 <= radius.^2;
           
           bckground_mask = circle_right;
           bckground_mask = logical(bckground_mask);
           
       case 'circleleftright'
           % Next create the circles in the image.
           radius  = (params.stim.gabor.img_sz_pix)/2;
           
           circle_left = (YY - params.stim.gabor.y0_pix(1)).^2 ...
               + (XX - params.stim.gabor.x0_pix(1)).^2 <= radius.^2;
           
           circle_right = (YY - params.stim.gabor.y0_pix(2)).^2 ...
               + (XX - params.stim.gabor.x0_pix(2)).^2 <= radius.^2;
           
           bckground_mask = circle_left + circle_right;
           bckground_mask = logical(bckground_mask);
           
       case 'singledots'
           
           % Get single dot radius
           radius = params.stim.dot.radius_pix;
           
           % Loop over 16 [x,y] locations
           dd = mod(pp-1,length(params.stim.dot.x0_pix))+1;
           dot_loc_x = params.stim.dot.x0_pix(dd) - xc; % pix coords are relative to upper right corner [0,0], but grid is relative to center
           dot_loc_y = params.stim.dot.y0_pix(dd) - yc;
           single_dot_im = (YY -  dot_loc_y ).^2 ...
               + (XX - dot_loc_x).^2 <= radius.^2;
           
           bckground_mask = single_dot_im;
           bckground_mask = logical(bckground_mask);
           
       case 'alldots'
           
           % Get single dot radius
           radius = params.stim.dot.radius_pix;
           
           % Loop over 16 [x,y] locations
           tmp = [];
           for dd = 1:length(params.stim.dot.x0_pix)
               dot_loc_x = params.stim.dot.x0_pix(dd) - xc; % pix coords are relative to upper right corner [0,0], but grid is relative to center
               dot_loc_y = params.stim.dot.y0_pix(dd) - yc;
               single_dot_im = (YY -  dot_loc_y ).^2 ...
                   + (XX - dot_loc_x).^2 <= radius.^2;
               if dd == 1
                   tmp = single_dot_im;
               else
                   tmp = tmp + single_dot_im;
               end
           end
           bckground_mask = tmp;
           bckground_mask = logical(bckground_mask);
           clear tmp
           
       case 'dotring'
           
           % Create the circles in the image.
           radius_inner  = (params.stim.dot.iso_eccen - (params.stim.dot.radius_deg)) *params.disp.ppd;
           radius_outer  = (params.stim.dot.iso_eccen + (params.stim.dot.radius_deg)) *params.disp.ppd;
           
           circle_inner = (YY - 0).^2 ...
               + (XX - 0).^2 <= radius_inner.^2;
           
           circle_outer = (YY - 0).^2 ...
               + (XX - 0).^2 <= radius_outer.^2;
           
           circMask = circle_outer-circle_inner;
           
           bckground_mask = logical(circMask);
           
       case 'square'
           
           % Create center square
           x_square = [1:(params.stim.ns.img_sz_pix)] - ((params.stim.ns.img_sz_pix)/2);
           bckground_mask(XX < min(x_square))=false;
           bckground_mask(XX > max(x_square))=false;
           bckground_mask(YY < min(x_square))=false;
           bckground_mask(YY > max(x_square))=false;
           
           bckground_mask = logical(bckground_mask);
   end
   
    % Invert contrast 
   bckground_mask = double(bckground_mask);
   
   % crop background mask to 1080x1080 screen height resolution
   bckground_mask_crop = placematrix(zebra_im0,bckground_mask,[]); 
   
   fig = figure;
   imagesc(bckground_mask_crop); axis image off 
   colormap gray
   set(gca,'CLim',[0 1]); colorbar off
   
   % PATH = '/opt/homebrew/bin/magick'
   % setenv PATH /usr/local/git/bin:/Users/kupers/bin:/opt/subversion/bin:/usr/local/bin:${PATH}
   mask0 = renderfigure(size(bckground_mask,1),2);  % double, decimal between 0 and 1
   close(fig);
   
   % crop zebra image with background mask (still hxh screen resolution)
   tmp_zebra_crop = tmp_zebra.*mask0;

   % Convert to [1 255] range
   tmp_zebra_crop = 128 + tmp_zebra_crop*254;
   
   % Duplicate for third channel
   zebra3D_im = repmat(tmp_zebra_crop, [1 1 3]);
   
   % Convert from double to uint8
   zebra3D_im = uint8(zebra3D_im);
   
   % Put in final array
   zebra_im(:,:,:,pp) = zebra3D_im;
   
   % save image if requested
   if store_imgs
       save_dir2 = fullfile(save_dir,'zebra_stim_pos',stimtype);
       if ~exist(save_dir2,'dir')
           mkdir(save_dir2);
       end
       if strcmp(stimtype,'singledots')
           imwrite(zebra3D_im, fullfile(save_dir2, sprintf('%03d_vcd_afloc_zebra_%s%02d.png', pp, stimtype, dd)));
       else
           imwrite(zebra3D_im, fullfile(save_dir2, sprintf('%03d_vcd_afloc_zebra_%s.png', pp, stimtype)));
       end
   end
   
   clear zebra_im3D tmp_zebra_crop tmp_zebra mask0 bckground_mask_crop bckground_mask_crop
end