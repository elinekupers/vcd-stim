%% s_preprocessObjects.m
% This script processes 'raw' object images, dealing with grayscale, position, size, luminance, and contrast.
%
% This script relies on KNK utils functions (github.com/cvnlab/knkutils).
%
% *** Order of operations ***
% We convert images to grayscale (using rgb2gray).
% We interpret images using a squaring luminance response.
% We center images with respect to each image's center of mass (calculated from the alpha mask).
%   The centering is done using integral pixels (so no resampling is performed).
%   The center of the images becomes 512.5 (for 1024 pixel images) using the 
%   convention that 1 and 1024 are the centers of the first and last pixels.
% Next, for each object, we compute the mean alpha mask across viewpoints. We
%   then determine the minimal square (centered at the center of the images) that
%   passes <squarethresh>.
% The determined square extent is used to scale the size of the object viewpoint images
%   such that <targetsize> is achieved. Cubic interpolation is used in this resizing.
%   Note that some actual object content CAN extend beyond <targetsize> because
%   the square merely covers most of the alpha masks.
% We next proceed to mean/std (luminance/contrast) normalization.
% We determine the additive offset to be applied to all viewpoint images such 
%   that the grand mean of luminance of the viewpoint images is <targetmnlum>.
%   (Each image's luminance is calculated as the mean within its alpha mask.)
% We then perform an iterative approach to determine the single scale factor
%   to be applied to all viewpoint images such that the grand mean of contrast of the
%   viewpoint images is within 1% of <targetsdlum>. (Each image's contrast is
%   calculated as the std dev within the full image. Note that we do not compute 
%   std dev only within the alpha mask.) When scaling each image, the scaling is 
%   performed using its mean luminance (within the mask) as the anchor point; 
%   hence, the scaling does not affect the mean luminance of the image.
%
% We write visualization figures. Review these carefully!
% We write some workspace variables to record.mat.
% We write the final prepared images as .png files.
%   Note that the interpretation of these files is intended for linear monitors;
%   specifically, we use integer range 1 through 255 and these are intended
%   to be shown on a linear display. If you look at the .png files on a
%   standard display, they will appear darker than they actually would be
%   on a linear display. 

%% DEAL WITH CONSTANTS

% Get display params
% Choose: '7TAS_BOLDSCREEN32'   - BOLD screen at the 7TAS MRI scanner
%         'KKOFFICE_AOSQ3277'   - external monitor in kendrick's CMRR office
%         'EKHOME_ASUSVE247'    - external monitor at Eline's home
%         'PPROOM_EIZOFLEXSCAN' - EizoFlexscan monitor @ CMRR's psychophysics room
%         'CCNYU_VIEWPIXX3D'    - ViewPixx monitor in Clay Curtis' psychophysics room 
dispname    = '7TAS_BOLDSCREEN32'; 
disp_params = vcd_getDisplayParams(dispname);

%% Define/Load stimulus params
% !!WARNING!! There is a randomization component involved in creating some
% stimuli (e.g., orientation of gabor stimuli or dot locations). If you
% don't want this, this leave the fifth argument:
% "overwrite_randomized_params" empty (default is set to false) or set to
% false.
%
% If you do want regenerate probabilistic params, set the fifth argument to
% true and some stimulus values will change.
load_stim_params                 = false; % if false, re-create params.
                                          % if true, we load stored mat file
store_stim_params                = false; % if false, we don't store params. 
                                          % if true, we store mat file in fullfile(vcd_rootPath,'workspaces','info')                                       
% Input 1: Display params struct (see vcd_getDisplayParams.m)
% Input 2: Load prior stored parameters or not. 
% Input 3: Store generated parameters or not. 
stim_params   = vcd_getStimParams('disp_name', dispname, ...
                                  'load_params', load_stim_params,...
                                   'store_params', store_stim_params); 

% Parameters (CHANGE ONLY THESE SETTINGS)
numrot       = 91;      	% number of image viewpoints
numobj       = 16;          % number of objects
ogsize       = 1024;        % original pixel size of object stimuli
targetsize   = stim_params.obj.og_res_stim_target_sz;         % number of pixels for one side of conformed square (354 for BOLD screen, 258 for PP room)
targetmnlum  = 0.5;         % desired mn-luminance for grand average of one object's viewpoints
targetsdlum  = 0.015;       % desired sd-luminance for grand average of one object's viewpoints
squaresizes  = 2:2:1000;    % square sizes (in pixels for one side) to evaluate
squarethresh = 0.99;        % previously 0.9; square size is chosen that includes at least 90% of the mass (of the mean mask)

% define folders and params
dir0         = fullfile(vcd_rootPath,'workspaces','stimuli','RAW','vcd_objects','all_to_process'); % where 'raw' images live
outputdir    = fullfile(vcd_rootPath,'workspaces','stimuli',dispname, 'vcd_objects_2degstep_lumcorrected');  % where final images are written to (this directory is REMOVED if it exists)
figuredir    = fullfile(vcd_rootPath,'workspaces','stimuli',dispname, ...
                sprintf('objects_preproc_lummn%s_sd%s_sz%dpct',strrep(num2str(targetmnlum),'.','pt'),strrep(num2str(targetsdlum),'.','pt'),squarethresh*100));

% pproom stimuli have wider histogram and smaller object size on 1024x1024.
% this is boldscreen uses more pixels than pproom (on the 1024x1024).
% to match contrast from boldscreen to pproom, let's use the strategy of
% embedding the final boldscreen image in a larger spatial extent when
% calculating contrast.
if strcmp(dispname,'PPROOM_EIZOFLEXSCAN')
    monitorfactor = 1;  % for PPROOM case
elseif strcmp(dispname,'7TAS_BOLDSCREEN32')
    % for BOLDscreen case (i.e. boldscreen ppd / pproom ppd)
    dp1 = vcd_getDisplayParams('7TAS_BOLDSCREEN32');
    dp2 = vcd_getDisplayParams('PPROOM_EIZOFLEXSCAN');
    monitorfactor = dp1.ppd/dp2.ppd; % monitorfactor = 88.189484478585300/63.902348145300280;  
    clear dp1 dp2
else
    error('[%s]: Welp, we don''t know what monitorfactor to use for this display setup!',mfilename)
end

assert(monitorfactor >= 1);  % we DO NOT want to reduce as we might risk truncation!!!!
fudgesize = round(ogsize * monitorfactor);

%% LOAD IMAGES

% find the images
files0 = matchfiles([dir0 '/*.png']);
assert(length(files0)==numrot*numobj);

% load all images
%   (images: convert to grayscale and then double. becomes 0-255.)
%   (alpha: convert to double. becomes 0-255.)
allimages = zeros(ogsize,ogsize,length(files0));
allalphas = zeros(ogsize,ogsize,length(files0));
for p=1:length(files0), p
  [im0,~,alpha0] = imread(files0{p});
  allimages(:,:,p) = double(rgb2gray(im0));
  allalphas(:,:,p) = double(alpha0);
end

% initialize final output
allimages2 = zeros(ogsize,ogsize,length(files0),'uint8');
allalphas2 = zeros(ogsize,ogsize,length(files0),'uint8');
validpct = zeros(1,length(files0));
OF = zeros(1,numobj);
SC = zeros(1,numobj);

% make dirs
rmdirquiet(outputdir);
mkdirquiet(outputdir);
mkdirquiet(figuredir);

%% PROCESS THE IMAGES

% loop over each object
for zz=1:numobj

  %% EXTRACT STUFF

  % extract images for the object and convert to luminance (assuming squaring)
  ims = (allimages(:,:,(zz-1)*numrot + (1:numrot))/255).^2;  % now in range 0-1
  
  % extract alphas and convert to fraction between 0-1
  als = allalphas(:,:,(zz-1)*numrot + (1:numrot))/255;       % now in range 0-1 (1 means object, 0 means background)
  
  %% CENTER BASED ON CENTER-OF-MASS
  
  % figure out center of each image and shift the image and alpha to be centered at 512.5.
  % note the rounding (hence, image pixel values don't actually get perturbed).
  for yy=1:size(als,3)
    com = centerofmass(als(:,:,yy),[1 2]);  % [row; col] in decimal matrix coordinates
    ims(:,:,yy) = circshift(ims(:,:,yy),round(((ogsize+1)/2)-com)');  % 512.5 be careful about the signs!!!
    als(:,:,yy) = circshift(als(:,:,yy),round(((ogsize+1)/2)-com)');
  end
  
  %% FIGURE OUT SQUARE SIZE
  
  % prepare mean mask and compute total mass
  meanmask = mean(als,3);
  tot0 = sum(flatten(meanmask));  % total mass
    
  % determine conformed square size
  rec = [];
  for ss=1:length(squaresizes)
    crop = ceil((ogsize+1)/2)-squaresizes(ss)/2 : floor((ogsize+1)/2)+squaresizes(ss)/2; % 513 vs 512
    rec(ss) = sum(flatten(meanmask(crop,crop)));
  end
  finalsqsz = squaresizes(firstel(find(rec/tot0 > squarethresh)));  % number of pixels on a side
  
  % visualize the results!!
  crop = ceil((ogsize+1)/2)-finalsqsz/2 : floor((ogsize+1)/2)+finalsqsz/2;   % 513 vs 512 pixels
  figureprep([100 100 800 800]); hold on;
  imagesc(meanmask); axis image tight;
  plotrectangle([crop(1)-0.5 crop(end)+0.5 crop(1)-0.5 crop(end)+0.5],'r-');
  set(gca,'YDir','reverse');
  figurewrite(sprintf('obj%02d',zz),[],[],figuredir);

  %% APPLY SCALING SUCH THAT SQUARE IS DESIRED SIZE
  
  % how much scale factor do we want?
  S = targetsize/finalsqsz;
  
  % apply scaling (note that rounding to integer pixels occurs and placement rounding occurs!)
  for yy=1:size(ims,3)
    im0 = normalizerange(imresize(ims(:,:,yy),S,'bicubic'),0,1,0,1);
    al0 = normalizerange(imresize(als(:,:,yy),S,'bicubic'),0,1,0,1);
    ims(:,:,yy) = placematrix(zeros(ogsize,ogsize),im0,[]);
    als(:,:,yy) = placematrix(zeros(ogsize,ogsize),al0,[]);
  end
  
  % check that we didn't enlarge too much
  assert(allzero(flatten(als([1 end],:,yy))));
  assert(allzero(flatten(als(:,[1 end],yy))));
  
 %% DEAL WITH MEAN AND SD OF LUMINANCE
  
  % compute mean of each image (within the alpha mask)
  mns = [];
  for yy=1:size(ims,3)
    mask0 = als(:,:,yy) > 0.5;  % note simple hack
    im0 = ims(:,:,yy);
    mns(yy) = mean(im0(mask0));
  end

  % figure out adjustment factor for mean
  OF(zz) = targetmnlum - mean(mns);  % what offset to add to ensure grand mean of mean-luminance meets the target      

  % after dealing with mean, compute sd of each image (use the whole image!!)
  sds = [];
  for yy=1:size(ims,3)

    % apply and mix
    newim0 = ims(:,:,yy) + OF(zz);
    newim0 = newim0 .* als(:,:,yy) + ims(:,:,yy) .* (1-als(:,:,yy));  % mix according to alpha (to avoid background going crazy)
    
    % make composite image
    comp0 = newim0 .* als(:,:,yy) + targetmnlum * (1-als(:,:,yy));

    % record
    sds(yy) = std(flatten(placematrix(targetmnlum*ones(fudgesize,fudgesize),comp0,[])));

  end

  % figure out initial guess for adjustment factor for std
  SC(zz) = targetsdlum / mean(sds);  % what scale to apply to ensure grand mean of std-luminance meets the target
  
  % figure out initial guess for adjustment factor for std
  SC(zz) = targetsdlum / mean(sds);  % what scale to apply to ensure grand mean of std-luminance meets the target
  
  % loop until we are happy
  while 1
  
    % apply adjustment factors and make RGB images again
    imsAFTER = zeros(size(ims));
    for yy=1:size(ims,3)
    
      % apply and mix
      newim0 = SC(zz) * (ims(:,:,yy) - mns(yy)) + mns(yy) + OF(zz);  % adjust scale and then adjust offset
      newim0 = newim0 .* als(:,:,yy) + ims(:,:,yy) .* (1-als(:,:,yy));  % mix according to alpha (to avoid background going crazy)
    
      % count percentage of OK pixels
      mask0 = als(:,:,yy) > 0.5;
      validpct((zz-1)*numrot + yy) = sum(newim0(mask0) >= 0 & newim0(mask0) <= 1) / sum(mask0(:)) * 100;
    
      % truncate
      newim0(newim0 < 0) = 0;
      newim0(newim0 > 1) = 1;
    
      % record
      imsAFTER(:,:,yy) = newim0;
      allimages2(:,:,(zz-1)*numrot + yy) = uint8(1+254*newim0);   %uint8(1+254*sqrt(newim0));   %uint8(255*sqrt(newim0));
      allalphas2(:,:,(zz-1)*numrot + yy) = uint8(255*als(:,:,yy));

    end

    % compute empirical std devs
    imtemp = allimages2(:,:,(zz-1)*numrot + (1:numrot));
    altemp = allalphas2(:,:,(zz-1)*numrot + (1:numrot));
    temp = ((double(imtemp)-1)/254) .* double(altemp)/255 + targetmnlum * (1-double(altemp)/255);  % NOTE that uint8 is already luminance (no squaring necessary)!
    estds = std(squish(placematrix(targetmnlum*ones(fudgesize,fudgesize,numrot),temp,[]),2),[],1);
    
    % if mean is within 1% of targetsdlum, we are done
    if abs(mean(estds)-targetsdlum) / targetsdlum < .01
      break;
    end

    % if we are off, apply guess for the update
    SC(zz) = SC(zz) * (targetsdlum / mean(estds));
    fprintf('**** ANOTHER (%d): sd lum target %1.3f vs mn empirical: %1.4f. SC(zz) is %.4f ****\n',zz, targetsdlum, mean(estds), SC(zz));

  end

  % visualize the histogram before and after
  figureprep([100 100 500 800]);
  subplot(2,1,1); hold on; hist(ims(als > 0.5),100); xlim([0 1]);
  subplot(2,1,2); hold on; hist(imsAFTER(als > 0.5),100); xlim([0 1]);
  figurewrite(sprintf('hist%02d',zz),[],[],figuredir);
    
end

% write images to disk
for ii=1:size(allimages2,3)
  imwrite(allimages2(:,:,ii),[outputdir '/' stripfile(files0{ii},1)],'Alpha',allalphas2(:,:,ii));
end

% visualize
figureprep;
bar(OF,1);
figurewrite('OF',[],[],figuredir);

% visualize
figureprep;
bar(SC,1);
figurewrite('SC',[],[],figuredir);

% visualize
figureprep([100 100 900 300]);
plot(validpct,'r.-');
figurewrite('validpct',[],[],figuredir);

% visualize
temp = ((double(allimages2)-1)/254) .* double(allalphas2)/255 + targetmnlum * (1-double(allalphas2)/255);  % NOTE that uint8 is already luminance (no squaring necessary)!
temp = placematrix(targetmnlum*ones(fudgesize,fudgesize,size(allimages2,3)),temp,[]);
emeans = mean(squish(temp,2),1);
estds = std(squish(temp,2),[],1);
figureprep([100 100 800 400]);
subplot(2,1,1); hold on;
plot(emeans);
straightline((1:numrot:numrot*numobj)-.5,'v','k-');
title('Empirical means');
subplot(2,1,2); hold on;
plot(estds);
straightline((1:numrot:numrot*numobj)-.5,'v','k-');
title('Empirical std devs');
figurewrite('empirical',[],[],figuredir);

% save workspace (except for some large variables)
saveexcept([figuredir '/record.mat'],{'allimages' 'allalphas' 'allimages2' 'allalphas2' 'ims' 'imsAFTER' 'als' 'temp' 'imtemp' 'altemp'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% ixtouse = [1:91];  % 
% 
%   sqim = zeros(1024,1024);
%   sqim(crop,crop) = 1;
% 
% %  [params,R2] = fitgaussian2d(mean(als,3),[],1);
% %  [xx,yy] = meshgrid(1:1024,1:1024);
% %  subplot(1,2,2); imagesc(evalgaussian2d(params,xx,yy)); axis image tight;
%   figure;
%   subplot(1,2,1); imagesc(meanmask); axis image tight;
%   subplot(1,2,2); imagesc(sqim); axis image tight;
% 
