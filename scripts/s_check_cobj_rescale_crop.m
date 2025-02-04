%% check_cobj_rescale_crop.m


% Create circle mask
imageSizeX = 1024;
imageSizeY = 1024;

[X, Y] = meshgrid(1:imageSizeX, 1:imageSizeY);

% Next create the circle in the image.
centerX = imageSizeX / 2;
centerY = imageSizeY / 2;
radius  = 400;

% translate
shift_x = 0; % (positive shifts make the circle go right, negative shifts make the circle go left)
shift_y = 0; % (positive shifts make the circle go down, negative shifts make the circle go up)
centerX = centerX + shift_x;
centerY = centerY + shift_y;

circlePixels = (Y - centerY).^2 ...
    + (X - centerX).^2 <= radius.^2;

figure;
imagesc(circlePixels) ;
title('circle mask');
axis square;


%% Load rescaled images

d = dir(fullfile(pwd,'*.png'));

clear nPixels ratio;

for ii = 1;%:length(d)
    fname = fullfile(d(ii).folder,d(ii).name);
    [im, map, alpha] = imread(fname);
        
    % Crop the image to the bounding box.
    circlePixelsRGB = bsxfun(@times, im, cast(circlePixels, class(im)));
    props = regionprops(circlePixels, 'BoundingBox');
    maskedImage = imcrop(circlePixelsRGB, props.BoundingBox);
    
    circleAlpha = bsxfun(@times, alpha, cast(circlePixels, class(alpha)));
    maskedAlpha = imcrop(circleAlpha,props.BoundingBox);
    
    
%     figure;
%     imagesc(maskedImage,'AlphaData',maskedAlpha);
    
    circleAlpha_bool = logical(circleAlpha);
    
    sz = size(maskedImage(circleAlpha_bool));
    nPixels(ii) = sz(1);
    ratio(ii) = sz(1)/sum(circlePixels(:));
    
    % get longest axis
    stats = regionprops(maskedAlpha,'Centroid',...
    'MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
    
    %% draw ellipse around longest angle
    phi = linspace(0,2*pi,50);
    cosphi = cos(phi); sinphi = sin(phi);
    figure; imagesc(maskedImage)
    hold on; 

    cmap = parula(10);

    for kk = 1:10

        xbar = stats(kk).Centroid(1);
        ybar = stats(kk).Centroid(2);
        a = stats(kk).MajorAxisLength/2;
        b = stats(kk).MinorAxisLength/2;
        theta = pi*stats(kk).Orientation/180;
        R = [ cos(theta) sin(theta)
            -sin(theta) cos(theta)];
        xy = [a*cosphi; b*sinphi];
        xy = R*xy;
        x = xy(1,:) + xbar;
        y = xy(2,:) + ybar;

        hold all;
        plot(stats(kk).Centroid(1),stats(kk).Centroid(2),'*','color',cmap(kk,:))
        plot(x,y,'color',cmap(kk,:),'LineWidth',2);
        title(kk)
        
%         waitforbuttonpress;
        
    end

    final_idx = [];
    im_majorAxisLength(ii) = stats(final_idx).MajorAxisLength;
end


%% Remove face image
nPixels(89) = [];
ratio(89) = [];

%% Plot the area 
num_objs = 22;
cat_name = {'bear','car','giraffe','parrot','banana', 'pizza','church', 'house', 'watertower','brush','drill','bus','suv'};
new_cat = 1:num_objs:length(ratio);

vline_x = (num_objs:num_objs:length(d));
vline_x = repmat(vline_x',1,2);

vline_y = repmat([0 0.5],1,length(vline_x));

for jj = 1:length(new_cat)
    within_cat_mn_ratio(jj) = mean(ratio(new_cat(jj):(new_cat(jj)+21)));
    within_cat_st_ratio(jj) = std(ratio(new_cat(jj):(new_cat(jj)+21)));
    within_cat_mn_abs(jj) = mean(nPixels(new_cat(jj):(new_cat(jj)+21)));
    within_cat_st_abs(jj) = std(nPixels(new_cat(jj):(new_cat(jj)+21)));
end



%% Plot object area ratio, and std
figure; stem(ratio); 
hold all;
plot(vline_x',vline_y','r-') 
plot([0,length(d)],[mean(ratio),mean(ratio)],'g-', 'linewidth',3)
plot([0,length(d)],[mean(ratio)-std(ratio),mean(ratio)-std(ratio)],'g:', 'linewidth',3)
plot([0,length(d)],[mean(ratio)+std(ratio),mean(ratio)+std(ratio)],'g:', 'linewidth',3)
ylabel('image to aperture size ratio (pixels)')
xlabel('image nr')
set(gca,'FontSize',14)

% yyaxis right
% plot([(num_objs/2):num_objs:length(ratio)],within_cat_st_ratio,'bo-');  
% ylabel('within im category std')

yyaxis right
plot([1:length(im_majorAxisLength)],im_majorAxisLength,'bo-');  
ylabel('Longest axis (pixels)')
