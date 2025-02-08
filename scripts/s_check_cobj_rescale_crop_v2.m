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
cat_name = {'bear','car','giraffe','parrot','damon','lisa','sophia', ...
            'banana', 'pizza','church', 'house', 'watertower','brush','drill','bus','suv'};

        
num_objs = 16;
num_views = 22;

% prep for ellipse plotting
phi    = linspace(0,2*pi,50);
cosphi = cos(phi); sinphi = sin(phi);

objstart = [1:(num_views/2):length(d)];

nPixels_object = [];
ratio = [];
idx = {};
avg_centroid_x = []; avg_centroid_y = [];
avg_majorAxis = []; avg_minorAxis = [];
avg_theta = []; avg_area = [];

object_vec = repelem(1:num_objs,2);
for cobj = 1:length(objstart)
    curr_cat = cat_name{object_vec(cobj)};
    
    figure; set(gcf,'Position',[1849,358,2282,412])
    file_idx = objstart(cobj)+[0:10];
    for vw = 1:length(file_idx)
        fname = fullfile(d(vw).folder,d(file_idx(vw)).name);
        [im, map, alpha] = imread(fname);
        
        % Crop the image to the bounding box.
        circlePixelsRGB = bsxfun(@times, im(1:1024,1:1024,:), cast(circlePixels, class(im)));
        
        propsIm = regionprops(rgb2gray(im(1:1024,1:1024,:)),'Area','Centroid',...
            'MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
        
        s = struct2cell(propsIm);
        centroids = reshape(cell2mat(s(2,:)),2,[])';
        p50 = prctile([centroids(:,1),centroids(:,2)],50); 
%         hold on; scatter(p50(1),p50(2),150,'y','filled')
        centerX_2 = p50(1);
        centerY_2 = p50(2);

        circlePixelsShifted = (Y - centerY_2).^2 ...
            + (X - centerX_2).^2 <= radius.^2;

        propsCircleMask = regionprops(circlePixelsShifted, 'BoundingBox');
        maskedImage = imcrop(circlePixelsRGB, propsCircleMask.BoundingBox);
        
        circleAlpha = bsxfun(@times, alpha(1:1024,1:1024), cast(circlePixelsShifted, class(alpha)));
        maskedAlpha = imcrop(circleAlpha,propsCircleMask.BoundingBox);

        circleAlpha_bool = logical(circleAlpha);
        
        sz = size(maskedImage(circleAlpha_bool));
        nPixels_object(cobj,vw) = sz(1);
        ratio(cobj,vw) = sz(1)/sum(circlePixels(:));
        
        % get cropped object stats:
        stats = regionprops(maskedAlpha,'Centroid','Area',...
            'MajorAxisLength','MinorAxisLength','Eccentricity','Orientation');
        
%         if ismember(objstart(cobj),[45,67,78,111,122,133,144,155,221,232,243,254, 265, 276, 309, 320])
            area_thresh = 180;
            majoraxis_thresh = 400;
%         else
%             area_thresh = 200;
%             majoraxis_thresh = 500;
%         end
            
        idx{cobj,vw} = find([stats.Area] > area_thresh & [stats.MajorAxisLength] > majoraxis_thresh);
        
        subplot(1,(num_views/2),vw);
        imagesc(maskedImage); axis square
        sgtitle(sprintf('%s (file nr %d-%d)', curr_cat, file_idx(vw)-10, file_idx(vw)));
        hold all;
        colors = lines(length(idx{cobj,vw}));
        
        clear centroid_x centroid_y a b theta
        for jj = 1:length(idx{cobj,vw})
            
            final_ellipse = idx{cobj,vw}(jj);
            
            centroid_x(jj) = stats(final_ellipse).Centroid(1);
            centroid_y(jj) = stats(final_ellipse).Centroid(2);
            a(jj) = stats(final_ellipse).MajorAxisLength/2;
            b(jj) = stats(final_ellipse).MinorAxisLength/2;
            theta(jj) = pi*stats(final_ellipse).Orientation/180;
            
            R = [ cos(theta(jj)) sin(theta(jj))
                -sin(theta(jj)) cos(theta(jj))];
            xy = [a(jj)*cosphi; b(jj)*sinphi];
            xy = R*xy;
            x = xy(1,:) + centroid_x(jj);
            y = xy(2,:) + centroid_y(jj);
            
            hold all;
            plot(stats(final_ellipse).Centroid(1),stats(final_ellipse).Centroid(2),'*','Color',colors(jj,:))
            plot(x,y,'Color',colors(jj,:),'LineWidth',2);
            
%             legendNames{jj} = sprintf('%i',idx(jj));
        end
        %     gg = gca;
        %     legend(gg.Children((length(gg.Children)-1):-2:1), legendNames)
        
        avg_centroid_x(cobj,vw) = mean(centroid_x);
        avg_centroid_y(cobj,vw) = mean(centroid_y);
        avg_majorAxis(cobj,vw) = mean(a);
        avg_minorAxis(cobj,vw) = mean(b);
        avg_theta(cobj,vw)     = circ_mean(theta,[],2);
        avg_area(cobj,vw)      = pi*mean(a)*mean(b);
        
        avg_R = [ cos(avg_theta(cobj,vw)) sin(avg_theta(cobj,vw))
            -sin(avg_theta(cobj,vw)) cos(avg_theta(cobj,vw))];
        avg_xy = [avg_majorAxis(cobj,vw)*cosphi; avg_minorAxis(cobj,vw)*sinphi];
        avg_xy = avg_R*avg_xy;
        avg_x = avg_xy(1,:) + avg_centroid_x(cobj,vw);
        avg_y = avg_xy(2,:) + avg_centroid_y(cobj,vw);
        
        hold all;
        plot(avg_centroid_x(cobj,vw),avg_centroid_y(cobj,vw),'*','Color','w')
        plot(avg_x,avg_y,'Color','w','LineWidth',4);
        title(vw);

    end
    
%     print(gcf,'-dpng','-r100',fullfile('~/Desktop/',sprintf('%s%d', curr_cat,file_idx(vw))))
end


%% Plot object specs

new_cat = 1:num_views:size(ratio,1);

vline_x = (num_views:num_views:length(d));
vline_x = repmat(vline_x',1,2);
within_cat_mn_ratio = mean(ratio,2);
within_cat_st_ratio = std(ratio,[],2);
within_cat_mn_abs   = mean(nPixels_object,2);
within_cat_st_abs  = std(nPixels_object,[],2);
within_cat_mn_area = mean(avg_area,2);
within_cat_st_area = std(avg_area,[],2);


%
figure(33); clf; set(gcf, 'Position', [1 1 1030 1277])
subplot(3,1,1);
stem(reshape(ratio',1,[]));
hold all;
vline_y = repmat([0 0.5],size(vline_x,1),1);
plot(vline_x',vline_y','r-')
plot([0,length(d)],[mean(ratio(:)),mean(ratio(:))],'g-', 'linewidth',3)
plot([0,length(d)],[mean(ratio(:))-std(ratio(:)),mean(ratio(:))-std(ratio(:))],'g:', 'linewidth',3)
plot([0,length(d)],[mean(ratio(:))+std(ratio(:)),mean(ratio(:))+std(ratio(:))],'g:', 'linewidth',3)
ylabel('image to aperture size ratio (pixels)')
xlabel('image nr')
set(gca,'FontSize',14)
xlim([0 length(ratio(:))])
set(gca,'XTick',[(num_views/2):num_views:length(ratio(:))], ...
        'XTickLabel',cat_name);
    
subplot(3,1,2); cla 
stem(reshape(avg_majorAxis',1,[]));
hold all;
vline_y = repmat([0 400],size(vline_x,1),1);
plot(vline_x',vline_y','r-')
plot([0,length(d)],[mean(avg_majorAxis(:)),mean(avg_majorAxis(:))],'g-', 'linewidth',3)
plot([0,length(d)],[mean(avg_majorAxis(:))-std(avg_majorAxis(:)),mean(avg_majorAxis(:))-std(avg_majorAxis(:))],'g:', 'linewidth',3)
plot([0,length(d)],[mean(avg_majorAxis(:))+std(avg_majorAxis(:)),mean(avg_majorAxis(:))+std(avg_majorAxis(:))],'g:', 'linewidth',3)

ylabel('Major axis (pixel radius)')
xlabel('image nr')
set(gca,'FontSize',14)
xlim([0 length(ratio(:))])
set(gca,'XTick',[(num_views/2):num_views:length(ratio(:))], ...
        'XTickLabel',cat_name);
    
subplot(3,1,3); cla 
stem(reshape(avg_area',1,[]));
hold all;
vline_y = repmat([0 10.^5.5],size(vline_x,1),1);
plot(vline_x',vline_y','r-')
plot([0,length(d)],[mean(avg_area(:)),mean(avg_area(:))],'g-', 'linewidth',3)
plot([0,length(d)],[mean(avg_area(:))-std(avg_area(:)),mean(avg_area(:))-std(avg_area(:))],'g:', 'linewidth',3)
plot([0,length(d)],[mean(avg_area(:))+std(avg_area(:)),mean(avg_area(:))+std(avg_area(:))],'g:', 'linewidth',3)

ylabel('Ellipse area (pixels)')
xlabel('image nr')
set(gca,'FontSize',14)
xlim([0 length(ratio(:))])
set(gca,'XTick',[(num_views/2):num_views:length(ratio(:))], ...
        'XTickLabel',cat_name);

