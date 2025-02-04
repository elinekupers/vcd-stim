load('stimulus_masks_NSD30_corrected.mat','masks')


im_per_cat = 6;
h_idx = 1:im_per_cat;
a_idx = im_per_cat+[1:im_per_cat];
f_idx = (2*im_per_cat)+[1:im_per_cat];
o_idx = (3*im_per_cat)+[1:im_per_cat];
p_idx = (4*im_per_cat)+[1:im_per_cat];

humans = squeeze(sum(masks(h_idx,:,:),1));
animals = squeeze(sum(masks(a_idx,:,:),1));
food = squeeze(sum(masks(f_idx,:,:),1));
objects = squeeze(sum(masks(o_idx,:,:),1));
places = squeeze(sum(masks(p_idx,:,:),1));

figure;
subplot(2,3,1)
imshow(humans,[0 6])
title('humans')
set(gca,'FontSize',20)
colorbar

subplot(2,3,2)
imshow(animals,[0 6])
title('animals')
set(gca,'FontSize',20)
colorbar

subplot(2,3,3)
imshow(food,[0 6])
title('food')
set(gca,'FontSize',20)
colorbar

subplot(2,3,4)
imshow(objects,[0 6])
title('objects')
set(gca,'FontSize',20)
colorbar

subplot(2,3,5)
imshow(places,[0 6])
title('places')
colorbar
set(gca,'FontSize',20)