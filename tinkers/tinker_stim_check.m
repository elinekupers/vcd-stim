% check stimuli

% Gabor
figure(1); 
for ii = 1:4 
    subplot(2,2,ii);
    for jj = 1:size(images.gabor.gabors{ii},4)
        for kk = 1:size(images.gabor.gabors{ii},5)
            imshow(images.gabor.gabors{ii}(:,:,:,jj,kk));
            title(sprintf('angle bin%d c%1.2f p%d',ii,jj,kk))
            drawnow;
        end
    end
end