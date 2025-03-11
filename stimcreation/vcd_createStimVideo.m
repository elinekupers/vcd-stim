function vcd_createStimVideo(frames, ifi, saveFolder, fname)

% Prepare the new file.
vidObj = VideoWriter(fullfile(saveFolder,[fname '.mp4']),'MPEG-4');
open(vidObj);

% Create an animation.
figure(101);
set(gca,'CLim', [0 255]);

if ndims(frames) == 3
    nframes = size(frames,3);
    colormap gray;
elseif ndims(frames) == 4
    nframes = size(frames,4);
else
    error('[%s]: Dimension mismatch, function expects ndims = 3 or 4',mfilename)
end

for k = 1:nframes
    
    cla
    if ndims(frames) == 3
        % Assume grayscale frame    
        imshow(frames(:,:,k)); title(sprintf('Frame %d - %3.3fs',k, k*ifi));
        
    elseif ndims(frames) == 4
        % Assume color frame
        imshow(frames(:,:,:,k)); title(sprintf('Frame %d - %3.3fs',k, k*ifi));
    end
    
    
    drawnow; axis equal tight off

    
    % Write each frame to the file.
    currFrame = getframe(101);
    writeVideo(vidObj,currFrame);
    
end

% Close the file.
close(vidObj);

end