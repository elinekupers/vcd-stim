function vcd_createStimVideo(frames, ifi, saveFolder, fname, printtitle)

% check inputs
if ~exist('printtitle','var')
   printtitle = false;
end

% Prepare the new file.
if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end
vidObj = VideoWriter(fullfile(saveFolder,[fname '.mp4']),'MPEG-4');
open(vidObj);

% Create an animation.
figure(101);
set(gcf,'Position',[500   500   700  700], 'Units','Pixels','Renderer','OpenGL','PaperUnits','normalized')
set(gca,'CLim', [1 255]);

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
        imshow(frames(:,:,k)); 
        axis off square tight; 
        if printtitle
            title(sprintf('Frame %d - %3.3fs',k, k*ifi));
        end
        
    elseif ndims(frames) == 4
        % Assume color frame
        imshow(frames(:,:,:,k)); 
        axis off square tight; 
        if printtitle
            title(sprintf('Frame %d - %3.3fs',k, k*ifi));
        end
    end

    drawnow;
    
    % Write each frame to the file.
    currFrame = getframe(gca);
    writeVideo(vidObj,currFrame);
    
end

% Close the file.
close(vidObj);

end