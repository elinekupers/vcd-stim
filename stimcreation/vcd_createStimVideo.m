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
vidObj.FrameRate = 1/ifi;
open(vidObj);

% Create an animation.
figure(101); set(gcf, 'Units','Pixels','Renderer','OpenGL','PaperUnits','normalized')
set(gcf,'Position', [0 0 size(frames,1), size(frames,2)]);
set(gca,'Units','normalized');
set(gca,'Position',[0 0 1 1]);
set(gca,'CLim', [1 255]);
axis off;
set(gcf, 'InvertHardCopy', 'off');
screenppd = get(0,'ScreenPixelsPerInch');

if ndims(frames) == 3
    nframes = size(frames,3);
    colormap gray;
elseif ndims(frames) == 4
    nframes = size(frames,4);
else
    error('[%s]: Dimension mismatch, function expects ndims = 3 or 4',mfilename)
end

for k = 1:nframes
    
    currFrame.cdata = [];
    currFrame.colormap = [];
    
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

    % write to temp.png
    files = printnice(gcf,[1 screenppd],'~/Desktop/temp');
    
    % read it in
    im = imread(files{1});
    currFrame.cdata = im;
    writeVideo(vidObj,currFrame);
    
end

% Close the file.
close(vidObj);

end