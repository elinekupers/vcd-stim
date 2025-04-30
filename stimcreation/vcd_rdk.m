function [rdks, masks, info] = vcd_rdk(params)
% VCD function:
%  [rdks, masks, info] = vcd_rdk(params)
%
% Purpose:
%   Create a random dot motion kinetograms for experimental display. 
%   For rdk stimulus parameters see vcd_getStimParams.m
% 
%   We have 24 core rdk images sorted as follows:
%   * im 25:32  8 motion directions for coherence level 1
%   * im 33:40  8 motion directions for coherence level 2
%   * im 41:48  8 motion directions for coherence level 3
% 
%   2 spatial locations (left/right) are evenly distributed across 24
%   core images.
% 
%   We also create 4 test images for WM task: -15, -5 +5 +15 deg
%   motion direction offsets for each core rdk motion direction, with 
%   unique image nrs 207-302.
% 
%   When params.stim.store_imgs = true, 
%   * We store individual RDK movie mat files for each orientation,
%   coherence, and delta (together with one transparency mask, one table
%   row, corresponding center positions for all 200 dots on each frame).
%   These individual mat files are stored in a separate folder: 
%   fullfile(vcd_rootPath,'workspaces', 'stimuli',params.disp.name, ... 
%   ['rdk_'datestr(now,'yyyymmdd')]).
%   * We store a single, big matfile with all the movie frames, as well as
%   all transparency masks, all center dot positions for each frame in each
%   video, and big info table. Single mat file is defined by
%   params.stim.rdk.stimfile: fullfile(vcd_rootPath,'workspaces','stimuli',
%   params.disp.name,'rdk_<params.disp.name>_<datestr(now,'yyyymmdd')>.mat')
%
%   When params.verbose = true, we create mp4 RDK movies stored in the
%   folder fullfile(vcd_rootPath, 'figs', ['rdk_'datestr(now,'yyyymmdd')]);
%
% INPUTS:
%   params                  : stim params struct (see vcd_setStimParams.m)
%    *** this function requires the following struct fields ***
%    bckgrnd_grayval        : (int) background gray value (128)
%    rdk.img_sz_deg         : (int) size of image support in pixels, needs
%                               to be even integer number!
%    rdk.dots_size          : (int) how big is a single dot (radius in 
%                              pixels)
%    rdk.duration           : (int) how many frames is each rdk video?
%    rdk.dots_interval      : (int) how many frames until do we update the dots?
%    rdk.dots_direction     : (double) motion direction of coherent moving 
%                               dots (deg), 0 deg = 12 o'clock
%    rdk.max_dots_per_frame : (int) how many dots do we want per frame?
%    rdk.dots_coherence     : (double) what percentage of dots will move 
%                               coherently in the same direction (fraction)
%    rdk.dots_speed         : (double) how fast do dots move (pixels/frames)
%    rdk.delta_from_ref     : (double) offset motion direction from core 
%                               rdk motion direction (deg)
%    rdk.dots_lifetime      : (double) how long will a single dot be on the
%                               screen and move before we kill it and draw
%                               a new dot? (frames)
%    rdk.unique_im_nrs_core      : (int) number for each unique core rdk
%    rdk.unique_im_nrs_wm_test   : (int) number for each unique WM test rdk
%
% OUTPUTS:
%   rdks         : (uint8) unique RDK images used for VCD experiment, 
%                   8x3x5 cell for each motion direction, coherence level and delta offset (0, -15, -5, +5, +15)
%                   Each cell contains an uint8 image with dimensions: height (pixels) x width (pixels) x 3 (rgb) x 60 frames
%   masks        : (uint8) unique alpha mask images used for VCD experiment,
%                   8x3x5 cell for each motion direction, coherence level and delta offset (0, -15, -5, +5, +15):
%                   Each cell contains an uint8 image with dimensions: height (pixels) x width (pixels)
%   info         : table with rdk stimulus information matching the gabor array
%      unique_im        : (double) unique image nr for each RDK movie: 
%                           range 25-48 for core RDKs, 207-302.   
%      stim_pos_i       : (double) stimulus position index. 1=left, 2=right
%      stim_pos         : (cell) stimulus position, same as stim_loc_i but
%                           human readable ({'left'} or {'right'})
%      dot_motdir_deg   : (double) motion direction in degrees (0 = 12
%                           o'clock).
%      dot_motdir_deg_i : (double) same as dot_motdir_deg but indexed 1
%                           (smallest: 22.5 deg) to 8 (largest: 337.5 deg)
%      dot_coh          : (double) coherence level (fraction of 1), 
%      dot_coh_i        : (double) same as dot_coh but indexed 1 (lowest: 0.064)
%                           2 (medium: 0.128) or 3 (highest: 0.512).
%      rel_motdir_deg   : (double) motion direction relative from
%                           corresponding core RDK -15, -5, +5, +15 (deg)
%      rel_motdir_deg_i : (double) same as rel_motdir_deg but indexed 1
%                           (-15 deg), 2 (-5 deg), 3 (+5 deg), or 4 (+15 deg).            
%      dot_pos          : {1×1 cell} 200x2xtime dot positions on each frame
%                           relative to the center of the aperture (pixels).
%      is_in_img_ltm    : (logical) whether the RDK movie is part of the
%                           subselected stimuli used in imagery and
%                           long-term memory task.
%
% Written by Eline Kupers 2024/12, updated 2025/04

%% Set Display & RDK parameters

% Define an elliptic aperture in screen
ap_center = [0,0]; % [x y] in pixels. center on zero for now

% Define dot aperture in pixels. 
% NOTE: we shave off 1 dot radius (3 pixels from each side) to avoid dots being plotting outside the aperture
ap_radius = [params.stim.rdk.img_sz_pix./2 params.stim.rdk.img_sz_pix./2]-(params.stim.rdk.dots_size); % [w h] in pixels

% When exporting the frame, the RDK image turns out to be 548.7 pixels for 7TAS
% BOLDscreen; This is ~1.5 times bigger than the 4 deg RDK aperture we specified. 
% We adjust for this by increasing the size of the gray RDK box such that 
% the inner dot aperture will have the 4 deg size we want.
rect_box = params.stim.rdk.img_sz_pix * (577./params.stim.rdk.img_sz_pix);

% Define number of frames within a single RDK video
num_frames = params.stim.rdk.duration/params.stim.rdk.dots_interval;

% Get number of dots within a frame
ndots = params.stim.rdk.max_dots_per_frame;

% Check if we have delta images for WM
if ~isempty(params.stim.rdk.delta_from_ref)
    rdk_motdir_ref = [0:length(params.stim.rdk.delta_from_ref)];
else
    rdk_motdir_ref = 0;
end

% Organize unique image numbers for core and WM test images
img_im_nrs = reshape(params.stim.rdk.unique_im_nrs_wm_test,length(params.stim.rdk.delta_from_ref),[]);
img_im_nrs2{1} = img_im_nrs(1,1:length(params.stim.rdk.dots_direction)); % dir 1:8, coherence 1, delta = nr+dd
img_im_nrs2{2} = img_im_nrs(1,(length(params.stim.rdk.dots_direction)+1):(length(params.stim.rdk.dots_direction)*2)); % dir 1:8, coherence 2, delta = nr+dd
img_im_nrs2{3} = img_im_nrs(1,(2*(length(params.stim.rdk.dots_direction))+1):(length(params.stim.rdk.dots_direction)*3)); % dir 1:8, coherence 3, delta = nr+dd

%% Preallocate space

% A cell array for each unique rdk video: 
% nr of directions x nr of coherence levels x nr of deltas (WM quiz images)
rdks     = cell(length(params.stim.rdk.dots_direction),length(params.stim.rdk.dots_coherence),length(rdk_motdir_ref)); 
dotlocs  = rdks;
masks    = rdks; 

% nr of unique videos = 72 : 8 directions x 3 coherence levels x 5 deltas (none + 4)
n_unique_videos = size(rdks,1)*size(rdks,2)*size(rdks,3);

% info table to log order of rdk videos
info = table(NaN(n_unique_videos,1), ... unique_im
             NaN(n_unique_videos,1), ... dot_motdir_deg
             NaN(n_unique_videos,1), ... dot_motdir_deg_i
             NaN(n_unique_videos,1), ... dot_coh
             NaN(n_unique_videos,1), ... dot_coh_i
             NaN(n_unique_videos,1), ... rel_motdir_deg
             NaN(n_unique_videos,1), ... rel_motdir_deg_i
             NaN(n_unique_videos,1), ... is_in_img_ltm
             cell(n_unique_videos,1)); % dot_pos
info.Properties.VariableNames = {'unique_im','dot_motdir_deg','dot_motdir_deg_i',...
    'dot_coh','dot_coh_i','rel_motdir_deg','rel_motdir_deg_i','dot_pos','is_in_img_ltm'};

%% RNG seed parameters
rseed = [1000 2010];
num_col = size(params.stim.rdk.dots_color,1);
counter = 1;

%% Folder to store RDK videos
tmpDir = fullfile(vcd_rootPath, 'workspaces','stimuli',params.disp.name,['rdk_' datestr(now,'yyyymmdd')]);
if ~exist(tmpDir,'dir'), mkdir(tmpDir); end
saveStimDir = fullfile(vcd_rootPath, 'figs', params.disp.name, ['rdk_' datestr(now,'yyyymmdd')]);
if ~exist(fullfile(saveStimDir,'export'),'dir'), mkdir(fullfile(saveStimDir,'export')); end

%% Create rdk images
fH = figure(1); clf; 
set(gcf,'Position',[0,0,2*params.stim.rdk.img_sz_pix,2*params.stim.rdk.img_sz_pix], ...
    'Units','Pixels','Renderer','OpenGL','PaperUnits','normalized')

for cc = 1:length(params.stim.rdk.dots_coherence)
    
    for bb = 1:length(params.stim.rdk.dots_direction)
    
        for dd = 0:length(params.stim.rdk.delta_from_ref)
        
            % Get motion direction in degrees and print for user
            if dd == 0
                curr_motdir_deg = params.stim.rdk.dots_direction(bb);
                if params.verbose,
                    fprintf('\n[%s]: Motion direction: %2.2f deg', mfilename, curr_motdir_deg);
                end
            else
                curr_motdir_deg = params.stim.rdk.dots_direction(bb)+params.stim.rdk.delta_from_ref(dd);
                if params.verbose, fprintf('\n[%s]: Motion direction: %2.2f: %2.2f with %2.2f delta deg', ...
                    mfilename, curr_motdir_deg, params.stim.rdk.dots_direction(bb), params.stim.gabor.delta_from_ref(dd)); 
                end
            end

            % TRICKY STUFF: Subtract 90 deg to ensure 0 deg desired motion direction is 12 o'clock in x,y-pixel space 
            curr_motdir_deg = curr_motdir_deg-90;

            % Initialize and reset rng for each RDK video
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',prod(rseed)));
            clock_seed = sum(100*clock);
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',clock_seed));
            rdk_seed(cc,bb,dd+1) = clock_seed;
            
            %% RDK MOVIE TIME!! 
            
            % Allocate space for new dots, their kill time and position in this video
            dots               = NaN(ndots,2);
            dot_kill_time      = NaN(ndots,1);
            stored_coh_dot_pos = NaN(ndots,2,n_unique_videos);

            % Calculate velocity of dots
            v = randn(ndots,2);                                             % sample from random *normal* distribution
            v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)));                       % normalize to unit length
            
            % Calculate direction of dots, for noise (incoherent) and signal (coherent)
            mot_dir = [cosd(curr_motdir_deg) -sind(curr_motdir_deg)];                  % convert angle (deg) to xy direction
            v_signal_pix_per_frames   = params.stim.rdk.dots_speed .* repmat(mot_dir,ndots,1); %[v_x; v_y] (in pixels/frames)
            v_noise_pix_per_frames    = params.stim.rdk.dots_speed .* v;                       %[v_x; v_y] (in pixels/frames)
            
            % Assign proportion of signal and noise dot velocities
            dot_vel_pix_per_frames = v_noise_pix_per_frames;
            dot_vel_pix_per_frames(1:round(params.stim.rdk.dots_coherence(cc) .* ndots),:) = v_signal_pix_per_frames(1:round(params.stim.rdk.dots_coherence(cc) .* ndots),:);

            % Assign dot index (we just count from 1)
            dot_idx = 1:ndots;
            
            % Select color (50:50 black:white)
            colori = ([0 randperm(ndots-1)]) < ceil(ndots/num_col);
            col = params.stim.rdk.dots_color(colori+1,:);
            
            % Create fresh batch of dots
            for ii = 1:ndots
    
                enoughSpace = 0;
                while ~enoughSpace
                    
                    % initialize dot spatial positions by assigning random
                    % [angle,eccen] coordinates
                    r       = sqrt(rand);
                    theta   = 2*pi*rand();
                    dots(ii,1) = r.*cos(theta).*ap_radius(1);
                    dots(ii,2) = r.*sin(theta).*ap_radius(2);
                    
                    % check spacing relative to the other dots
                    mydist = sqrt((dots(:,1) - dots(ii,1)).^2 + (dots(:,2) - dots(ii,2)).^2);
                    mydist(ii) = 0;
                    
                    if (min(mydist(:)) >= 0)
                        enoughSpace = 1;
                    end
                end
                
                % Initialize dot kill time for all dots (starting from frame 1)
                dot_kill_time(ii) = 1 + rand().*params.stim.rdk.dots_lifetime; 
            end

            
            % Reset video frames
            frames = [];
            
            % loop over time frames
            for curr_frame = 1:num_frames

                % loop over dots
                for ii = 1:ndots

                    % update spatial position
                    dots(ii,:) = dots(ii,:) + dot_vel_pix_per_frames(ii,:) .* 1; %(1/params.stim.presentationrate_hz);

                    % Check whether the dot has reached its 'kill time'
                    % based on frames
                    if curr_frame > dot_kill_time(ii)
                        
                        % if so, we create dot with a new kill time and a
                        % new random x and y coordinate
                        enoughSpace = 0;
                        while ~enoughSpace

                            % Update dot spatial positions by assigning new
                            % random [angle,eccen] coordinate
                            r = sqrt(rand);
                            theta = 2*pi*rand();
                            
                            dots(ii,1) = r.*cos(theta).*ap_radius(1);
                            dots(ii,2) = r.*sin(theta).*ap_radius(2);
                            
                            % check spacing relative to the other dots
                            mydist = sqrt((dots(:,1) - dots(ii,1)).^2 + (dots(:,2) - dots(ii,2)).^2);
                            mydist(ii) = 0;
                            
                            if (min(mydist(:)) >= 0)
                                enoughSpace = 1;
                            end
                        end

                        dot_kill_time(ii) = curr_frame + params.stim.rdk.dots_lifetime;

                    end

                    % Check if dot has moved outside of the aperture and update and flip contrast if needed
                    inside = isInsideAperture(dots(ii,:), ap_center, ap_radius);
                    if ~inside
                        % move to point symmetric location
                        dots(ii,:) = -dots(ii,:); 
                        col(ii,:) = setdiff(params.stim.rdk.dots_color,col(ii,:),'rows'); % flip contrast
                    end
                end % ndots loop

                
                % Create frame with gray background
                figure(fH);
                clf; hold all;
                ax = gca;
                ax.Units = 'pixels';
                r=rectangle(ax,'Position', [(ap_center -(rect_box/2)), ...
                    rect_box, rect_box], ...
                    'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); %
                colormap gray; axis off square tight; 
                set(gca, 'CLim',[0 1]);
                for mm = 1:size(dots,1)
                    drawcircle('Parent',ax,'Center',dots(mm,:),'Radius',params.stim.rdk.dots_size,...
                        'Color',col(mm,:), 'InteractionsAllowed', 'none', 'FaceAlpha', 1, 'LineWidth', 1);

                end

                % store location in case we want to check it later
                stored_coh_dot_pos(:,:,curr_frame) = dots;
                
                offsetRect = ax.Children(end).Parent.Position;
                if mod(round(offsetRect(3)),2)==1
                    scf = floor(offsetRect(3))/offsetRect(3);
                    offsetRect = [0,0,scf*offsetRect(3),scf*offsetRect(3)];
                end
                f = getframe(ax,offsetRect);
                im = frame2im(f);
                
                fn = fullfile(saveStimDir,'export',sprintf('vcd_rdk_coh%02d_dir%02d_delta%02d_t%02d.png',cc,bb,dd,curr_frame));
                imwrite(im,fn);%,'Location',[pbRect(1:2)],'ScreenSize',[pbRect(3:4)]);

                % store frame
                frames = cat(4,frames,im);
                
                % clean up
                clear f im
                
            end % frame loop
            
            % Store rdk video and dot positions in cell array
            rdks{bb,cc,dd+1} = frames;
            dotlocs{bb,cc,dd+1} = stored_coh_dot_pos;
            
            % Log info 
            info.dot_motdir_deg(counter) = curr_motdir_deg + 90; % shift 90 deg back
            info.dot_coh(counter) = params.stim.rdk.dots_coherence(cc);
            info.dot_motdir_deg_i(counter) = bb;
            info.dot_coh_i(counter) = cc;
            info.dot_pos{counter} = {stored_coh_dot_pos};
            
            if dd == 0
                info.rel_motdir_deg(counter) = 0;
                info.rel_motdir_deg_i(counter) = 0;
                info.unique_im(counter) = params.stim.rdk.unique_im_nrs_core(bb + ((cc-1)*length(params.stim.rdk.dots_direction)));
            else
                info.rel_motdir_deg(counter) = params.stim.gabor.delta_from_ref(dd);
                info.rel_motdir_deg_i(counter) = dd;
                info.unique_im(counter) = img_im_nrs2{cc}(bb)+(dd-1);
            end

            info.is_in_img_ltm(counter) = (info.unique_im(counter)==params.stim.rdk.unique_im_nrs_specialcore);
            
            % Create binary circular alpha mask for rdk frames
            [XX,YY]       = meshgrid((1:size(frames,1))-(size(frames,1)/2),(1:size(frames,1))-(size(frames,1)/2)); 
            mask_radius   = ap_radius(1)+(4*params.stim.rdk.dots_size); % add 4 x dot pixel radius to avoid cutting off dots
            circlemask    = (YY - ap_center(1)).^2 + (XX - ap_center(2)).^2 <= mask_radius.^2;
            mask          = double(circlemask);
            mask(mask==1) = 255;
            mask          = uint8(mask);
    
            masks{bb,cc,dd+1} = mask;
    
            if params.store_imgs
                % save intermediate stage in case matlab crashes 
                rdk_info = info(counter,:);
                im_name  = bb + ((cc-1)*length(params.stim.rdk.dots_direction));
                save(fullfile(tmpDir, sprintf('%d_rdk_ori%d_coh%d_delta%d.mat', im_name, bb,cc,dd)),'frames','rdk_info','mask','stored_coh_dot_pos','-v7.3');
            end
            
            clear frames
            counter = counter +1;
        end
    end
end

% add stimulus location for each rdk orientation and coherence level
% (uneven numbers are left, even numbers are right).
stim_pos    = repmat(repelem({'left','right'},1+length(params.stim.gabor.delta_from_ref)), 1,length(params.stim.rdk.dots_direction)/2*length(params.stim.rdk.dots_coherence))';
stim_pos_i  = repmat(repelem([1,2],1+length(params.stim.gabor.delta_from_ref)), 1, length(params.stim.rdk.dots_direction)/2*length(params.stim.rdk.dots_coherence))';
info.stim_pos   = stim_pos;
info.stim_pos_i = stim_pos_i;

if params.store_imgs
    % in addition to storing separate conditions, we also store larger rdk file 
    saveDir = fileparts(fullfile(params.stim.rdk.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',params.stim.rdk.stimfile,datestr(now,30))),'rdks','info','masks','dotlocs','rdk_seed','-v7.3');
    
    saveDir = fileparts(fullfile(params.stim.rdk.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',params.stim.rdk.infofile,datestr(now,30))))
end

%% RDK

if params.verbose
    
    if params.store_imgs
        saveFigDir1 = fullfile(vcd_rootPath,'figs',params.disp.name,'rdk','visual_checks');
        saveFigDir2 = fullfile(vcd_rootPath,'figs',params.disp.name,'rdk','export');
        if ~exist(saveFigDir1,'dir'); mkdir(saveFigDir1); end
        if ~exist(saveFigDir2,'dir'); mkdir(saveFigDir2); end
    end
    
    fH1 = figure(1); 
    set(fH1, 'Position', [0 0 1024 1080], 'color','w')
    deltalabels = cellfun(@num2str, (num2cell(params.stim.rdk.delta_from_ref)),'UniformOutput', false);
    for ll = 1:length(deltalabels), 
        if strfind(deltalabels{ll},'-')
            deltalabels{ll} = strrep(deltalabels{ll},'-','min'); 
        else
            deltalabels{ll} = ['plus' deltalabels{ll}]; 
        end
    end
    deltalabels = [{'00'}, deltalabels(:)'];
    
    cohlabels = cellfun(@num2str, (num2cell(params.stim.rdk.dots_coherence*100)),'UniformOutput', false);
    cohlabels = cellfun(@(x) strrep(x,'.','pt'), cohlabels,'UniformOutput',false);
    
    % Loop over deltas
    for dd = 1:size(rdks,3)
        
        % Loop over coherence levels
        for cc = 1:size(rdks,1)
            
            % Loop over motion directions
            for bb = 1:size(rdks,1)
                
                im_idx = (info.dot_motdir_deg_i == params.stim.rdk.dots_direction(dd) & ...
                          info.dot_coh_i == params.stim.rdk.dots_coherence(cc) & ...
                          info.dot_motdir_deg_i == rdk_motdir_ref);
                im_nr  = info.unique_im(im_idx); 
                
                % Create RDK movie (mp4, with compression)
                vcd_createStimVideo(rdks{bb,cc,dd}, 1/params.stim.presentationrate_hz, ...
                    fullfile(vcd_rootPath,'figs',params.disp.name,'rdk'),sprintf('%03d_vcd_rdk_coh%02d_dir%02d_delta%02d',im_nr,cc,bb,dd-1));

                % Plot dot position & motion vector
                dot_pos = dotlocs{bb, cc, dd};
                dxdy    = dot_pos(:,:,2:end)-dot_pos(:,:,1:end-1);
                dxdy    = cat(3,zeros(size(dot_pos,1),size(dot_pos,2)),dxdy);

                for tt = 1:size(dot_pos,3) % loop over time
                    figure(fH1); cla;
                    xlim([-250, 250]); ylim([-250, 250]); box off;
                    plot(dot_pos(~isnan(dot_pos(:,1,tt)),1,tt),dot_pos(~isnan(dot_pos(:,2,tt)),2,tt),'o'); hold on;
                    quiver(dot_pos(~isnan(dot_pos(:,1,tt)),1,tt),dot_pos(~isnan(dot_pos(:,2,tt)),2,tt), ...
                        dxdy(~isnan(dot_pos(:,1,tt)),1,tt),dxdy(~isnan(dot_pos(:,2,tt)),2,tt));
                    title(sprintf('frame %02d: coh: %3.2f%% dir:%02d delta: %s',...
                            tt, params.stim.rdk.dots_coherence(cc)*100,params.stim.rdk.dots_direction(dd), deltalabels{dd})); 
                    axis square;
                    
                    if params.store_imgs
                        % store motion vector figure
                        if ~exist(fullfile(saveFigDir,sprintf('motdir%02d',dd)), 'dir')
                            mkdir(fullfile(saveFigDir,sprintf('motdir%02d',dd))); end
                        print(gcf,'-dpng','-r300',fullfile(saveFigDir1,sprintf('motdir%02d',dd),...
                            sprintf('%02d_vcd_rdk_coh%03d_dir%02d_delta%02d_motionvec%d.png',im_nr,cc,bb,dd,ii)));
                        
                        % store PNG of movie frame
                        imwrite(rdks{bb,cc,dd}, fullfile(saveFigDir2, ...
                            sprintf('%03d_vcd_rdk_coh%02d_dir%02d_delta%02d',im_nr,cc,bb,dd-1)));
                    end
                end
                
            end
        end
    end
end
return

