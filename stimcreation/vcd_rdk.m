function [rdks, masks, info] = vcd_rdk(params, verbose, store_imgs)
% VCD function:
%  [rdks, masks, info] = vcd_rdk(params, verbose, store_imgs)
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
%   We also create 4 test images for WM task: -20, -10 +10 +20 deg
%   motion direction offsets for each core rdk motion direction, with
%   unique image nrs 143-174.
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
%   When verbose = true, we create mp4 RDK movies stored in the
%   folder fullfile(vcd_rootPath, 'figs', ['rdk_'datestr(now,'yyyymmdd')]);
%
% INPUTS:
%  params                  : stim params struct (see vcd_setStimParams.m)
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
%  verbose           : (logical) show debug figures
%  store_imgs        : (logical) store stimuli and debug figures as pngs
%
% OUTPUTS:
%  rdks         : (uint8) unique RDK images used for VCD experiment,
%                   8x3x5 cell for each motion direction, coherence level
%                   and delta offset (0, -20, -10, +10, +20)
%                   Each cell contains an uint8 image with dimensions:
%                   height (pixels) x width (pixels) x 3 (rgb) x 60 frames.
%                   Note that delta offset movies are only made for the
%                   highest coherence level.
%  masks        : (uint8) unique alpha mask images used for VCD experiment,
%                   8x3x5 cell for each motion direction, coherence level
%                   and delta offset (0, -20, -10, +10, +20):
%                   Each cell contains an uint8 image with dimensions:
%                   height (pixels) x width (pixels).
%                   Note that delta offset masks are only made for the
%                   highest coherence level.
%  info         : table with rdk stimulus information
%      unique_im        : (double) unique image nr for each RDK movie:
%                           range 25-48 for core RDKs, 143-174 for wm test movies.
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
%                           corresponding core RDK 0, -20, -10, +10, +20 (deg)
%      rel_motdir_deg_i : (double) same as rel_motdir_deg but indexed 1
%                           (-20 deg), 2 (-10 deg), 3 (+10 deg), or 4 (+20 deg).
%      dot_pos          : {1×1 cell} 200x2xtime dot positions on each frame
%                           relative to the center of the aperture (pixels).
%      is_specialcore    : (logical) whether the RDK movie is part of the
%                           subselected stimuli used in imagery and
%                           long-term memory task.
%
% Written by Eline Kupers 2024/12, updated 2025/04

%% Set spatial and temporal RDK movie parameters

% The total size of the RDK frame is ~1.5x larger than the support expected
% for a 4 deg RDK aperture. The BOLDscreen 4 deg aperture is 352 x 352 pixels
% and frame support is 544 x 544 pixels. The Eizoflexscan 4 deg aperture is
% 256 x 256 pixels and frame support is 396 x 396 pixels. Note that this
% support scale factor is not exactly 1.5 because we want an even number of
% pixels for the total frame size.
if strcmp(params.disp.name, '7TAS_BOLDSCREEN32')
    scf = 544/params.stim.rdk.img_sz_pix; % scf = 1.545454545454545; convert square support image from [92.52, 78.44, 545.6, 573.76] into [0 0 544 544];
    rect_box_sz = 546;
elseif strcmp(params.disp.name, 'PPROOM_EIZOFLEXSCAN')
    scf = 396/params.stim.rdk.img_sz_pix; % scf = 1.546875000000000; convert square support image from [67.56, 57.32, 396.8, 417.28] into [0 0 396 396];
    rect_box_sz = 396;
else
    scf = 1.545;
    rect_box_sz = scf*params.stim.rdk.img_sz_pix; % infer the rectangular support of the background by comparing the pixels per degree to the ppd of the BOLD screen.
end

% Define an elliptic aperture in screen
stim_center = [0,0]; % [x y] in pixels. center on zero for now

% Define smaller dot aperture in pixels. (=actual stimulus)
% NOTE ap_radius is the 2-degree radius in pixels for the particular display,
% minus the size of a pixel radius from horizontal and vertical sides.
% * BOLDscreen:   352/2 = 176 --> 176 - 3 pixels from each side = 173 pixels
% * Eizoflexscan: 256/2 = 128 --> 128 - 2 pixels from each side = 126 pixels
stim_radius = [params.stim.rdk.img_sz_pix params.stim.rdk.img_sz_pix]./2 - params.stim.rdk.dots_size_pix;     % [w h] in pixels

% Define larger dot aperture in pixels. (=support stimulus)
% We then define a larger circular aperture where we draw individual dots
% where the diameter of the large circle is 3 x Full movement length
% (i.e., nr of pixels traveled by params.stim.rdk.dots_speed) so we have
% plenty of padding. We then mask this larger dot frame with a smaller
% circle (radius of 2 deg - radius of one dot to avoid dots being placed
% outside the smaller circular aperture).
% Dots that are outside the smaller aperture will be discarded.
% This means a particular dot can initially exist OUTSIDE the 4-diameter
% stimulus aperture (and thus will be ignored) and travel into the aperture
% during their lifetime (and thus will be present in the movie).
% Similarly, a particular dot can initially exist INSIDE the 4-diameter
% stimulus aperture (present) and travel outside the aperture during their
% lifetime (ignored).
support_radius    = [rect_box_sz rect_box_sz]./2;     % [w h] in pixels

% how many dots do we need to preserve the overall density.
support_area_deg2 = (pi*((support_radius(1)/params.disp.ppd).^2));  % 30.2 deg2 for PProom

% We calculate the new number of dots that you'll need
% Get number of dots within a frame
ndots_support = round(params.stim.rdk.dots_density * support_area_deg2); % 480 dots for PProom
ndots_stim    = params.stim.rdk.max_dots_per_frame;

% Define number of frames within a single RDK video
num_frames = params.stim.rdk.duration/params.stim.rdk.dots_interval;

% Check if we have delta images for WM
if ~isempty(params.stim.rdk.delta_from_ref)
    rdk_motdir_ref = [0:length(params.stim.rdk.delta_from_ref)]; %#ok<NBRAK>
else
    rdk_motdir_ref = 0;
end

% Organize unique image numbers for core and WM test images
wm_im_nrs = reshape(params.stim.rdk.unique_im_nrs_wm_test,length(params.stim.rdk.delta_from_ref),[]);

%% Preallocate space

% A cell array for each unique rdk video:
% nr of directions x nr of coherence levels x nr of deltas (WM quiz images)
rdks     = cell(length(params.stim.rdk.dots_direction),length(params.stim.rdk.dots_coherence),length(rdk_motdir_ref));
dotlocs  = rdks;
masks    = rdks;

% actual nr of unique videos = 56:
% 8 directions x 3 coherence levels x delta = 0 deg
% + % 8 directions x 3 coherence levels x 4 deltas
% But here we overestimate and then trim down
n_unique_videos = size(rdks,1)*size(rdks,2)*size(rdks,3);

% info table to log order of rdk videos
info = table(NaN(n_unique_videos,1), ... unique_im
    NaN(n_unique_videos,1), ... dot_motdir_deg
    NaN(n_unique_videos,1), ... dot_motdir_deg_i
    NaN(n_unique_videos,1), ... dot_coh
    NaN(n_unique_videos,1), ... dot_coh_i
    NaN(n_unique_videos,1), ... rel_motdir_deg
    NaN(n_unique_videos,1), ... rel_motdir_deg_i
    NaN(n_unique_videos,1), ... is_special_core
    cell(n_unique_videos,1)); % dot_pos
info.Properties.VariableNames = {'unique_im','dot_motdir_deg','dot_motdir_deg_i',...
    'dot_coh','dot_coh_i','rel_motdir_deg','rel_motdir_deg_i','is_special_core','dot_pos'};

%% RNG seed parameters
rseed    = [1000 2010];
num_col  = size(params.stim.rdk.dots_color,1);
rdk_seed = NaN(size(rdks,1),size(rdks,2),size(rdks,3));

%% Folder to store RDK videos
tmpDir = fullfile(vcd_rootPath, 'workspaces','stimuli',params.disp.name,'rdk',['rdk_' datestr(now,'yyyymmdd')]);
if ~exist(tmpDir,'dir'), mkdir(tmpDir); end
saveStimDir = fullfile(vcd_rootPath, 'figs', params.disp.name, 'rdk',['rdk_' datestr(now,'yyyymmdd')]);
if ~exist(fullfile(saveStimDir),'dir'), mkdir(fullfile(saveStimDir)); end

%% Create rdk images
fH = figure(1); clf;
set(gcf,'Position',[1 1 2*params.stim.rdk.img_sz_pix, 2*params.stim.rdk.img_sz_pix], ...
    'Units','Pixels','Renderer','OpenGL','PaperUnits','normalized')

counter = 1;
for cc = 1:length(params.stim.rdk.dots_coherence)
    if verbose
        fprintf('\n[%s]: Dot coherence: %2.2f %%', mfilename, 100*params.stim.rdk.dots_coherence(cc));
    end
    for bb = 1:length(params.stim.rdk.dots_direction)
        
        for dd = 0:length(params.stim.rdk.delta_from_ref)
            
            % Get motion direction in degrees and print for user
            if dd == 0
                curr_motdir_deg = params.stim.rdk.dots_direction(bb);
                if verbose
                    fprintf('\n[%s]: Motion direction: %2.2f deg', mfilename, curr_motdir_deg);
                end
            else
                if dd > 0 && cc == 3
                    curr_motdir_deg = params.stim.rdk.dots_direction(bb)+params.stim.rdk.delta_from_ref(dd);
                    if verbose, fprintf('\n[%s]: Motion direction: %2.2f: %2.2f with %2.2f delta deg', ...
                            mfilename, curr_motdir_deg, params.stim.rdk.dots_direction(bb), params.stim.rdk.delta_from_ref(dd));
                    end
                else
                    curr_motdir_deg = NaN;
                end
            end
            
            
            if ~isnan(curr_motdir_deg)
                % TRICKY STUFF: Subtract 90 deg to ensure 0 deg desired motion direction is 12 o'clock in x,y-pixel space
                curr_motdir_deg = curr_motdir_deg-90;
                
                % Initialize and reset rng for each RDK video (will apply to
                % rand, randn, and randi)
                RandStream.setGlobalStream(RandStream('mt19937ar','seed',prod(rseed)));
                clock_seed = sum(100*clock);
                RandStream.setGlobalStream(RandStream('mt19937ar','seed',clock_seed));
                rdk_seed(cc,bb,dd+1) = clock_seed;
                
                %% RDK MOVIE TIME!!
                
                % Allocate space for new dots for the entire support,
                % as well as their kill time and position in this rdk movie
                all_dots            = NaN(ndots_support,2);
                all_dots_kill_time  = NaN(ndots_support,1);
                
                stored_coh_dot_pos  = NaN(ndots_support,2,num_frames); % only the ones that are in the stimulus aperture
                
                % Calculate direction of dots, for noise (incoherent) and signal (coherent)
                v_coh = [cos(deg2rad(curr_motdir_deg)) -sin(deg2rad(curr_motdir_deg))];    % convert angle from deg to radians, and then into dxdy direction (is already unit length)
                
                % Get direction vectors with random orientation (unit length)
                v_incoh = rand(ndots_support,1)*(2*pi);                                 % sample from random *uniform* motion directions (in radians)
                v_incoh = [cos(v_incoh) -sin(v_incoh)];
                v_incoh = bsxfun(@rdivide,v_incoh,sqrt(sum(v_incoh.^2,2)));     % normalize to unit length (EK: Is this needed?)
                
                % Create coherence (single motion direction) and incoherent (random motion direction) velocities
                v_signal_pix_per_frames   = params.stim.rdk.dots_speed .* repmat(v_coh,ndots_support,1); %[v_x; v_y] (in pixels/frames)
                v_noise_pix_per_frames    = params.stim.rdk.dots_speed .* v_incoh;               %[v_x; v_y] (in pixels/frames)
                
                % Assign proportion of signal and noise dot velocities
                dot_vel_pix_per_frames = v_noise_pix_per_frames;
                dot_vel_pix_per_frames(1:round(params.stim.rdk.dots_coherence(cc) .* ndots_support),:) = ...
                    v_signal_pix_per_frames(1:round(params.stim.rdk.dots_coherence(cc) .* ndots_support),:);
                
                % Select color (50:50 black:white)
                colori = ([0 randperm(ndots_support-1)]) < ceil(ndots_support/num_col);
                col    = params.stim.rdk.dots_color(colori+1,:);
                assert(isequal(sum(colori==0),sum(colori==1)));
                
                % Create fresh batch of dots in larger circle
                for ii = 1:ndots_support
                    
                    enoughSpace = 0;
                    while ~enoughSpace
                        
                        % initialize dot spatial positions by assigning random
                        % [r,theta] for [eccen,angle] coordinates
                        r          = sqrt(rand)*support_radius(1); % eccentricity in pixels relative to aperture center [0,0]
                        theta      = (2*pi)*rand(); % angle in radians
                        all_dots(ii,1) = (r.*cos(theta));
                        all_dots(ii,2) = (r.*sin(theta));
                        
                        % check spacing relative to the other dots
                        dotdist = sqrt((all_dots(~isnan(all_dots(~isnan(all_dots(:,1)),1)),1) - all_dots(ii,1)).^2 + ...
                            (all_dots(~isnan(all_dots(:,2)),2) - all_dots(ii,2)).^2);
                        dotdist(ii) = 0; % distance of dot to itself is 0
                        
                        if (min(dotdist(:)) >= 0)
                            enoughSpace = 1;
                        end
                    end
                    
                    % Initialize dot kill time for all dots (starting from frame 1)
                    all_dots_kill_time(ii) = 1 + rand().*params.stim.rdk.dots_lifetime;
                end
                
                
                % Reset video frames
                frames = [];
                
                % loop over time frames
                for curr_frame = 1:num_frames
                    
                    % Check whether the dot has reached its 'kill time'
                    % based on frames.  If it's time to die 
                    % get a new random position.
                    % Even if center is outside the stimulus circle,
                    % it will still "live" off screen in case it may
                    % move inside the circle later..
                    for ii = 1:ndots_support
                        if curr_frame > all_dots_kill_time(ii)
                            
                            % if so, we create dot with a new kill time and a
                            % new random x and y coordinate
                            enoughSpace = 0;
                            while ~enoughSpace
                                
                                % Update dot spatial positions by assigning new
                                % random [angle,eccen] coordinate
                                % (anywhere in the large circle). keep the same color
                                r     = sqrt(rand)*support_radius(1);
                                theta = (2*pi)*rand();
                                
                                all_dots(ii,1) = (r.*cos(theta));
                                all_dots(ii,2) = (r.*sin(theta));
                                
                                % check spacing relative to the other dots
                                dotdist = sqrt((all_dots(~isnan(all_dots(:,1)),1) - all_dots(ii,1)).^2 + ...
                                    (all_dots(~isnan(all_dots(:,2)),2) - all_dots(ii,2)).^2);
                                dotdist(ii) = 0; % distance of dot to itself is 0
                                
                                if (min(dotdist(:)) >= 0)
                                    enoughSpace = 1;
                                end
                            end
                        
                            % update dot_kill_time
                            all_dots_kill_time(ii) = curr_frame + params.stim.rdk.dots_lifetime -1;
                        end
                    end
                    
                    % update spatial position
                    all_dots = all_dots + dot_vel_pix_per_frames;
                    
                    % check which dots are inside the stimulus circle
                    inside = isInsideAperture(all_dots, stim_center, stim_radius);
                    
                    % Now grab stim dots (if dots aren't inside stim
                    % circle, then don't bother drawing)
                    stim_dots = all_dots(inside,:);
                    stim_col  = col(inside,:);
                    
                    
                    % Create frame with gray background
                    figure(fH);
                    clf; hold all;
                    ax = gca;
                    ax.Units = 'pixels';
                    
                    % Draw gray background as a square patch
                    rct = rectangle('Position',[stim_center(1)-(rect_box_sz/2),stim_center(2)-(rect_box_sz/2), ...
                        rect_box_sz,rect_box_sz],'FaceColor',[ones(1,3).*params.stim.bckgrnd_grayval]./255, 'EdgeColor', 'none');
                    
                    for mm = 1:size(stim_dots,1)
                        set(drawellipse(stim_dots(mm,1),stim_dots(mm,2),0, ...
                            params.stim.rdk.dots_size_pix,params.stim.rdk.dots_size_pix, ...
                            [],[],{stim_col(mm,:)./255},180),'EdgeColor','none');
                    end
                    
                    colormap gray; axis off square tight;
                    set(gca, 'CLim',[0 1]);
                    
                    % Store location of each individual dot in case we want to check it later
                    stored_coh_dot_pos(inside,:,curr_frame) = all_dots(inside,:);
                    
                    % % debug visualization thing
                    % scatter(0,0,'rx');
                    
                    % set figure window position
                    setfigurepos([0 0 scf*params.stim.rdk.img_sz_pix scf*params.stim.rdk.img_sz_pix]);
                    
                    % set axis position
                    set(gca,'Units','normalized');
                    set(gca,'Position',[0 0 1 1]);
                    axis off;
                    set(gcf, 'InvertHardCopy', 'off');
                    
                    % write to temp.png
                    screenppd = get(0,'ScreenPixelsPerInch');
                    files = printnice(fH,[1 screenppd],'~/Desktop/temp');
                    
                    % read it in
                    im = imread(files{1});
                    
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
                    info.rel_motdir_deg(counter)   = params.stim.rdk.delta_from_ref(dd);
                    info.rel_motdir_deg_i(counter) = dd;
                    info.unique_im(counter) = wm_im_nrs(dd,bb);
                end
                
                special_core_idx = find(info.unique_im(counter)==params.stim.rdk.unique_im_nrs_specialcore); %#ok<EFIND>
                if ~isempty(special_core_idx)
                    info.is_specialcore(counter) = true;
                else
                    info.is_specialcore(counter) = false;
                end
                
                % Create binary circular alpha mask for rdk frames
                [XX,YY]       = meshgrid((1:size(frames,1))-(size(frames,1)/2),(1:size(frames,1))-(size(frames,1)/2));
                mask_radius   = stim_radius(1)+(4*params.stim.rdk.dots_size_pix); % add 4 x dot pixel radius to avoid cutting off dots
                circlemask    = (YY - stim_center(1)).^2 + (XX - stim_center(2)).^2 <= mask_radius.^2;
                mask          = double(circlemask);
                mask(mask==1) = 255;
                mask          = uint8(mask);
                
                masks{bb,cc,dd+1} = mask;
                
                if store_imgs
                    % save intermediate stage in case matlab crashes
                    rdk_info = info(counter,:);
                    rdk_rng_seed = clock_seed;
                    save(fullfile(tmpDir, sprintf('%04d_vcd_rdk_ori%02d_coh%02d_delta%02d.mat', info.unique_im(counter), bb,cc,dd)),'frames','rdk_info','mask','rdk_rng_seed','stored_coh_dot_pos','-v7.3');
                    
                    % Create RDK movie (mp4, with compression)
                    vcd_createStimVideo(rdks{bb,cc,dd+1}, 1/params.stim.presentationrate_hz, ...
                        saveStimDir,sprintf('%04d_vcd_rdk_coh%02d_dir%02d_delta%02d',info.unique_im(counter),cc,bb,dd), []);
                    
                end
                
                clear frames
                counter = counter +1;
                
            end % nan check motion direction
        end % delta loop
    end % direction loop
end % coherence loop

% add stimulus location for each rdk orientation and coherence level
% (uneven numbers are left, even numbers are right).
stim_pos    = repmat(repelem({'left','right'},1+length(params.stim.rdk.delta_from_ref)), 1,length(params.stim.rdk.dots_direction)/2*length(params.stim.rdk.dots_coherence))';
stim_pos_i  = repmat(repelem([1,2],1+length(params.stim.rdk.delta_from_ref)), 1, length(params.stim.rdk.dots_direction)/2*length(params.stim.rdk.dots_coherence))';
info.stim_pos   = stim_pos;
info.stim_pos_i = stim_pos_i;
info = info(:,[1, 10, 11, 2:9]); % reorder info table columns
info = info(~isnan(info.unique_im),:);

if store_imgs
    % in addition to storing individual rdk movies, we also store larger rdk file
    saveDir = fileparts(fullfile(params.stim.rdk.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',params.stim.rdk.stimfile,datestr(now,30))),'rdks','info','masks','dotlocs','rdk_seed','-v7.3');
    
    saveDir = fileparts(fullfile(params.stim.rdk.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',params.stim.rdk.infofile,datestr(now,30))))
end

return


%% Debug figures

% if verbose
%
%     if store_imgs
%         saveFigDir1 = fullfile(vcd_rootPath,'figs',params.disp.name,'rdk','visual_checks');
% %         saveFigDir2 = fullfile(vcd_rootPath,'figs',params.disp.name,'rdk','export');
%         if ~exist(saveFigDir1,'dir'); mkdir(saveFigDir1); end
% %         if ~exist(saveFigDir2,'dir'); mkdir(saveFigDir2); end
%     end

% % Plot mean of all dots across frames for a single movie (times some scale factor)
% nrMovies = length(params.stim.rdk.dots_direction)*length(params.stim.rdk.dots_coherence);
% mn = cellfun(@(x) mean(x,4), reshape(rdks(:,:,1), nrMovies, 1), 'UniformOutput', false);
% for tt = 1:nrMovies,
%     f2 = mn{tt}; imwrite(uint8(f2(:,:,1)*4-300),fullfile(saveFigDir1,sprintf('mean_frames%02d.png',tt)));
% end
%
%
% % Plot the motion vector of each frame
%     fH1 = figure(1);
%     set(fH1, 'Position', [0 0 1024 1080], 'color','w')
%     deltalabels = cellfun(@num2str, (num2cell(params.stim.rdk.delta_from_ref)),'UniformOutput', false);
%     for ll = 1:length(deltalabels),
%         if strfind(deltalabels{ll},'-')
%             deltalabels{ll} = strrep(deltalabels{ll},'-','min');
%         else
%             deltalabels{ll} = ['plus' deltalabels{ll}];
%         end
%     end
%     deltalabels = [{'00'}, deltalabels(:)'];
%
%     cohlabels = cellfun(@num2str, (num2cell(params.stim.rdk.dots_coherence*100)),'UniformOutput', false);
%     cohlabels = cellfun(@(x) strrep(x,'.','pt'), cohlabels,'UniformOutput',false);
%
%    % Loop over deltas
%     for dd = 1:size(rdks,3)
%
%         % Loop over coherence levels
%         for cc = 1:size(rdks,2)
%
%             % Loop over motion directions
%             for bb = 1:size(rdks,1)
%
%                 im_idx = ismember(info.dot_motdir_deg_i, bb) & ...
%                          ismember(info.dot_coh,params.stim.rdk.dots_coherence(cc)) & ...
%                          ismember(info.rel_motdir_deg_i,rdk_motdir_ref(dd));
%                 im_nr  = info.unique_im(im_idx);
%
%                 % Plot dot position & motion vector
%                 dot_pos = dotlocs{bb, cc, dd};
%                 dxdy    = dot_pos(:,:,2:end)-dot_pos(:,:,1:end-1);
%                 dxdy    = cat(3,zeros(size(dot_pos,1),size(dot_pos,2)),dxdy);
%
%                 for tt = 1:size(dot_pos,3) % loop over time
%                     figure(fH1); cla;
%                     xlim([-250, 250]); ylim([-250, 250]); box off;
%                     plot(dot_pos(~isnan(dot_pos(:,1,tt)),1,tt),dot_pos(~isnan(dot_pos(:,2,tt)),2,tt),'o'); hold on;
%                     quiver(dot_pos(~isnan(dot_pos(:,1,tt)),1,tt),dot_pos(~isnan(dot_pos(:,2,tt)),2,tt), ...
%                         dxdy(~isnan(dot_pos(:,1,tt)),1,tt),dxdy(~isnan(dot_pos(:,2,tt)),2,tt));
%                     title(sprintf('frame %02d: coh: %3.2f%% dir:%3.2f delta: %s',...
%                             tt, params.stim.rdk.dots_coherence(cc)*100,params.stim.rdk.dots_direction(dd), deltalabels{dd}));
%                     axis square;
%
%                     if store_imgs
%                         % store motion vector figure
%                         if ~exist(fullfile(saveFigDir1,sprintf('motdir%02d',dd)), 'dir')
%                             mkdir(fullfile(saveFigDir1,sprintf('motdir%02d',dd))); end
%                         print(gcf,'-dpng','-r300',fullfile(saveFigDir1,sprintf('motdir%02d',dd),...
%                             sprintf('%02d_vcd_rdk_coh%03d_dir%02d_delta%02d_motionvec%d.png',im_nr,cc,bb,dd,tt)));
%
% %                         % store PNG of movie frame
% %                         imwrite(rdks{bb,cc,dd}(:,:,3,tt), fullfile(saveFigDir2, ...
% %                             sprintf('%03d_vcd_rdk_coh%02d_dir%02d_delta%02d.png',im_nr,cc,bb,dd-1)));
%                     end
%                 end
%
%             end
%         end
%     end
% end
% return

