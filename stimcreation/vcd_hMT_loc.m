function [mot_fields, dotlocs, rng_seed] = vcd_hMT_loc(params, verbose, store_imgs)

% Contracting, expanding or stationary moving-dot fields for human MT+
% (hMT+) localizer
%
% STIMULUS PARAMETERS:
% * Moving dots traveled toward and away from fixation.
% * Speed: (8°/sec)
% * Aperture size: 12° diameter circular aperture.
% * Movie duration: 4 s (240 frames; assuming 60 Hz rate);
% * Dot color: 50/50 white/black dots on a mean luminance square gray background
% * Dot diameter: 0.25°.
% * Dot motion coherence: 100% (all dots move eccentric in either inward, outward or stationary position)
% * Dot lifetime: 0.1 s (6 frames; assuming 60 Hz rate)
% * Dot update rate: 1 (every frame)

% * Nr of movies:


% SPATIAL PARAMETERS
ap_diam_deg   = floor(params.disp.h_deg); % 12 deg visual angle
ap_diam_pix   = round((ap_diam_deg*params.disp.ppd)/2)*2;
dot_diam_pix  = 2*params.stim.rdk.dots_size_pix; % 6 pixels (0.068 deg) // x2 diameter of RDK stim (0.25 deg was used for Huk et al. 2002 J Neurosci)
rect_box_sz   = params.disp.h_pix; % height of the display (1080 pixels for BOLDscreen)

% Define an circular aperture in screen
stim_center = [0,0]; % [x y] in pixels. center on zero for now

% Define smaller dot aperture in pixels. (=actual stimulus)
% NOTE ap_radius is the 6-degree radius (in pixels) for the particular display,
% minus the size of a pixel radius from horizontal and vertical sides.
% * BOLDscreen:   2116/2 = 1058 --> 1058 - 6 pixels from each side = total diameter of 1046 pixels
stim_radius = [ap_diam_pix ap_diam_pix]./2 - dot_diam_pix;     % [w h] in pixels

% Define larger dot aperture in pixels. (=support stimulus)
% We then define a larger circular aperture where we draw individual dots
% where the diameter of the large circle is 3 x Full movement length
% (i.e., nr of pixels traveled by params.stim.rdk.dots_speed) so we have
% plenty of padding. We then mask this larger dot frame with a smaller
% circle (radius of 6 deg - radius of one dot to avoid dots being placed
% outside the smaller circular aperture).
% Dots that are outside the smaller aperture will be discarded.
% This means a particular dot can initially exist OUTSIDE the 12-diameter
% stimulus aperture (and thus will be ignored) and travel into the aperture
% during their lifetime (and thus will be present in the movie).
% Similarly, a particular dot can initially exist INSIDE the 12-diameter
% stimulus aperture (present) and travel outside the aperture during their
% lifetime (ignored).
support_radius_pix = [ap_diam_pix ap_diam_pix]./2;        % [w h] in pixels

% how many dots do we need to preserve the overall density.
support_area_deg2 = (pi*((rect_box_sz/params.disp.ppd).^2));            % 447.04 deg^2 for BOLD screen
dot_area_deg2     = (pi*(((0.5*dot_diam_pix)/params.disp.ppd).^2));     % 0.0036 deg^2 for BOLD screen

% We calculate the density of dots within the stimulus aperture.
% We ensure a round number such that we have an equal nr of white and blank
% dots. To get dot density, aim for same as 2% in Toottel et al. 1995
% (https://doi.org/10.1523/JNEUROSCI.15-04-03215.1995). It reads like Huk
% et al. 2002 JNeuro also used this density percentage.
ndot_density  = 0.0025; % percentages of total area within circular aperture (dots/deg^2) // for reference: RDK stim density = 15.9 dots / deg^2;
ndots_support = round(((ndot_density * support_area_deg2)/dot_area_deg2) /2)*2; % XX dots for BOLDscreen

% Define number of frames within a single RDK video
dur_frames   = 4*params.stim.presentationrate_hz;             % duration of individual dot movie in nr of presentation frames (each frame = 16.67 ms)
dot_interval = params.stim.rdk.dots_interval;
num_frames   = dur_frames/dot_interval;

% Define speed of contraction/expansion
speed_deg_s     = 8; % speed in deg per frame (currently set to 8 deg/s).  Same as rdk dots_speed (and Huk et al. 2002 MT loc)
speed_pix_frame = (speed_deg_s*params.disp.ppd)/params.stim.presentationrate_hz; % speed in pixels per frame

dot_interval   = 1;                                            % update dots every frame (60 Hz)
dot_lifetime   = 0.1 * params.stim.presentationrate_hz;        % 6 frames (0.1 seconds)

nr_iterations  = 1;

%% Preallocate space
% A cell array for each unique rdk video:
% 3 directions (in vs out vs stationary)
mot_fields     = cell(3,1);
dotlocs  = mot_fields;

%% RNG seed parameters
rseed    = [1000 2010];
num_col  = size(params.stim.rdk.dots_color,1);
rng_seed = NaN(size(mot_fields,1),1);

%% Folder to store RDK videos
saveStimDir = fullfile(vcd_rootPath, 'workspaces','stimuli',params.disp.name,'localizer',['hMT_' datestr(now,'yyyymmdd')]);
if ~exist(saveStimDir,'dir'), mkdir(saveStimDir); end
saveFigsDir = fullfile(vcd_rootPath, 'figs', params.disp.name, 'localizer',['hMT_' datestr(now,'yyyymmdd')]);
if ~exist(fullfile(saveFigsDir),'dir'), mkdir(fullfile(saveFigsDir)); end

%% Create localizer images
fH = figure(1); clf;
set(gcf,'Position',[0 0 rect_box_sz, rect_box_sz], ...
    'Units','Pixels','Renderer','OpenGL','PaperUnits','normalized')

dot_dir = [-1 1 0]; % inward, outward, stationary
dot_dir_name = {'inward','outward','stationary'};

for iter = 1:nr_iterations
    
    for bb = 3:length(dot_dir)
        
        % Initialize and reset rng for each RDK video (will apply to
        % rand, randn, and randi)
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',prod(rseed)));
        clock_seed = sum(100*clock);
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',clock_seed));
        rng_seed(bb) = clock_seed;
        
        % Get motion direction in degrees and print for user
        curr_motdir = dot_dir(bb);
        if verbose
            fprintf('\n[%s]: Motion direction: %s', mfilename, dot_dir_name{bb});
        end
        
        %% CREATE DOT MOVIE FRAMES
        
        % Allocate space for new dots for the entire support,
        % as well as their kill time and position in this rdk movie
        all_dots            = NaN(ndots_support,2);
        all_dots_kill_time  = NaN(ndots_support,1);
        stored_dot_pos  = NaN(ndots_support,2,num_frames); % only the ones that are in the stimulus aperture
        
        % Select color (50:50 black:white)
        colori = ([0 randperm(ndots_support-1)]) < ceil(ndots_support/num_col);
        col    = params.stim.rdk.dots_color(colori+1,:);
        assert(isequal(sum(colori==0),sum(colori==1)));
        
        % Create fresh batch of dots in larger circle
        for ii = 1:ndots_support
            
            enoughSpace = 0;
            while ~enoughSpace
                
                % initialize dot spatial positions by assigning random
                % [r,theta] for [eccen,angle] coordinates for the entire
                % support frame (later we will trim this)
                r          = sqrt(rand)*support_radius_pix(1); % eccentricity in pix relative to aperture center [0,0]
                theta      = (2*pi)*rand(); % angle in radians
                
                all_dots(ii,1) = (r.*cos(theta)); % cartesian coords (in pix)
                all_dots(ii,2) = (r.*sin(theta)); % cartesian coords (in pix)
                
                % check spacing relative to the other dots
                dotdist = sqrt((all_dots(~isnan(all_dots(~isnan(all_dots(:,1)),1)),1) - all_dots(ii,1)).^2 + ...
                    (all_dots(~isnan(all_dots(:,2)),2) - all_dots(ii,2)).^2);
                dotdist(ii) = 0; % distance of dot to itself is 0
                
                % check distance to center
                r0 = sqrt(all_dots(ii,1).^2 + all_dots(ii,2).^2);
                
                % If dots do not overlap and are not in the center, then we
                % continue
                if (min(dotdist(:)) >= 0) && (r0 > 0)
                    enoughSpace = 1;
                end
            end
            
            % Initialize dot kill time for all dots (starting from frame 1)
            if curr_motdir == 0
                num_frames = 1;
                all_dots_kill_time(ii) = num_frames+1;
            else
                all_dots_kill_time(ii) = 1 + rand().*dot_lifetime;
            end
        end
        
        
        % Reset video frames
        frames = [];
        
        % loop over time frames
        for curr_frame = 1:num_frames
            
            % Check whether the dot has reached its 'kill time'
            % based on frames.  If it's time to die
            % get a new random position.
            % Even if center is outside the stimulus circle
            % it will still "live" off screen in case it may
            % move inside the circle later..
            for ii = 1:ndots_support
                
                % convert to polar coords
                [pa, ecc] = cart2pol(all_dots(ii,1),all_dots(ii,2)); % [x,y] in pix to polar angle (radians) and eccentricity in pixels
                
                % Update eccentricity: If curr_motdir = -1 = contraction, +1 = expansion, 0 = stationary
                new_r = ecc + curr_motdir*speed_pix_frame;
                
                % If dot is at center, generate new dot location (leave dot kill time as is?)
                % If dot moves outside the aperture, we don't care? (or deal
                % with it later?)
                if (new_r <= 0)
                    new_r = sqrt(rand)*support_radius_pix(1);
                end
                
                % convert dot locs back to cartesian
                [x,y] = pol2cart(pa, new_r); % convert polar angle in rad and ecc in pix and updated eccentricity to [x,y] in pix
                all_dots(ii,:) = [x,y]; %
                
                if curr_frame > all_dots_kill_time(ii)
                    
                    % if so, we create dot with a new kill time and a
                    % new random x and y coordinate
                    enoughSpace = 0;
                    while ~enoughSpace
                        
                        % Update dot spatial positions by assigning new
                        % random [angle,eccen] coordinate
                        % (anywhere in the large circle). keep the same color
                        r     = sqrt(rand)*support_radius_pix(1);
                        theta = (2*pi)*rand();
                        
                        % cartesian coords (in deg vis angle)
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
                    all_dots_kill_time(ii) = curr_frame + dot_lifetime -1;
                end
            end
            
            
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
                    dot_diam_pix,dot_diam_pix, ...
                    [],[],{stim_col(mm,:)./255},180),'EdgeColor','none');
            end
            
            colormap gray; axis off square tight;
            set(gca, 'CLim',[0 1]);
            
            % Store location of each individual dot in case we want to check it later
            stored_dot_pos(inside,:,curr_frame) = all_dots(inside,:);
            
            % % debug visualization thing
            %         scatter(0,0,'rx');
            
            % set axis position
            axis off
            set(gca,'Units','normalized', 'Position',[0 0 1 1]);
            set(gcf,'InvertHardCopy', 'off');
            
            % set figure window position
            set(gcf,'Position',[0 0 884 884]) %setfigurepos([0 0 rect_box_sz rect_box_sz])
            
            % write to temp.png
            screenppd = get(0,'ScreenPixelsPerInch');
            files = printnice(fH,[1 88],'~/Desktop/temp');
            
            % read it in
            im0 = imread(files{1});
            
            % Crop
            assert(isequal(size(im0,1),size(im0,2)))
            if size(im0,1)~=1080
                im = imresize(im0,rect_box_sz/size(im0,1));
            else
                im = im0;
            end
            assert(isequal(size(im),[rect_box_sz rect_box_sz 3]))
            clear im0
            
            % store frame
            frames = cat(4,frames,im);
            
            % clean up
            clear f im
            
        end % frame loop
        
        % Store rdk video and dot positions in cell array
        mot_fields{bb} = frames;
        dotlocs{bb}    = stored_dot_pos;
        
        if store_imgs
            % save intermediate stage in case matlab crashes
            rdk_rng_seed = clock_seed;
            save(fullfile(saveStimDir, sprintf('%03d_vcd_hMT_%s.mat', iter, dot_dir_name{bb})),'frames','rdk_rng_seed','stored_dot_pos','-v7.3');
            
            % Create RDK movie (mp4, with compression)
            vcd_createStimVideo(frames, 1/params.stim.presentationrate_hz, ...
                saveFigsDir,sprintf('%03d_vcd_hMT_%s',iter, dot_dir_name{bb}), []);
            close(101)
        end
        
        clear frames
        
    end % motion direction
    
    if store_imgs
        % save iteration set
        save(fullfile(saveStimDir, sprintf('%03d_vcd_hMT_%s.mat', iter, datestr(now,30))),'mot_fields','dotlocs','rng_seed','-v7.3');
    end
end

return