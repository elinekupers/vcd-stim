function [rdk, mask, info, p] = vcd_rdk(p)
% VCD function:
%  [rdk, mask, info, p] = vcd_rdk(p)
%
% Purpose:
%   Create a random dot motion kinetograms for  experimental display. 
%   For rdk stimulus parameters see vcd_getStimParams.m
% 
% INPUTS:
%  p            : rdk params
%   p.disp         : monitor display params (struct) pixels per degree in field of view (pix)
%
% OUTPUTS:
%   rdk     : RDK images, 8x3x5 cell for each motion direction, coherence
%               level and delta offset (0, -15, -5, +5, +15)
%               Each cell contains w (pixels) x h (pixels) x 3 x frames
%   info    : table with rdk stimulus information (motion direction, coherence,
%               location of dots for each individual frame), matching the rdk array
%   p       : updated params struct
%%
% Written by Eline Kupers 2024/12
%
% For reference: unique rdks images
%   uniqe_im   stim_loc   mot_dir    coh
%     1.0000    1.0000   18.0000    0.0640
%     2.0000    2.0000   62.0000    0.0640
%     3.0000    1.0000   98.0000    0.0640
%     4.0000    2.0000  152.0000    0.0640
%     5.0000    1.0000  192.0000    0.0640
%     6.0000    2.0000  236.0000    0.0640
%     7.0000    1.0000  282.0000    0.0640
%     8.0000    2.0000  326.0000    0.0640
%     9.0000    1.0000   18.0000    0.1280
%    10.0000    2.0000   62.0000    0.1280
%    11.0000    1.0000   98.0000    0.1280
%    12.0000    2.0000  152.0000    0.1280
%    13.0000    1.0000  192.0000    0.1280
%    14.0000    2.0000  236.0000    0.1280
%    15.0000    1.0000  282.0000    0.1280
%    16.0000    2.0000  326.0000    0.1280
%    17.0000    1.0000   18.0000    0.5120
%    18.0000    2.0000   62.0000    0.5120
%    19.0000    1.0000   98.0000    0.5120
%    20.0000    2.0000  152.0000    0.5120
%    21.0000    1.0000  192.0000    0.5120
%    22.0000    2.0000  236.0000    0.5120
%    23.0000    1.0000  282.0000    0.5120
%    24.0000    2.0000  326.0000    0.5120

%% Set Display & RDK parameters

% Define dot apeture in pixels
p.stim.rdk.dots_aperture = floor([0 0 p.stim.rdk.img_sz_deg./2 p.stim.rdk.img_sz_deg./2].*p.disp.ppd); % [x y w h] in pixels

% Define an elliptic aperture in screen
ap_center = p.stim.rdk.dots_aperture(1:2);
ap_radius = p.stim.rdk.dots_aperture(3:4)-p.stim.rdk.dots_size; % shave off 3 pixels (one dot) to avoid dot falling outside the aperture

% Define number of frames within a single RDK video
num_frames = p.stim.rdk.duration/p.stim.rdk.dots_interval;

% Get number of dots within a frame
ndots = min(p.stim.rdk.max_dots_per_frame, ...
    round(p.stim.rdk.dots_density .* (p.stim.rdk.dots_aperture(:,3).*p.stim.rdk.dots_aperture(:,4)) / p.stim.framedur_s));

% Check if we have delta images for WM
if ~isempty(p.stim.rdk.delta_from_ref)
    rdk_motdir_ref = [0:length(p.stim.rdk.delta_from_ref)];
else
    rdk_motdir_ref = 0;
end

%% Preallocate space

% A cell array for each unique rdk video: 
% nr of directions x nr of coherence levels x nr of deltas (WM quiz images)
rdk = cell(length(p.stim.rdk.dots_direction),length(p.stim.rdk.dots_coherence),length(rdk_motdir_ref)); 

% nr of unique videos = 72 : 8 directions x 3 coherence levels x 5 deltas (none + 4)
n_unique_videos = size(rdk,1)*size(rdk,2)*size(rdk,3);

% info table to log order of rdk videos
info = table(NaN(n_unique_videos,1), NaN(n_unique_videos,1), cell(n_unique_videos,1),NaN(n_unique_videos,1),NaN(n_unique_videos,1));
info.Properties.VariableNames = {'dot_dir','dot_coh','dot_pos','motdir_deg_ref','unique_im'};

%% RNG seed parameters
rseed = [1000 2010];
num_col = size(p.stim.rdk.dots_color,1);
counter = 1;

%% Folder to store RDK videos
tmpDir = fullfile(vcd_rootPath, 'workspaces','stimuli',p.disp.name,['rdk_' datestr(now,'yyyymmdd')]);
if ~exist('tmpDir','dir'), mkdir(tmpDir); end
saveStimDir = fullfile(vcd_rootPath, 'figs', ['rdk_' datestr(now,'yyyymmdd')]);

%% Create rdk images
fH = figure(1); clf; 
set(gcf,'Position',[0,0,p.stim.rdk.img_sz_pix,p.stim.rdk.img_sz_pix], 'Units','Pixels','Renderer','zbuffer')

for cc = 1:length(p.stim.rdk.dots_coherence)
    
    for bb = 1:length(p.stim.rdk.dots_direction)
    
        for dd = 0:length(p.stim.rdk.delta_from_ref)
        
            if dd == 0
                curr_motdir_deg = p.stim.rdk.dots_direction(bb);
                fprintf('\nMotion direction: %d deg', curr_motdir_deg)
            else
                curr_motdir_deg = p.stim.rdk.dots_direction(bb)+p.stim.rdk.delta_from_ref(dd);
                fprintf('\nMotion direction: %d: %d with %d delta deg', ...
                    curr_motdir_deg, p.stim.rdk.dots_direction(bb), p.stim.gabor.delta_from_ref(dd))
            end

        
            % Initialize and reset rng for each RDK video
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',prod(rseed)));
            RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
           
            
            % Get new dots for trial
            dots = [];
            dot_kill_time = [];
            distance_before_kill = [];

            v = randn(ndots,2); % sample from random *normal* distribution
            v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2))); % Normalize to unit length
            
            mot_dir = [cosd(curr_motdir_deg) -sind(curr_motdir_deg)];
            v_signal_deg_per_sec   = p.stim.rdk.dots_speed .* repmat(mot_dir,ndots,1); % v_x, v_y          % experiment speeds per eye (in deg/frames)
            v_noise_deg_per_sec    = p.stim.rdk.dots_speed .* v;
            
            % Assign proportion of signal and noise dot velocities
            dot_vel_deg_per_sec = v_noise_deg_per_sec;
            dot_vel_deg_per_sec(1:round(p.stim.rdk.dots_coherence(cc) .* ndots),:) = v_signal_deg_per_sec(1:round(p.stim.rdk.dots_coherence(cc) .* ndots),:);
            
            stored_coh_dot_pos = NaN(ndots,2,n_unique_videos);
            
            % Assign the x and y dot positions
            dot_idx = 1:ndots;
            
            % Select color (50:50 black:white)
            colori = ([0 randperm(ndots-1)]) < ceil(ndots/num_col);
            col = p.stim.rdk.dots_color(colori+1,:);
            
            % Create fresh batch of dots
            for d = 1:ndots
    
                enoughSpace = 0;
                while ~enoughSpace
                    
                    % create spatial position
                    r = sqrt(rand);
                    theta = 2*pi*rand();
                    
                    dots(dot_idx(d),1) = r.*cos(theta).*ap_radius(1);
                    dots(dot_idx(d),2) = r.*sin(theta).*ap_radius(1);
                    
                    % check spacing relative to the other dots
                    mydist = sqrt((dots(:,1) - dots(dot_idx(d),1)).^2 + (dots(:,2) - dots(dot_idx(d),2)).^2);
                    mydist(dot_idx(d)) = 0;
                    
                    if (min(mydist(:)) >= 0)
                        enoughSpace = 1;
                    end
                end
                
                % Initialize dot kill time for all dots (starting from
                % frame 1)
                dot_kill_time(dot_idx(d)) = 1 + rand().*p.stim.rdk.dots_lifetime; 
            end

            
            % Reset video frames
            frames = [];
            
            % loop over time frames
            for curr_frame = 1:num_frames

                % loop over dots
                for ii = 1:ndots

                    % update spatial position
                    dots(ii,:) = dots(ii,:) + dot_vel_deg_per_sec(ii,:) .* (1/p.stim.presentationrate_hz);

                    % Check whether the dot has reached its 'kill time'
                    % based on frames
                    if curr_frame > dot_kill_time(ii)
                        
                        % if so, we create dot with a new kill time and a
                        % new random x and y coordinate
                        enoughSpace = 0;
                        while ~enoughSpace

                            r = sqrt(rand);
                            theta = 2*pi*rand();
                            
                            dots(ii,1) = r.*cos(theta).*ap_radius(1);
                            dots(ii,2) = r.*sin(theta).*ap_radius(1);
                            
                            % check spacing relative to the other dots
                            mydist = sqrt((dots(:,1) - dots(ii,1)).^2 + (dots(:,2) - dots(ii,2)).^2);
                            mydist(ii) = 0;
                            
                            if (min(mydist(:)) >= 0)
                                enoughSpace = 1;
                            end
                        end

                        dot_kill_time(ii) = curr_frame + p.stim.rdk.dots_lifetime; % why not rand*lifetime?

                    end

                    % Check if dot has moved outside of the aperture and update and flip contrast if needed
                    if sqrt(dots(ii,1)^2+dots(ii,2)^2) > ap_radius(1)
                        % move to point symmetric location
                        dots(ii,:) = -dots(ii,:); 
                        col(ii,:) = setdiff(p.stim.rdk.dots_color,col(ii,:),'rows'); % flip contrast
                    end
                end % ndots loop

                
                % Create frame with gray background
                figure(fH);
                clf; hold all;
                ax = gca;
                ax.Units = 'pixels';
                r = rectangle(ax,'Position', [(ap_center -1.*ap_radius), ...
                    2.*ap_radius], ...
                    'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none'); %
                colormap gray; axis square; axis off; axis tight; axis manual; axis image
                set(gca, 'CLim',[0 1]);
                for mm = 1:size(dots,1)
                    drawcircle('Parent',ax,'Center',dots(mm,:),'Radius',p.stim.rdk.dots_size,...
                        'Color',col(mm,:), 'InteractionsAllowed', 'none', 'FaceAlpha', 1, 'LineWidth', 1);

                end

                %store location in case we want to use it
                stored_coh_dot_pos(:,:,curr_frame) = dots;

                f = getframe(ax);
                im = frame2im(f);

                frames = cat(4,frames,im);
                clear f im
            end % frame loop

            % debug video
            vcd_createStimVideo(frames, 1/p.stim.presentationrate_hz, ...
                fullfile(saveStimDir),sprintf('vcd_rdk_coh%03d_dir%02d_delta%02d',...
                p.stim.rdk.dots_coherence(cc)*100,bb,dd));           

            rdk{bb,cc,dd+1} = frames;
            dotlocs{bb,cc,dd+1} = stored_coh_dot_pos;
            
            info.dot_dir(counter) = curr_motdir_deg;
            info.dot_coh(counter) = p.stim.rdk.dots_coherence(cc);
            info.dot_dir_i(counter) = bb;
            info.dot_coh_i(counter) = cc;
            info.dot_pos{counter} = {stored_coh_dot_pos};
            
            if dd == 0
                info.motdir_deg_ref(counter) = 0;
                info.motdir_deg_ref_i(counter) = 0;
                info.unique_im(counter) = bb + ((cc-1)*length(p.stim.rdk.dots_direction));
            else
                info.motdir_deg_ref(counter) = p.stim.gabor.delta_from_ref(dd);
                info.motdir_deg_ref_i(counter) = dd;
                info.unique_im(counter) = NaN;
            end
            
            im_name = bb + ((cc-1)*length(p.stim.rdk.dots_direction));
            
            % Create binary circular alpha mask for rdk frames
            [XX,YY]       = meshgrid((1:size(frames,1))-(size(frames,1)/2),(1:size(frames,1))-(size(frames,1)/2)); % EK: 549 happens to be size of the RDK image.
            mask_radius   = ap_radius(1);
            circlemask    = (YY - ap_center(1)).^2 + (XX - ap_center(2)).^2 <= mask_radius.^2;
            mask          = double(circlemask);
            mask(mask==1) = 255;
            mask          = uint8(mask);
            
            % save intermediate stage in case matlab crashes 
            rdk_info = info(counter,:);
            save(fullfile(tmpDir, sprintf('%d_rdk_ori%d_coh%d_delta%d.mat', im_name, bb,cc,dd)),'frames','rdk_info','mask','stored_coh_dot_pos','-v7.3');

             clear frames
            counter = counter +1;
        end
    end
end


if p.store_imgs
    % in addition to storing separate conditions, also store ginormous file 
    saveDir = fileparts(fullfile(p.stim.rdk.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    save(fullfile(sprintf('%s_%s.mat',p.stim.rdk.stimfile,datestr(now,30))),'rdk','info','mask','dotlocs','-v7.3');
    
    saveDir = fileparts(fullfile(p.stim.rdk.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',p.stim.rdk.infofile,datestr(now,30))))
end

return

