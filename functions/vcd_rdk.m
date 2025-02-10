function [rdk,info,p] = vcd_rdk(p, disp)
% VCD function:
%  [rdk,info,p] = vcd_rdk(p, disp)
%
% Purpose:
%   Create a random dot motion kinetograms for  experimental display. 
%   For rdk stimulus parameters see vcd_getStimParams.m
% 
% INPUTS:
%  p            : rdk params
%  disp         : monitor display params (struct) pixels per degree in field of view (pix)
%
% OUTPUTS:
%   rdk     : RDK images, 8x3 cell for each motion direction and coherence
%               where each cell contains w (pixels) x h (pixels) x 3 x frames
%   info    : table with rdk stimulus information (motion direction, coherence,
%               location of dots for each individual frame), matching the rdk array
%   p       : updated params struct
%%
% Written by Eline Kupers 2024/12
%

%% Set Display & RDK parameters

p.stim.rdk.dots_aperture = floor([0 0 p.stim.rdk.img_sz_deg./2 p.stim.rdk.img_sz_deg./2].*disp.ppd); % [x y w h] in pixels
p.stim.rdk.dots_angle = pi*p.stim.rdk.dots_direction/180;             % convert deg 2 rad

%define an elliptic aperture in screen
ap_center = p.stim.rdk.dots_aperture(1:2);
ap_radius = p.stim.rdk.dots_aperture(3:4);

ndots = min(p.stim.rdk.max_dots_per_frame, ...
    round(p.stim.rdk.dots_density .* (p.stim.rdk.dots_aperture(:,3).*p.stim.rdk.dots_aperture(:,4)) / disp.refresh_hz));

rdk = cell(length(p.stim.rdk.dots_direction),length(p.stim.rdk.dots_coherence));
n_unique_conds = size(rdk,1)*size(rdk,2);
info = table(NaN(n_unique_conds,1), NaN(n_unique_conds,1), cell(n_unique_conds,1));
info.Properties.VariableNames = {'dot_dir','dot_coh','dot_pos'};

% RNG seed
rseed = [1000 2010];
num_col = size(p.stim.rdk.dots_color,1);
counter = 1;

for cc = 1:length(p.stim.rdk.dots_coherence)
    
    for dd = 1:length(p.stim.rdk.dots_direction)
        %find the xy displacement of coherent
        dxdy = repmat(p.stim.rdk.dots_speed * p.stim.rdk.dots_interval/ disp.refresh_hz *...
            [cos(p.stim.rdk.dots_direction(dd)) -sin(p.stim.rdk.dots_direction(dd))], ndots, 1) * disp.ppd;
        d_ppd = repmat(ap_radius, ndots, 1);
        dot_pos = (rand(ndots,2,p.stim.rdk.dots_interval)-0.5)*2;
        
        % Reset rng
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',prod(rseed)));
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
        
        for jj = 1 : p.stim.rdk.dots_interval
            dot_pos(:,:,jj) = dot_pos(:,:,jj) .* d_ppd;
        end
        
        % reset frames
        frames = [];
        
        % store dot position
        all_pos = zeros(size(dot_pos,1),size(dot_pos,2),p.stim.rdk.duration*p.stim.rdk.dots_interval);
        
        % Draw dots
        max_t = p.stim.rdk.duration; %p.duration/disp.refresh_hz; % start_t = GetSecs;
        t = 0;
        
        while 1
            %     cur_t = GetSecs - start_t;
            t = t + 1;
            
            if t >= max_t
                break;
            else
                
                colori = ([0 randperm(ndots-1)]) < ceil(ndots/num_col);
                dot_pos_col = p.stim.rdk.dots_color(colori+1,:);
                
                for loopi = 1: p.stim.rdk.dots_interval % OG code states length(scr_rfsh/ dots_interval)
                    
                    % update dots positions and draw them on the
                    % find the index of coherently moving dots in this
                    L = rand(ndots,1) < p.stim.rdk.dots_coherence(cc);
                    
                    % move the coherent
                    dot_pos(L,:,loopi) = dot_pos(L,:,loopi) + dxdy(L,:);
                    
                    % replace the other
                    dot_pos(~L,:,loopi) = (rand(sum(~L),2)-0.5)*2 .* d_ppd(~L,:);
                    
                    % wrap
                    L = dot_pos(:,1,loopi) > d_ppd(:,1);
                    dot_pos(L,1,loopi) = dot_pos(L,1,loopi) - 2*d_ppd(L,1);
                    L = dot_pos(:,1,loopi) < -d_ppd(:,1);
                    dot_pos(L,1,loopi) = 2*d_ppd(L,1) - dot_pos(L,1,loopi);
                    L = dot_pos(:,2,loopi) > d_ppd(:,2);
                    dot_pos(L,2,loopi) = dot_pos(L,2,loopi) - 2*d_ppd(L,2);
                    L = dot_pos(:,2,loopi) < -d_ppd(:,2);
                    dot_pos(L,2,loopi) = 2*d_ppd(L,2) - dot_pos(L,2,loopi);
                    
                    % Find the dots that will be shown in the aperture. note that
                    %is calculated relative to the center of the aperture.
                    L = isInsideAperture(dot_pos(:,:,loopi), ap_center, ap_radius-1);
                    
                    pos_idx = find(L);
                    
                    fH = figure(1); clf; hold all;
                    r = rectangle('Position',  [(ap_center -1.*ap_radius), ...
                        2.*ap_radius], ...
                        'FaceColor', [127 127 127]./255, 'EdgeColor', 'none');
%                     xlim(ap_center(1) + [-1,1].*ap_radius)
%                     ylim(ap_center(2) + [-1,1].*ap_radius)
                    colormap gray; axis square; axis off;
                    
                    %round dot_pos and transpose it because Screen wants positions in row
                    pos = round(dot_pos(L,:,loopi));
                    col = dot_pos_col(L,:);
                    
                    % Create a blank image for each frame
                    %             img = zeros(size(dot_pos,1),size(dot_pos,2));
                    
                    for ii = 1:size(pos,1)
                        h(ii) = drawcircle('Center',pos(ii,:),'Radius',p.stim.rdk.dots_size,...
                            'Color',col(ii,:), 'InteractionsAllowed', 'none', 'FaceAlpha', 1, 'LineWidth', 1);
                        all_pos(pos_idx(ii),:,t) = pos(ii,:);
                    end
                    g = gca;
                    f = getframe(g);
                    frames = cat(4,frames,f.cdata);
                    clear f
                    %             %draw on the
                    %             if any(isnan(prod(pos,1))==0)
                    %                 Screen('DrawDots', screen_struct.cur_window, pos, dots_struct.dot_size, dots_struct.dot_color', AP.spec.center);
                    %
                    %                 %update the loop
                    %                 loopi = loopi + 1;
                    %                 if loopi > dots_struct.interval
                    %                     loopi = 1;
                    %                 end
                    %
                    %
                    %                 %flip the screen to make things
                    %                 Screen('Flip', screen_struct.cur_window);
                    %             end
                    
                    %             %flip the screen to clean it
                    %             Screen('Flip', screen_struct.cur_window);
                    
                    %         end
                end
            end
        end
        rdk{dd,cc} = frames;

        info.dot_dir(counter) = p.stim.rdk.dots_direction(dd);
        info.dot_coh(counter) = p.stim.rdk.dots_coherence(cc);
        info.dot_pos{counter} = {all_pos};
        clear frames all_pos
        counter = counter +1;
    end
end


if p.store_imgs
    saveDir = fileparts(fullfile(p.stim.rdk.stimfile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    tmp = strsplit(p.stim.rdk.stimfile,'.mat');
    save(fullfile(sprintf('%s_%s.mat',tmp{1},datestr(now,30))),'rdk','info','-v7.3');
    
    saveDir = fileparts(fullfile(p.stim.rdk.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    tmp = strsplit(p.stim.rdk.infofile,'.csv');
    writetable(info, fullfile(sprintf('%s_%s.csv',tmp{1},datestr(now,30))))
end

return

