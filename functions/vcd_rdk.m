function [rdk,info,p] = vcd_rdk(p)
% VCD function:
%  [rdk,info,p] = vcd_rdk(p)
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
%   rdk     : RDK images, 8x3 cell for each motion direction and coherence
%               where each cell contains w (pixels) x h (pixels) x 3 x frames
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

p.stim.rdk.dots_aperture = floor([0 0 p.stim.rdk.img_sz_deg./2 p.stim.rdk.img_sz_deg./2].*p.disp.ppd); % [x y w h] in pixels
p.stim.rdk.dots_angle = pi*p.stim.rdk.dots_direction/180;             % convert deg 2 rad

%define an elliptic aperture in screen
ap_center = p.stim.rdk.dots_aperture(1:2);
ap_radius = p.stim.rdk.dots_aperture(3:4);

ndots = min(p.stim.rdk.max_dots_per_frame, ...
    round(p.stim.rdk.dots_density .* (p.stim.rdk.dots_aperture(:,3).*p.stim.rdk.dots_aperture(:,4)) / p.disp.refresh_hz));

if ~isempty(p.stim.rdk.delta_from_ref)
    rdk_motdir_ref = [0:length(p.stim.rdk.delta_from_ref)];
else
    rdk_motdir_ref = 0;
end

% Preallocate space
rdk = cell(length(p.stim.rdk.dots_direction),length(p.stim.rdk.dots_coherence),length(rdk_motdir_ref));
n_unique_conds = size(rdk,1)*size(rdk,2)*size(rdk,3);

% Preallocate info table
info = table(NaN(n_unique_conds,1), NaN(n_unique_conds,1), cell(n_unique_conds,1),NaN(n_unique_conds,1),NaN(n_unique_conds,1));
info.Properties.VariableNames = {'dot_dir','dot_coh','dot_pos','motdir_deg_ref','unique_im'};

% RNG seed
rseed = [1000 2010];
num_col = size(p.stim.rdk.dots_color,1);
counter = 1;

tmpDir = fullfile(vcd_rootPath, 'workspaces','stimuli',['rdk_' datestr(now,'yyyymmdd')]);
if ~exist('tmpDir','dir'), mkdir(tmpDir); end

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
        
            if curr_motdir_deg < 0
                curr_motdir_deg = 360+curr_motdir_deg;
            end
        
            %find the xy displacement of coherent
            dxdy = repmat(p.stim.rdk.dots_speed * p.stim.rdk.dots_interval/ p.disp.refresh_hz *... %% 
                [cos(curr_motdir_deg) -sin(curr_motdir_deg)], ndots, 1) * p.disp.ppd;
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
            max_t = p.stim.rdk.duration; %p.duration/p.disp.refresh_hz; % start_t = GetSecs;
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
                        
                        clf; hold all;
                        ax = gca;
                        ax.Units = 'pixels';
                        r = rectangle(ax,'Position', [(ap_center -1.*ap_radius), ...
                            2.*ap_radius], ...
                            'FaceColor', [repmat(p.stim.bckgrnd_grayval,1,3)]./255, 'EdgeColor', 'none');
%                         xl = ap_center(1) + [-1,1].*ap_radius + [-1 1];
%                         yl = ap_center(2) + [-1,1].*ap_radius + [-1 1];
%                         xlim(xl)
%                         ylim(yl)
                        colormap gray; axis square; axis off; axis tight; axis manual; axis image
                        
                        %round dot_pos and transpose it because Screen wants positions in row
                        pos = round(dot_pos(L,:,loopi));
                        col = dot_pos_col(L,:);
                        
                        % Create a blank image for each frame
                        %             img = zeros(size(dot_pos,1),size(dot_pos,2));
                        
                        
%                         rect = [-axti(1), -axti(2), axpos(3)+axti(1)+axti(3), axpos(4)+axti(2)+axti(4)];

                        
                        for ii = 1:size(pos,1)
                            drawcircle('Parent',ax,'Center',pos(ii,:),'Radius',p.stim.rdk.dots_size,...
                                'Color',col(ii,:), 'InteractionsAllowed', 'none', 'FaceAlpha', 1, 'LineWidth', 1);
                            all_pos(pos_idx(ii),:,t) = pos(ii,:);
                        end
                        
%                         f2 = hardcopy(fH,'-dzbuffer','-r0');
                        
                        f = getframe(ax);
                        im = frame2im(f);
                        f_bw = rgb2gray(im);
%                         scale_factor = (2*ap_radius)./size(f_bw);
%                         f_bw_rz = imresize(f_bw,scale_factor(1));
                        
                        frames = cat(3,frames,f_bw);
                        clear f im f_bw scale_factor f_bw_rz
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
%             rdk{bb,cc,dd+1} = frames;
            
            info.dot_dir(counter) = curr_motdir_deg;
            info.dot_coh(counter) = p.stim.rdk.dots_coherence(cc);
            info.dot_dir_i(counter) = bb;
            info.dot_coh_i(counter) = cc;
            info.dot_pos{counter} = {all_pos};
            
            if dd == 0,
                info.motdir_deg_ref(counter) = 0;
                info.motdir_deg_ref_i(counter) = 0;
                info.unique_im(counter) = bb + ((cc-1)*length(p.stim.rdk.dots_direction));
            else
                info.motdir_deg_ref(counter) = p.stim.gabor.delta_from_ref(dd);
                info.motdir_deg_ref_i(counter) = dd;
                info.unique_im(counter) = NaN;
            end
            
            im_name = bb + ((cc-1)*length(p.stim.rdk.dots_direction));
            % save intermediate stage in case matlab crashes 
            rdk_info = info(counter,:);
            save(fullfile(tmpDir, sprintf('%d_rdk_ori%d_coh%d_delta%d.mat', im_name, bb,cc,dd)),'frames','rdk_info','-v7.3');

            clear frames all_pos
            counter = counter +1;

        end
    end
end





if p.store_imgs
    % we store separate conditions, otherwise the file gets ginormous
%     saveDir = fileparts(fullfile(p.stim.rdk.stimfile));
%     if ~exist(saveDir,'dir'), mkdir(saveDir); end
%     save(fullfile(sprintf('%s_%s.mat',p.stim.rdk.stimfile,datestr(now,30))),'rdk','info','-v7.3');
    
    saveDir = fileparts(fullfile(p.stim.rdk.infofile));
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',p.stim.rdk.infofile,datestr(now,30))))
end

return

