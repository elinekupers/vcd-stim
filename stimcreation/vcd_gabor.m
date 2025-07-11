function [gabors, masks, info] = vcd_gabor(params, verbose, store_imgs)
% VCD function to create gabor images for experimental display.
% 
%  [gabors, masks, info] = vcd_gabor(params, verbose, store_imgs)
%
% Purpose:
%   Make gabors images and corresponding alpha transparency masks that crop
%   the top and bottom corners when presenting images.
%
%   We have 24 core gabor images sorted as follows:
%   * im 1:8   8 orientations for contrast level 1
%   * im 9:16  8 orientations for contrast level 2
%   * im 17:24 8 orientations for contrast level 3
%
%   4 quadrature phases (0,90,180,270 deg) will be evenly distributed
%   across these 24 core images. 2 spatial locations (left/right) will be
%   evenly distributed across 24 core images: uneven image nrs are left,
%   even image nrs are right.
%
%   We also create 4 test images for WM task: -16, -8 +8 +16 deg
%   orientation tilt offsets from core image orientation, with unique image
%   nrs 111-142. These WM test images inherit the spatial location,
%   and phase from their reference gabor, but all use the highest contrast 
%   level.
%
%   When params.stim.store_imgs = true, this function will store generated
%   object images as a single mat file in params.stim.obj.stimfile (e.g.:
%   fullfile(vcd_rootPath,'workspaces','stimuli',<disp_name>, ...
%   'gabor_<disp_name>_<date>.mat').
%
% INPUTS:
%  params                  : stim params struct (see vcd_setStimParams.m)
%    *** this function requires the following struct fields ***
%    bckgrnd_grayval       : (int) background gray value (128)
%    gabor.img_sz_pix      : (int) size of image support in pixels, needs
%                               to be even integer number!
%    gabor.gauss_std_pix   : (double) std of gaussian window
%    gabor.cycles_per_pix  : (double) spatial frequency (cycles per pixel)
%    gabor.contrast        : (double) gabor contrasts (fraction) (Michelson)
%    gabor.ori_deg         : (double) gabor orientation (deg), 0 deg = 12 o'clock
%    gabor.ph_deg          : (double) gabor phase (deg)
%    gabor.delta_from_ref  : (double) offset gabor orientation from core
%                               gabor orientation (deg)
%    gabor.unique_im_nrs_core    : (int) number for each unique core gabor
%    gabor.unique_im_nrs_wm_test : (int) number for each unique WM test gabor
%  verbose                 : (logical) show debug figures
%  store_imgs              : (logical) store debug figures as pngs 
% 
% OUTPUTS:
%  gabors      : (uint8) Gabor images used for VCD experiment, 5-dim array:
%                   height (BOLDscreen: 352 pixels) 
%                   x width (BOLDscreen: 352 pixels) 
%                   x 3 (rgb) 
%                   x 8 orientations (11.25, 33.75, 56.25, 78.75, 101.25,
%                   123.75, 146.25, 168.75 deg)
%                   x 3 contrasts (low: 5%, medium: 20%, high: 80%)
%                   x 5 deltas (1 core orientation (no delta), 
%                               + 4 wm deltas: [-16, -8, 8, 16] deg). 
%                   Note that dims gabors(:,:,:,:,[1,2],:) are all zeros
%                   as all wm test stimuli use the highest contrast level.
%  masks       : (uint8) alpha mask images used for VCD experiment, 4-dim array:
%                   height (BOLDscreen: 352 pixels; PProom: 256 pixels) 
%                   x width (BOLDscreen: 352 pixels; PProom: 256 pixels) 
%                   x 8 orientations (11.25, 33.75, 56.25, 78.75, 101.25,
%                   123.75, 146.25, 168.75 deg)
%                   x 3 contrasts (low: 5%, medium: 20%, high: 80%)
%                   x 5 deltas (1 core orientation (no delta) 
%                               + 4 wm deltas: [-16, -8, 8, 16] deg)
%                   Note that dims masks(:,:,:,[1,2],:) are all zeros
%                   as all wm test stimuli use the highest contrast level.
%  info        : 56x12 table with gabor stimulus information matching the 
%                   gabor array
%      unique_im     : (double) unique image nr for each Gabor: range 1-24,
%                       111-142 for wm test images.
%      stim_pos_i    : (double) stimulus position index. 1=left, 2=right
%      stim_pos      : (cell) stimulus position, same as stim_pos_i but
%                       human readable ({'left'} or {'right'})
%      orient_deg    : (double) tilt in degrees (0 = 12 o'clock).
%      orient_i      : (double) same as orient_deg but indexed 1
%                       (smallest: 11.25 deg) to 8 (largest: 168.75 deg)
%      contrast      : (double) contrast level (fraction of 1), ranging
%                       from low to high: 0.05, 0.20, 0.80.
%      contrast_i    : (double) same as 'contrast' but indexed 1 (lowest: 0.05)
%                       2 (medium: 0.2) or 3 (highest: 0.8).
%      phase_deg     : (double) gabor phase (deg), one of four quadrature
%                       phases (0, 90, 180, 270).
%      phase_i       : (double) same as phase_deg but indexed 1 (lowest: 0
%                       deg) to 4 (highest: 270).
%      delta_deg     : (double) Gabor orientation relative from
%                       corresponding core Gabor -16, -8, +8, +16 (deg)
%      delta_i       : (double) same as delta_deg but indexed 1
%                       (-16 deg), 2 (-8 deg), 3 (+8 deg), or 4 (+16 deg).
%      is_specialcore : (logical) whether the Gabor stimulus is part of the
%                       subselected stimuli used in imagery and long-term
%                       memory task.
%
% Written by Eline Kupers 2024/12
% 2025/4: updated 
% 2025/6: cleaned up header info


%% Check inputs

% Make sure the image has an even number of pixels, so we have center pix
if mod(params.stim.gabor.img_sz_pix,2)~=0
    error('[%s]: image support size does not have an even nr of pixels!', mfilename)
end

% If we have delta refs, add a "zero" orientation for the baseline (no delta change)
if ~isempty(params.stim.gabor.delta_from_ref)
    gbr_deltas = [0, params.stim.gabor.delta_from_ref];
else
    gbr_deltas = 0;
end

% Define all phases:
all_phases_deg = repelem(params.stim.gabor.ph_deg, length(params.stim.gabor.contrast)+length(params.stim.gabor.delta_from_ref))';
all_phases_deg = repmat(all_phases_deg,length(params.stim.gabor.ori_deg)/length(params.stim.gabor.ph_deg),1);
[~,phase_idx]  = ismember(all_phases_deg,params.stim.gabor.ph_deg);

% Define all orientations:
all_orient_deg = repelem(params.stim.gabor.ori_deg, length(params.stim.gabor.contrast)+length(params.stim.gabor.delta_from_ref))';
[~,orient_idx] = ismember(all_orient_deg, params.stim.gabor.ori_deg);

% Define all contrasts:
all_contrasts = [repmat(params.stim.gabor.contrast, 1, 1)'; ... 24 core stimuli (3 contrast levels * 8 orientations)
                 repmat(params.stim.gabor.contrast(end), 1, length(params.stim.gabor.delta_from_ref))'];     
all_contrasts = repmat(all_contrasts, length(params.stim.gabor.ori_deg),1);
[~,contrast_idx]  = ismember(all_contrasts,params.stim.gabor.contrast);

% Define all stim locations (we use NaN for deltas as they can be used for either left or right stimulus locations):
stimloc_names = {'left','right'};
all_stim_loc  = cell(size(contrast_idx));
stimloc_idx   = NaN(size(contrast_idx));
stimloc_idx0  = repmat([1,2],1,0.5*length(params.stim.gabor.ori_deg)*length(params.stim.gabor.contrast))';
st_idx        = [1:length(params.stim.gabor.contrast)]'+ (length(params.stim.gabor.ori_deg)-1)*[0:length(params.stim.gabor.ori_deg)-1];
st_idx        = st_idx(:);
stimloc_idx(st_idx(:)) = stimloc_idx0;
all_stim_loc(~isnan(stimloc_idx)) = stimloc_names(stimloc_idx0);

% Define all deltas:
all_deltas    = repmat([zeros(length(params.stim.gabor.contrast),1); params.stim.gabor.delta_from_ref'],length(params.stim.gabor.ori_deg),1);
[~,delta_idx] = ismember(all_deltas, gbr_deltas);
delta_idx     = delta_idx-1;

% get unique image nrs
img_im_nrs1  = reshape(params.stim.gabor.unique_im_nrs_core,8,[])';
img_im_nrs2  = reshape(params.stim.gabor.unique_im_nrs_wm_test,length(params.stim.gabor.delta_from_ref),[]);
img_im_nrs   = cat(1,img_im_nrs1,img_im_nrs2); % ori 1, contrast 1:3, delta = nr+dd, ori 2 ... etc
img_im_nrs   = img_im_nrs(:);
specialcore_bool = ismember(img_im_nrs,params.stim.gabor.unique_im_nrs_specialcore);
assert(isequal(sum(specialcore_bool),length(params.stim.gabor.unique_im_nrs_specialcore)));

% Create info table
info = table(img_im_nrs(:), ...
    all_stim_loc(:), ...
    stimloc_idx(:),...
    all_orient_deg(:), ...
    orient_idx(:), ...
    all_contrasts(:),...
    contrast_idx(:), ...
    all_phases_deg(:), ...
    phase_idx, ...
    all_deltas(:), ...
    delta_idx, ...
    specialcore_bool);

% add column names
info.Properties.VariableNames = {'unique_im','stim_pos','stim_pos_i','orient_deg',...
    'orient_i', 'contrast','contrast_i','phase_deg','phase_i','delta_deg','delta_i', 'is_specialcore'};

% Preallocate space
gabors = uint8(zeros(params.stim.gabor.img_sz_pix,params.stim.gabor.img_sz_pix,3, ...
    length(params.stim.gabor.ori_deg), ...
    length(params.stim.gabor.contrast), ...
    length(gbr_deltas)));

masks = uint8(zeros(params.stim.gabor.img_sz_pix,params.stim.gabor.img_sz_pix, ...
    length(params.stim.gabor.ori_deg), ...
    length(params.stim.gabor.contrast), ...
    length(gbr_deltas)));

% Loop over orientation
for bb = 1:length(orient_idx)
    
    curr_angle_deg = all_orient_deg(bb);
    curr_delta = all_deltas(bb);
    curr_phase = all_phases_deg(bb);
    curr_contrast = all_contrasts(bb);
    
    final_angle_deg = curr_angle_deg+curr_delta;
    
    % Print corresponding gabor params
    fprintf('\nImage idx',bb)
    if curr_delta == 0
        fprintf('\nAngle: %2.2f deg', curr_angle_deg)
    else
        fprintf('\nAngle: %2.2f + %2.2f delta deg', ...
            curr_angle_deg, all_deltas(bb))
    end
    fprintf('\nPhase %d deg', curr_phase)
    fprintf('\nContrast %1.2f', curr_contrast)
    
    % TRICKY STUFF: DO we subtract 90 deg to ensure 0 deg desired
    % orientation is 12 o'clock in x,y-pixel space? NO?
    %     curr_angle_deg = curr_angle_deg-90;
    
    % Wrap around 360
    final_angle_deg = mod(final_angle_deg,360);

    
    % Create gabor!
    img_c = vcd_create_gabor(...
        params.stim.gabor.img_sz_pix, ...
        params.stim.gabor.gauss_std_pix, ...
        params.stim.gabor.cycles_per_pix,...
        final_angle_deg,...
        curr_phase,...
        curr_contrast,...
        params.stim.bckgrnd_grayval);
    
    % create transparency mask to crop out corners of rectangular support
    support_x  = params.stim.gabor.img_sz_pix;
    alpha_mask = makecircleimage(support_x, params.stim.gabor.img_sz_pix/2);
    alpha_mask = alpha_mask*255;
    
    % Store image
    gabors(:,:,:,orient_idx(bb),contrast_idx(bb),delta_idx(bb)+1) = uint8(repmat(img_c,[1 1 3]));
    masks(:,:,orient_idx(bb),contrast_idx(bb),delta_idx(bb)+1)    = uint8(alpha_mask);
    
    if store_imgs
        saveFigDir = fullfile(params.saveFigsFolder,'gabor');
        if ~exist(saveFigDir,'dir'); mkdir(saveFigDir); end
        filename = sprintf('%04d_vcd_gabor_ori%02d_c%02d_delta%02d.png',img_im_nrs(bb),orient_idx(bb),contrast_idx(bb),delta_idx(bb));
        imwrite(gabors(:,:,:,orient_idx(bb),contrast_idx(bb),delta_idx(bb)+1), fullfile(saveFigDir,filename));
    end
    
    % Clean up
    clear img_c
    
end % num orientations

fprintf('Done!\n')

%% Store images and info table
if store_imgs
    fprintf('[%s]:Storing images and info..\n',mfilename)
    saveStimMatFileDir = fileparts(fullfile(params.stim.gabor.stimfile));
    saveStimMatFileDir = fullfile(saveStimMatFileDir, 'gabor');
    if ~exist(saveStimMatFileDir,'dir'), mkdir(saveStimMatFileDir); end
    save(fullfile(sprintf('%s_%s.mat',params.stim.gabor.stimfile,datestr(now,30))),'gabors','masks','info','-v7.3');
    
    saveStimMatFileDir = fileparts(fullfile(params.stim.gabor.infofile));
    if ~exist(saveStimMatFileDir,'dir'), mkdir(saveStimMatFileDir); end
    writetable(info, fullfile(sprintf('%s_%s.csv',params.stim.gabor.infofile,datestr(now,30))))
end

%% Visualize Gabors if requested

if verbose
    makeprettyfigures;
    
    % create folders for visual checks
    if store_imgs
        saveFigDir2 = fullfile(saveFigDir,'visual_checks','gabor_hist');
        if ~exist(saveFigDir2,'dir'); mkdir(saveFigDir2); end
    end
    
    fH1 = figure(1); set(fH1, 'Position', [1 400 750 578], 'color','w'); % histogram
    
    counter = 1;
    
    % loop over deltas
    for dd = 1:size(gabors,6)
        if dd==1
            dlta = 0;  
        else
            dlta = params.stim.gabor.delta_from_ref(dd-1); 
        end
        
        % loop over contrasts
        for ii = 1:size(gabors,5)
            
            if dd > 1
                ii = size(gabors,5);
            end
            
            % loop over orientations
            for jj=1:size(gabors,4)
                
                % check image nr
                im_idx = (info.orient_deg==params.stim.gabor.ori_deg(jj) & ...
                          info.contrast==params.stim.gabor.contrast(ii) & ...
                          info.delta_deg==dlta);
                
                im_nr = info.unique_im(im_idx);
                
                % Plot pix histogram
                clf
                histogram(gabors(:,:,:,jj,ii,dd), 'NumBins', 30);
                title(sprintf('ori:%3.2f deg c:%1.2f ph:%3.0f delta:%02d', ...
                    params.stim.gabor.ori_deg(jj),...
                    params.stim.gabor.contrast(ii),...
                    params.stim.gabor.ph_deg(mod(ii-1,4)+1),...
                    dlta));
                box off; axis square; xlim([-10 260]); ylim([0,3.5].*10^5)
                drawnow;
                axis square
                set(gca,'XTick',[0 64 128 190 255])
                set(gca,'YTick',[1:3].*10^5)
                if store_imgs
                    filename = sprintf('%04d_vcd_gabor_ori%02d_c%02d_delta%02d_hist.png',im_nr,jj,ii,dd-1);
                    print(fH1,'-dpng','-r300',fullfile(saveFigDir2,filename));
                end
                counter = counter+1;
            end
        end
    end
end

return

