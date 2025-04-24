%% s_createStim.m
%
% Stand-alone script to create and store the stimuli shown in VCD core 
% experiment.

%% %%%%%%%%%%%%%%%%%%%
%%%%%% PARAMETERS %%%% 
%%%%%%%%%%%%%%%%%%%%%%

p.verbose        = true; % visualize stimuli or not
p.store_imgs   = true; % store visualization figures

% Get display params
dispname = '7TAS_BOLDSCREEN32'; % Choose from: '7TAS_BOLDSCREEN32';'KKOFFICE_AOCQ3277';'EKHOME_ASUSVE247';'PPROOM_EIZOFLEXSCAN'
p.disp   = vcd_getDisplayParams(dispname);


saveFigsFolder = fullfile(vcd_rootPath,'figs',dispname); % where to store visualization figures
if ~exist(saveFigsFolder,'dir'); mkdir(saveFigsFolder); end

% Get stimulus parameters
p.load_params                 = true; % if false, re-create params.
                                      % if true, we load stored mat file
p.store_params                = true; % if false, we don't store params. 
                                      % if true, we store mat file in fullfile(vcd_rootPath,'workspaces','info')
p.overwrite_randomized_params = false; % if false, we will use now hardcoded (once probabilistic) stimulus/experimental design params. Use this for reproducible results.
                                       % if true, we do redefine all stimulus/experimental design params, including probabilistic params

% Reset random number generator with arbitrary number (based on system clock) 
rand('seed', sum(100*clock));
p.rng.rand_seed  = rng; % store 
randn('seed', sum(100*clock)); % early versions of Matlab use different generators for rand and randn.
p.rng.randn_seed  = rng; % store 

%% Define/Load stimulus params 
% !!WARNING!! There is a randomization component involved in creating some
% stimuli (e.g., orientation of gabor stimuli or dot locations). If you
% don't want this, this leave the fifth argument:
% "overwrite_randomized_params" empty (default is set to false) or set to
% false.
%
% If you do want regenerate probabilistic params, set the fifth argument to
% true and some stimulus values will change.
%
% Input 1: Stimulus class, choose from 'gabor','rdk','dot','obj','ns','all' (default is 'all')
% Input 2: Display params struct (see vcd_getDisplayParams.m)
% Input 3: Load prior stored parameters or not. 
% Input 4: Store generated parameters or not. 
% Input 5: Overwrite stored parameters and regenerate probabilistic params
p.stim   = vcd_getStimParams('stim_class','all', ...
                             'disp_name',p.disp.name, ...
                             'load_params', p.load_params,...
                             'store_params', p.store_params, ...
                             'overwrite_randomized_params', p.overwrite_randomized_params); 

%% %%%%%%%%%%%%%%%%
%%%%%% STIMULI %%%% 
%%%%%%%%%%%%%%%%%%%

%% Create background
% Input 1: params struct 
% Input 2: type: 'puzzle', 'dotring', or 'comb'
% Input 3: rim width, choose from: 'skinny' (no blank space between noise
%           and stim) or 'fat' (2 dva space between noise and stim) 
% Input 4: number of unique noise images
%
% OUTPUT: bckgrnd_im uint8 images:
%   height (1080 pixels) x width (1920 pixels) x 3 (rbg)

gaptype     = 'comb';
borderwidth = 'fat';
num         = 1;
bckgrnd_im  = vcd_pinknoisebackground(p, ...
                                     'gaptype',gaptype, ...
                                     'borderwidth', borderwidth,...
                                     'num',num); 

%% Fixation dot images
% fix_im (uint8) is a 5D array containing unique fixation dot images: 
%   height (24 pixels) x width (24 pixels) x 3 (rgb) x 5 luminance levels x 5 dot rims types (thin, thick, thick-left, thick-right, thick-both)
% fix_mask (uint8) is a 4D array with alpha transparency masks: 
%   height (24 pixels) x width (24 pixels) x 2 dot rims types (thin, thick)
[fix_im, fix_mask, fix_info] = vcd_fixationDot(p);

%% Gabor images
% gabor (uint8) is a uint8 6D array containing unique gabor images: 
%   height (354 pixels) x width (354 pixels) x 3 (rgb) x 24 unique images (8 ori x 3 contrasts) x 5 tilt offsets (0 + -15, -5, +5, +15 deg)
% gbr_masks (uint8) is a uint8 4D array with alpha transparency masks: 
%   height (354 pixels) x width (354 pixels) x 24 unique images x 5 tilt offsets
[gabor, gbr_masks, gbr_info] = vcd_gabor(p);

%% RDK movies
% rdk (uint8) is a 8 x 3 x 5 cell array: 
%   8 motion directions x 3 coherence levels x 5 motion direction offsets (0, -15, -5, +5, +15 deg). 
% Each cell contains a movie (uint8):
%   height (549 pixels) x width (549 pixels) x 3 (rgb) x 30 frames (total of 2 s, 33 ms per frame).
% rdk_masks (uint8) is a uint8 matrix with a single alpha transparency mask: 
%   height (549 pixels) x width (549 pixels)
[rdk, rdk_masks, rdk_info] = vcd_rdk(p);

%% Single dot
% * single_dot (uint8) is a single matrix: 
%   height (94 pixels) x width (94 pixels) x 3 (rgb).
% * dot_masks (uint8) is a alpha transparency mask:
%   height (94 pixels) x width (94 pixels)
[single_dot, dot_masks, dot_info] = vcd_singledot(p);

%% Objects
% Output images:
% * objects (uint8) is a 5D array:
%   height (1024 pixels) x width (1024 pixels) x 3 (rgb) x 16 object categories (subordinate level) x 4 rotation offsets (0, -8, -4, +4, +8 deg). 
% * cobj_masks (uint8) is a 4D array containing alpha transparency mask:
%   height (1024 pixels) x width (1024 pixels) x 16 object categories (subordinate level) x 5 rotation offsets (0, -8, -4, +4, +8 deg). 
[objects, obj_masks, ~, obj_info] = vcd_objects(p);

%% Natural scenes
% Output images:
% * scenes (uint8) is a 6D array: 
%   height (BOLDscreen 743 pixels) x width (BOLDscreen 743 pixels) x 3 (rgb) x 5 superordinate categories x 2 scene locations (indoor/outdoor) x 3 obj locations (left/middle/right)
% * lure_im (uint8) is a 7D array:
%   height (BOLDscreen 743 pixels) x width (BOLDscreen 743 pixels) x 3 (rgb) x 5 superordinate categories x 2 scene locations (indoor/outdoor) x 3 obj locations (left/middle/right) x 4 lure types:
%       1:very similar/difficult ... 4:least similar/easy
% * cblind_im (uint8) is a 7D array:
%   height (743 pixels) x width (743 pixels) x 3 (rgb) x 5 superordinate categories x 2 scene locations (indoor/outdoor) x 3 obj locations (left/middle/right) x 4 change types:
%       1:easy add -- scene is altered by adding something big/obvious 
%       2:hard add -- scene is altered by adding something small/subtle
%       3:easy remove -- scene is altered by removing something big/obvious
%       4:hard remove -- scene is altered by removing something small/subtle
% Note: There are no alpha masks for this stimulus class.
[scenes, lure_im, cblind_im, ~, ns_info] = vcd_naturalscenes(p);


%% %%%%%%%%%%%%%%%%
%%%% Visualize %%%% 
%%%%%%%%%%%%%%%%%%%

if p.verbose
    makeprettyfigures;
    
    %% Background
    saveDir = fullfile(vcd_rootPath,'figs',dispname,'background','visual_check');
    if ~exist(saveDir,'dir'); mkdir(saveDir); end
    
    fH = figure(100);
    set(fH, 'Position', [0 0 p.disp.w_pix p.disp.h_pix], 'color','w')

    for jj = 1:size(bckgrnd_im,4)
        clf;
        imagesc(bckgrnd_im(:,:,jj)); colormap gray
        axis image
        title(sprintf('Background im %s %s %03d',gaptype, borderwidth, jj))
        set(gca, 'TickDir','out', 'LineWidth',2,'FontSize',20)
        drawnow;
        if p.store_imgs        
            filename = sprintf('vcd_background_%s%s%03d.png',gaptype, borderwidth, jj);
            imwrite(bckgrnd_im(:,:,:,jj), fullfile(vcd_rootPath,'figs',dispname,'background',filename));
            print(fH,'-dpng','-r300',fullfile(saveDir,filename));
        end
    end

    %% Fixation dot
    saveDir = fullfile(vcd_rootPath,'figs',dispname,'fix','visual_check');
    if ~exist(saveDir,'dir'); mkdir(saveDir); end
 
    fH = figure(101); clf;
    set(fH, 'Position', [1   400   750   578], 'color','w')
    counter = 1;
    for ii = 1:size(fix_im,4)
        for jj = 1:size(fix_im,5)
            subplot(5,5,counter)
            imagesc(squeeze(fix_im(:,:,:,ii,jj)), 'AlphaData',p.stim.fix.dotopacity)
            axis square
            axis off
            set(gca,'CLim',[1 255])
            title(sprintf('%01d, %01d',ii, jj))
            counter = counter+1;
            if p.store_img
            imwrite(fix_im(:,:,:,ii,jj), fullfile(vcd_rootPath,'figs',dispname,'fix',sprintf('vcd_fixdots_w_alphamask_%01d_%01d.png', ii, jj)));
            end
        end
    end
    if p.store_imgs
        filename = sprintf('vcd_fixdots_w_alphamask.png');
        print(fH,'-dpng','-r300',fullfile(saveDir,filename));
    end
    
    fH = figure(101); clf;
    set(fH, 'Position', [1   400   750   578], 'color','w')
    counter = 1;
    for ii = 1:size(fix_im,4)
        for jj = 1:size(fix_im,5)
            subplot(5,5,counter)
            imagesc(squeeze(fix_im(:,:,:,ii,jj)))
            axis square
            axis off
            set(gca,'CLim',[1 255])
            title(sprintf('%01d, %01d',ii, jj))
            counter = counter+1;
            if p.store_img
                imwrite(fix_im(:,:,:,ii,jj), fullfile(vcd_rootPath,'figs',dispname,'fix',sprintf('vcd_fixdots_%01d_%01d.png', ii, jj)));
            end
        end
    end
    if p.store_imgs
        filename = sprintf('vcd_fixdots.png');
        print(fH,'-dpng','-r300',fullfile(saveDir,filename));
    end
    
    %% Gabor
    saveDir1 = fullfile(vcd_rootPath,'figs',dispname,'gabor','visual_checks','gabor_hist');
    saveDir2 = fullfile(vcd_rootPath,'figs',dispname,'gabor','visual_checks','gabor_im');
    if ~exist(saveDir1,'dir'); mkdir(saveDir1); end
    if ~exist(saveDir2,'dir'); mkdir(saveDir2); end
    
    fH1 = figure(1); fH2 = figure(2); 
    set(fH1, 'Position', [0 0 1024 1080], 'color','w')
    set(fH2, 'Position', [1   400   750   578], 'color','w')
    if p.verbose
        counter = 1;
        for dd = 1:size(gabor,6)
            if dd==1, dlta = 0; else, dlta = p.stim.gabor.delta_from_ref(dd-1); end
            for ii = 1:size(gabor,5)
                for jj=1:size(gabor,4)
                    figure(fH1); clf;
                    imagesc(gabor(:,:,:,jj,ii,dd));
                    
                    colormap gray; set(gca, 'CLim',[1 255])
                    title(sprintf('ori:%3.2f deg c:%1.2f ph:%3.0f delta:%02d', ...
                        p.stim.gabor.ori_deg(jj),...
                        p.stim.gabor.contrast(ii),...
                        p.stim.gabor.ph_deg(mod(ii-1,4)+1),...
                        dlta));
                    cb = colorbar; cb.Ticks = [1, 50, 100, 150, 200, 250];
                    set(gca, 'FontSize',20, 'TickDir','out', 'LineWidth', 2); axis image;
                    drawnow;
                    axis image
                    if p.store_imgs
                        filename = sprintf('vcd_gabor_ori%02d_c%02d_delta%02d.png',jj,ii,dd);
                        print(fH1,'-dpng','-r300',fullfile(saveDir1,filename));
                        imwrite(gabor(:,:,:,jj,ii,dd), fullfile(vcd_rootPath,'figs',dispname,'gabor',filename));
                    end
                    
                    % Plot pix histogram
                    figure(fH2); clf
                    histogram(gabor(:,:,:,ii,dd), 'NumBins', 30);
                    title(sprintf('ori:%3.2f deg c:%1.2f ph:%3.0f delta:%02d', ...
                        p.stim.gabor.ori_deg(jj),...
                        p.stim.gabor.contrast(ii),...
                        p.stim.gabor.ph_deg(mod(ii-1,4)+1),...
                        dlta));
                    box off; axis square; xlim([-10 260]); ylim([0,3.5].*10^5)
                    drawnow;
                    axis square
                    set(gca,'XTick',[0 64 128 190 255])
                    set(gca,'YTick',[1:3].*10^5)
                    if p.store_imgs
                        filename = sprintf('vcd_gabor_ori%02d_c%02d_delta%02d_hist.png',jj,ii,dd);
                        print(fH2,'-dpng','-r300',fullfile(saveDir2,filename));
                    end
                    counter = counter+1;
                end
            end
        end
    end

    %% RDK
    saveDir = fullfile(vcd_rootPath,'figs',dispname,'rdk','visual_checks');
    if ~exist(saveDir,'dir'); mkdir(saveDir); end
    
    fH1 = figure(1); 
    set(fH1, 'Position', [0 0 1024 1080], 'color','w')
    deltasToPlot = [0, p.stim.rdk.delta_from_ref];
    deltalabels = cellfun(@num2str, (num2cell(p.stim.rdk.delta_from_ref)),'UniformOutput', false);
    for ll = 1:length(deltalabels), 
        if strfind(deltalabels{ll},'-')
            deltalabels{ll} = strrep(deltalabels{ll},'-','min'); 
        else
            deltalabels{ll} = ['plus' deltalabels{ll}]; 
        end
    end
    deltalabels = [{'00'}, deltalabels(:)'];
    
    cohlabels = cellfun(@num2str, (num2cell(p.stim.rdk.dots_coherence*100)),'UniformOutput', false);
    cohlabels = cellfun(@(x) strrep(x,'.','pt'), cohlabels,'UniformOutput',false);
    for dl = 1:length(deltalabels)
        for cc = 1:length(p.stim.rdk.dots_coherence)
            for dd = 1:length(p.stim.rdk.dots_direction)
                
%                 % create movie
                vcd_createStimVideo(rdk{dd,cc,dl}, 1/p.stim.presentationrate_hz, ...
                    fullfile(vcd_rootPath,'figs',dispname,'rdk'), ...
                    sprintf('vcd_rdk_coh%s_dir%02d_delta%s',...
                    cohlabels{cc},p.stim.rdk.dots_direction(dd),deltalabels{dl}));

                % plot dot position & motion vector
                dot_pos = dotlocs{dd,cc, dl};
                dxdy    = dot_pos(:,:,2:end)-dot_pos(:,:,1:end-1);
                dxdy    = cat(3,zeros(size(dot_pos,1),size(dot_pos,2)),dxdy);
                if p.store_imgs
                    saveDir = fullfile(vcd_rootPath,'figs',dispname,'rdk','motion_vecs',sprintf('dir%02d',dd));
                    if ~exist(saveDir,'dir'); mkdir(saveDir); end
                end
                for ii = 1:size(dot_pos,3) % loop over time
                    cla;
                    xlim([-250, 250]); ylim([-250, 250]); box off;
                    plot(dot_pos(~isnan(dot_pos(:,1,ii)),1,ii),dot_pos(~isnan(dot_pos(:,2,ii)),2,ii),'o'); hold on;
                    quiver(dot_pos(~isnan(dot_pos(:,1,ii)),1,ii),dot_pos(~isnan(dot_pos(:,2,ii)),2,ii), ...
                        dxdy(~isnan(dot_pos(:,1,ii)),1,ii),dxdy(~isnan(dot_pos(:,2,ii)),2,ii));
                    title(sprintf('frame %02d: coh: %3.2f%% dir:%02d delta: %s',...
                    ii, p.stim.rdk.dots_coherence(cc)*100,p.stim.rdk.dots_direction(dd), deltalabels{dl})); axis square;
                    filename = sprintf('vcd_rdk_coh%03d_dir%02d_delta%02d_motionvec%d.png',p.stim.rdk.dots_coherence(cc)*100,dd,dl,ii);
                    if p.store_imgs
                        print(gcf,'-dpng','-r300',fullfile(saveDir,filename));
                    end
                end
            end
        end
    end
    
    % ek start here, check dot post vs coherence..
    figure(100);
    for ii = 1:60; clf;
        
        

    end
    
    %% single dot: Main dot locations per hemifield
    saveDir = fullfile(vcd_rootPath,'figs',dispname,'dot','visual_checks');
    if ~exist(saveDir,'dir'); mkdir(saveDir); end
    
    bckground = uint8(ones(p.disp.h_pix,p.disp.w_pix))*p.stim.bckgrnd_grayval;
    
    im1 = repmat(bckground,[1 1 3]);
    
    for ang = p.stim.dot.ang_deg
        
        angle = deg2rad(ang);
        [x_shift,y_shift] = pol2cart(angle,p.stim.dot.iso_eccen);
    
        ys = p.disp.yc + round(y_shift*p.disp.ppd);
        xs = p.disp.xc + round(x_shift*p.disp.ppd);
        dot_halfsz = (size(single_dot,1)/2);
        dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz -1);
        dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz -1);
    
        im1( dot_coords_y, dot_coords_x,:) = single_dot;
    end
    
    figure;
    imshow(im1,[1 255]);
    hold on;
    h0 = drawcircle('Center',[p.disp.xc,p.disp.yc],'Radius',11,'color', [1 1 1],'LineWidth',6);
    h0.InteractionsAllowed = 'none';
    title('single dot locations - both hemifields');
    axis image;
    if p.store_imgs
        set(gcf, 'InvertHardCopy', 'off');
        print(fullfile(saveDir,'singledot_all_loc'),'-dpng','-r150');
    end
    %% single dot: visualize individual dot locations

    saveDir = fullfile(vcd_rootPath,'figs',dispname,'dot', 'visual_checks');
    if ~exist(saveDir,'dir'); mkdir(saveDir); end

    bckground = double(ones(p.disp.h_pix,p.disp.w_pix))*0.5;
    
    im1 = bckground;
    dot_halfsz = (size(single_dot,1)/2)-0.5;
    cmap = varysat(parula(16),4);
    xpos_dots_to_plot = cat(1,p.stim.dot.x0_pix, p.stim.dot.x0_pix_delta);
    ypos_dots_to_plot = cat(1,p.stim.dot.y0_pix, p.stim.dot.y0_pix_delta);

    dot_ref_locs = [0, p.stim.dot.delta_from_ref];
    for aa = 1:size(xpos_dots_to_plot,2)
        for bb = 1:size(xpos_dots_to_plot,1)
            figure(99); clf;
            imshow(im1,[0 1]);
            hold all;
            h0 = drawcircle('Center',[p.disp.xc,p.disp.yc],'Radius',11,'color', [1 1 1],'LineWidth',6);
            h0.InteractionsAllowed = 'none';
            h1 = drawcircle('Center',[xpos_dots_to_plot(bb,aa),ypos_dots_to_plot(bb,aa)],'Radius',dot_halfsz,'color',[1 1 1],'EdgeAlpha',0);
            h1.InteractionsAllowed = 'none';
            h1.FaceAlpha=1;
            
            title(sprintf('Angle: %02d; Delta: %02d',p.stim.dot.ang_deg(aa),dot_ref_locs(bb)), 'FontSize',20);
            axis image;
            set(gcf, 'InvertHardCopy', 'off');
            print(fullfile(saveDir,sprintf('%02d_%02d_singledot', aa, bb)),'-dpng','-r150');
        end
    end

    imwrite(single_dot, fullfile(vcd_rootPath,'figs',dispname,'dot',sprintf('singledot.png')));

    %% Complex objects
    obj_masks = masks;
    counter = 1;
    saveDir = fullfile(vcd_rootPath,'figs',dispname, 'objects','visual_checks');
    if ~exist(saveDir,'dir'); mkdir(saveDir); end
    figure(99); set(gcf, 'Position', [300   584   868   753]);
    for objectNr = 1:size(objects,4)
        for rot = 1:size(objects,5)
            clf;
            I = imshow(objects(:,:,:,objectNr,rot),[1 255]);
            I.AlphaData = obj_masks(:,:,objectNr,rot);
            axis image
            
            title(sprintf('Object:%02d Rot:%02d Delta:%02d',objectNr,im_order.abs_rot(counter),im_order.rel_rot(counter)), 'FontSize',20);
            
            if rot == 1
                dd = 0;
            else
                dd = rot;
            end
            print(fullfile(saveDir,sprintf('object%02d_rot%02d_delta%02d', objectNr, rot, dd)),'-dpng','-r150');
            imwrite(objects(:,:,:,objectNr,rot), fullfile(vcd_rootPath,'figs',dispname,'objects',sprintf('object%02d_rot%02d_delta%02d.png', objectNr, rot, dd)),'Alpha',obj_masks(:,:,objectNr,rot));
            counter = counter+1;
        end
    end
    
    %% Natural scenes
    saveDir = fullfile(vcd_rootPath,'figs',dispname,'ns','visual_checks','resized_and_squared');
    if ~exist(saveDir,'dir'); mkdir(saveDir); end
    figure; set(gcf,'Position', [156    91   881   706],'color','w');
    counter = 1;
    for ss = 1:size(scenes,4)
        for bb = 1:size(scenes,5)
            for cc = 1:size(scenes,6)
                clf;
                imagesc(scenes(:,:,:,ss,bb,cc));
                title(sprintf('Im %02d resized & squared',counter), 'FontSize',20);
                axis image; box off; axis off
                set(gca,'CLim',[1 255]);
                if p.stim.store_imgs
                    saveDir = fullfile(vcd_rootPath,'figs','ns','resized_and_squared');
                    if ~exist(saveDir,'dir'), mkdir(saveDir); end
                    print(fullfile(saveDir, sprintf('ns_%02d', counter)),'-dpng','-r150');
                    imwrite(scenes(:,:,:,ss,bb,cc), fullfile(vcd_rootPath,'figs',dispname,'ns',sprintf('ns_%02d.png', counter)));
                end
        
                for ll = 1:size(lures,7)
                    clf;
                    imagesc(lures(:,:,:,ss,bb,cc,ll));
                    title(sprintf('Im %02d, lure %02d resized & squared',counter,ll), 'FontSize',20);
                    axis image; box off; axis off
                    set(gca,'CLim',[1 255]);
                    if p.stim.store_imgs
                        print(fullfile(saveDir, sprintf('ns_%02d_lure%02d', counter,ll)),'-dpng','-r150');
                        imwrite(lures(:,:,:,ss,bb,cc,ll), fullfile(vcd_rootPath,'figs',dispname,'ns',sprintf('ns_%02d_lure%02d.png', counter,ll)));
                    end
                    
                    clf;
                    imagesc(cblind(:,:,:,ss,bb,cc,ll));
                    title(sprintf('Im %02d, cblindness %02d resized & squared',counter,ll), 'FontSize',20);
                    axis image; box off; axis off
                    set(gca,'CLim',[1 255]);
                    
                    if p.stim.store_imgs
                        print(fullfile(saveDir, sprintf('ns_%02d_cblind%02d', counter,ll)),'-dpng','-r150');
                        imwrite(cblind(:,:,:,ss,bb,cc,ll), fullfile(vcd_rootPath,'figs',dispname,'ns',sprintf('ns_%02d_lcblind%02d.png', counter,ll)));
                    end
                end
                counter = counter + 1;
            end
        end
    end
end