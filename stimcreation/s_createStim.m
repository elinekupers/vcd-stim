%% s_createStim.m
%
% Stand-alone script to create and store the stimuli shown in VCD core 
% experiment.

%% %%%%%%%%%%%%%%%%%%%
%%%%%% PARAMETERS %%%% 
%%%%%%%%%%%%%%%%%%%%%%

verbose        = true; % visualize stimuli or not
p.store_imgs   = true; % store visualization figures
saveFigsFolder = fullfile(vcd_rootPath,'figs'); % where to store visualization figures

% Get display params
dispname = '7TAS_BOLDSCREEN32'; % Choose from: '7TAS_BOLDSCREEN32';'KKOFFICE_AOCQ3277';'EKHOME_ASUSVE247';'PPROOM_EIZOFLEXSCAN'
p.disp   = vcd_getDisplayParams(dispname);

% Get stimulus parameters
p.load_params                 = false; % if false, re-create params.
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
gaptype     = 'comb';
borderwidth = 'fat';
num         = 1;
bckgrnd_im  = vcd_pinknoisebackground(p, ...
                                     'gaptype',gaptype, ...
                                     'borderwidth', borderwidth,...
                                     'num',num); 

%% Fixation dot images
% fix_im (uint8) dimensions are w (pix) x h (pix) x 3 x 5 luminance levels
% x 5 dot rims types (thin, thick, thick-left, thick-right, thick-both)
[fix_im, mask, fix_info, p] = vcd_fixationDot(p);

%% Gabor images
% gabor (uint8) dimensions are w (pix) x h (pix) x 3 x 24 unique images (8
% ori x 3 contrasts) x 5 delta levels (0 + 4 deltas)
[gabor, gbr_masks, gbr_info, p] = vcd_gabor(p);

%% RDK images
% Output images:
% * rdk (uint8) is a 8x3x5 cell for each 8 motion directions x 3
%   coherence levels x 5 delta offsets (0, -15, -5, +5, +15 deg). Each cell
%   contains w (pixels) x h (pixels) x 3 x frames.
[rdk, rdk_masks, rdk_info, p] = vcd_rdk(p);

%% Simple dot
% Output images:
% * simple_dot (uint8) is a single matrix: w (pixels) by h (pixels).
[simple_dot, dot_masks, dot_info, p] = vcd_simpledot(p);

%% Objects
% Output images:
% * objects (uint8) dimensions are: width (pixels) x height (pixels) x 3
%   (rgb) x 16 object categories (subordinate level) x 22 rotations
% * alpha transparency mask: cobj_masks (uint8) dimensions are: width
%   (pixels) x height (pixels) x 16 object categories (subordinate level) x
%   22 rotations
[objects, cobj_masks, cobj_order, cobj_info, p] = vcd_objects(p);

%% Natural scenes
% Output images:
% * scenes (uint8) dimensions are width x height x 3 (rgb) x 5 superordinate
% categories x 2 scene locations (indoor/outdoor) x 3 obj locations (left/middle/right)
% * lure_im (uint8) dimensions are width x height x 3 (rgb) x 5 superordinate
% categories x 2 scene locations (indoor/outdoor) x 3 obj locations
% (left/middle/right) x 4 lure types:
%   1:very similar/difficult ... 4:least similar/easy
% * cblind_im (uint8) dimensions are width x height x 3 (rgb) x 5 superordinate
% categories x 2 scene locations (indoor/outdoor) x 3 obj locations
% (left/middle/right) x 4 change types:
%   1:easy add -- scene is altered by adding something big/obvious 
%   2:hard add -- scene is altered by adding something small/subtle
%   3:easy remove -- scene is altered by removing something big/obvious
%   4:hard remove -- scene is altered by removing something small/subtle
[scenes, lure_im, cblind_im, ns_order, ns_info, p] = vcd_naturalscenes(p);


%% %%%%%%%%%%%%%%%%
%%%% Visualize %%%% 
%%%%%%%%%%%%%%%%%%%

if verbose
    makeprettyfigures;
    
    %% Background
    fH = figure(100);
    set(fH, 'Position', [0 0 1920 1080], 'color','w')

    for jj = 1:num
        clf;
        imagesc(bckgrnd_im(:,:,jj)); colormap gray
        axis image
        title(sprintf('Background im %s %s %03d',gaptype, borderwidth, jj))
        set(gca, 'TickDir','out', 'LineWidth',2,'FontSize',20)
        drawnow;
        if p.store_imgs
            filename = sprintf('vcd_background_%s%s%03d.png',gaptype, borderwidth, jj);
            print(fH,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
    end

    %% Fixation dot
    fH = figure(101);
    set(fH, 'Position', [1   400   750   578], 'color','w')
    counter = 1;
    for ii = 1:size(fix_im,5)
        for jj = 1:size(fix_im,5)
            subplot(5,5,counter)
            imshow(squeeze(fix_im(:,:,1,ii,jj)),[1 255])
            colormap gray; axis square
            title(sprintf('%01d, %01d',ii, jj))
            counter = counter+1;
        end
    end
    if p.store_imgs
        filename = sprintf('vcd_fixdots.png');
        print(fH,'-dpng','-r300',fullfile(saveFigsFolder,filename));
    end
    
    %% Gabor
    fH1 = figure(1); fH2 = figure(2); 
    set(fH1, 'Position', [0 0 1024 1080], 'color','w')
    set(fH2, 'Position', [1   400   750   578], 'color','w')
    if verbose
        for dd = 1:size(gbr_im,5)
            if dd==1, dlta = 0; else, dlta = p.stim.gabor.delta_from_ref(dd-1); end
            for ii = 1:size(gbr_im,4)
%                 figure(fH1); clf;
%                 imagesc(gbr_im(:,:,:,ii,dd)); colormap gray; set(gca, 'CLim',[1 255])
%                 title(sprintf('ori:%03d deg c:%1.2f ph:%03d delta:%02d', ...
%                     gbr_info(gbr_info.unique_im==ii,:).ori_deg,...
%                     gbr_info(gbr_info.unique_im==ii,:).contrast,...
%                     gbr_info(gbr_info.unique_im==ii,:).phase_deg,...
%                     gbr_info(gbr_info.unique_im==ii,:).delta_deg));
%                 cb = colorbar; cb.Ticks = [1, 50, 100, 150, 200, 250];
%                 set(gca, 'FontSize',20, 'TickDir','out', 'LineWidth', 2); axis image;
%                 drawnow;
%                 if p.store_imgs
%                     filename = sprintf('vcd_gabor_im%02d_delta%02d.png',ii,dd);
%                     print(fH1,'-dpng','-r300',fullfile(saveFigsFolder,filename));
%                 end
                
                % Plot pix histogram
                figure(fH2); clf
                histogram(gbr_im(:,:,:,ii,dd), 'NumBins', 30); 
                title(sprintf('ori:%03d deg c:%1.2f ph:%03d delta:%02d', ...
                    gbr_info(gbr_info.unique_im==ii,:).ori_deg,...
                    gbr_info(gbr_info.unique_im==ii,:).contrast,...
                    gbr_info(gbr_info.unique_im==ii,:).phase_deg,...
                    gbr_info(gbr_info.unique_im==ii,:).delta_deg));
                box off; axis square; xlim([-10 260]); ylim([0,3.5].*10^5)
                drawnow;
                set(gca,'XTick',[0 64 128 190 255])
                set(gca,'YTick',[1:3].*10^5)
                if p.store_imgs
                    filename = sprintf('vcd_gabor_im%02d_delta%02d_hist.png',ii,dd);
                    print(fH2,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                end
            end
        end
    end

    %% RDK
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
    for dl = length(deltalabels)
        for cc = 1:length(p.stim.rdk.dots_coherence)
            for dd = 1:length(p.stim.rdk.dots_direction)
                vcd_createStimVideo(rdk{dd,cc,dl}, 1/p.stim.presentationrate_hz, ...
                    fullfile(saveFigsFolder), ...
                    sprintf('vcd_rdk_coh%s_dir%02d_delta%s',...
                    cohlabels{cc},p.stim.rdk.dots_direction(dd),deltalabels{dl}));
            end
        end
    end
    
    
    %% Simple dot: Main dot locations per hemifield
    bckground = uint8(ones(p.disp.h_pix,p.disp.w_pix))*p.stim.bckgrnd_grayval;
    
    im1 = bckground;
    
    for ang = p.stim.dot.loc_deg
        
        angle = deg2rad(ang);
        [x_shift,y_shift] = pol2cart(angle,p.stim.dot.iso_eccen);
    
        ys = p.disp.yc + round(y_shift*p.disp.ppd);
        xs = p.disp.xc + round(x_shift*p.disp.ppd);
        dot_halfsz = (size(simple_dot,1)/2);
        dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz -1);
        dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz -1);
    
        im1( dot_coords_y, dot_coords_x) = simple_dot;
    end
    
    figure;
    imshow(im1,[1 255]);
    hold on;
    h0 = drawcircle('Center',[p.disp.xc,p.disp.yc],'Radius',16,'color', [1 1 1]);
    h0.InteractionsAllowed = 'none';
    title('simple dot locations - both hemifields');
    axis image;
    if p.store_imgs
        set(gcf, 'InvertHardCopy', 'off');
        print(fullfile(vcd_rootPath,'figs','simple_dot','simpledot_all_loc'),'-dpng','-r150');
    end
    %% Simple dot: visualize individual dot locations
    if ~isfield(p,'disp')
        p.disp = vcd_getDisplayParams;
    end
    bckground = double(ones(p.disp.h_pix,p.disp.w_pix))*0.5;
    
    im1 = bckground;
    dot_halfsz = (size(simple_dot,1)/2)-0.5;
    cmap = varysat(parula(16),4);
    xpos_dots_to_plot = cat(1,p.stim.dot.x0_pix, p.stim.dot.x0_pix_delta);
    ypos_dots_to_plot = cat(1,p.stim.dot.y0_pix, p.stim.dot.y0_pix_delta);

    dot_ref_locs = [0, p.stim.dot.delta_from_ref];
    for aa = 1:size(xpos_dots_to_plot,2)
        for bb = 1:size(xpos_dots_to_plot,1)
            figure(99); clf;
            imshow(im1,[0 1]);
            hold all;
            h0 = drawcircle('Center',[p.disp.xc,p.disp.yc],'Radius',16,'color', [1 1 1]);
            h0.InteractionsAllowed = 'none';
            h1 = drawcircle('Center',[xpos_dots_to_plot(bb,aa),ypos_dots_to_plot(bb,aa)],'Radius',dot_halfsz,'color',[1 1 1],'EdgeAlpha',0);
            h1.InteractionsAllowed = 'none';
            h1.FaceAlpha=1;
            
            title(sprintf('Angle: %02d; Delta: %02d',p.stim.dot.ang_deg(aa),dot_ref_locs(bb)), 'FontSize',20);
            axis image;
            set(gcf, 'InvertHardCopy', 'off');
            print(fullfile(vcd_rootPath,'figs','simple_dot',sprintf('%02d_%02d_simpledot', aa, bb)),'-dpng','-r150');
        end
    end

    %% Complex objects
    
    for objectNr = 1:size(objects,3)
        figure(objectNr); clf;
        for ii = 1:size(objects,4)
            subplot(2,11,ii); hold all;
            imshow(objects(:,:,objectNr,ii),[1 255]);
        end
    end
    
    %% Natural scenes
    
end