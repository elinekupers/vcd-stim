%% s_createStim.m

verbose = true; % plot stimuli or not

% Get display params
dispname = '7TAS_BOLDSCREEN32'; % or 'KKOFFICEQ3277' or psychophys room??
p.disp = vcd_getDisplayParams(dispname);

%% Get stimulus parameters
p.load_params  = true;
p.store_params = true;
p.store_imgs   = true;
saveFigsFolder = fullfile(vcd_rootPath,'figs');

p.stim   = vcd_getStimParams('all',p.disp,p.load_params,p.store_params); % Choose from <'gabor'> <'rdk'> <'dot'> <'cobj'> <'ns'> <'all'>
p.exp    = vcd_getSessionParams(p,p.load_params,p.store_params);
%%
p.trials = vcd_makeTrials(p,p.load_params,p.store_params);
subject_sessions = vcd_createSessions(p,p.load_params,p.store_params);


%% background

bckgrnd_im = vcd_pinknoisebackground(p, 'comb', 'fat', 1); % 'skinny' or 'fat', type: 'puzzle', 'dotring', or 'comb'

%% fixation dot
fix_im = vcd_fixationDot(p);

%% Gabor
gbr_im = vcd_gabor(p);

%% RDK
rdk_im = vcd_rdk(p);

%% Simple dot
img_dot  = vcd_simpledot(p);

%% Complex object
img_cobj = vcd_complexobjects(p);


%% Natural scenes
img_ns  = vcd_naturalscenes(p);


%% Visualize stimuli

if verbose
    
    %% Gabor
    figure(1); clf;
    set(gcf, 'Position', [0 0 1024 1080], 'color','w')
    if verbose
        for dd = 1:size(gbr_im,6)
            if dd==1, dlta = 0; else, dlta = p.stim.gabor.delta_from_ref(dd-1); end
            for cont = 1:size(gbr_im,5)
                for pp = 1:size(gbr_im,4)
                    for ori = 1:size(gbr_im,3)
                        clf;
                        imshow(gbr_im(:,:,ori,pp,cont,dd),[0 256]);
                        
                        title(sprintf('%d:deg c:%1.2f ph:%d dlta:%d', ...
                            p.stim.gabor.ori_deg(ori),...
                            p.stim.gabor.contrast(cont),...
                            p.stim.gabor.ph_deg(pp),...
                            dlta))
                        drawnow;
                        if p.store_imgs
                            filename = sprintf('vcd_gabor_o%dc%dph%dd%d.png',ori,cont,pp,dd);
                            print(fH,'-dpng','-r300',fullfile(saveFigsFolder,filename));
                        end
                    end
                end
            end
        end
    end

    %% RDK
    
    for cc = 1:length(p.stim.rdk.dots_coherence)
        for dd = 1:length(p.stim.rdk.dots_direction)
            vcd_createStimVideo(rdk_im{dd,3}, p.stim.rdk.dots_interval/display.refresh_hz, ...
                fullfile('~/Desktop'), sprintf('vcd_rdk_coh%1.3f_dir%02d',p.stim.rdk.dots_coherence(cc),p.stim.rdk.dots_direction(dd)));
        end
    end
    
    
    
    %% Simple dot
    bckground = uint8(ones(disp.h_pix,disp.w_pix))*p.stim.bckgrnd_grayval;
    
    im1 = bckground;
    
    for ang = p.stim.dot.loc_deg
        angle = deg2rad(ang-90);
        [x_shift,y_shift] = pol2cart(angle,p.stim.dot.iso_eccen);
    
        ys = display.yc + round(y_shift*display.ppd);
        xs = display.xc + round(x_shift*display.ppd);
        dot_halfsz = (size(img_dot,1)/2)-0.5;
        dot_coords_x = (xs - dot_halfsz) : (xs + dot_halfsz);
        dot_coords_y = (ys - dot_halfsz) : (ys + dot_halfsz);
    
        im1( dot_coords_y, dot_coords_x) = img_dot;
    end
    
    figure;
    imshow(im1,[1 256]);
    hold on;
    plot(display.xc+[-5:5], display.yc.*ones(1,11),'k',display.xc*ones(1,11),display.yc+[-5:5],'k')
    title('simple dot - one hemifield');
    axis image;

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