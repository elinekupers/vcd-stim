%% s_createStim.m

verbose = true; % plot stimuli or not 

dispname = '7TASBOLDSCREEN32'; % or 'KKOFFICEQ3277'
display = vcd_getDisplayParams(dispname);

%% Get stimulus parameters
p.store_params = true;
p.store_imgs   = true;
saveFigsFolder = fullfile(vcd_rootPath,'figs');

p.stim   = vcd_getStimParams('all',dispname,p.store_params); % Choose from <'gabor'> <'rdk'> <'dot'> <'cobj'> <'ns'> <'all'>
p.task   = vcd_getTaskParams;
trials   = vcd_makeTrials(p);


%% Gabor
[gbr_im,gbr_info,p] = vcd_gabor(p);

%%
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
[rdk_im,rdk_info,p] = vcd_rdk(p, display);

if verbose
    
    for cc = 1:length(p.stim.rdk.dots_coherence)
        for dd = 1:length(p.stim.rdk.dots_direction)
            vcd_createStimVideo(rdk_im{dd,3}, p.stim.rdk.dots_interval/display.refresh_hz, ...
                fullfile('~/Desktop'), sprintf('vcd_rdk_coh%1.3f_dir%02d',p.stim.rdk.dots_coherence(cc),p.stim.rdk.dots_direction(dd)));
        end
    end
end
%% Simple dot
params  = vcd_getStimParams('dot');
[img_dot,coord,params.dot]  = vcd_simpledot(params);

%% Complex object
params  = vcd_getStimParams('cobj');
img_obj  = vcd_complexobjects(params);

%% Natural scenes
params  = vcd_getStimParams('ns');
img_ns  = vcd_naturalscenes(params);

