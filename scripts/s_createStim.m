%% s_createStim.m

verbose = true; % plot stimuli or not 

dispname = '7TASBOLDSCREEN32'; % or 'KKOFFICEQ3277'
display = vcd_getDisplayParams(dispname);

%% Get stimulus parameters

p.stim   = vcd_getStimParams('all',dispname); % Choose from 'gabor', 'rdk','dot','cobj', or 'ns'
p.task   = vcd_getTaskParams;
trials   = vcd_makeTrials(p);

p.store_imgs = true;

%% Gabor
[gbr_im,gbr_info,p] = vcd_gabor(p);

%%
figure(1); clf;
set(gcf, 'Position', [0 0 1024 1080], 'color','w')
if verbose
    for delta = 1:size(gbr_im,6) 
        for contrast = 1:size(gbr_im,5) 
            for pp = 1:size(gbr_im,4) 
                for ori = 1:size(gbr_im,3) 
                    fH =figure(1); clf;
                    imshow(gbr_im(:,:,ori,pp,contrast,delta),[0 256]);
                    title(sprintf('%d deg Contrast%d Phase%d Delta%d',p.stim.gabor. ori,contrast,pp,delta))
                    drawnow;
                    filename = sprintf('vcd_gabor_o%dc%dph%dd%d.png',ori,contrast,pp,delta); 
                    print(fH,'-dpng','-r300',fullfile('~/Desktop',filename));  
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

