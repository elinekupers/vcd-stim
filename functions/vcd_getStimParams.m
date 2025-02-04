function stim = vcd_getStimParams(type, dispname, store_params)

if ~exist('type','var') || isempty(type)
    type = {'all'};
end

if ischar(type)
    type = {type};
end

if any(strcmp(type{:},'all'))
    type = {'gabor','rdk','dot','cobj','ns'};
end

if ~exist('store_params','var') || isempty(store_params)
    store_params = false;
end

% Get display params
disp = vcd_getDisplayParams(dispname);

% Setup struct
stim = struct();

%% GENERAL PARAMS for Nova1x32 (some will get inherited below) 

stim.store_imgs                 = true;
stim.frame_bin                  = 2;                                       % we update disp every 10 VBL monitor refreshes
stim.frame_dur_s                = stim.frame_bin * (1/disp.refresh_hz);    % duration of frame bin (seconds)
stim.iscolor                    = true;                                % use color or not
stim.bckgrnd_grayval            = uint8(127);                          % background color (rgb)

% Default params to inherit (or overwrite)
duration                   = 120;                                       % nr of frame bins
x0_deg                     = [-3.5 0 3.5];                             % [Left Center Right] (degrees) 
y0_deg                     = [0 0 0];                                  % [Left Center Right] (degrees) 
x0_pix                     = round(x0_deg.*disp.ppd);                  % [Left Center Right] (pixels)  
y0_pix                     = round(y0_deg.*disp.ppd);                  % [Left Center Right] (pixels)  
parafov_circle_diam_deg    = 4;                                        % desired parafoveal circular diameter aperture (degrees) 
ctr_square_deg             = 8.4;                                      % desired center square side length (degrees) 
parafov_circle_diam_pix    = round(parafov_circle_diam_deg.*disp.ppd); % desired parafoveal circular diameter aperture (pixels) 
ctr_square_pix             = round(ctr_square_deg.*disp.ppd);          % desired center square side length in pixels (pixels) 

%% FIXATION DOT

% TEMPORAL
stim.fix.dotmeanchange          = 3;                       % dot changes occur with this average interval (in seconds)
stim.fix.dotchangeplusminus     = 2;                       % plus or minus this amount (in seconds)

% SPATIAL
stim.fix.dotsize_deg            = 0.2;                     % dot diameter in deg
stim.fix.dotsizeborder_pix      = 3;                       % pixel-width for dot border
stim.fix.lumstep                = 10;                      
stim.fix.dotcol_rgba            = {uint8(repmat([0:15:255]',[1 3])) 0.5};  % dot color, 16 gray levels (incl. blank, white)  50% opacity
stim.fix.dotdiam_pix            = [2*round(stim.fix.dotsize_deg * disp.ppd) stim.fix.dotsizeborder_pix];  % fixation diameter in pixels (EVEN)

fprintf('*** FIXATION MARK: fixation [diam rim] = [%d %d] pixels ***\n',stim.fix.dotdiam_pix);

for ii = 1:length(type)
   
    p = [];
        
    switch type{ii}
        
        case {'gabor',1}
            
            % Where to store stimulus images?
            p.stimfile = fullfile(vcd_rootPath,'workspaces','gabors.mat');
            p.infofile = fullfile(vcd_rootPath,'workspaces','gabors_info.csv');
            % TEMPORAL
            % Fixed params
            p.duration        	  = duration;                               % frames (nr of monitor refreshes)
            
            % Manipulated params
            p.t_gausswin_N        = duration;                               % number of timepoints for gaussian time window (contrast ramp)
            p.t_gausswin_std      = 8;                                      % standard devation of gaussian window in time (ms)

            % SPATIAL
            % Fixed params
            p.img_sz_deg      = parafov_circle_diam_deg;                    % height (or width) of stimulus support (deg)
            p.img_sz_pix      = parafov_circle_diam_pix;                    % height (or width) of square stimulus support (pix)
            p.gauss_std_deg   = 0.5;                                        % standard deviation of gaussian window (deg)
            p.gauss_std_pix   = round(p.gauss_std_deg.*disp.ppd);
            
            p.sf_cpd          = 4;                                          % spatial frequency (cycles/deg)
            p.cycles_per_pix  = (p.sf_cpd/p.img_sz_deg)/disp.ppd;           % nr of cycles per image (pix)
            
            p.x0_deg          = x0_deg;                                     % x-center loc in deg (translation from 0,0)
            p.y0_deg          = y0_deg;                                     % y-center loc in deg (translation from 0,0)
            p.x0_pix          = x0_pix;                                     % x-center loc in pix (translation from 0,0)
            p.y0_pix          = y0_pix;                                     % y-center loc in pix (translation from 0,0)
            
            p.ph_deg          = [0:(180/4):179];                            % 4 quadrature Gabor phases
            
            % Manipulated params
            p.contrast        = [0.05, 0.10, 1];                            % Michelson contrasts [0-1] (fraction)
            p.ori_jitter_sd   = 2;                                          % std of normal distribution to sample orientation jitter 
            p.ori_jitter_mu   = 1;                                          % mean of normal distribution to sample orientation jitter 
            p.num_ori         = 8;                                          % orientation "bins" from which we create final gabor orientations (deg), 0 = 12 o'clock
            p.ori_bins        = [10:(180/p.num_ori):170];                   % rotate 10 deg away from vertical, to avoid ill-defined response options
            p.ori_jitter      = p.ori_jitter_mu + (p.ori_jitter_sd.*randn(1,p.num_ori)); % add a small amount of jitter around the orientation
            p.ori_deg         = round(p.ori_bins + p.ori_jitter);           % final gabor orientations
            
            p.delta_from_ref  = [-15, -5, 5, 15];                           % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                                                                            % the bigger the delta, the easier the trial
            % Add params to struct
            stim.gabor = p;

            
        case {'rdk',2}
            % Where to store stimulus images?
            p.stimfile = fullfile(vcd_rootPath,'workspaces','rdk.mat');
            p.infofile = fullfile(vcd_rootPath,'workspaces','rdk_info.csv');

            % TEMPORAL
            p.duration        = duration;                                   % frames (nr of monitor refreshes)
            p.dots_coherence  = [0.064, 0.128, 0.512];                      % Kiani lab uses usually one of these [0 0.032 0.064 0.128 0.256 0.512]
            p.dots_speed      = 5;                                          % pix/ms? Kiani lab uses usually 5 to 10
            p.dots_interval   = 4;                                          % frames  to approx match Kiani 3 x 75 hz frames
            
            % SPATIAL
            p.img_sz_deg      = parafov_circle_diam_deg;                    % stimulus aperture diameter (deg)
            p.img_sz_pix      = parafov_circle_diam_pix;                    % stimulus aperture diameter (pix)
            
            p.num_mot_dir     = 8; 
            p.motdir_bins      = [0:(360/p.num_mot_dir):359]+10;             % sample direction of coherent motion from [0-359] in deg (0 deg is aligned with 12 o'clock)
            p.motdir_jitter_sd = 2;                                          % std of normal distribution to sample orientation jitter 
            p.motdir_jitter_mu = 1;                                          % mean of normal distribution to sample orientation jitter 
            p.motdir_jitter    = p.motdir_jitter_mu + (p.motdir_jitter_sd.*randn(1,p.num_mot_dir)); % add a small amount of jitter around the orientation
            p.dots_direction   = round(p.motdir_bins + p.motdir_jitter);     % final gabor orientations

            % RDK specific                                                        
            p.dots_density     = 16.7;                                      % dots/deg??
            p.dots_size        = 3;                                         % radius in pixels??
            p.dots_color       = [255 255 255;0 0 0]./255;                  % 50:50 white:black, color in RGB [0-1]
            p.max_dots_per_frame = 200;                                     % from Kiani lab (roughly matches to nr of pixels in aperture)
            p.dots_contrast    = 1;                                         % Michelson [0-1] (fraction)
            
            p.x0_deg          = x0_deg;                                     % x-center loc in deg (translation from 0,0)
            p.y0_deg          = x0_deg;                                     % y-center loc in deg (translation from 0,0)
            p.x0_pix          = x0_pix;                                     % x-center loc in pix (translation from 0,0)
            p.y0_pix          = y0_pix;                                     % y-center loc in pix (translation from 0,0)
            
            p.delta_from_ref  = [-15, -5, 5, 15];                           % how much should stim iso-eccen loc deviate from reference (WM: double epochs)

            % Add params to struct
            stim.rdk = p;
            
        case {'dot',3}
            % Where to store stimulus images?
            p.stimfile = fullfile(vcd_rootPath,'workspaces','dots.mat');
            
            % TEMPORAL
            p.duration        = duration;                                   % frames (nr of monitor refreshes)
            
            % SPATIAL
            p.img_sz_deg      = 1;                                          % radius in deg
            p.img_sz_pix      = p.img_sz_deg * disp.ppd;                    % radius in pix
            p.color           = [255 255 255]./255;                         % white, color in RGB [0-1]
            p.contrast        = 1;                                          % Michelson [0-1] (fraction)
            p.x0_deg          = x0_deg;                                     % x-center loc in deg (translation from 0,0)
            p.y0_deg          = x0_deg;                                     % y-center loc in deg (translation from 0,0)
            p.x0_pix          = x0_pix;                                     % x-center loc in pix (translation from 0,0)
            p.y0_pix          = y0_pix;                                     % y-center loc in pix (translation from 0,0)
            
            p.iso_eccen       = 5.5;                                        % eccentricity of dots from center (degrees)
            p.num_loc         = 16;                                          % orientation "bins" from which we create final gabor orientations (deg), 0 = 12 o'clock
            p.loc_jitter_sd   = 1;                                          % std of normal distribution to sample orientation jitter
            p.loc_jitter_mu   = 1;                                          % mean of normal distribution to sample orientation jitter 
            p.loc_bins        = [3:(180/p.num_loc):176];                   % rotate 10 deg away from vertical, to avoid ill-defined response options
            p.loc_jitter      = p.loc_jitter_sd + (p.loc_jitter_mu.*randn(1,p.num_loc)); % add a small amount of jitter around the orientation
            p.loc_deg         = round(p.loc_bins + p.loc_jitter);           % final gabor orientations
            
            p.delta_from_ref  = [-15, -5, 5, 15];                           % how much should stim iso-eccen loc deviate from reference (WM: double epochs)
                                                                            % the bigger the delta, the easier the trial
            % Add params to struct
            stim.dot = p;
            
        case {'cobj',4}
            % Where to store stimulus images?
            p.stimfile = fullfile(vcd_rootPath,'workspaces','objects.mat');
            
            % Where to find stimulus info?
            p.infofile = fullfile(vcd_rootPath,'workspaces','objects_info.csv');
            
            % TEMPORAL
            p.duration      = duration;                                     % frames (nr of monitor refreshes)
            
            % SPATIAL
            p.img_sz_deg    = parafov_circle_diam_deg;                      % height (or width) of square stimulus support (deg)
            p.img_sz_pix    = parafov_circle_diam_pix;                      % height (or width) of square stimulus support (pix)
            p.contrast      = 1;                                            % Michelson [0-1] (fraction)
            p.x0_deg        = x0_deg;                                       % x-center loc in deg (translation from 0,0)
            p.y0_deg        = y0_deg;                                       % y-center loc in deg (translation from 0,0)
            p.x0_pix        = x0_pix;                                       % x-center loc in pix (translation from 0,0)
            p.y0_pix        = y0_pix;                                       % y-center loc in pix (translation from 0,0)
            
            p.num_unique_objects = 16;                                      % orientation "bins" from which we create final gabor orientations (deg), 0 = 12 o'clock
            p.facing_dir_deg        = repmat(90 + [-34, 34],1,p.num_unique_objects/2); % rotate 10 deg away from canonical view

            p.super_cat    = {'human','animal','object','place'};           % 4 superordinate categories
            
            p.basic_cat{1} = {'facemale','facefemale','facefemale'};                             
            p.basic_cat{2} = {'small','small','big','big'};                               
            p.basic_cat{3} = {'tool','tool','food','food','vehicle','vehicle'};                      
            p.basic_cat{4} = {'building','building','building'};                               
            
            p.sub_cat{1} = {'damon','lisa','sophia'};                    
            p.sub_cat{2} = {'parrot','cat','bear','giraffe'};
            p.sub_cat{3} = {'drill','brush','pizza','banana','bus','car'};      
            p.sub_cat{4} = {'church','house','watertower'};

            p.delta_from_ref = [-15, -5, 5, 15];                            % how much should stim pose rotate from reference (WM: for double epochs)
                                                                            % the bigger the delta, the easier the trial. Negative is counter-clockwise, positive is clockwise
            % Add params to struct
            stim.cobj = p;
            
        case {'ns',5}
            % Where to store stimulus images?
            p.stimfile = fullfile(vcd_rootPath,'workspaces','scenes.mat');
            
            % Where to find stimulus info?
            p.infofile = fullfile(vcd_rootPath,'workspaces','scenes_info.csv');
            
            % TEMPORAL
            p.duration    = duration;                                       % frames (nr of monitor refreshes)
            
            % SPATIAL
            p.img_sz_deg  = ctr_square_deg;                                 % height (or width) of square stimulus support (deg)
            p.img_sz_pix  = ctr_square_pix;                                 % height (or width) of square stimulus support (pix)
                                                            
            p.og_res_stim = 425;                                            % original resolution of NSD stimuli
            p.rz_res_stim = 714;                                            % resized resolution of NSD stimuli
            p.dres = ((ctr_square_deg/disp.h_deg * disp.h_pix) / p.og_res_stim);  % scale factor to apply
%             assert(isequal(stim.ns.rz_res_stim, round(stim.ns.img_sz_pix * -stim.dres))); % check: slightly larger than OG NSD stim size?? = 744 pixels 7TAS BOLDscreen Nova1x32
            
            p.x0_deg        = x0_deg;                                       % x-center loc in deg (translation from 0,0)
            p.y0_deg        = y0_deg;                                       % y-center loc in deg (translation from 0,0)
            p.x0_pix        = x0_pix;                                       % x-center loc in pix (translation from 0,0)
            p.y0_pix        = y0_pix;                                       % y-center loc in pix (translation from 0,0)
            
            p.super_cat    = {'human','animal','food','objects','places'};  % 5 superordinate categories
            
            p.basic_cat{1} = {'indoor','outdoor'};
            p.basic_cat{2} = {'indoor','outdoor'};
            p.basic_cat{3} = {'indoor','outdoor'};
            p.basic_cat{4} = {'indoor','outdoor'};
            p.basic_cat{5} = {'indoor','outdoor'};

            p.sub_cat{1,1} = {'face1_left','face2_center','face3_right'};                  
            p.sub_cat{1,2} = {'face1_left','face2_center','face3_right'};                  
            
            p.sub_cat{2,1} = {'cat1_left','cat2_center','cat3_right'};                  
            p.sub_cat{2,2} = {'giraffe1_left','giraffe2_center','giraffe3_right'};                                                
           
            p.sub_cat{3,1} = {'donut1_left','donut2_center','donut3_right'};
            p.sub_cat{3,2} = {'banana1_left','banana2_center','banana3_right'};              
            
            p.sub_cat{4,1} = {'vase1_left','vase2_center','vase3_right'};
            p.sub_cat{4,2} = {'bus1_left','bus2_center','bus3_right'};                  

            p.sub_cat{5,1} = {'bathroom1_left','bathroom2_center','bathroom3_right'};     
            p.sub_cat{5,2} = {'building1_left','building2_center','building3_right'};     

            p.change_im = {'easy','hard'};
%             p.delta_from_ref  = [-5 15 5 15];                               % how much should scene pov rotate L/R for WM
                                                                            % the bigger the delta, the easier the trial 
            
            % Add params to struct
            stim.ns = p;
            
            fprintf('*** %s: dres (scale factor) = %.4f ***\n',upper(type{ii}),stim.ns.dres);
    end
    
    fprintf('*** %s: Using a stimulus size of %2.2f deg (%3.2f pixels). ***\n', upper(type{ii}), p.img_sz_deg, p.img_sz_pix);
    fprintf('*** %s: Using a stimulus duration of %d frame bins (%3.2f seconds), where 1 bin is %d frames of 1/%d Hz ***\n', ...
        upper(type{ii}), p.duration, p.duration*stim.frame_dur_s, stim.frame_bin, disp.refresh_hz);
end

if store_params
    save(fullfile(vcd_rootPath,'workspaces','p_stim.mat'),'stim') 
end

return