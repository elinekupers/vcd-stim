function img_quiz_dots = vcd_getImageryQuizDots(params, varargin)
% Function to create uint8 quiz dot images, quiz dot locations, and load
% text prompts to probe subject's mental image of special core Gabors,
% RDKs, dot, objects, or naturalistic scenes.
%
%  img_quiz_dots = vcd_getImageryQuizDots(params, [update_info_file], [verbose])
%
% Purpose:
%   Create 20 specific quiz dot images that are called after the
%   delay period, in the second stimulus presentation window of the imagery
%   trial. The goal of the subject is to say how the quiz dots relate to
%   the mental image, according to the text prompt. For example, do the
%   quiz dots overlap the giraffe in the imagined giraffe scene? Yes/No.
%   Every imagined special core stimulus has 10 "yes" test images and 10 
%   "no" test images. Dot images are only made for BOLDscreen.
%
% Some notes on imagery quiz dot images: 
% - The imagery quiz dot are 0.5° diameter and white (uint 255).
% - Imagery is done only for the special core.
% - The general trial structure for imagery is: spatial cue, small text
%   phrase at center of display (stim1) [allow eye movements], 8-s delay,
%   two-dot test image in both left and right hemifields for classic 
%   stimuli (stim2).
% - Exact text prompts can be found here:
%       ./workspaces/instructions/img_image_prompts/ (pngs)
%       ./workspaces/instructions/img_text_prompts/  (txt files)
% - Creation of imagery quiz dots for Gabors/RDKs/Dots are automated.
% - Creation of imagery quiz dots for Objects and Scenes are manually drawn 
%   dots that are then postprocessed. Postprocessing will deal with  
%   stimulus size, detect x-y-positions of the two dots, etc.
% - Gabor: two dots (within the 4° aperture) perfectly aligned or off
% - RDK: two dots (within the 4° aperture) perfectly aligned or off
% - Single dot: two dots (on the 4° iso-ecc circle) straddling or not. 
% - Object: two dots either both fully within object mask or not [within 4° circle]
% - NS: two dots either both fully within some fixed object or not
% - For object and NS, the yes’s are similar in difficulty, whereas the 
%   no’s vary a bit in difficulty.
% - Parameters for a specific display environment can be derived by calling
%   vcd_getStimParams.m.
%
% Quiz image support:
% * Gabor and RDK support area: circular aperture with 4-degree diameter 
%      (352 pixels for BOLDscreen). Support location: same as the imagined
%      special core stimulus [-4,0] or [4,0] degrees relative from the
%      center of the display ([0,0]). (608 or 1312 pixels for BOLDscreen).
% * Single Dot support area: 4-deg iso-eccentric ring (BOLDscreen: 352
%      pixels). Support location: ring is centered on the center of the display.
% * Object support area: circular aperture with 5-degree diameter (XX
%     pixels for BOLDscreen). Support location: centered on corresponding 
%     core object that subjects imagine, either [-4,0] or [4,0] degrees 
%     relative from the center of the display ([0,0]). (608 or 1312 pixels 
%     for BOLDscreen).
% * NS support area: square aperture width/height: 8.4 degrees.
%   (741 pixels for BOLDscreen). Support location: centered on display:
%   [0,0] degrees or [960,540] pixels.
%
%
% INPUTS:
%   params             : (struct) parameter struct should contain a field for 
%                         display params (params.disp), and stimulus parameters 
%                         for each stimulus class:
%                           - params.stim.bckgrnd_grayval - uint8 value used for mean luminance (between [0 255])
%                           - params.stim.(stim_class).img_sz_pix - special core stimulus size in pixels (diameter or width)
%                           - params.stim.(stim_class).imagery_sz_deg - visual field location size in degrees where imagery quiz dots can be placed
%                           - params.stim.(stim_class).imagery_sz_pix - visual field location size in pixels where imagery quiz dots can be placed
%                           - params.stim.(stim_class).quiz_dot_diam_pix - diameter of imagery quiz dot in pixels
%                           - params.stim.(stim_class).imagery_quiz_images - list of match (1) and no match (2) imagery quiz dot images
%                           -
%                           params.stim.(gabor/rdk/dot).imagery_min_dot_dist_deg - minimum distance between two quiz dots on one side of the fixation circle (in degrees visual angle)
%                           - params.stim.(gabor/rdk).imagery_minmax_anglemismatch - minimum and maximum angular distance between the aligned tilt/motion direction (in angular degrees)
%                           - params.stim.(stim_class).stimfile
%                           - params.stim.(stim_class).unique_im_nrs_img_test
%                           - params.stim.(stim_class).unique_im_nrs_specialcore
%   [update_info_file] : [optional] (logical) if true, update info table in stimulus mat file and corresponding csv info file (default is false)
%   [verbose]          : [optional] (logical) if true, plot debug figures (default is false).
%
% OUTPUTS:
%   all_quiz_dots.(stimClass) : struct with the following fields for each stimulus class: 
%       * im                  : quiz dot images: height (pixels) by width (pixels) x 3 (rgb)
%       * mask                : alpha masks for each quiz dot image to crop out image edges
%       * xy_coords_deg       : (double) is a 4D array with [angle,eccentricity] coordinates in deg 
%                               (relative to center of imagined core stimulus [0,0]):
%                               Dimensions of array: 
%                               nr special core images (8 for Gabor/RDK/DOT/OBJ, 14 for NS)
%                               x 20 unique quiz images (1-10 = "yes", 11-20 = "no")
%                               x quiz dot x-pos (deg) 
%                               x quiz dot y-pos (deg)
%       * xy_coords_pix       : (double) is a 4D array with [x,y] coordinates in pixels 
%                               (relative to center of imagined core stimulus [0,0])
%                               Dimensions of array: 
%                               nr special core images (8 for Gabor/RDK/DOT/OBJ, 14 for NS)
%                               x 20 unique quiz images (1-10 = "yes", 11-20 = "no")
%                               x quiz dot x-pos (deg) 
%                               x quiz dot y-pos (deg)
%
% EXAMPLE:
% params.disp = vcd_getDisplayParams('7TAS_BOLDSCREEN32');
% params.exp  = vcd_getSessionParams('disp_name','7TAS_BOLDSCREEN32');
% params.stim = vcd_getStimParams('disp_name','7TAS_BOLDSCREEN32');
% all_quiz_dots = vcd_getImageryQuizDots(params, true, true);

% Check inputs
if nargin == 1
    update_info_file = false;
    verbose          = false;
elseif nargin == 2
    update_info_file = varargin{1};
    verbose = false;
elseif nargin == 3
    update_info_file = varargin{1};
    verbose = varargin{2};
else
    error('[%s]: Wrong number of inputs!',mfilename)
end

% Create output struct
img_quiz_dots = struct();

%%%%%% CREATE DOT IMAGE %%%%%%

% Create spatial support
x = 0:(params.stim.img.img_sz_pix - 1); % spatial support width (pixels)
y = 0:(params.stim.img.img_sz_pix - 1); % spatial support height (pixels)
x = x - x(end) / 2;
y = y - y(end) / 2;
[X, Y] = meshgrid(x, y);

% clean up
clear x y

% Center at zero for now (stimpresentation code will deal with x,y offset)
centerY = 0;
centerX = 0;

% Create quiz dot image (0.5-degree white dot on gray background)
quiz_dot              = double( (Y - centerY).^2 + (X - centerX).^2 <= (params.stim.img.quiz_dot_diam_pix/2).^2 ); % range is [0 1]
quiz_dot(quiz_dot==0) = 0.5;                        % insert gray background
quiz_dot              = 1+(quiz_dot*254);           % range is [1 255]
quiz_dot              = repmat(quiz_dot, [1,1,3]);  % replicate for rgb channel
quiz_dot              = uint8(quiz_dot);            % convert double to uint8
assert(isequal(unique(quiz_dot(:))',uint8([params.stim.bckgrnd_grayval, 255])));

% Alpha mask is same size as dot image but has slightly larger radius (+2
% pixels) than quiz dot image.
mask          = double( (Y - centerY).^2 + (X - centerX).^2 <= (2+ params.stim.img.quiz_dot_diam_pix/2).^2 ); % range is [0 1]
mask(mask==1) = 255;            % everything inside the alpha mask us visible, outside mask is invisible.
mask          = uint8(mask);    % convert double to uint8
assert(isequal(unique(mask(:))',uint8([0, 255])));

img_quiz_dots.im    = quiz_dot;
img_quiz_dots.mask  = mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sc = 1:length(params.exp.stimclassnames)
    
    % reset variables
    quiz_dot_xypos_pix  = []; % (8 or 14) nr of special core stim  x 20 quiz dot images (1-10 = yes, 11-20 = no) x 2 dots x 2 center coords ([x;y] pixels relative to center of stim loc)
    quiz_dot_xypos_deg  = []; % (8 or 14) nr of special core stim  x 20 quiz dot images (1-10 = yes, 11-20 = no) x 2 dots x 2 center coords ([x;y] deg relative to center of stim loc)
    quiz_dot_orient_deg = []; % (8 or 14) nr of special core stim  x 20 quiz dot images (1-10 = yes, 11-20 = no) x 2 dots x 1 angle (deg), where 0 = 12 o'clock.
    quiz_dot_overlap    = [];  % (8 or 14) nr of special core stim  x 20 quiz dot images (1 = yes, 2 = no);
    stimClass = params.exp.stimclassnames{sc};
    
    img_quiz_dots.(stimClass).special_core_stim_nr  = params.stim.(stimClass).unique_im_nrs_specialcore';
    img_quiz_dots.(stimClass).quiz_dot_stim_nr      = reshape(params.stim.(stimClass).unique_im_nrs_img_test, [], length(img_quiz_dots.(stimClass).special_core_stim_nr))'; % nr special core stim x 20 quiz dot image
    
    n_match    = sum(params.stim.(stimClass).imagery_quiz_images==1); % quiz dots overlap (1) or not (2)
    n_mismatch = sum(params.stim.(stimClass).imagery_quiz_images==2);
    
    % For visualization: red for mismatch, green for matching quiz dots
    cols = cmapturbo(n_match+n_mismatch);
   
    % load stim info file.
    if update_info_file
        d = dir(sprintf('%s*.mat',params.stim.(stimClass).stimfile));
        stim_mat_file = fullfile(d(end).folder,d(end).name);
        a = load(stim_mat_file,'info');
        q_info = a.info([],:); % create empty table with the same columns
        
        % add new columns
        q_info.img_quiz_dots_overlap = zeros(0);
        q_info.img_quiz_dot1_x_pix = zeros(0);
        q_info.img_quiz_dot1_y_pix = zeros(0);
        q_info.img_quiz_dot2_x_pix = zeros(0);
        q_info.img_quiz_dot2_y_pix = zeros(0);
        q_info.img_quiz_dot1_x_deg = zeros(0);
        q_info.img_quiz_dot1_y_deg = zeros(0);
        q_info.img_quiz_dot2_x_deg = zeros(0);
        q_info.img_quiz_dot2_y_deg = zeros(0);
        q_info.img_quiz_dot1_orient_deg = zeros(0);
        q_info.img_quiz_dot2_orient_deg = zeros(0);
        q_info.img_quiz_dot_specialcore_stim_nr = zeros(0);
        q_info.img_quiz_dot_filename            = {};

    end
    
% Create or extract imagery quiz dot locations for each stimulus class 
switch stimClass

    case {'gabor','rdk'}
        
        % What is the minimum distance between two quiz dots?
        min_dot_dist = params.disp.ppd * params.stim.(stimClass).imagery_min_dot_dist_deg;
        
        % Get the angles of the core stimuli that are imagined by the subject
        if strcmp(stimClass,'gabor')
            orient_deg = params.stim.gabor.ori_deg; % 8 tilts (deg)
        elseif  strcmp(stimClass,'rdk')
            orient_deg = params.stim.rdk.dots_direction; % 8 motion directions (deg)
        end
        max_radius_pix     = params.stim.(stimClass).imagery_sz_pix/2; % pixels
        dt = 0.25; % granularity of possible quiz dot tilts/motion directions

        % Whole dots are inside gabor, so we first create aperture slightly smaller than size of gabor/rdk (4-0.5 deg diameter):
        ap = params.disp.ppd * (params.stim.(stimClass).img_sz_deg - params.stim.(stimClass).imagery_aperture_buffer_deg); %
        
        % Define mismatch angles of two quiz dots (deviation from Gabor tilt/RDK motion direction/location of single dot)
        all_no_angles = cat(2, -1* (params.stim.(stimClass).imagery_minmax_anglemismatch(1):dt:params.stim.(stimClass).imagery_minmax_anglemismatch(2)), ...
                                      params.stim.(stimClass).imagery_minmax_anglemismatch(1):dt:params.stim.(stimClass).imagery_minmax_anglemismatch(2));
        
        % preallocate space
        quiz_dot_xypos_pix = NaN(length(orient_deg),length(params.stim.(stimClass).imagery_quiz_images),2,2);

        for tt = 1:length(orient_deg)

            % start with NaNs
            angles = NaN(length(params.stim.(stimClass).imagery_quiz_images),1);

            % add stimulus tilt orientation
            angles(params.stim.(stimClass).imagery_quiz_images==1,1) = orient_deg(tt); % in deg, og tilt
            
            % add mismatched quiz dot angles from stimulus tilt orientation
            % (randomly sample angles from a predefined list WITHOUT replacement)
            angles(params.stim.(stimClass).imagery_quiz_images==2,1) = mod(orient_deg(tt) + datasample(all_no_angles,n_mismatch,'Replace',false),360);

            quiz_dot_orient_deg(tt,:,:) = [angles, angles+180]; % special core stim nr x 20 quiz dot test images x 2 dot angles
            quiz_dot_overlap(tt,:)      = params.stim.(stimClass).imagery_quiz_images'; % aligned with orientation/motion direction (1) or not (2)
            
            % get xy-pos for imagined special core dot
            [x00,y00] = pol2cart(ones(1,2).*deg2rad(orient_deg(tt)-90),(params.stim.(stimClass).img_sz_pix/2).*[-1 1]); % note -90 to align 0 deg with 12 o'clock
            
            if verbose
                % For visualization: plot stimulus aperture and gabor tilt orientation / rdk motion direction
                figure(99); clf; xlim(max_radius_pix.*[-1.5 1.5]); ylim(max_radius_pix.*[-1.5 1.5])
                set(gca,'YDir','reverse')
                axis square;
                h = drawcircle('Center',[0,0],'Radius',max_radius_pix,'Color','black');
                h.InteractionsAllowed = 'none';
                h.HandleVisibility = 'off';
                h.Deletable = false;
                h.FaceAlpha = 0;
                hold on; plot(x00, y00,'bx-'); 
            end
            
            for qq = 1:length(params.stim.(stimClass).imagery_quiz_images)
                
                while 1
                    % For dot 1, choose first dot randomly positioned within the 4-0.5 deg circle;
                    r1     = sqrt(rand) * ap - (ap/2); % in pixels [-1.75 1.75 deg]
                    theta1 = rand * 2*pi; % in radians
                    x1     = r1.*cos(theta1); % in pixels
                    y1     = r1.*sin(theta1); % in pixels                   

                    % For dot 2, choose a random distance uniform between [mingap, 3.5];
                    r2 =  sqrt(rand) * (ap-min_dot_dist) + min_dot_dist;
                
                    theta2 = angles(qq) - 90; % note -90 to align 0 deg with 12 o'clock
                    
                    % calculate [x,y] coords for dot 2 (in pixels
                    x2  = x1 + r2 * cos(deg2rad(theta2));
                    y2  = y1 + r2 * sin(deg2rad(theta2));
                    
                    % see if second dot position is valid (within 4°-0.5° circle); 
                    % if not valid, repeat from beginning
                    if isInsideAperture([x2,y2], [0,0], ap/2)
                        
                        if verbose
                            % visualize quiz dots
                            scatter(x1,y1,70,cols(qq,:),'x','LineWidth',3)
                            plot([x1, x2],[y1, y2],'o-','color',cols(qq,:));
                            drawnow;
                        end
                        
                        break;
                    end
                end

                quiz_dot_xypos_pix(tt,qq,1,:) = round([x1 y1]); % quiz dot 1, [x,y] in pixels centered on 0, round to nearest pixel
                quiz_dot_xypos_pix(tt,qq,2,:) = round([x2 y2]); % quiz dot 2, [x,y] in pixels centered on 0, round to nearest pixel
                
                [th1,rho1] = cart2pol(x1,y1);
                [th2,rho2] = cart2pol(x2,y2);
                quiz_dot_xypos_deg(tt,qq,1,:) = [rad2deg(th1), rho1] ./params.disp.ppd; % quiz dot 1, [x,y] in deg where 0 = upper left corner (based on rounded pix values)
                quiz_dot_xypos_deg(tt,qq,2,:) = [rad2deg(th2), rho2] ./params.disp.ppd; % quiz dot 2, [x,y] in deg where 0 = upper left corner (based on rounded pix values)
                
                % Add quiz image properties to info file
                if update_info_file
                    idx = find(a.info.unique_im==img_quiz_dots.(stimClass).special_core_stim_nr(tt));
                    assert(a.info.is_specialcore(idx)==1)
                    tmp_info = a.info(idx,:); % copy stim info
                    
                    % update unique image nr
                    tmp_info.unique_im = img_quiz_dots.(stimClass).quiz_dot_stim_nr(tt,qq);
                    
                    tmp_info.img_quiz_dots_overlap = quiz_dot_overlap(tt,qq);
                    tmp_info.img_quiz_dot1_x_pix = quiz_dot_xypos_pix(tt,qq,1,1);
                    tmp_info.img_quiz_dot1_y_pix = quiz_dot_xypos_pix(tt,qq,1,2);
                    tmp_info.img_quiz_dot2_x_pix = quiz_dot_xypos_pix(tt,qq,2,1);
                    tmp_info.img_quiz_dot2_y_pix = quiz_dot_xypos_pix(tt,qq,2,2);
                    tmp_info.img_quiz_dot1_x_deg = quiz_dot_xypos_deg(tt,qq,1,1);
                    tmp_info.img_quiz_dot1_y_deg = quiz_dot_xypos_deg(tt,qq,1,2);
                    tmp_info.img_quiz_dot2_x_deg = quiz_dot_xypos_deg(tt,qq,2,1);
                    tmp_info.img_quiz_dot2_y_deg = quiz_dot_xypos_deg(tt,qq,2,2);
                    tmp_info.img_quiz_dot1_orient_deg = quiz_dot_orient_deg(tt,qq,1);
                    tmp_info.img_quiz_dot2_orient_deg = quiz_dot_orient_deg(tt,qq,2);
                    tmp_info.img_quiz_dot_specialcore_stim_nr = img_quiz_dots.(stimClass).special_core_stim_nr(tt);
                    tmp_info.img_quiz_dot_filename            = {NaN};
                    
                    q_info = cat(1,q_info,tmp_info);
                end
            end
        end
         
    case {'dot'}
        % Center of imagined dot is the midpoint of the arc
        % Arc length uniform random (-/+ [5-20] deg).
        % Dots can never cross vertical meridian
        
        % Get the angles of the special core stimuli that are imagined by the subject
        orient_deg = params.stim.dot.ang_deg(ismember(params.stim.dot.unique_im_nrs_core,params.stim.dot.unique_im_nrs_specialcore)); % 8 dot locations on iso-eccentric ring (deg)
        orient_deg_aligned = orient_deg-90; % -90 deg to align 0 deg with 12 o'clock
        eccen_deg  = params.stim.dot.imagery_sz_deg;
        eccen_pix  = params.stim.dot.imagery_sz_deg*params.disp.ppd;
        
        % Get the eccentricity of the imagined dot
        dot_radius_deg = params.stim.dot.imagery_dot_radius_deg;
        dot_radius_pix = params.stim.dot.imagery_dot_radius_deg*params.disp.ppd;
        
        dt = 0.25; % granularity of possible quiz dot angles
        
        % preallocate space
        quiz_dot_xypos_pix = NaN(length(orient_deg),length(params.stim.dot.imagery_quiz_images),2,2);
        
        for tt = 1:length(orient_deg)
        
            if verbose
                % For visualization: plot stimulus aperture and gabor tilt orientation / rdk motion direction
                figure(99); clf; xlim(eccen_pix.*[-1.5 1.5]); ylim(eccen_pix.*[-1.5 1.5])
                set(gca,'YDir','reverse')
                axis square;
                h = drawcircle('Center',[0,0],'Radius',eccen_pix,'Color','black');
                h.InteractionsAllowed = 'none';
                h.HandleVisibility = 'off';
                h.Deletable = false;
                h.FaceAlpha = 0;
                [x_core,y_core] = pol2cart(deg2rad(orient_deg_aligned(tt)),eccen_pix);  % get xypos in pixels for current special core stimulus
                hold on; plot(x_core, y_core,'k*','MarkerSize',10,'linewidth',3);
            end
            
            if orient_deg(tt) < 360 && orient_deg(tt) > 180 % left hemifield
                all_angles = 180:dt:360;
            elseif orient_deg(tt) < 180 && orient_deg(tt) > 0 % right hemifield
                all_angles = 0:dt:180;
            end
            quiz_dot_overlap(tt,:) = params.stim.dot.imagery_quiz_images;
            
            % Define straddling angles of two quiz dots (absolute angles in deg)
            minmax_angle_yes_ccw  = orient_deg(tt) - params.stim.dot.imagery_minmax_anglemismatch; % no mod(X,360) --> angles beyond 360 and less than 0 will be removed in the next line
            minmax_angle_yes_cw   = orient_deg(tt) + params.stim.dot.imagery_minmax_anglemismatch; % no mod(X,360) --> angles beyond 360 and less than 0 will be removed in the next line

            all_yes_angles_dot1 = all_angles(all_angles >= minmax_angle_yes_ccw(2) & all_angles <=  minmax_angle_yes_ccw(1)); % quiz dot 1 will always be ccw
            all_yes_angles_dot2 = all_angles(all_angles >= minmax_angle_yes_cw(1) & all_angles <=  minmax_angle_yes_cw(2)); % quiz dot 2 will always be cw
            
            % Define non-straddling angles of two quiz dots (absolute angles in deg)
            minmax_angle_no  = [orient_deg(tt) - params.stim.dot.imagery_minmax_anglemismatch(1), ...
                                orient_deg(tt) + params.stim.dot.imagery_minmax_anglemismatch(1)];
            
            % quiz dot 1 2 can be both ccw or cw
            all_no_angles_dot1 = setdiff(all_angles, minmax_angle_no(1):dt:minmax_angle_no(2));
            all_no_angles_dot2 = all_no_angles_dot1; 
            
            % %%% YES DOTS %%%
            while 1
                % sample straddling ("yes") angles for quiz dots [5:0.25:20]
                yes_angles = cat(2,datasample(all_yes_angles_dot1, n_match,'Replace',false)',  ... % ccw from special core stimulus loc
                                   datasample(all_yes_angles_dot2, n_match,'Replace',false)'); % cw from special core stimulus loc
                
               % Check 1) we leave enough distance between quiz dots 1 and 2
               %       2) we not sample the same angle for quiz dots 1 and 2
                if all(abs(circulardiff(yes_angles(:,1),yes_angles(:,2),360)) > params.stim.dot.imagery_min_dot_dist_deg) && ...  
                    all(~isequal(yes_angles(:,1),yes_angles(:,2))) 

                    % Check 3): quiz dots 1 and 2 are not in other hemifield
                    if orient_deg(tt) < 360 && orient_deg(tt) > 180 % left hemifield
                        if all(yes_angles(:,1) < 360) && all(yes_angles(:,1) > 180)
                            break;
                        end
                    elseif orient_deg(tt) < 180 && orient_deg(tt) > 0 % right hemifield
                        if all(yes_angles(:,1) < 180) && all(yes_angles(:,1) > 0)
                            break;
                        end
                    end
                end
            end
            
            % %%% NO DOTS %%%
            while 1
                % sample non-straddling ("no") angles of two quiz dots
                no_angles1 = datasample(all_no_angles_dot1, n_mismatch,'Replace',false);
                
                cw_no  = (no_angles1 > orient_deg(tt));
                ccw_no = (no_angles1 < orient_deg(tt));
                
                no_angles2         = zeros(size(no_angles1));
                no_angles2(cw_no)  = datasample(all_no_angles_dot2(all_no_angles_dot2 > orient_deg(tt)),sum(cw_no));
                no_angles2(ccw_no) = datasample(all_no_angles_dot2(all_no_angles_dot2 < orient_deg(tt)),sum(ccw_no));
                
                % combine angles of quiz dot 1 and 2
                no_angles = cat(1,no_angles1,no_angles2)';
                no_angles = sort(no_angles,2);
                
                % Check 1) we leave enough distance between quiz dots 1 and 2
                %       2) we not sample the same angle for quiz dots 1 and 2
                if all(abs(circulardiff(no_angles(:,1),no_angles(:,2),360)) > params.stim.dot.imagery_min_dot_dist_deg) && ...
                        all(~isequal(no_angles(:,1),no_angles(:,2)))
                    
                    % Check 3): quiz dots 1 and 2 are not in other hemifield
                    if orient_deg(tt) < 360 && orient_deg(tt) > 180 % left hemifield
                        if all(no_angles(:,1) < 360) && all(no_angles(:,1) > 180)
                            break;
                        end
                    elseif orient_deg(tt) < 180 && orient_deg(tt) > 0 % right hemifield
                        if all(no_angles(:,1) < 180) && all(no_angles(:,1) > 0)
                            break;
                        end
                    end
                end
            end
                
            % %%% COMBINE YES/NO QUIZ DOTS, CONVERT TO XY-POS in PIXELS %%%
            
            % start with 20 NaNs
            angles = NaN(length(params.stim.dot.imagery_quiz_images),2);
            
            % add 10 randomly sampled left/right quiz dot locations that straddle the core imagined dot
            angles(params.stim.dot.imagery_quiz_images==1,:) = yes_angles; % in deg, 0 deg = 12 o'clock
            
            % add 10 randomly sampled left/right quiz dot locations that do NOT straddle the core imagined dot
            angles(params.stim.dot.imagery_quiz_images==2,:) = no_angles; % in deg, 0 deg = 12 o'clock
            
            % record angles
            quiz_dot_orient_deg(tt,:,:) = round(angles); % in deg, location on iso-eccentric circle. dims: nr special core stim x 20 quiz dot stim x 2 dots
            
            % convert angles/eccen to xy-pos in pixel space
            [x00,y00] = pol2cart(deg2rad(angles-90),params.disp.ppd*params.stim.dot.imagery_sz_deg.*ones(size(angles)));  % in pixels
            
            x1 = x00(:,1);
            y1 = y00(:,1);
            x2 = x00(:,2);
            y2 = y00(:,2);
            
            % final check if second dot position is not crossing the vertical meridian
            if orient_deg(tt) > 180 && orient_deg(tt) < 360
                % if core dot is in the left hemisphere
                assert(all(x1 < 0));
            elseif orient_deg(tt) > 0 && orient_deg(tt) < 180
                % if core dot is in the right hemisphere
                assert(all(x1 > 0));
            else
                error('[%s]: Check xy-pos of imagery quiz dots for SINGLE DOT',mfilename)
            end

            % record dot location
            quiz_dot_xypos_pix(tt,1:length(params.stim.dot.imagery_quiz_images),1,:) = round([x1 y1]); % quiz dot 1, [x,y] center location in pixels; (rounded to nearest pixel)
            quiz_dot_xypos_pix(tt,1:length(params.stim.dot.imagery_quiz_images),2,:) = round([x2 y2]); % quiz dot 2, [x,y] center location in pixels; (rounded to nearest pixel)
            
            [th1,rho1] = cart2pol(x1,y1);
            [th2,rho2] = cart2pol(x2,y2);
            quiz_dot_xypos_deg(tt,1:length(params.stim.dot.imagery_quiz_images),1,:) = [rad2deg(th1), rho1]./params.disp.ppd;
            quiz_dot_xypos_deg(tt,1:length(params.stim.dot.imagery_quiz_images),2,:) = [rad2deg(th2), rho2]./params.disp.ppd;
            
            % Add quiz image properties to info file
            if update_info_file
                idx = find(a.info.unique_im==img_quiz_dots.(stimClass).special_core_stim_nr(tt));
                assert(a.info.is_specialcore(idx)==1)
                tmp_info = repmat(a.info(idx,:), length(params.stim.dot.imagery_quiz_images), 1); % copy stim info
                
                % update unique image nr
                tmp_info.unique_im(1:length(params.stim.dot.imagery_quiz_images)) = img_quiz_dots.(stimClass).quiz_dot_stim_nr(tt,:);
                
                tmp_info.img_quiz_dots_overlap(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_overlap(tt,:)';
                tmp_info.img_quiz_dot1_x_pix(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_xypos_pix(tt,:,1,1)';
                tmp_info.img_quiz_dot1_y_pix(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_xypos_pix(tt,:,1,2)';
                tmp_info.img_quiz_dot2_x_pix(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_xypos_pix(tt,:,2,1)';
                tmp_info.img_quiz_dot2_y_pix(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_xypos_pix(tt,:,2,2)';
                tmp_info.img_quiz_dot1_x_deg(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_xypos_deg(tt,:,1,1)';
                tmp_info.img_quiz_dot1_y_deg(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_xypos_deg(tt,:,1,2)';
                tmp_info.img_quiz_dot2_x_deg(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_xypos_deg(tt,:,2,1)';
                tmp_info.img_quiz_dot2_y_deg(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_xypos_deg(tt,:,2,2)';
                tmp_info.img_quiz_dot1_orient_deg(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_orient_deg(tt,:,1)';
                tmp_info.img_quiz_dot2_orient_deg(1:length(params.stim.dot.imagery_quiz_images)) = quiz_dot_orient_deg(tt,:,2)';
                tmp_info.img_quiz_dot_specialcore_stim_nr(1:length(params.stim.dot.imagery_quiz_images)) = repmat(img_quiz_dots.(stimClass).special_core_stim_nr(tt),length(params.stim.dot.imagery_quiz_images),1);
                tmp_info.img_quiz_dot_filename(1:length(params.stim.dot.imagery_quiz_images))            = repmat({NaN},length(params.stim.dot.imagery_quiz_images),1);
                q_info = cat(1,q_info,tmp_info);
            end
            
            if verbose
                % visualize quiz dots
                scatter(x1,y1,70,cols,'x','LineWidth',3);
                scatter(x2,y2,70,cols,'o','LineWidth',3);
                drawnow; pause(0.5);
            end
        end
        
    case {'obj','ns'}
        
        if strcmp(stimClass,'ns')
            [filenames, unique_im_nrs] = vcd_getNSfilenames;
            img_filenames = filenames(ismember(unique_im_nrs, params.stim.(stimClass).unique_im_nrs_img_test));
        end
        
        d = dir(params.stim.(stimClass).imagery_raw_dot_folder);

        if isempty(d)
                error('[%s]: can''t find imagery quiz dot test image file',mfilename)
        end
        assert(isequal(length(d),length(params.stim.(stimClass).unique_im_nrs_specialcore)))
        
        for tt = 1:length(d)
            fnames_yes = dir(fullfile(d(tt).folder,d(tt).name,'*yes*.png'));
            fnames_no  = dir(fullfile(d(tt).folder,d(tt).name,'*no*.png'));
            assert(isequal(length(fnames_yes),sum(params.stim.(stimClass).imagery_quiz_images==1)))
            assert(isequal(length(fnames_no),sum(params.stim.(stimClass).imagery_quiz_images==2)))
            assert(isequal(length(fnames_no),length(fnames_yes)));

            fnames_yesno = cat(1,fnames_yes,fnames_no);
            quiz_dot_overlap(tt,:) = params.stim.(stimClass).imagery_quiz_images;
            
            for qq = 1:length(fnames_yesno)
                
                % Read in images
                im = imread(fullfile(fnames_yesno(qq).folder,fnames_yesno(qq).name));
                im2 = double(im(:,:,1)); % take first channel, convert uint8 to double
                im3 = (im2 > 255/2);  % dots are 0 on background 1
                
                if strcmp(stimClass,'ns') % resize image from OG resolution (425 x 425 pixels) to 7TAS resolution matching 8.4 x 8.4 deg (741 x 741 pixels)
                    if size(im3,1) && (size(im3,1) == size(im3,2))
                        im3 = im3(1:425, 1:425,:); % crop if photoshop editing added a single row/column of pixels.
                    end
                    im3_rz = imresize(im3,params.stim.ns.dres);
                    
                    % check for clipping
                    assert( min(im3_rz(:))>=0 && max(im3_rz(:))<=1)
                else
                    im3_rz = im2;
                end

                idx = find(im3_rz==0);
                rowii = mod2(idx,size(im3_rz,1));  % 1 through 425, from top to bottom
                colii = ceil(idx/size(im3_rz,1));  % 1 through 425, from left to right
                
                [~,centers] = runkmeans([rowii(:) colii(:)],2); % [y; x] (or [row; col]) in decimal matrix coordinates relative to upper left pixel
                
                assert(isequal(size(im3_rz,1),params.stim.(stimClass).imagery_sz_pix)); % check if image size is the same as expected quiz dot image size in params 
                centers_rel                   = [centers(:,2), centers(:,1)] - params.stim.(stimClass).imagery_sz_pix/2; % [x, y] relative to center of image 
                quiz_dot_xypos_pix(tt,qq,:,:) = centers_rel;
                quiz_dot_xypos_deg(tt,qq,:,:) = centers_rel ./ params.disp.ppd;
                
                % calculate angle of line that connects quiz dots and vertical line
                ref_vec = [0, params.stim.(stimClass).imagery_sz_pix/2];
                im_vec  = abs([diff(centers(:,2)), diff(centers(:,1))]);
                
                u1 = ref_vec/norm(ref_vec);
                u2 = im_vec/norm(im_vec);
                angle = rad2deg(acos(dot(u1,u2)));
                
                if (centers_rel(1,1) < centers_rel(2,1)) && (centers_rel(1,2) < centers_rel(2,2)) || ...
                    (centers_rel(1,1) > centers_rel(2,1)) && (centers_rel(1,2) > centers_rel(2,2))
                    angle = -angle; % flip sign
                    angle = mod(angle, 360);
                end
                
                if verbose
                    % plot centers for debugging purposes
                    figure(101); clf; hold on;
                    imagesc(im3_rz);
                    set(gca,'YDir','reverse');
                    scatter(centers(:,2),centers(:,1),'cx');
                    axis image;
                    
                    [x00,y00] = pol2cart(deg2rad([angle - 90, angle + 90]),[1 1]);
                    x00 = x00*params.stim.(stimClass).imagery_sz_pix/2 + (params.stim.(stimClass).imagery_sz_pix/2);  % in pixels
                    y00 = y00*params.stim.(stimClass).imagery_sz_pix/2 + (params.stim.(stimClass).imagery_sz_pix/2);  % in pixels
                    hold on; plot(x00,y00,'b*-','LineWidth',3);
                end
                
                % record angles
                quiz_dot_orient_deg(tt,qq,:,:) = [angle, angle+180]; % in deg, angle relative to center of scene. dims: nr special core stim x 20 quiz dot stim x 2 dots (same angle for both dots on left or right side).
                
                % Add quiz image properties to info file
                if update_info_file
                    idx = find(a.info.unique_im==img_quiz_dots.(stimClass).special_core_stim_nr(tt));
                    assert(a.info.is_specialcore(idx)==1)
                    tmp_info = a.info(idx,:); % copy stim info
                    
                    % update unique image nr
                    tmp_info.unique_im = img_quiz_dots.(stimClass).quiz_dot_stim_nr(tt,qq);
                    
                    tmp_info.img_quiz_dots_overlap = quiz_dot_overlap(tt,qq);
                    tmp_info.img_quiz_dot1_x_pix = quiz_dot_xypos_pix(tt,qq,1,1);
                    tmp_info.img_quiz_dot1_y_pix = quiz_dot_xypos_pix(tt,qq,1,2);
                    tmp_info.img_quiz_dot2_x_pix = quiz_dot_xypos_pix(tt,qq,2,1);
                    tmp_info.img_quiz_dot2_y_pix = quiz_dot_xypos_pix(tt,qq,2,2);
                    tmp_info.img_quiz_dot1_x_deg = quiz_dot_xypos_deg(tt,qq,1,1);
                    tmp_info.img_quiz_dot1_y_deg = quiz_dot_xypos_deg(tt,qq,1,2);
                    tmp_info.img_quiz_dot2_x_deg = quiz_dot_xypos_deg(tt,qq,2,1);
                    tmp_info.img_quiz_dot2_y_deg = quiz_dot_xypos_deg(tt,qq,2,2);
                    tmp_info.img_quiz_dot1_orient_deg = quiz_dot_orient_deg(tt,qq,1);
                    tmp_info.img_quiz_dot2_orient_deg = quiz_dot_orient_deg(tt,qq,2);
                    tmp_info.img_quiz_dot_specialcore_stim_nr = img_quiz_dots.(stimClass).special_core_stim_nr(tt);
                    if strcmp(stimClass,'ns')
                        tmp_info.img_quiz_dot_filename            = {img_filenames(qq + ((tt-1)*20))};
                    else
                        tmp_info.img_quiz_dot_filename            = {NaN};
                    end
                    q_info = cat(1,q_info,tmp_info);
                end
                

                
            end
        end
end
    

img_quiz_dots.(stimClass).xy_coords_pix       = quiz_dot_xypos_pix;
img_quiz_dots.(stimClass).xy_coords_deg       = quiz_dot_xypos_deg;
img_quiz_dots.(stimClass).quiz_dot_orient_deg = quiz_dot_orient_deg;
img_quiz_dots.(stimClass).quiz_dot_overlap    = quiz_dot_overlap;

if update_info_file
    % Concatenate new columns to old info table
    info = q_info([],:);
    info(1:size(a.info,1),1:size(a.info,2))     = a.info;
    info{:,(size(a.info,2)+1):(size(info,2)-1)} = NaN; % replace zeros with NaNs
    info{:, size(info,2)} = NaN(size(info,1),1);
    info = cat(1,info, q_info);
    
    % Store images and info table into a new file
    fprintf('[%s]:Storing updated info table and csv file..\n',mfilename)
    new_file = fullfile(sprintf('%s_%s.mat',params.stim.(stimClass).stimfile,datestr(now,30)));
    copyfile(stim_mat_file,new_file);
    save(new_file,'info','-append');
    writetable(info, fullfile(sprintf('%s_%s.csv',params.stim.(stimClass).infofile,datestr(now,30))))
end


end

% Store images and info table into a new file
fprintf('[%s]: Storing imagery quiz dot info..\n',mfilename)
fname = fullfile(fileparts(params.stim.gabor.infofile), sprintf('imagery_quiz_dot_info_%s.mat',datestr(now,30)));
save(fname,'img_quiz_dots');



