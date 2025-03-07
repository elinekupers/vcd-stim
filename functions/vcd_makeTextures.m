function [im_tex, im_rect, txt_tex, txt_rect, fix_tex, fix_rect, framecolor] = ...
    vcd_makeTextures(params, win, rect, flipfun, scan, timing, taskscript)

%% Create text blocks & rects



%% Create background and fixation textures prior to exp onset (as we need them throughout the experiment)
% bckground_rect = CenterRect([0 0 round(size(scan.bckground,1)) round(size(scan.bckground,2))],rect);
bckground_rect    = rect;
bckrgound_texture = Screen('MakeTexture', win, feval(flipfun,scan.bckground));

% make fixation dot texture
fix_texture_thin_full   = {};
fix_texture_thick_full  = {};
fix_texture_thick_left  = {};
fix_texture_thick_right = {};
fix_rect_thin = CenterRect([0 0 round(size(scan.fix_im,1)) round(size(scan.fix_im,2))],rect);
fix_rect_thick = CenterRect([0 0 round(size(scan.fix_im,1)) round(size(scan.fix_im,2))],rect);
for ll = 1:size(scan.fix_im,4) % loop over luminance values
    fix_texture_thin_full{ll} = Screen('MakeTexture',win,feval(flipfun,  cat(3, scan.fix_im(:,:,:,ll,1), scan.fix_alpha_mask.*params.stim.fix.dotopacity)));
    fix_texture_thick_full{ll} = Screen('MakeTexture',win,feval(flipfun,  cat(3, scan.fix_im(:,:,:,ll,2), scan.fix_alpha_mask.*params.stim.fix.dotopacity)));
    fix_texture_thick_left{ll} = Screen('MakeTexture',win,feval(flipfun,  cat(3, scan.fix_im(:,:,:,ll,3), scan.fix_alpha_mask.*params.stim.fix.dotopacity)));
    fix_texture_thick_right{ll} = Screen('MakeTexture',win,feval(flipfun,  cat(3, scan.fix_im(:,:,:,ll,4), scan.fix_alpha_mask.*params.stim.fix.dotopacity)));
end

%% create single vector for fix textures

% Preallocate space
im_tex  = cell(length(timing.trig_stim),1);
im_rect = im_tex;
txt_tex = im_tex;
txt_rect = im_tex;
fix_tex = {};
fix_rect = {};
framecolor = im_tex;
% im_global_alpha = im_tex;

for frame = 1:length(timing.trig_stim)
    
    framecolor{frame} = 255*ones(1,3); % <framecolor> can also be size(<frameorder>,2) x 1 with values in [0,1] indicating an alpha change.

    blockID = timing.trig_block(frame);
    
    % set up fixation dot textures
    opacity_idx = timing.trig_fix(frame);
    
    if isnan(timing.seq_spatial_cue{frame})
        fix_tex{frame} = fix_texture_thin_full{opacity_idx};
        fix_rect{frame} = fix_rect_thin;
    else
        switch timing.trig_spatial_cue(frame)
            case 1
                fix_tex{frame} = fix_texture_thick_left{opacity_idx};
                fix_rect{frame} = fix_rect_thick;
            case 2
                fix_tex{frame} = fix_texture_thick_right{opacity_idx};
                fix_rect{frame} = fix_rect_thick;
            case 0
                fix_tex{frame} = fix_texture_thick_full{opacity_idx};
                fix_rect{frame} = fix_rect_thick;
        end
    end
    
    switch timing.trig_stim(frame,1)
        
        % 0  : pre/post blank
        % 93 : exp_session.miniblock.response_ID
        % 94 : exp_session.miniblock.trial_start_ID
        % 95 : exp_session.miniblock.spatial_cue_ID
        % 96 : exp_session.miniblock.delay_ID
        % 97 : exp_session.miniblock.task_cue_ID
        % 98 : exp_session.miniblock.ITI_ID
        % 99 : exp_session.miniblock.IBI_ID
   
        % Draw background + thin fix dot on top
        case {0, 93, 94, 95, 96, 98, 99} 
            
             % DrawTextures 
             % * TexturePointers  need to be: n vector (where n is the number of textures)
             % * DestinationRects need to be: 4 row x n columns (where n is the number of textures)            
            im_tex{frame,:} = cat(1, bckrgound_texture, fix_tex);
            im_rect{frame,:} = cat(1, bckground_rect, fix_rect);
%             im_global_alpha{frame} = 1;
            
        case 97 % task_cue_ID  
            
            script = taskscript{~cellfun(@isempty, regexp(taskscript,sprintf('%02d',blockID),'match'))};
            [task_instr, task_rect] = vcd_getInstructionText(params, script, rect);
            
            im_tex{frame,:} = cat(1, bckrgound_texture, fix_texture_thick_full{opacity_idx});
            im_rect{frame,:} = cat(1, bckground_rect, fix_rect_thick);
            
            txt_tex{frame,:} = task_instr;
            txt_rect{frame,:} = task_rect;

            % Draw stimulus textures
        % 1-30 = all 2 peripheral stimulus aperture stim-task crossings:
        %   01-xx = Gabors
        %   xx-xx = RDKs
        %   xx-xx = Simple dot
        %   xx-30 = Complex objects
    
    case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30} 
            if blockID <= 30
                % trig_seq_exp_im_w_cd is a cell with dims: frames x 1, where each cell has 1 or 2 sides (1:l, 2:r)
                for side = 1:length(find(~cellfun(@isempty, timing.trig_seq_exp_im_w_cd{frame})))
                    
                    % add mask if we have one
                    if ~isempty(timing.trig_seq_exp_im_mask{frame}) || ~isempty(timing.trig_seq_exp_im_mask{frame}{1})
                        if  length(timing.trig_seq_exp_im_mask{frame})==1 && ~isequalwithequalnans(timing.trig_seq_exp_im_mask{frame}{1},NaN)
                            
                            stim_tex = feval(flipfun, cat(3, timing.trig_seq_exp_im_w_cd{frame}{side}, timing.trig_seq_exp_im_mask{frame}{1}));
                            stim_rect = scan.rects{frame,side};
                        else
                            stim_tex = feval(flipfun, timing.trig_seq_exp_im_w_cd{frame}{side});
                            stim_rect = scan.rects{frame,side};
                        end
                        
                    else
                        stim_tex = feval(flipfun, timing.trig_seq_exp_im_w_cd{frame}{side});
                        stim_rect = scan.rects{frame,side};
                    end
                    
                    stim_texture = Screen('MakeTexture',win, txttemp);
                    
                    
                    im_tex{frame,:} = cat(1, bckrgound_texture, stim_texture, fix_tex{frame});
                    im_rect{frame,:} = cat(1, bckground_rect, stim_rect, fix_rect{frame});
                                        
                end
            
            else % 31-39 = Natural Scene stim-task crossings ({31,32,33,34,35,36,37,38,39})
                txttemp = feval(flipfun,timing.trig_seq_exp_im_w_cd{frame}{1});
                stim_rect = scan.rects{frame,1};
                stim_texture = Screen('MakeTexture',win, txttemp);
                
                im_tex{frame,:} = cat(1, bckrgound_texture, stim_texture, fix_tex{frame});
                im_rect{frame,:} = cat(1, bckground_rect, stim_rect, fix_rect{frame});
            end
    end  
end

return




