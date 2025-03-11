function timing = vcd_getImageTiming_framelocked30Hz(params, subj_run,im_seq_order, exp_im, fix_im, exp_im_mask) 

% Fixation order and fixation
fixsoafun = @() round(params.stim.fix.dotmeanchange + (params.stim.fix.dotchangeplusminus*(2*(rand-.5))))*params.stim.framedur_s;

% Contrast decrement gaussian time window onset
cdsoafun = @() round(params.stim.cd.meanchange + params.stim.cd.changeplusminus*(2*(rand-.5)))*params.stim.framedur_s;


%% TIMING

seq_stim = {}; 
seq_timing = []; 
spatial_cue = []; 
seq_block = [];

seq_exp_im = {};
seq_exp_im_mask = {};

% 6 fields (name, ID, within_session_repeat, trial, trial_type, timing)
% 8 blocks: 1:run, 2:block, 3:stimtaskID, 4:unique_im, 5:spatial_cue, 6:onset_time, 7:event_dur, 8:run_time
cellblock = squeeze(struct2cell(subj_run.block));

for bb = 1:size(cellblock,2)
    
    block_ID = cellblock{2,bb};
    
    tmp_timing = cellblock{6,bb};
    im_nr = tmp_timing.unique_im;
    if ~iscell(im_nr)
        % convert to cell
        im_nr = num2cell(im_nr);
    end

    % check for nans
    if sum(cell2mat(cellfun(@isempty, cellfun(@isnan, im_nr,'UniformOutput',false),'UniformOutput',false)))>0
        for xi = 1:length(im_nr)
            tmp_im = im_nr{xi};
            if isnan(tmp_im)
                tmp_im(isnan(tmp_im))=0;
                im_nr{xi} = tmp_im;
            end
            clear tmp_im
        end  
    end
    
    if bb == 1
        cumultime = 0;
    end
    % add second stim if it is a double epoch trial
    for tt = 1:length(im_nr)
 
        
        if im_nr{tt}==0 % blockstart
            seq_stim = cat(1, seq_stim, im_nr(tt));
            seq_timing = cat(1,seq_timing, cumultime);
            seq_block = cat(1,seq_block,0);
            cumultime = cumultime + tmp_timing.event_dur(tt);
            seq_exp_im = cat(1,seq_exp_im, {0});
        
        elseif im_nr{tt}==97 % 97: task cue
            seq_stim = cat(1, seq_stim, im_nr(tt));
            seq_timing = cat(1,seq_timing, cumultime);
            seq_block = cat(1,seq_block,block_ID);
            cumultime = cumultime + tmp_timing.event_dur(tt);
            seq_exp_im = cat(1,seq_exp_im, {0});

        elseif any(im_nr{tt}==([98,99])) % 98: iti or 99: ibi
            seq_stim = cat(1, seq_stim, im_nr(tt));
            seq_timing = cat(1,seq_timing, cumultime);

            cumultime = cumultime + tmp_timing.event_dur(tt);
            
            if im_nr{tt}==98
                seq_block = cat(1,seq_block,block_ID);
            elseif im_nr{tt}==99
                seq_block = cat(1,seq_block,0);
            end
            seq_exp_im = cat(1,seq_exp_im, {0});

        % 93: response ID // 94: trial_start_ID // 95: spatial_cue_ID // 96:delay_ID = 96    
        else 
            idx = tt/2;
            assert(isequal(cell2mat(squeeze(im_seq_order(bb,idx,:,1)))', im_nr{tt}))
            
            trial_start_im       = {params.exp.miniblock.trial_start_ID}; % 96: trial start 
            trial_start_run_time = cumultime;
            cumultime            = cumultime + params.exp.trial.start_cue_dur;
            
            spatial_cue_im       = {params.exp.miniblock.spatial_cue_ID};
            spatial_cue_run_time = cumultime;
            cumultime            = cumultime + params.exp.trial.spatial_cue_dur;
            
            stim_im              = im_nr(tt);
            stim_run_time        = cumultime; 
            cumultime = cumultime + params.exp.trial.stim_array_dur;

            
            if cellblock{5,bb} == 2
                delay_im            = {params.exp.miniblock.delay_ID};
                delay_run_time      = cumultime; 
                cumultime           = cumultime + params.exp.trial.delay_dur; 
                
                query_im            = {catcell(2,squeeze(im_seq_order(bb,idx,:,2)))};
                query_run_time      = cumultime; 
                cumultime           = cumultime + params.exp.trial.stim_array_dur;
                
                response_im         = {params.exp.miniblock.response_ID};
                response_run_time   = cumultime; 
                cumultime           = cumultime + params.exp.trial.response_win_dur;
                
                seq_stim = cat(1, seq_stim, ...
                    trial_start_im,...
                    spatial_cue_im, ...
                    stim_im, ...
                    delay_im,...
                    query_im, ...
                    response_im);

                seq_timing = cat(1,seq_timing, ...
                    trial_start_run_time, ...
                    spatial_cue_run_time, ...
                    stim_run_time, ...
                    delay_run_time, ...
                    query_run_time, ...
                    response_run_time);
                
                seq_block = cat(1,seq_block,repmat(block_ID,6,1));
                
                % if there is a mask, Concatenate stim and alpha mask 
                if ~isempty(squeeze(exp_im_mask(bb,idx,1,1))) 
                    
                    for jj = 1:size(exp_im,4)
                        tmp_im = squeeze(exp_im(bb,idx,:,jj));
                        tmp_mask = squeeze(exp_im_mask(bb,idx,:,jj));
                        numImages = sum(~cellfun(@isempty, tmp_im));
                        alphaIm = tmp_im;
                        for nn = 1:numImages
                            if ndims(tmp_im{nn})==4
                                alphaIm{nn,jj} = cat(3, tmp_im{nn}, repmat(tmp_mask{nn}, [1 1 1 size(tmp_im{nn},4)]));
                            else
                                alphaIm{nn,jj} = cat(3, tmp_im{nn}, tmp_mask{nn});
                            end
                        end
                    end
                    

                    seq_exp_im = cat(1,seq_exp_im, repmat({0},2,1), ... % trial start, spatial cue
                                                   {alphaIm(:,1)'}, ... % stim interval 1
                                                   {0}, ... % delay period
                                                   {alphaIm(:,2)'}, ... % stim interval 2
                                                   {0}); % response 
                % Otherwise we just use the images                           
                else                           
                    seq_exp_im = cat(1,seq_exp_im, repmat({0},2,1), ... % trial start, spatial cue
                                               {squeeze(seq_exp_im(bb,idx,:,1))'}, ... % stim interval 1
                                               {0}, ... % delay period
                                               {squeeze(seq_exp_im(bb,idx,:,2))'}, ... % stim interval 2
                                               {0}); % response 
                end
                
            else
                response_im = {params.exp.miniblock.response_ID};
                response_run_time = cumultime;
                cumultime = cumultime + params.exp.trial.response_win_dur;
                
                seq_stim = cat(1, seq_stim, ...
                    trial_start_im,...
                    spatial_cue_im, ...
                    stim_im, ...
                    response_im);

                seq_timing = cat(1,seq_timing, ...
                    trial_start_run_time, ...
                    spatial_cue_run_time, ...
                    stim_run_time, ...
                    response_run_time);
                
                seq_block = cat(1,seq_block,repmat(block_ID,4,1));

                % if there is a mask, Concatenate stim and alpha mask 
                if ~isempty(squeeze(exp_im_mask(bb,idx,1,1)))
                    tmp_im = squeeze(exp_im(bb,idx,:,1));
                    tmp_mask = squeeze(exp_im_mask(bb,idx,:,1));
                    numImages = sum(~cellfun(@isempty, tmp_im));
                    alphaIm = tmp_im;
                    for nn = 1:numImages
                        if ndims(tmp_im{nn})==4
                            alphaIm{nn} = cat(3, tmp_im{nn}, repmat(tmp_mask{nn}, [1 1 1 size(tmp_im{nn},4)]));
                        else
                            alphaIm{nn} = cat(3, tmp_im{nn}, tmp_mask{nn});
                        end
                    end
                    seq_exp_im = cat(1,seq_exp_im, repmat({0},2,1), {alphaIm'}, {0});
                
                % Otherwise we just use the images in stim interval 1    
                else
                    seq_exp_im = cat(1,seq_exp_im, repmat({0},2,1), {squeeze(exp_im(bb,idx,:,1))'},{0});
                end

            end
            spatial_cue = cat(1,spatial_cue, [spatial_cue_im, tmp_timing.spatial_cue(tt)]); % same as params.exp.trial.spatial_cue_dur?;
        end
       
    end
end

seq_stim        = cat(1, seq_stim(2:end), 0);
seq_timing      = cat(1, seq_timing(2:end), cumultime);
seq_block       = cat(1, seq_block(2:end), 0);
seq_exp_im      = cat(1, seq_exp_im(2:end), {0});

timing.seq_stim        = seq_stim;
timing.seq_spatial_cue = spatial_cue;
timing.seq_timing      = seq_timing;
timing.seq_block       = seq_block;
timing.seq_exp_im      = seq_exp_im;

%% Convert sequence of events from seconds into frames

trig_timing = [0:params.stim.framedur_s:seq_timing(end)]'; % seconds

trig_stim   = zeros(size(trig_timing,1),2);
trig_block  = zeros(size(trig_timing,1),1);
trig_spatial_cue = zeros(size(trig_timing,1),1);
trig_seq_exp_im = cell(size(trig_timing,1),1);

spc_idx = 1;
for tt = 1:length(seq_timing)
    event_time = seq_timing(tt);
    if tt == length(seq_timing)
        next_event = event_time;
    else
        next_event = seq_timing(tt+1);
    end
    
    [t_diff,t_idx] = min(abs(trig_timing-event_time));
    if isempty(t_idx) || ~nearZero(t_diff)
        error('[%s]: Event_time doesn''t match monitor refresh rate', mfilename)
    end
    
    event_dur = next_event-event_time;
    if event_dur > 0
        event_dur = round(event_dur/params.stim.framedur_s);
        t_idx_total = t_idx:(t_idx+event_dur-1);
    else
        t_idx_total = t_idx;
    end
    
    if length(seq_stim{tt})==2
        trig_stim(t_idx_total,:) = repmat(seq_stim{tt},length(t_idx_total),1);
    else
        trig_stim(t_idx_total,:) = repmat(seq_stim{tt},length(t_idx_total),2);
    end
    
%     disp([t_idx_total(1),t_idx_total(end)])
    trig_block(t_idx_total) = repmat(seq_block(tt), length(t_idx_total),1);
    if any(seq_stim{tt}==95)
        trig_spatial_cue(t_idx_total,:) = repmat(spatial_cue{spc_idx,2}, length(t_idx_total),1);
        spc_idx = spc_idx+1;
    end

    if iscell(seq_exp_im{tt}) && ndims(seq_exp_im{tt}{1}) == 4 % we are dealing with rdk, which has time dim
        rdk_images = [squeeze(mat2cell(seq_exp_im{tt}{1}, size(seq_exp_im{tt}{1},1),size(seq_exp_im{tt}{1},2),size(seq_exp_im{tt}{1},3),ones(1,size(seq_exp_im{tt}{1},4)))), ...
            squeeze(mat2cell(uint8(seq_exp_im{tt}{2}), size(seq_exp_im{tt}{2},1),size(seq_exp_im{tt}{2},2),size(seq_exp_im{tt}{2},3), ones(1,size(seq_exp_im{tt}{2},4))))];
        for oob = 1:size(rdk_images,1)
            rdk_images2{oob} = [rdk_images(oob,1),rdk_images(oob,2)];
        end
        if size(rdk_images2,1) < size(rdk_images2,2)
            rdk_images2 = rdk_images2';
        end
        trig_seq_exp_im(t_idx_total) = rdk_images2(1:length(t_idx_total),:);
    else
        
%         if isempty(seq_exp_im{tt}) || isequalwithequalnans(seq_exp_im_mask{tt},NaN)
            trig_seq_exp_im(t_idx_total) = repmat(seq_exp_im(tt), length(t_idx_total),1);
%         end
    end
end


timing.trig_stim            = trig_stim;
timing.trig_timing          = trig_timing;
timing.trig_block           = trig_block;
timing.trig_spatial_cue     = trig_spatial_cue;
timing.trig_seq_exp_im      = trig_seq_exp_im;

%% FIXATION DOT SEQUENCE
fudge   = 5; % shave off 2 second at the end to avoid going over run dur

fix_tmp    = [0]; 
fix_timing = 0; 
fix_seq = [randi(length(params.stim.fix.dotlum),1)];

while fix_timing < (trig_timing(end)-fudge)
    fix_tmp = [fix_tmp, fixsoafun()];
    fix_timing = cumsum(fix_tmp);
    fix_seq = [fix_seq, randi(length(params.stim.fix.dotlum),1)];
end

trig_fix_im = cell(size(trig_timing,1),1);
trig_fix = zeros(size(trig_timing,1),1);

for tt = 1:length(fix_timing)
    event_type = trig_stim(tt,1);
    event_time = fix_timing(tt);
    if tt==length(fix_timing)
        next_event = trig_timing(end);
    else
        next_event = fix_timing(tt+1);
    end
    
    [t_diff,t_idx] = min(abs(trig_timing-event_time));
    if isempty(t_idx) || ~nearZero(t_diff)
        error('[%s]: Event_time doesn''t match monitor refresh rate', mfilename)
    end
    
    event_dur = next_event-event_time;
    event_dur = round(event_dur/params.stim.framedur_s);
    if  tt==length(fix_timing)
        event_dur = event_dur+1;
    end
    
    t_idx_total = t_idx:(t_idx+event_dur-1);
    
    trig_fix(t_idx_total,:)  = repmat(fix_seq(tt),length(t_idx_total),1);
    
    if any(event_type==[0,98,99]) % pre/post-block, iti, ibi 
        trig_fix_im(t_idx_total) = repmat({fix_im(:,:,:,fix_seq(tt), 1)},length(t_idx_total),1); % thin rim
    elseif any(event_type==[93,94,96,97]) % trial start
        trig_fix_im(t_idx_total) = repmat({fix_im(:,:,:,fix_seq(tt), 2)},length(t_idx_total),1); % thick rim
    elseif event_type==95 && trig_spatial_cue(tt)==1
        trig_fix_im(t_idx_total) = repmat({fix_im(:,:,:,fix_seq(tt), 3)},length(t_idx_total),1); % thick rim left
    elseif event_type==95 && trig_spatial_cue(tt)==2
        trig_fix_im(t_idx_total) = repmat({fix_im(:,:,:,fix_seq(tt), 4)},length(t_idx_total),1); % thick rim right
    end
end


timing.seq_fix      = fix_seq;
timing.trig_fix     = trig_fix;
timing.trig_fix_im  = trig_fix_im;


%% Contrast change event sequence
cdID = [];
for nn = 1:length(subj_run.block)
    if ~isempty(regexp(subj_run.block(nn).name,'cd-*','ONCE'))
        cdID(nn) = 1;
    else
        cdID(nn) = 0;
    end
    
end

cdID = find(cdID);


cd_timing = ones(size(trig_timing,1),2);
cd_seq    = [];
trig_seq_exp_im_w_cd = timing.trig_seq_exp_im;
if ~isempty(cdID)
    for ii = cdID
        
        blockOnset = find(trig_block(:,1)==subj_run.block(ii).ID);
        blockOnset = blockOnset(1);
        foo = (trig_block(blockOnset:end,1)~=trig_block(blockOnset,1));
        foo = find(foo); blockOffset = blockOnset+foo(1)-1; clear foo
        
        stimOnset  = (trig_stim(blockOnset:blockOffset,1)>0 & trig_stim(blockOnset:blockOffset,1)<90);
        stimOnset  = blockOnset+find(stimOnset)-1;
        stimOnset_all  = [stimOnset(1); stimOnset((find(abs(diff(stimOnset))>1)+1))];
        time_table = subj_run.block(ii).timing;
        stim0 = time_table.onset_time(2:2:end-1);
        
        for ss = 1:length(stim0)
            trig_stim_time0 = trig_timing(stimOnset_all(ss));
            assert(nearZero((stim0(ss)+params.exp.trial.start_cue_dur+params.exp.trial.spatial_cue_dur) -trig_stim_time0))
            for nn = 1:size(subj_run.block(ii).trial,2)
                onset_decrement_s = trig_timing(stimOnset_all(ss)) + cdsoafun();
                onset_decrement_f = round(onset_decrement_s/params.stim.framedur_s);
                t_idx = onset_decrement_f:(onset_decrement_f+length(params.stim.cd.t_gauss)-1);
                cd_timing(t_idx,nn) = params.stim.cd.t_gauss';
                cd_seq = [cd_seq, onset_decrement_s];
                
                for tt = 1:length(t_idx)
                    clear tmp_im_c
                    if strcmp(subj_run.block(ii).name,'cd-ns') % range from [0 254]
                        tmp_im = (timing.trig_seq_exp_im{t_idx(tt)}{1});
                        sz0 = size(tmp_im);
                        if ndims(tmp_im)==3 && size(tmp_im,3)==4
                            has_mask = true;
                            tmp_im_mask = tmp_im(:,:,4);
                            tmp_im = tmp_im(:,:,1:3);
                        end
                        
                        tmp_im_g = rgb2gray(tmp_im);  % rgb to gray 
                        tmp_im_g_norm  = (double(tmp_im_g)./255).^2; % range [0-1];   
                        tmp_im_norm = (double(tmp_im)./255).^2; 
                        mn_g = mean(tmp_im_g_norm(:));
                        
                        % subtract the mean luminance of this scene
                        tmp_im_c = ((tmp_im_norm-mn_g).*params.stim.cd.t_gauss(tt)) + mn_g;
                        tmp_im_c = uint8(255.*sqrt(tmp_im_c)); % bring back to 0-255
                        
                        if has_mask
                            tmp_im_c = cat(3,tmp_im_c,tmp_im_mask);
                            has_mask = false;
                        end
                        
                        trig_seq_exp_im_w_cd{t_idx(tt)}{nn} = tmp_im_c;
                    elseif strcmp(subj_run.block(ii).name,'cd-cobj')
                        tmp_im = double(timing.trig_seq_exp_im{t_idx(tt)}{nn});
                        sz0 = size(tmp_im);

                        if ndims(tmp_im)==3 && size(tmp_im,3)==4
                            has_mask = true;
                            tmp_im_mask = tmp_im(:,:,4);
                            tmp_im = tmp_im(:,:,1:3);
                        end
                        
                        tmp_im_c = ((tmp_im-1)./254).^2; % range [0 1]
                        tmp_im_c = tmp_im_c-0.5; % center around 0, min/max range [-0.5 0.5]
                        tmp_im_c = tmp_im_c.*params.stim.cd.t_gauss(tt); % scale
                        tmp_im_c = uint8(254.*sqrt(tmp_im_c))+1; % bring back to 1-255
                        
                        if has_mask
                            tmp_im_c = cat(3,tmp_im_c,tmp_im_mask);
                            has_mask = false;
                        end
                        trig_seq_exp_im_w_cd{t_idx(tt)}{nn} = reshape(tmp_im_c,sz0);
                    else
                        tmp_im = double(timing.trig_seq_exp_im{t_idx(tt)}{nn});
                        sz0 = size(tmp_im);
                        if ndims(tmp_im)==3 && size(tmp_im,3)==4
                            has_mask = true;
                            tmp_im_mask = tmp_im(:,:,4);
                            tmp_im = tmp_im(:,:,1:3);
                        end
                        tmp_im_c = (tmp_im./255)-0.5; % center around 0, range [-0.5 0.5]
                        tmp_im_c = tmp_im_c.*params.stim.cd.t_gauss(tt); % scale
                        tmp_im_c = uint8( bsxfun(@plus, (255.*tmp_im_c), double(params.stim.bckgrnd_grayval))); % bring back to 1-255
                        
                        if has_mask
                            tmp_im_c = cat(3,tmp_im_c,tmp_im_mask);
                            has_mask = false;
                        end  
                        
                        trig_seq_exp_im_w_cd{t_idx(tt)}{nn} = reshape(tmp_im_c,sz0);
                    end
                end
            end
        end
    end
end

timing.seq_cd  = cd_seq;
timing.trig_cd = cd_timing;
timing.trig_seq_exp_im_w_cd = trig_seq_exp_im_w_cd;

% repmat(Expand(params.taskColor,1,2), length(params.seq)/4,1);
