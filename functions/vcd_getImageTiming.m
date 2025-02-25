function timing = vcd_getImageTiming(params, subj_run) 

% Fixation order and fixation
fixsoafun = @() round(params.stim.fix.dotmeanchange + (params.stim.fix.dotchangeplusminus*(2*(rand-.5))))*params.stim.fps;

% Contrast decrement gaussian time window onset
cdsoafun = @() round(params.stim.cd.meanchange + params.stim.cd.changeplusminus*(2*(rand-.5)))*params.stim.fps;


%% TIMING

seq_stim = {}; 
seq_timing = []; 
spatial_cue = []; 

% 6 fields (name, ID, within_session_repeat, trial, trial_type, timing)
% 8 blocks: 1:run, 2:block, 3:stimtaskID, 4:unique_im, 5:spatial_cue, 6:onset_time, 7:event_dur, 8:run_time
cellblock = squeeze(struct2cell(subj_run.block));

for bb = 1:size(cellblock,2)

    tmp_timing = cellblock{6,bb};
    im_nr = tmp_timing.unique_im;
    
    if ~iscell(im_nr)
        % convert to cell
        im_nr = num2cell(im_nr);
    end

    % check for nans
    for xi = 1:length(im_nr)
        tmp_im = im_nr{xi};
        if isnan(tmp_im)
            tmp_im(isnan(tmp_im))=0;
            im_nr{xi} = tmp_im;
        end
        clear tmp_im
    end  
    seq_stim = cat(1, seq_stim, im_nr);

    seq_timing = cat(1,seq_timing, tmp_timing.run_time);
    spatial_cue = cat(1,spatial_cue, tmp_timing.spatial_cue);
end

timing.seq_stim        = seq_stim;
timing.seq_spatial_cue = spatial_cue;
timing.seq_timing      = seq_timing;


%% Convert sequence of events from seconds into frames

trig_timing = [0:params.stim.fps:seq_timing(end)]'; % seconds

trig_stim = zeros(size(trig_timing,1),2);
for tt = 1:length(seq_timing)
    event_time = seq_timing(tt);
    
    [t_diff,t_idx] = min(abs(trig_timing-event_time));
    if isempty(t_idx) || ~nearZero(t_diff)
        error('[%s]: Event_time doesn''t match monitor refresh rate', mfilename)
    end
    
    if length(seq_stim{tt})==2
        trig_stim(t_idx,:) = seq_stim{tt};
    else
        trig_stim(t_idx,:) = repmat(seq_stim{tt},1,2);
    end

end

timing.trig_stim    = trig_stim;
timing.trig_timing  = trig_timing;

%% FIXATION DOT SEQUENCE
fudge   = 5; % shave off 2 second at the end to avoid going over run dur

fix_tmp = [0]; fix_timing = 0; fix_seq = [randi(length(params.stim.fix.dotlum),1)];
while fix_timing < (trig_timing(end)-fudge)
    fix_tmp = [fix_tmp, fixsoafun()];
    fix_timing = cumsum(fix_tmp);
    fix_seq = [fix_seq, randi(length(params.stim.fix.dotlum),1)];
end

trig_fix = zeros(size(trig_timing,1),1);
for tt = 1:length(fix_timing)
    event_time = fix_timing(tt);
    
    [t_diff,t_idx] = min(abs(trig_timing-event_time));
    if isempty(t_idx) || ~nearZero(t_diff)
        error('[%s]: Event_time doesn''t match monitor refresh rate', mfilename)
    end
    
    trig_fix(t_idx,:) = fix_seq(tt);
end


timing.seq_fix      = fix_seq;
timing.trig_fix     = trig_fix;

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
if ~isempty(cdID)
    for ii = cdID
        
        blockOnset = find(trig_stim(:,1)==subj_run.block(ii).ID);
        time_table = subj_run.block(ii).timing;
        stim0 = time_table.onset_time(2:2:end);
        for ss = 1:length(stim0)
            for nn = 1:size(subj_run.block(ii).trial,2)
                onset_decrement_s = stim0(ss) + cdsoafun();
                onset_decrement_f = round(onset_decrement_s/params.stim.fps);
                cd_timing(onset_decrement_f:(onset_decrement_f+length(params.stim.cd.t_gauss)-1),nn) = params.stim.cd.t_gauss';
                cd_seq = [cd_seq, onset_decrement_s];
            end
        end
    end
end

timing.seq_cd  = cd_seq;
timing.trig_cd = cd_timing;



% repmat(Expand(params.taskColor,1,2), length(params.seq)/4,1);
