function fix_matrix = vcd_createFixationSequence(params,fixsoafun,run_dur,blank_onset,blank_offset)
% VCD function to create a continuous sequence of fixation circle luminance
% changes for a given run duration, given a particular
% stimulus-onset-asynchrony (SOA).
%
% INPUTS:
%  * WRITE ME
%
% OUTPUTS:
%  * fix_matrix     : (double) Nx4 matrix where N is the number of time
%                       points. Column 1: time points (in frames), column
%                       2: absolute luminance values, 3: relative change in
%                       luminance values compared to the previous time
%                       point, 4: correct button press associated with the
%                       relative luminance change (1=brighter, 2=dimmer).

% %%%% Get fixation timing
fix_seq = []; f_fix = 0;
while f_fix(end) < run_dur
    fix_seq = [fix_seq, fixsoafun()];
    f_fix = cumsum(fix_seq);
end

% trim in case we accidentally went overtime
f_fix = f_fix(f_fix<run_dur);

% freeze fixation lum during black periods
delete_me = [];
for ii = 1:length(blank_onset)
    delete_me = cat(2,delete_me, find(f_fix > blank_onset(ii) & f_fix < blank_offset(ii)));
end

f_fix(delete_me) = [];

% Update luminance of fixation dot randomly (WITHOUT replacement)
lum_shuffled_idx = [];
for kk = 1:(1+ceil(length(f_fix)/length(params.stim.fix.dotlum)))
    lum_shuffled_idx = cat(2,lum_shuffled_idx, datasample(double(params.stim.fix.dotlum),length(params.stim.fix.dotlum),'Replace',false));
end

lum_shuffled_idx = lum_shuffled_idx(1:(length(f_fix)+1));

% figure out when dot is brighter, dimmer or the same
% relative to previous time point
dimmer   = find(diff([0, lum_shuffled_idx])<0);
nochange = find(diff([0, lum_shuffled_idx])==0);
brighter = find(diff([0, lum_shuffled_idx])>0);

[r_sort,r_ai] = sort([dimmer, nochange,brighter]);
response_fix_vex = [-1.*ones(size(dimmer)), zeros(size(nochange)),1.*ones(size(brighter))];
response_fix_vex = response_fix_vex(r_ai);
response_fix_vex(1) = 0; % set first lum to no response;

% Add pre and post fixation sequence duration
dur_fix_frames = [f_fix(1), diff(f_fix), (run_dur-f_fix(end))];

fix_timing = []; fix_abs_lum = []; fix_rel_lum = [];
for ff = 1:length(dur_fix_frames)
    fix_timing = cat(1, fix_timing, Expand(ff, 1, dur_fix_frames(ff)));
    fix_abs_lum = cat(1, fix_abs_lum, Expand(lum_shuffled_idx(ff), 1, dur_fix_frames(ff)));
    fix_rel_lum = cat(1,fix_rel_lum, [response_fix_vex(ff); zeros(dur_fix_frames(ff)-1,1)]);
end
% translate changes in luminance into button responses
button_response_fix = fix_rel_lum;
button_response_fix(button_response_fix==1)  = 1; % brighter
button_response_fix(button_response_fix==-1) = 2; % dimmer

fix_matrix = [fix_timing,fix_abs_lum,fix_rel_lum,button_response_fix];

end