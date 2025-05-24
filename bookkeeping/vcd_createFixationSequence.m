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
% fix_seq = [0]; f_fix = 0;
% while f_fix(end) < run_dur
%     fix_seq = [fix_seq, fixsoafun()];
%     f_fix = cumsum(fix_seq);
% end

% % trim in case we accidentally went overtime
% f_fix = f_fix(f_fix<run_dur);
% 
% % fill up fixation sequence until run dur
% if f_fix(end)~=run_dur
%     f_fix(end+1) = run_dur;
%     fix_seq(end+1) = f_fix(end) - f_fix(end-1);
% end


% Determine frozen periods
blank_onset(1)    = 1; % add 1 because we can't zero index 
blank_offset(end) = blank_offset(end)+1;

freeze_me = zeros(run_dur,1);
for ii = 1:length(blank_onset)
    freeze_me(blank_onset(ii):(blank_offset(ii)-1)) = 1;
end

freeze_me(end) = 1;

% Ignore 128 lum value (that's for IBI and black periods only)
lum_vals = setdiff(params.stim.fix.dotlum,128);

fix_interval = fixsoafun()-1;

while 1
    % create sequence of shuffled luminance values (so randomly selecting WITHOUT replacement)
    lum_vec = shuffle_concat(lum_vals,(ceil(run_dur/fix_interval)/length(lum_vals)));
    repeat_lum = find(diff(lum_vec)==0);
    lum0 = lum_vec(repeat_lum);
    
    % if we happen to sample the same luminance in a row
    for ii = 1:length(lum0)
        % we check neighoring values
        lum_next0  = lum_vec(repeat_lum(ii)-1);
        lum_next2  = lum_vec(repeat_lum(ii)+2);
        lum_next3  = lum_vec(repeat_lum(ii)+3);
        % and swap them if we can
        if lum0(ii) ~= lum_next0
            tmp1 = lum0(ii);
            tmp2 = lum_next0;
            lum_vec(repeat_lum(ii)) = tmp2;
            lum_vec(repeat_lum(ii)-1) = tmp1;
        elseif lum0(ii) ~= lum_next2
            tmp1 = lum0(ii);
            tmp2 = lum_next2;
            lum_vec(repeat_lum(ii)) = tmp2;
            lum_vec(repeat_lum(ii)+2) = tmp1;
        elseif lum0(ii) ~= lum_next3
            tmp1 = lum0(ii);
            tmp2 = lum_next3;
            lum_vec(repeat_lum(ii)) = tmp2;
            lum_vec(repeat_lum(ii)+3) = tmp1;
        end
    end
    

    % if we are good, we break out of the loop
    if sum(diff(lum_vec)==0)==0
        break;
    end
    % otherwise we start over
end

fix_abs_lum = NaN(run_dur,1);
for kk = 1:run_dur
    
    if freeze_me(kk) == 1 % we freeze luminance value at mid-gray level (128)
       fix_abs_lum(kk) = 128;
       counter = 0;
    elseif freeze_me(kk) == 0 % add luminance value
       if kk > 1 && (fix_abs_lum(kk-1) == 128)
           % if we just had a frozen time frame
           % then we want a new luminance sample
           lum_sample      = lum_vec(1);
           fix_abs_lum(kk) = lum_sample;
           lum_vec(1)      = []; % and toss it out 
       elseif counter==fix_interval
           % if reached our update rate, we move on the the luminance sample
           % and reset the counter
           counter          = 0;
           lum_sample       = lum_vec(1);
           lum_vec(1)       = [];
           fix_abs_lum(kk)  = lum_sample;
       else % if we are in a stable state, we continue with the same 
           % luminance sample
           fix_abs_lum(kk)  = lum_sample;
       end
       counter = counter + 1;
    end
end

% figure out when dot is brighter, dimmer or the same
% relative to previous time point
fix_timing  = [0:(run_dur-1)]';
fix_rel_lum = [0; diff(fix_abs_lum)];

% translate changes in luminance into button responses
button_response_fix = NaN(run_dur,1);
button_response_fix(fix_rel_lum>0) = 1; % brighter
button_response_fix(fix_rel_lum<0) = 2; % dimmer

fix_matrix = [fix_timing,fix_abs_lum,fix_rel_lum,button_response_fix];
fix_matrix = fix_matrix(1:run_dur,:);

end