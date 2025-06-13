function fix_matrix = vcd_createFixationSequence(params,fixsoafun,run_dur,blank_onset,blank_offset)
% VCD function to create a continuous sequence of fixation circle luminance
% changes for a given run duration, given a particular
% stimulus-onset-asynchrony (SOA).
%
% INPUTS:
%  * params         : (struct) parameter struct
%  * fixsoafun      : (function handle) function to determine the rate of 
%                       the fixation circle changing luminance during
%                       stimulus blocks.
%  * run_dur        : (double) duration of single run in time frames
%                       (expected to be the same across all runs in a session)
%  * blank_onset    : (double) first time frames related to onset of blank
%                        or "non-stimulus block" period. During these
%                        periods the fixation circle luminance will be
%                        frozen to mean luminance gray. This can be the 
%                        pre-run rest period, post-run rest period, or
%                        interblock interval (IBI).
%  * blank_offset   : (double) first time frames related to offset of blank
%                        or "non-stimulus block" period. 
%
% OUTPUTS:
%  * fix_matrix     : (double) Nx4 matrix where N is the number of time
%                       points. Column 1: time points (in frames), column
%                       2: absolute luminance values, 3: relative change in
%                       luminance values compared to the previous time
%                       point, 4: correct button press associated with the
%                       relative luminance change (1=brighter, 2=dimmer).
%
% Written by E Kupers @ UMN 2025/06

%% Determine frozen fixation periods

% Convert time frames to indices
blank_onset_idx  = blank_onset + 1;
blank_offset_idx = blank_offset + 1;

% Get frozen time periods
freeze_me = zeros(run_dur,1);
for ii = 1:length(blank_onset_idx)
    freeze_me(blank_onset_idx(ii):(blank_offset_idx(ii)-1)) = 1;
end

% Get twinkle luminance values (i.e., ignore 128 lum value as that's for IBI and blank periods only)
lum_vals = setdiff(params.stim.fix.dotlum,params.stim.bckgrnd_grayval);

% Assuming fixed fixation update interval
fix_interval = fixsoafun();

while 1
    % create sequence of shuffled luminance values (so randomly selecting WITHOUT replacement)
    lum_vec    = shuffle_concat(lum_vals,(ceil(run_dur/fix_interval)/length(lum_vals)));
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
for rr = 1:run_dur
    
    if freeze_me(rr) == 1 % we freeze luminance value at mid-gray level (128)
       fix_abs_lum(rr) = 128;
       counter = 0;
    elseif freeze_me(rr) == 0 % add luminance value
       if rr > 1 && (fix_abs_lum(rr-1) == 128)
           % if we just had a frozen time frame
           % then we want a new luminance sample
           lum_sample      = lum_vec(1);
           fix_abs_lum(rr) = lum_sample;
           lum_vec(1)      = []; % and toss it out 
           counter         = counter + 1; % update counter
       elseif counter==fix_interval
           % if reached our update rate, we move on the the luminance sample
           % and reset the counter
           lum_sample       = lum_vec(1);
           lum_vec(1)       = [];
           fix_abs_lum(rr)  = lum_sample;
           counter         = 1; % reset counter
       else % if we are in a stable state, we continue with the same 
           % luminance sample
           fix_abs_lum(rr) = lum_sample;
           counter         = counter + 1; % update counter
       end
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

%% Clean up button response
% the first transition from frozen to one of the 5 luminance levels is
% considered moot.

% Get frames corresponding to button presses.
button_press_frame  = find(~isnan(button_response_fix));

% All button press frames should correspond to a change in relative luminance.
assert(all(abs(fix_rel_lum(button_press_frame))>0)); 

% Find those corresponding to mean luminance gray.
button_press_mean_lum = button_press_frame(find(fix_abs_lum(button_press_frame)==params.stim.bckgrnd_grayval)); % indices corresponding to button presses for first transition to mean luminance gray 

% make sure the registered button press is linked to mean gray luminance
assert(isequal(unique(fix_abs_lum(button_press_mean_lum)),params.stim.bckgrnd_grayval)) 

% make sure this is indeed registered as a button press
assert(all(button_response_fix(button_press_mean_lum)~=0))        

% make sure this is indeed the first button press after a period of nothing
assert(all(isnan(button_response_fix(button_press_mean_lum-1)))); 

% check if there is a first button press related to the onset of that fixation block
blank_onset_mn_lum = blank_onset_idx(ismember(blank_onset_idx,button_press_frame));
blank_offset_mn_lum = blank_offset_idx(ismember(blank_offset_idx,button_press_frame));

delete_me_onset  = fix_abs_lum(blank_onset_mn_lum)==128;
delete_me_offset = fix_abs_lum(blank_offset_mn_lum)~=128;

if sum(delete_me_onset)>0
    button_response_fix(blank_onset_mn_lum(delete_me_onset)) = NaN;
end

if sum(delete_me_offset)>0
    button_response_fix(blank_offset_mn_lum(delete_me_onset)) = NaN;
end

% Record timing, absolute, relative fixation luminance, and button response
fix_matrix = [fix_timing,fix_abs_lum,fix_rel_lum,button_response_fix];
fix_matrix = fix_matrix(1:run_dur,:); % ensure duration is as we expect

end