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

% fill up fixation sequence until run dur
if f_fix(end)~=run_dur
    f_fix(end+1) = run_dur;
    fix_seq(end+1) = f_fix(end) - f_fix(end-1);
end

f_fix = [0, f_fix];
blank_onset(1) = 0;

% Ignore 128 lum value (that's for IBI and black periods only)
lum_vals = setdiff(params.stim.fix.dotlum,128);

% Update luminance of fixation dot randomly (WITHOUT replacement)
lum_shuffled_vec = [];
for kk = 1:(1+ceil(length(f_fix)/length(lum_vals)))
    lum_sample = datasample(double(lum_vals),length(lum_vals),'Replace',false);
    if kk > 1 && (lum_sample(1) == lum_shuffled_vec(end)) % if we happen to sample the same luminance for the first of this series and the last of previous series, we swap first and second lum values
        lum_sample([1,2]) = lum_sample([2,1]);
    end
    lum_shuffled_vec = cat(2,lum_shuffled_vec, lum_sample);
end
lum_shuffled_vec = lum_shuffled_vec(1:length(f_fix));

assert(isequal(sort(unique(lum_shuffled_vec)),lum_vals)); % no mid gray level!
assert(all(diff(lum_shuffled_vec)~=0)) % no repeats!

blank_offset(end) =  f_fix(end);

% Expand luminance values to time frames
fix_abs_lum = [];
for ff = 1:length(lum_shuffled_vec)
    curr_time_frames = f_fix(ff):(f_fix(ff)+fix_seq(ff)-1);
    clear ia0 ib0
    ia = []; ib = [];
    for ii = 1:length(blank_onset)
        freeze_me = blank_onset(ii):blank_offset(ii);
        [~,ia0] = intersect(curr_time_frames,freeze_me);
        [~,ib0] = setdiff(curr_time_frames,freeze_me);
        ia = cat(1,ia,ia0);
        ib = cat(1,ib,ib0);
    end
    
    curr_fix = [];
    if length(fix_abs_lum)+fix_seq(ff) >= run_dur
       fix_abs_lum = cat(1,fix_abs_lum,128*ones(run_dur-length(fix_abs_lum),1));
    else
        if ~isempty(ia)
            % freeze fixation lum during IBI/pre/post blank periods
            curr_fix(ia) = 128*ones(length(ia),1); % mean luminance gray

            if ia(1) == 1 && (ia(end) < length(curr_time_frames)) % if we start with frozen fix, but not finish with it, we continue counting the same luminance
                curr_fix(ib(1):(ib(1)+length(curr_time_frames)-1)) = Expand(lum_shuffled_vec(ff), 1, length(curr_time_frames));

                % and push the onset of all the other luminance changes forward
                f_fix((ff+1):end) = f_fix((ff+1):end)+ib(1);

                future_time_frames =  f_fix(ff+1):(f_fix(ff+1)+fix_seq(ff+1)-1);
                [~,ic] = intersect(future_time_frames,freeze_me);
                assert(isempty(ic));

            elseif ia(1) > 1 && (ia(end) == length(curr_time_frames)) % if we finish on frozen fix
                curr_fix(ib) = Expand(lum_shuffled_vec(ff), 1, length(ib));
            end
        else
            curr_fix = Expand(lum_shuffled_vec(ff), 1, length(curr_time_frames));
        end
        if size(curr_fix,1) < size(curr_fix,2)
            curr_fix = curr_fix';
        end
        fix_abs_lum = cat(1, fix_abs_lum, curr_fix);
        if length(fix_abs_lum) > run_dur
            fix_abs_lum = fix_abs_lum(1:run_dur); % trim the end
            break;
        end
    end
end

% figure out when dot is brighter, dimmer or the same
% relative to previous time point
dimmer   = find(diff([0, fix_abs_lum'])<0);
nochange = find(diff([0, fix_abs_lum'])==0);
brighter = find(diff([0, fix_abs_lum'])>0);

[r_sort,r_ai] = sort([dimmer, nochange,brighter]);
response_fix_vex = [-1.*ones(size(dimmer)), zeros(size(nochange)),1.*ones(size(brighter))];
response_fix_vex = response_fix_vex(r_ai);
response_fix_vex(1) = 0; % set first lum to no response;

fix_timing  = [0:(run_dur-1)]';
fix_rel_lum = [0; diff(fix_abs_lum)];

% translate changes in luminance into button responses
button_response_fix = response_fix_vex';
button_response_fix(response_fix_vex==1)  = 1; % brighter
button_response_fix(response_fix_vex==-1) = 2; % dimmer

fix_matrix = [fix_timing,fix_abs_lum,fix_rel_lum,button_response_fix];
fix_matrix = fix_matrix(1:run_dur,:);

end