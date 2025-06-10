% tinker_optimalBlockAllocation.m

pres_rate = 60; % hz
s_range = 0:7;
d_range = 0:7;

single_block_dur = 2670; % time frames // 44.5 s (includes task_cue_dur)
double_block_dur = 3630; % time frames // 60.5 s (includes task_cue_dur)
frames2min = 3600;
min_IBI = 300; % time frames // 5 sec
max_IBI = 540; % time frames // 9 sec
mn_IBI  = 420; % time frames // 7 sec

task_cue_dur = 5*pres_rate; % time frames // 5 sec
lower_limit  = 5*60*pres_rate;  % time frames // 4 min
upper_limit  = 388.8*pres_rate; % time frames // ~6.5 min (6 min & 29 s)

et_block_dur = 18.5*pres_rate; % time frames // 18.5 sec
prepost_dur  = (12+4)*pres_rate; % time frames // 16 sec
overhead_dur = et_block_dur + prepost_dur;

clear min_block_durs max_block_durs mean_block_durs min_rundur max_rundur mean_rundur
for ii = 1:length(s_range)
    
    for jj = 1:length(d_range)
        
        nr_blocks = s_range(ii)+d_range(jj);
        min_block_durs(ii,jj)   = (single_block_dur*s_range(ii) + double_block_dur*d_range(jj)) + (choose((nr_blocks-1)<0, 0, nr_blocks-1)*min_IBI);
        max_block_durs(ii,jj)   = (single_block_dur*s_range(ii) + double_block_dur*d_range(jj)) + (choose((nr_blocks-1)<0, 0, nr_blocks-1)*max_IBI);
        mean_block_durs(ii,jj)  = (single_block_dur*s_range(ii) + double_block_dur*d_range(jj)) + (choose((nr_blocks-1)<0, 0, nr_blocks-1)*mn_IBI);
        
        min_rundur(ii,jj) = min_block_durs(ii,jj) + overhead_dur; 
        max_rundur(ii,jj) = max_block_durs(ii,jj) + overhead_dur; 
        mean_rundur(ii,jj) = mean_block_durs(ii,jj) + overhead_dur; 
    end
end 

% find potential block combinations
[potential_blocks_min_d,potential_blocks_min_s] = find( (min_rundur >= lower_limit) & (min_rundur <= upper_limit));
[potential_blocks_max_d,potential_blocks_max_s] = find( (max_rundur >= lower_limit) & (max_rundur <= upper_limit));
[potential_blocks_mean_d,potential_blocks_mean_s] = find((mean_rundur >= lower_limit) & (mean_rundur <= upper_limit));

all_potential_blocks = unique([find( (min_rundur >= lower_limit) & (min_rundur <= upper_limit));find( (max_rundur >= lower_limit) & (max_rundur <= upper_limit))]);
[all_potential_s_blocks,all_potential_d_blocks] = ind2sub([length(s_range), length(d_range)],all_potential_blocks);

figure(1); clf; 
subplot(131); imagesc([min_rundur/frames2min])
hold on; plot(potential_blocks_min_s,potential_blocks_min_d,'bo','linewidth',3)
hold on; plot(potential_blocks_mean_s,potential_blocks_mean_d,'rx','linewidth',5)
title('Min run duration (min)')
ylabel('# single blocks')
xlabel('# double blocks')
colorbar
axis square image
set(gca,'XTick',[1:8],'XTickLabel',s_range, 'YTick',[1:8],'YTickLabel', d_range, 'CLim', [0 15])
subplot(132); imagesc(max_rundur/frames2min)
hold on; plot(potential_blocks_max_s,potential_blocks_max_d,'bo','linewidth',3)
hold on; plot(potential_blocks_mean_s,potential_blocks_mean_d,'rx','linewidth',5)
title('Max run duration (min)')
ylabel('# single blocks')
xlabel('# double blocks')
colorbar
set(gca,'XTick',[1:8],'XTickLabel',s_range, 'YTick',[1:8],'YTickLabel', d_range, 'CLim', [0 15])
axis square image
subplot(133); imagesc(mean_rundur/frames2min)
hold on; plot(potential_blocks_mean_s,potential_blocks_mean_d,'rx','linewidth',5)
title('Mean run duration (min)')
ylabel('# single blocks')
xlabel('# double blocks')
colorbar
set(gca,'XTick',[1:8],'XTickLabel',s_range, 'YTick',[1:8],'YTickLabel', d_range, 'CLim', [0 15])
axis square image

% total_run_dur_frames = repmat((et_block_dur+prepost_dur),length(double_b_sub),1)+sum(block_durs([single_b_sub,double_b_sub]),2);
% total_run_dur_sec    = total_run_dur_frames/pres_rate;

tmp = [s_range(all_potential_s_blocks)',d_range(all_potential_d_blocks)', ...
    min_rundur(all_potential_blocks)/frames2min, ...
   max_rundur(all_potential_blocks)/frames2min, ...
   mean_rundur(all_potential_blocks)/frames2min, ...
        (upper_limit - max_rundur(all_potential_blocks))/60];

disp(tmp)

%% Run dur = 6 min & 28.8 s
%% single    double   min min   max min   mean min  unaccount time (s) from max run time
% 6.0000         0    5.4417    5.7750    5.6083   42.1000
% 7.0000         0    6.2667    6.6667    6.4667  -11.4000
% 4.0000    1.0000    4.8833    5.1500    5.0167   79.6000
% 5.0000    1.0000    5.7083    6.0417    5.8750   26.1000
% 3.0000    2.0000    5.1500    5.4167    5.2833   63.6000
% 4.0000    2.0000    5.9750    6.3083    6.1417   10.1000
% 2.0000    3.0000    5.4167    5.6833    5.5500   47.6000
% 3.0000    3.0000    6.2417    6.5750    6.4083   -5.9000
%      0    4.0000    4.8583    5.0583    4.9583   85.1000
% 1.0000    4.0000    5.6833    5.9500    5.8167   31.6000
%      0    5.0000    5.9500    6.2167    6.0833   15.6000


 
 

 
  
 


  



  
  
  
  
  
  
  