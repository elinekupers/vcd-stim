%% Just a test run
all_images = struct;
[data,all_images] = runme_vcdcore(1, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', true, ...
      'all_images',all_images);


%% KK's automatic data generation stuff [FOR DEVELOPMENT ONLY]
all_images = struct;
for p=1:12
  [data,all_images] = runme_vcdcore(1, 1, 1, p, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', false, ...
      'all_images',all_images,'wantdatabypass',true,'exp_env',2, ...
      'timetable_file','/Users/psphuser/Desktop/cvnlab/VCD/vcd-stim/data/BEHAVIOR/vcd_subj001/vcd_subj001_time_table_master_PPROOM_EIZOFLEXSCAN_20250618T081915.mat');
end
