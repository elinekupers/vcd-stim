%% a run

all_images = struct;
[data,all_images] = runme_vcdcore(4, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', false, ...
'all_images',all_images, 'exp_env',2);

%% Run demo

demo_session = 1; % can be 1, 2, or 3.
all_images = struct;
[data,all_images] = runme_vcdcore(999, demo_session, 1, 1, 'PPROOM_EIZOFLEXSCAN', ...
    'wanteyetracking', false, 'all_images',all_images, ...
    'exp_env',2, 'is_demo', true);

%% rerun behavioral performance analysis;
subj_folder = 'vcd_subj999_ses01';
mat_file    = 'behavior_20250623124130_vcd_subj999_ses01_A_run01.mat';
load(fullfile('data','BEHAVIOR',subj_folder,mat_file));
performance = vcdbehavioralanalysis(fullfile(params.savedatafolder,params.behaviorfile));
%% KK's automatic data generation stuff [FOR DEVELOPMENT ONLY]
all_images = struct;
for p=1:12
  [data,all_images] = runme_vcdcore(1, 1, 1, p, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', false, ...
      'all_images',all_images,'wantdatabypass',true,'exp_env',2, ...
      'timetable_file','/Users/psphuser/Desktop/cvnlab/VCD/vcd-stim/data/BEHAVIOR/vcd_subj001/vcd_subj001_time_table_master_PPROOM_EIZOFLEXSCAN_20250624T090728.mat');
end

%% making videos
all_images = struct;
[data,all_images] = runme_vcdcore(3, 1, 1, 3, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', false, ...
      'wantsynctest', false, 'all_images',all_images, 'exp_env',2);  
  
%% Rsync stimuli
%#history | grep rsync
%#rsync -avrh kupers@surly.cmrr.umn.edu:"/var/fairstate-ext2/GoogleDrive/VCD/experimental_design/stimuli/final_stimuli/PPROOM_EIZOFLEXSCAN/mat_files/*.mat" workspaces/stimuli/PPROOM_EIZOFLEXSCAN/

%% Rsync data

%#rsync -avrh /Users/psphuser/Desktop/cvnlab/VCD/vcd-stim/data/BEHAVIOR/vcd_subj003/vcd_subj003_time_table_master_PPROOM_EIZOFLEXSCAN_20250624T134224.mat poole163@surly.cmrr.umn.edu:"/var/fairstate-ext2/GoogleDrive/VCD/vcdalpha/20250624_adrian/"
