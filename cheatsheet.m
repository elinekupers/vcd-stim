%% Official vcdpp01 experiment (SUBJECTID is an integer up to 3 digits, RUNNUMBER is 1-12)
all_images = struct;
[data,all_images] = runme_vcdcore(SUBJECTID, 1, 1, RUNNUMBER, 'PPROOM_EIZOFLEXSCAN', ...
  'wanteyetracking', true, 'all_images', all_images, 'exp_env', 2);

%% Official vcdwide01 experiment (SUBJECTID is an integer up to 3 digits, RUNNUMBER is 1-10)
% SESS_TYPE is 1 or 2 where 1=A, 2=B
all_images = struct;
[data,all_images] = runme_vcdcore(SUBJECTID, 1, SESSTYPE, RUNNUMBER, '7TAS_BOLDSCREEN32', ...
  'is_wide',true,'wanteyetracking', true, 'all_images', all_images, 'exp_env', 1);

%% Official vcddeep08-45 experiment (SUBJECTID is an integer up to 3 digits, RUNNUMBER is 1-10)
% For SESNR 8-45 --> SESSTYPE MUST BE 1, where 1=A
% Flag 'is_wide' = false by default.
% 
% To prepare a time table file for a subject (you only need to do this once, 
% before the first run of deep08):
% Step 1: Move all previous subject mat files to a new folder with a different name in ../data/MRI/
% Step 2: Run the following to generate new subject time table:
%   all_images = struct;
%   [data,all_images] = runme_vcdcore(SUBJECTID, 8, 1, 1, '7TAS_BOLDSCREEN32', ...
%     'wanteyetracking', true, 'all_images', all_images, 'exp_env', 1);
%
% Step 3: Run the other runs and sessions where you define the time table:

SUBJECTID = 999; % is an integer up to 3 digits
SESNR     = 8; % integer between 8 through 46 (inclusive)
SESSTYPE  = 1; % SES_NR 8-45 -> ALWAYS use 1=A

data_dir = fullfile(vcd_rootPath, 'data','MRI',sprintf('vcd_subj%03d',SUBJECTID));
tt_file  = dir(fullfile(data_dir, sprintf('vcd_subj%03d_time_table_master_deep_7TAS_BOLDSCREEN32_2026*.mat',SUBJECTID))); % we will always use the latest time table file

all_images = struct;
for RUNNUMBER = 1:10
    [data,all_images] = runme_vcdcore(SUBJECTID, SESNR, SESSTYPE, RUNNUMBER, '7TAS_BOLDSCREEN32', ...
    'wanteyetracking', true, 'all_images', all_images, 'exp_env', 1, ...
    'timetable_file',fullfile(tt_file(end).folder, tt_file(end).name));
end

%% For Official vcddeep46 --> SESSTYPE is 1 or 2 where 1=A, 2=B

SUBJECTID = 999; % is an integer up to 3 digits
SESNR     = 46; % integer between 8 through 46 (inclusive)
SESSTYPE  = [];  % use 1=A or 2=B

data_dir = fullfile(vcd_rootPath, 'data','MRI',sprintf('vcd_subj%03d',SUBJECTID));
tt_file  = dir(fullfile(data_dir, sprintf('vcd_subj%03d_time_table_master_deep_7TAS_BOLDSCREEN32_2026*.mat',SUBJECTID))); % we will always use the latest time table file

all_images = struct;
for RUNNUMBER = 1:10
    [data,all_images] = runme_vcdcore(SUBJECTID, SESNR, SESSTYPE, RUNNUMBER, '7TAS_BOLDSCREEN32', ...
    'wanteyetracking', true, 'all_images', all_images, 'exp_env', 1, ...
    'timetable_file',fullfile(tt_file(end).folder, tt_file(end).name));
end

                
%% A run for PP room

all_images = struct;
[data,all_images] = runme_vcdcore(4, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', false, ...
'all_images',all_images, 'exp_env',2);

%% Run demo
% Demo Session 1: FIX-GBR, WM-RDK, CD-DOT, HOW-OBJ, WHAT-NS, SCC-ALL (2x)
% Demo Session 2: WM-GBR, PC-RDK, FIX-DOT, PC-DOT, PC-OBJ, CD-NS, WHERE-NS
% Demo session 3: PC-GBR, CD-RDK, WM-DOT, WHAT-OBJ, FIX-NS, WM-NS, HOW-NS                                           
demo_session = 1; % can be 1, 2, or 3.
all_images = struct;
[data,all_images] = runme_vcdcore(999, demo_session, 1, 1, 'PPROOM_EIZOFLEXSCAN', ...
    'wanteyetracking', false, 'all_images',all_images, ...
    'exp_env',2, 'is_demo', true, 'is_wide', true);

%% Run DEEP demo
% Demo Session 1-4 A = LTM-only REGULAR TEST runs (10 runs per session).
% Demo Session 1-4 B = LTM-only STUDY trials (10 runs per session).
% Demo Session 5-8 A = IMG-only REGULAR TEST runs (10 runs per session). 
% Demo Session 5-8 B = IMG-only PERCEPTION runs (10 runs per session).
demo_session = 1; % can be 1 through 8
demo_session_type = 1; % can be 1 (A - regular) or 2 (B - LTM study/ IMG perception)
all_images = struct;
[data,all_images] = runme_vcdcore(999, demo_session, 1, 1, 'PPROOM_EIZOFLEXSCAN', ...
    'wanteyetracking', false, 'all_images',all_images, ...
    'exp_env',2, 'is_demo', true, 'is_wide', false);

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
