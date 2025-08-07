% EK testbed

%% DRY-RUN 3 SESSIONS OF DEMO (BEHAVIORAL) EXPERIMENT
subj_nr = 998;
data_dir = fullfile(vcd_rootPath, 'data','BEHAVIOR',sprintf('vcd_subj%03d',subj_nr));

all_images = struct;
[data,all_images] = runme_vcdcore(subj_nr, 1, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', false, ...
'all_images',all_images,'wantdatabypass',true,'exp_env',4, 'is_demo',true); 

tt_file    = dir(fullfile(data_dir, sprintf('vcd_subj%03d_time_table_master_demo_PPROOM_EIZOFLEXSCAN_*.mat',subj_nr)));
all_images = struct;
for p = 1:3
[data,all_images] = runme_vcdcore(subj_nr, p, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', false, ...
'all_images',all_images,'wantdatabypass',true,'exp_env',4, 'is_demo',true,'timetable_file',fullfile(tt_file(end).folder, tt_file(end).name)); 
end


load(fullfile(tt_file(end).folder, tt_file(end).name));
writetable(time_table_master,'~/Desktop/time_table_master_demo.csv');
writetable(all_run_frames,'~/Desktop/all_run_frames_demo.csv');


for p = 1:3
    subj_folder = sprintf('vcd_subj%03d_ses%02d',subj_nr,p);
    mat_file    = sprintf('behavior_*_vcd_subj%03d_ses%02d_A_run01.mat',subj_nr,p);
    dd = dir(fullfile(vcd_rootPath,'data','BEHAVIOR',subj_folder,mat_file));
    load(fullfile(dd(end).folder,dd(end).name));
    performance = vcdbehavioralanalysis(fullfile(params.savedatafolder,params.behaviorfile));
end



%% DRY-RUN 12 RUNS OF BEHAVIORAL EXPERIMENT
subj_nr = 999;
ses_nr = 1;
data_dir = fullfile(vcd_rootPath, 'data','BEHAVIOR',sprintf('vcd_subj%03d',subj_nr));

all_images = struct;
[data,all_images] = runme_vcdcore(subj_nr, ses_nr, 1, 1, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', false, ...
'all_images',all_images,'wantdatabypass',false,'wantptbbypass',true,'exp_env',4); 

tt_file    = dir(fullfile(data_dir, sprintf('vcd_subj%03d_time_table_master_PPROOM_EIZOFLEXSCAN_*.mat',subj_nr)));
load(fullfile(tt_file(end).folder, tt_file(end).name));
writetable(time_table_master,'~/Desktop/time_table_master.csv');
writetable(all_run_frames,'~/Desktop/all_run_frames.csv');

all_images = struct;
for rr = 1:12
    [data,all_images] = runme_vcdcore(subj_nr, ses_nr, 1, rr, 'PPROOM_EIZOFLEXSCAN', 'wanteyetracking', false, ...
        'all_images',all_images,'wantdatabypass',true,'exp_env',4, 'timetable_file',fullfile(tt_file(end).folder, tt_file(end).name));
end

clear behresults
for rr = 1:12
    subj_folder = sprintf('vcd_subj%03d_ses%02d',subj_nr,ses_nr);
    mat_file    = sprintf('behavior_*_vcd_subj%03d_ses%02d_A_run%02d.mat',subj_nr,ses_nr,rr);
    dd = dir(fullfile(vcd_rootPath,'data','BEHAVIOR',subj_folder,mat_file));
    load(fullfile(dd(end).folder,dd(end).name));
    performance = vcdbehavioralanalysis(fullfile(params.savedatafolder,params.behaviorfile));

    writetable(run_frames,sprintf('~/Desktop/run_frames%02d.csv',rr));
    writetable(run_table,sprintf('~/Desktop/run_table%02d.csv',rr));
    behresults(rr) = performance;
end

vcd_checkTimeTable(time_table_master,behresults)

%% STANDALONE BEHAVIORAL RESULTS CHECK
subj_nr = 999;
ses_nr  = 1;

clear behresults
for rr = 1:12
    subj_folder = sprintf('vcd_subj%03d_ses%02d',subj_nr,ses_nr);
    mat_file    = sprintf('behavior_*_vcd_subj%03d_ses%02d_A_run%02d.mat',subj_nr,ses_nr,rr);
    dd = dir(fullfile(vcd_rootPath,'data','BEHAVIOR',subj_folder,mat_file));
    load(fullfile(dd(end).folder,dd(end).name));
    performance = vcdbehavioralanalysis(fullfile(dd(end).folder,dd(end).name));

    writetable(run_frames,sprintf('~/Desktop/run_frames%02d.csv',rr));
    writetable(run_table,sprintf('~/Desktop/run_table%02d.csv',rr));
    behresults(rr) = performance;
end

tt_file    = dir(fullfile(vcd_rootPath,'data','BEHAVIOR',sprintf('vcd_subj%03d',subj_nr), ...
    sprintf('vcd_subj%03d_time_table_master_PPROOM_EIZOFLEXSCAN_*.mat',subj_nr)));
load(fullfile(tt_file(end).folder, tt_file(end).name));
vcd_checkTimeTable(time_table_master,behresults)


%% DRY-RUN 10 RUNS OF WIDE MRI EXPERIMENT
subj_nr = 000;
ses_nr = 1;
ses_type = [1,2];

data_dir = fullfile(vcd_rootPath, 'data','MRI',sprintf('vcd_subj%03d',subj_nr));

for st = ses_type
    all_images = struct;
    [data,all_images] = runme_vcdcore(subj_nr, ses_nr, st, 1, '7TAS_BOLDSCREEN32', 'wanteyetracking', false, ...
        'all_images',all_images,'wantdatabypass',true,'exp_env',4);
    
    tt_file    = dir(fullfile(data_dir, sprintf('vcd_subj%03d_time_table_master_wide_7TAS_BOLDSCREEN32_*.mat',subj_nr)));
    load(fullfile(tt_file(end).folder, tt_file(end).name));
    writetable(time_table_master,'~/Desktop/time_table_master_wide.csv');
    writetable(all_run_frames,'~/Desktop/all_run_frames_wide.csv');
    
    all_images = struct;
    for rr = 1:10
        [data,all_images] = runme_vcdcore(subj_nr, ses_nr, 1, rr, '7TAS_BOLDSCREEN32', 'wanteyetracking', false, ...
            'all_images',all_images,'wantdatabypass',true,'exp_env',4, 'timetable_file',fullfile(tt_file(end).folder, tt_file(end).name));
    end
end

is_wide = true;
for st = ses_type
    clear behresults
    for rr = 1:10
        subj_folder = sprintf('vcd_subj%03d_ses%02d',subj_nr,ses_nr);
        mat_file    = sprintf('behavior_*_vcd_subj%03d_ses%02d_A_run%02d.mat',subj_nr,ses_nr,rr);
        dd = dir(fullfile(vcd_rootPath,'data','MRI',subj_folder,mat_file));
        load(fullfile(dd(end).folder,dd(end).name));
        performance = vcdbehavioralanalysis(fullfile(params.savedatafolder,params.behaviorfile));
        
        writetable(run_frames,sprintf('~/Desktop/run_frames%02d_wide.csv',rr));
        writetable(run_table,sprintf('~/Desktop/run_table%02d_wide.csv',rr));
        behresults(rr) = performance;
    end
    
    
    vcd_checkTimeTable(time_table_master,'behresults',behresults,'is_wide',is_wide)
end