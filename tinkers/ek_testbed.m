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
'all_images',all_images,'wantdatabypass',true,'exp_env',4); 

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
subj_nr = 090;
ses_nr = 1;
ses_type = 1; %[1,2];

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

%%
ses_nr   = 1;
subj_nr  = [1,2];
ses_type = [1,2];
is_wide  = true;
results  = {};
ddir     = pwd; %fullfile(vcd_rootPath,'data','MRI');
for st = ses_type
    clear behresults
    for rr = 1:10
        subj_folder = sprintf('testwide1%s/mat_files_from_scan',choose(ses_type(st)==1,'A','B')); %sprintf('vcd_subj%03d_ses%02d',subj_nr,ses_nr);
        mat_file    = sprintf('behavior_*_vcd_subj%03d_ses%02d_%s_run%02d.mat',subj_nr(st),ses_nr, choose(ses_type(st)==1,'A','B'), rr);
        dd = dir(fullfile(ddir,subj_folder,mat_file));
        performance = vcdbehavioralanalysis(fullfile(dd(end).folder,dd(end).name));

        tt_file = dir(fullfile(ddir, subj_folder, sprintf('vcd_subj%03d_time_table_master_wide_7TAS_BOLDSCREEN32_*.mat',subj_nr(st))));
        load(fullfile(tt_file(end).folder, tt_file(end).name));
        
%         writetable(run_frames,fullfile(ddir,sprintf('/run_frames%02d_wide.csv',rr)));
%         writetable(run_table,fullfile(ddir,sprintf('run_table%02d_wide.csv',rr)));
        behresults(rr) = performance;
    end
    
    
    results{st} = vcd_checkTimeTable(time_table_master,'behresults',behresults,'is_wide',is_wide);
end


%%

%% DRY-RUN 10 RUNS OF deep MRI EXPERIMENT
subj_nr = 999; %996;
ses_nrs = 46;
ses_types = 2; %[1,2];

data_dir = fullfile(vcd_rootPath, 'data','MRI',sprintf('vcd_subj%03d',subj_nr));

for ses_nr = ses_nrs
    for st = ses_types

%         if ses_nr == 8
%              all_images = struct;
%             [data,all_images] = runme_vcdcore(subj_nr, ses_nr, st, 1, '7TAS_BOLDSCREEN32', 'wanteyetracking', false, ...
%             'all_images',all_images,'wantdatabypass',true,'exp_env',4, 'is_wide', false);
%             
%           tt_file = '/Users/kupers/projects/git/toolboxes/vcd-stim/data/MRI/vcd_subj999/vcd_subj999_condition_master_deep_7TAS_BOLDSCREEN32_20260219T140250.mat';
%             tt_file    = dir(fullfile(data_dir, sprintf('vcd_subj%03d_time_table_master_deep_7TAS_BOLDSCREEN32_2026*.mat',subj_nr)));
%             load(fullfile(tt_file(end).folder, tt_file(end).name));
%             writetable(time_table_master,'~/Desktop/time_table_master_deep.csv');
%             writetable(all_run_frames,'~/Desktop/all_run_frames_deep.csv');
%         end
        
        all_images = struct;
        for rr = 1:10
%             if rr == 1 && ses_nr == 1 && st == 1
%                 [data,all_images] = runme_vcdcore(subj_nr, ses_nr, st, rr, '7TAS_BOLDSCREEN32', 'wanteyetracking', false, 'is_wide',false, ...
%                     'all_images',all_images,'wantdatabypass',true,'exp_env',4);
%             else
                tt_file = dir(fullfile(data_dir, sprintf('vcd_subj%03d_time_table_master_deep_7TAS_BOLDSCREEN32_2026*.mat',subj_nr)));

                [data,all_images] = runme_vcdcore(subj_nr, ses_nr, st, rr, '7TAS_BOLDSCREEN32', 'wanteyetracking', false, 'is_wide',false, ...
                    'all_images',all_images,'wantdatabypass',true,'exp_env',4, 'timetable_file',fullfile(tt_file(end).folder, tt_file(end).name));
%             end
        end
    end
end

%%
ses_nr   = 1;
subj_nr  = 996;
ses_type = 1;
is_wide  = false;
results  = {};
ddir     = fullfile(vcd_rootPath,'data','MRI');

clear behresults
for ses = 2:7
    for st = 1

        for rr = 1:10
            subj_folder = sprintf('vcd_subj%03d_ses%02d',subj_nr,ses);
            mat_file    = sprintf('behavior_*_vcd_subj%03d_ses%02d_%s_run%02d.mat',subj_nr(st),ses,choose(st==1,'A','B'), rr);
            dd = dir(fullfile(ddir,subj_folder,mat_file));
            performance = vcdbehavioralanalysis(fullfile(dd(end).folder,dd(end).name));
            behresults(ses,st,rr) = performance;
            
%             tt_file = dir(fullfile(ddir, sprintf('vcd_subj%03d',subj_nr), sprintf('vcd_subj%03d_time_table_master_deep_7TAS_BOLDSCREEN32_*.mat',subj_nr)));
%             load(fullfile(tt_file(end).folder, tt_file(end).name));
%             run_frames = all_run_frames(all_run_frames.session_nr==ses & all_run_frames.run_nr==rr & all_run_frames.session_type==st,:);
%             run_table = time_table_master(time_table_master.session_nr==ses & time_table_master.run_nr==rr & time_table_master.session_type==st,:);
%             writetable(run_frames,fullfile(ddir,sprintf('/run_frames%02d_deep_session%02d.csv',rr,ses)));
%             writetable(run_table,fullfile(ddir,sprintf('run_table%02d_deep_session%02d.csv',rr,ses)));
            
        end
    end
    
    % move files to folder called 'mat_files_from_scan' as this is the
    % place where analysis_behavior looks for the files
    olddir = fullfile(pwd,'data/MRI',sprintf('vcd_subj%03d_ses%02d/',subj_nr,ses_nr));
    newdir = fullfile(pwd,'data/MRI',sprintf('vcd_subj%03d_ses%02d',subj_nr,ses_nr),'mat_files_from_scan/');
    mkdirquiet(newdir)
    movefile(sprintf('%s/*',olddir), newdir)
    
    analysis_behavior(fullfile(pwd,'data/MRI',sprintf('vcd_subj%03d_ses%02d',subj_nr,ses_nr)),fullfile(pwd,'results/MRI',sprintf('vcd_subj%03d_ses%02d',subj_nr,ses_nr)))
    
end
results{st} = vcd_checkTimeTable(time_table_master,'behresults',behresults,'is_wide',is_wide);


