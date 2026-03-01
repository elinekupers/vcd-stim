% s_swapGlobalTrialNumbersConditionMaster.m

% Step 1: load condition master
% Step 2: load/get params: 
% % params.disp = vcd_getDisplayParams('7TAS_BOLDSCREEN32');
% % params.stim = vcd_getStimParams('disp_name','7TAS_BOLDSCREEN32');
% % params.exp = vcd_getSessionParams;
% Step 3: run code to randomize blocks:
% % condition_master_shuffled = ...
% %     vcd_randomizeBlocksAndRunsWithinSession(params, ...
% %     condition_master, env_type, 'subj_id', subj_id, 'saveDir',saveDir);

st = 1;

% Step 4: Copy suggested global trial numbers into A1 and A2
A1 = 12975;
A2 = 13182;

A1_idx = condition_master.global_trial_nr==A1; % get condition master table index
A2_idx = condition_master.global_trial_nr==A2; % get condition master table index

% Step 5: Define trials to swap
% NOTE: We start from column 9 and leave column 1:8 with session_nr/session_type/block_nr/trial_nr/etc. 
% We only want to swap the stimulus-related trial info.
tmp1 = condition_master(A1_idx & condition_master.session_type==st,9:end); 
tmp2 = condition_master(A2_idx & condition_master.session_type==st,9:end);

% Step 6: Display trials to swap, and the trials in the blocks. Confirm
% there are indeed 3 or more of the same stimuli in one block
tmp1
tmp2
condition_master(condition_master.global_block_nr==condition_master.global_block_nr(A1_idx & condition_master.session_type==st) & condition_master.session_type==st,:)
condition_master(condition_master.global_block_nr==condition_master.global_block_nr(A2_idx & condition_master.session_type==st) & condition_master.session_type==st,:)

%% Step 7: Do the swap
condition_master(A2_idx & condition_master.session_type==st,9:end) = tmp1; 
condition_master(A1_idx & condition_master.session_type==st,9:end) = tmp2;

%% Step 8: Resave large condition_master
params.is_wide = false;
params.is_demo = false;
fname = sprintf('condition_master_%s%s%s_%s.mat',choose(params.is_wide,'wide_','deep_'), choose(params.is_demo,'demo_',''),params.disp.name,datestr(now,30));
saveDir = fullfile(vcd_rootPath, 'workspaces','info');
save(fullfile(saveDir,fname),'condition_master','all_unique_im','all_cond');

%% Step 9: Repeat
subj_id = '999';
env_type = 'MRI';

condition_master_shuffled = ...
     vcd_randomizeBlocksAndRunsWithinSession(params, ...
     condition_master, env_type, 'subj_id', subj_id, 'saveDir',saveDir);
 
 
 