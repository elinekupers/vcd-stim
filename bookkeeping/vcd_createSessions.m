function [params,time_table_master, all_subj_run_frames] = vcd_createSessions(params,varargin)
% VCD function to create image and event order for individual subject 
% sessions and runs
% 
%  [params,time_table_master,all_subj_run_frames] = vcd_createSessions(params,condition_master,...
%         'load_params',[load_params],'store_params',[store_params])
%  

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'           , @isstruct);
p0.addParameter('condition_master', [], @istable);
p0.addParameter('load_params'     , true, @islogical);
p0.addParameter('store_params'    , true, @islogical);
p0.addParameter('session_env'     , 'MRI', @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));
p0.addParameter('randomization_file_pth', [], @ischar);
p0.addParameter('prefix'          , 'subjXXX_', @ischar);

% Parse inputs
p0.parse(params,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

%% Load params if requested and we can find the file
if load_params
   
    % load trial info
    d = dir(fullfile(vcd_rootPath,'workspaces','info','condition_master*.mat'));
    if isempty(d)
        error('[%s]: Can''t find trial.mat files! Please check or run vcd_makeTrials.m\n', mfilename);
    elseif ~isempty(d(end).name)
        if length(d) > 1
            warning('[%s]: Multiple trial.mat files! Will pick the most recent one.\n', mfilename);
        end
        load(fullfile(d(end).folder,d(end).name),'condition_master');
    end
    
    
%     d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('time_table_master*%s*.mat',params.disp.name)));
%     if isempty(d)
%         error('[%s]: Can''t find time table file with subject session!\n', mfilename)
%     elseif ~isempty(d(end).name)
%         if length(d) > 1
%             warning('[%s]: Multiple trial .mat files! Will pick the most recent one.\n', mfilename);
%         end
%         load(fullfile(d(end).folder,d(end).name),'time_table_master');
%         params.time_table = time_table_master;   
%     end
    
else
    % Create subject sessions
    
    % check if trial struct is already defined and load it if needed
    if ~exist('condition_master','var') || isempty(condition_master)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info',['condition_master*' session_env '*.mat']));
        
        if ~isempty(d(end).name)
            if length(d) > 1
                warning('[%s]: Multiple condition_master .mat files! Will pick the most recent one', mfilename);
            end
            load(fullfile(d(end).folder,d(end).name),'condition_master');
        else
            error('[%s]: Can''t find condition_master .mat files! Please check or run vcd_createConditions.m', mfilename);
        end
    end
    
    % check if experimental parameters are already defined and load it if needed
    if ~isfield(params,'exp') || isempty(params.exp)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info','exp*.mat'));
        
        if length(d) == 1
            load(fullfile(d(end).folder,d(end).name),'exp');
            params.exp = exp;
        elseif length(d) > 1
            warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            load(fullfile(d(end).folder,d(end).name),'exp');
            params.exp = exp;
        else
            warning('[%s]: Can''t find exp session .mat files! Will run vcd_getSessionParams.m', mfilename);
            params.exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
        end
    end
    
    % check if stimulus parameters are already defined and load it if needed
    if ~isfield(params,'stim') || isempty(params.stim)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info','stim*.mat'));
        
        if ~isempty(d)
            if length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            end
            load(fullfile(d(end).folder,d(end).name),'stim');
            params.stim = stim;
        else
            warning('[%s]: Can''t find stim session .mat files! Will run vcd_getStimParams.m', mfilename);
            params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);        
        end
    end
    
    if ~exist('randomization_file_pth','var') || isempty(randomization_file_pth)
        
        %% 1. We randomize the blocks and runs within a session and store a file that records this order
        [condition_master_shuffled, ~, ~, ~] = ...
            vcd_randomizeBlocksAndRunsWithinSession(params, condition_master, session_env, 'prefix', prefix);
    else
        % load provided randomized condition master file
        load(fullfile(randomization_file_pth,sprintf('%scondition_master_shuffled_%s_%s_*.mat',prefix, params.disp.name, session_env)),'condition_master_shuffled')
    end
    
    %% 2. We expand the condition master table and add all trial events 
    % (in units of presentationrate_hz frames), and we now call it
    % "time_table_master"
    % We also shuffle blocks within a run for each subject session
    time_table_master = vcd_createRunTimeTables(params,condition_master_shuffled,session_env);
    
    %% 3. We expand the "time_table_master" with the fixation sequence 
    % and onset of contrast dip, and correct button presses for FIX and CD
    % task-crossings
    [time_table_master,all_subj_run_frames] = vcd_addFIXandCDtoTimeTableMaster(params,time_table_master,session_env);
    
    
end


return
 
