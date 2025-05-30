function [params,time_table_master, all_subj_run_frames] = vcd_createSessions(params,varargin)
% VCD function to create image and event order for individual subject 
% sessions and runs
% 
%  [params,time_table_master,all_subj_run_frames] = vcd_createSessions(params,condition_master,...
%         'load_params',[load_params],'store_params',[store_params])
%  

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'             , @isstruct);
p0.addParameter('load_params'  , true, @islogical);
p0.addParameter('store_params' , true, @islogical);
p0.addParameter('session_type' , 'MRI', @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));

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
    % check if trial struct is already defined and load it if needed
    if ~isfield(params,'trials') || isempty(params.trials)
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info','trials*.mat'));
        if isempty(d)
            error('[%s]: Can''t find trial.mat files! Please check or run vcd_makeTrials.m\n', mfilename);
        elseif ~isempty(d(end).name)
            if length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one.\n', mfilename);
            end
            load(fullfile(d(end).folder,d(end).name),'condition_master');
        end
    end
    
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('time_table_master*%s*.mat',params.disp.name)));
    if isempty(d)
        error('[%s]: Can''t find time table file with subject session!\n', mfilename)
    elseif ~isempty(d(end).name)
        if length(d) > 1
            warning('[%s]: Multiple trial .mat files! Will pick the most recent one.\n', mfilename);
        end
        load(fullfile(d(end).folder,d(end).name),'time_table_master');
        params.time_table = time_table_master;   
    end
    
else
    % Create subject sessions
    
    % check if trial struct is already defined and load it if needed
    if ~exist('condition_master','var') || isempty(condition_master)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info',['trials*' session_type '*.mat']));
        
        if ~isempty(d(end).name)
            if length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            end
            load(fullfile(d(end).folder,d(end).name),'condition_master');
        else
            error('[%s]: Can''t find trial.mat files! Please check or run vcd_makeTrials.m', mfilename);
        end
    end
    
    % check if session struct is already defined and load it if needed
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
            error('[%s]: Can''t find exp session .mat files! Please check or run vcd_getSessionParams.m', mfilename);
        end
    end

    
    %% 1. Create condition master, defining all unique trials and repeats of trials.
    % Per task/stim crossing: get nr of trials per block
    % Different trial order per block (already accomplished in vcd_makeTrials.m)
    % Across a single session, each subject will experience the same
    % blocks and unique images
    condition_master = vcd_allocateBlocksToRuns(params,condition_master,session_type);

    %% 2. We expand the condition master table and add all trial events 
    % (in units of presentationrate_hz frames), and we now call it
    % "time_table_master"
    % We also shuffle blocks within a run for each subject session
    time_table_master = vcd_createRunTimeTables(params,condition_master,session_type);
    
    %% 3. We expand the "time_table_master" with the fixation sequence 
    % and onset of contrast dip, and correct button presses for FIX and CD
    % task-crossings
    [time_table_master,all_subj_run_frames] = vcd_addFIXandCDtoTimeTableMaster(params,time_table_master,session_type);
    
    
end


return
 
