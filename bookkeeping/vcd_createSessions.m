function [params,condition_master_shuffled,time_table_master_shuffled, all_subj_run_frames] = vcd_createSessions(params,varargin)
% VCD function to create image and event order for individual subject
% sessions and runs.
%
%  [params,condition_master_shuffled,time_table_master_shuffled, all_subj_run_frames] = ...
%      vcd_createSessions(params,condition_master,...
%         'load_params',[load_params],'store_params',[store_params])
%
% [WRITE ME]
%
%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'                 , @isstruct);
p0.addParameter('condition_master'      , [], @istable);
p0.addParameter('load_params'           , false, @islogical);
p0.addParameter('store_params'          , true, @islogical);
p0.addParameter('store_imgs'            , false, @islogical);
p0.addParameter('verbose'               , false, @islogical);
p0.addParameter('session_env'           , 'MRI', @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));
p0.addParameter('randomization_file_pth', [], @ischar);
p0.addParameter('subj_id'               , 'vcd_subj000', @ischar);
p0.addParameter('saveDir'               , [], @ischar)

% Parse inputs
p0.parse(params,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0


%% Load subject files if we can
if ~isempty(randomization_file_pth)
    % load provided randomized condition master file
    load(fullfile(randomization_file_pth,sprintf('%scondition_master_shuffled_%s_%s_*.mat',[subj_id '_'], params.disp.name, session_env)),'condition_master_shuffled');
    load(fullfile(randomization_file_pth,sprintf('%stime_table_master_shuffled_%s_%s_*.mat',[subj_id '_'], params.disp.name, session_env)),'time_table_master_shuffled','all_subj_run_frames');
else
    
    if load_params % load condition_master and params if requested
        
        % load condition master
        d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('condition_master*%s*.mat',params.disp.name)));
        if isempty(d)
            error('[%s]: Can''t find condition master .mat files! Please check or run vcd_createConditions.m\n', mfilename);
        elseif ~isempty(d(end).name)
            if length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one.\n', mfilename);
            end
            a1 = load(fullfile(d(end).folder,d(end).name));
            condition_master = a1.condition_master; clear a1
        end
        
        d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('stim_%s*.mat',params.disp.name)));
        if ~isempty(d)
            if verbose
                fprintf('[%s]: Found %d stim params .mat file(s)\n',mfilename,length(d));
                if length(d) > 1
                    warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
                end
                fprintf('[%s]: Loading stim params .mat file: %s\n', mfilename, d(end).name);
            end
            load(fullfile(d(end).folder,d(end).name),'stim');
            params.stim = stim; clear stim;
        else
            error('[%s]: Can''t find stim params file!\n', mfilename)
        end
        
        d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('exp_%s*.mat',params.disp.name)));
        if ~isempty(d)
            if verbose
                fprintf('[%s]: Found %d exp params .mat file(s)\n',mfilename,length(d));
                if length(d) > 1
                    warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
                end
                fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
            end
            load(fullfile(d(end).folder,d(end).name),'exp');
            params.exp = exp; clear exp;
        else
            error('[%s]: Can''t find stim params file!\n', mfilename)
        end
        
    else
        
        % check if experimental parameters are already defined and load it if needed
        if ~isfield(params,'exp') || isempty(params.exp)
            warning('[%s]: Can''t find exp field in params! Will run vcd_getSessionParams.m', mfilename);
            params.exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
        end
        
        % check if stimulus parameters are already defined and load it if needed
        if ~isfield(params,'stim') || isempty(params.stim)
            warning('[%s]: Can''t find stim field in params! Will run vcd_getStimParams.m', mfilename);
            params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose', false);
        end
    end  
        
    
    %% 0. set saveDir if we haven't already
    if store_params && isempty(saveDir)
        saveDir = fullfile(vcd_rootPath,'data',session_env, subj_id);
    end
    
    %% 1. We shuffle the block order within a session (as well as trial order within a block)
    % and allocate blocks randomly (under some constraints) to each run.
    % We then store a file that records this new block order
    [condition_master_shuffled, ~, ~, ~, ~] = ... % other outputs are: fname, condition_master_shuffle_idx, session_crossing_matrix, session_block_matrix
        vcd_randomizeBlocksAndRunsWithinSession(params, ...
        condition_master, session_env, 'subj_id', subj_id, 'saveDir',saveDir);
    
    %% 2. We expand the condition master table and add all trial events
    % (in units of presentationrate_hz frames), and we now call it
    % "time_table_master".
    
    [time_table_master_shuffled,all_subj_run_frames] = ...
        vcd_createRunTimeTables(params, ...
        'load_params', load_params, ...
        'store_params', store_params, ...
        'verbose', verbose, ...
        'condition_master',condition_master_shuffled,...
        'session_env',session_env, ...
        'saveDir',saveDir, ...
        'subj_id',subj_id);
    
end



return

