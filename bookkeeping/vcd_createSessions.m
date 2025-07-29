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
% Example:
% { ...
% params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);
% params.exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
%    vcd_createSessions(params);
% }
%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'                 , @isstruct);
p0.addParameter('condition_master'      , [], @istable);
p0.addParameter('load_params'           , false, @islogical);
p0.addParameter('store_params'          , true, @islogical);
p0.addParameter('store_imgs'            , false, @islogical);
p0.addParameter('verbose'               , false, @islogical);
p0.addParameter('env_type'              , '', @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));
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


% Infer environment if we haven't set this parameter  
if ~exist('env_type','var') || isempty(env_type)
    if strcmp(params.disp.name,'7TAS_BOLDSCREEN32')
        env_type = 'MRI';
    elseif ismember(params.disp.name,{'CCNYU_VIEWPIXX3D','PPROOM_EIZOFLEXSCAN'}) 
        env_type = 'BEHAVIOR';
    elseif ismember(params.disp.name,{'KKOFFICE_AOCQ3277','EKHOME_ASUSVE247'}) 
        env_type = 'BEHAVIOR';
    end
end

%% Load subject files if we can
if ~isempty(randomization_file_pth)
    % load provided randomized time table master file
    a1 = load(fullfile(randomization_file_pth),'time_table_master','all_run_frames');
    % allocate variables to output arguments
    time_table_master_shuffled  = a1.time_table_master;
    all_subj_run_frames         = a1.all_run_frames;
    
    % find subject-specific condition_master file (assume it lives in the
    % same folder as subject-specific time table master file.
    if ~isfield(params,'is_demo'), params.is_demo = false; end
    if ~isfield(params,'is_wide'), params.is_wide = false; end
    
    subj_data_folder = fileparts(randomization_file_pth);
    fname = sprintf('%s_condition_master_wide_%s%s%s_%s_*.mat',subj_id,choose(params.is_wide,'wide_',''), choose(params.is_demo,'demo_',''), params.disp.name);
    a1 = load(fullfile(subj_data_folder,fname),'condition_master_shuffled');
    condition_master_shuffled   = a1.condition_master_shuffled;
else
    
    if load_params % load condition_master and params if requested
        if ~isfield(params,'is_demo'), params.is_demo = false; end
        if ~isfield(params,'is_wide'), params.is_wide = false; end
        % load condition master
        fname = sprintf('condition_master_%s%s%s*.mat',choose(params.is_wide,'wide_',''),choose(params.is_demo,'demo_',''),params.disp.name);
            
        d = dir(fullfile(vcd_rootPath,'workspaces','info',fname));
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
        saveDir = fullfile(vcd_rootPath,'data',env_type, subj_id);
    end

    
    %% 1. We shuffle the block order within a session (as well as trial order within a block)
    % and allocate blocks randomly (under some constraints) to each run.
    % We then store a file that records this new block order
    [condition_master_shuffled, ~, ~, ~, ~] = ... % other outputs are: fname, condition_master_shuffle_idx, session_crossing_matrix, session_block_matrix
        vcd_randomizeBlocksAndRunsWithinSession(params, ...
        condition_master, env_type, 'subj_id', subj_id, 'saveDir',saveDir);
    
    %% 2. We expand the condition master table and add all trial events
    % (in units of presentationrate_hz frames), and we now call it
    % "time_table_master".
    
    [time_table_master_shuffled,all_subj_run_frames] = ...
        vcd_createRunTimeTables(params, ...
        'load_params', load_params, ...
        'store_params', store_params, ...
        'verbose', verbose, ...
        'condition_master',condition_master_shuffled,...
        'env_type',env_type, ...
        'saveDir',saveDir, ...
        'subj_id',subj_id);
    
end



return

