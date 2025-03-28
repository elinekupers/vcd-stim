function [params,time_table_master] = vcd_createSessions(params,varargin)
% VCD function to create image and event order for individual subject 
% sessions and runs
% 
%  [params,time_table_master] = vcd_createSessions(params,...
%         'load_params',[load_params],'store_params',[store_params])
%  

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'             , @isstruct);
p0.addParameter('load_params'  , true, @islogical);
p0.addParameter('store_params' , true, @islogical);

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
            error('[%s]: Can''t find trial.mat files! Please check or run vcd_makeTrials.m', mfilename);
        elseif ~isempty(d(end).name)
            if length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            end
            load(fullfile(d(end).folder,d(end).name),'condition_master');
            params.trials = condition_master;
        end
    end
    
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('time_table*.mat')));
    if isempty(d)
        error('[%s]: Can''t find time table file with subject session!', mfilename)
    elseif ~isempty(d(end).name)
        if length(d) > 1
            warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
        end
        load(fullfile(d(end).folder,d(end).name),'time_table_master');
        params.time_table = time_table_master;   
    end
    
else
    % Create subject sessions
    
    % check if trial struct is already defined and load it if needed
    if ~isfield(params,'trials') || isempty(params.trials)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info','trials*.mat'));
        
        if ~isempty(d(end).name)
            if length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            end
            load(fullfile(d(end).folder,d(end).name),'condition_master');
            params.trials = condition_master;
        else
            error('[%s]: Can''t find trial.mat files! Please check or run vcd_makeTrials.m', mfilename);
        end
    end
    
    % check if session struct is already defined and load it if needed
    if ~isfield(params,'exp') || isempty(params.exp)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info','exp_session*.mat'));
        
        if length(d) == 1
            load(fullfile(d(end).folder,d(end).name),'exp_session');
            params.exp = exp_session;
        elseif length(d) > 1
            warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            load(fullfile(d(end).folder,d(end).name),'exp_session');
            params.exp = exp_session;
        else
            error('[%s]: Can''t find exp session .mat files! Please check or run vcd_getSessionParams.m', mfilename);
        end
    end

    
    %% 1. Create condition master, defining all unique trials and repeats of trials.
    % Per task/stim crossing: get nr of trials per block
    % Different trial order per block (already accomplished in vcd_makeTrials.m)
    % Across a single session, each subject will experience the same
    % blocks and unique images
    if isempty(find(strcmp(params.trials.Properties.VariableNames,'session_nr')))
        params.trials = vcd_allocateBlocksToRuns(params);
    end
    %% We expand the condition master table and add all trial events (in units of presentationrate_hz frames), 
    % We also shuffle blocks within a run for each subject session
    time_table_master = vcd_createRunTimeTables(params);
    
    params.time_table = time_table_master;
    
    %% 
    
    
end

% if p.verbose
%     for sj = 1:length(p.exp.total_subjects)
%         
%         for ses = 1:length(p.exp.n_sessions)
%             
%             figure(ses); clf; set(gcf,'Position',[1,1,1920,1080])
%             for rr = 1:p.exp.n_runs_per_session
%                 
%                 %                 run_data = subj_master_time_table.subj
%                 
%                 %                 timepoints = 0:
%                 
%                 subplot(ceil(p.exp.n_runs_per_session/2),2,rr)
%                 stem(timepoints,conditions); hold all;
%                 stem(imID(:,1),imID(:,2),'r');
%                 xlabel('time (s)')
%                 ylabel('stim-task crossing')
%                 set(gca,'YTick',[0:length(p.exp.stimTaskLabels)+2],'YTickLabel',[{'blank'}; p.exp.stimTaskLabels;{'taskcue'};{'IBI'}])
%                 set(gca,'TickDir', 'out')
%                 title(sprintf('run %d',ii))
%             end
%             %         if p.store_imgs
%             %             print(gcf,'-dpng','-r300',fullfile(vcd_rootPath,'figs',sprintf('ses%02d_run_order_master',ses)));
%             %         end
%         end
%     end
% end

return
 
