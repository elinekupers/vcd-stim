function [p,subj_master_time_table] = vcd_createSessions(p,load_params,store_params)

if ~exist('store_params','var') || isempty(store_params)
    store_params = true;
end

if ~exist('load_params','var') || isempty(load_params)
    load_params = true;
end

if load_params
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('subject_sessions*.mat')));
    if ~isempty(d)
        load(fullfile(d(end).folder,d(end).name),'subj_master_time_table');
    else
        error('[%s]: Can''t find subject sessions file!', mfilename)
    end
    
else
    % Create subject sessions
    
    % check if trial struct is already defined and load it if needed
    if ~isfield(p,'trials') || isempty(p.trials)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info','trials*.mat'));
        
        if ~isempty(d(end).name)
            if length(d) > 1
                warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            end
            load(fullfile(d(end).folder,d(end).name),'all_trials');
            p.trials = all_trials;
        else
            error('[%s]: Can''t find trial.mat files! Please check or run vcd_makeTrials.m', mfilename);
        end
    end
    
    % check if session struct is already defined and load it if needed
    if ~isfield(p,'exp') || isempty(p.exp)
        
        % load trial info
        d = dir(fullfile(vcd_rootPath,'workspaces','info','exp_session*.mat'));
        
        if length(d) == 1
            load(fullfile(d(end).folder,d(end).name),'exp_session');
            p.exp = exp_session;
        elseif length(d) > 1
            warning('[%s]: Multiple trial.mat files! Will pick the most recent one', mfilename);
            load(fullfile(d(end).folder,d(end).name),'exp_session');
            p.exp = exp_session;
        else
            error('[%s]: Can''t find exp session .mat files! Please check or run vcd_getSessionParams.m', mfilename);
        end
    end
    
    %% Per task/stim crossing: get nr of trials per miniblock
    % Different trial order per block (already accomplished in vcd_makeTrials.m)
    % Across a single session, each subject will experience the same
    % miniblocks and unique images
    [p, ~] = vcd_allocateMiniblocksToRuns(p);
    
    
    %% Now we expand the condition master table and add all trial events, 
    % as well as shuffling blocks within a run for each subject session
    subj_master_time_table = vcd_createRunTimeTables(p);
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
 
