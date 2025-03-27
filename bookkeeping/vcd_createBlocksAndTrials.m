function [p, condition_master, all_unique_im, all_cond] = vcd_makeMiniblocksAndTrials(p,varargin)
% VCD function to define, shuffle, and organize trials into miniblocks:
%
%   [condition_master, all_unique_im, all_cond] = ...
%       vcd_makeMiniblocksAndTrials(p,'load_params',<load_params>,'store_params',<store_params>)
%
% Purpose:
%   1. Define all possible stimulus-task combinations depending on session,
%   run, and subject.
% For every stimulus-task combination,
%   2. Define unique images (all images within a stimulus class, where
%       stimulus features of interest are fully-crossed.
%   3. Pseudo-randomly shuffle unique images according to prioritize
%       stimulus feature of interest.
%   4. Merge unique stimuli with left/right parafoveal stimulus locations
%       into a single trial, pseudo-randomly assign covert spatial
%       attention cue and fixation dot rim thickening direction.
%   5. Allocate trials to miniblocks
%
% INPUTS:
%  p              : (required) Struct with stimulus and session params
%  [load_params]  : (optional) Load prior stored parameters or not. Will
%                     search for a file name starting with trials_<dispname>*.mat
%                     in fullfile(vcd_rootPath,'workspaces','info').
%                     Default: true
%  [store_params] : (optional) Store generated parameters or not.
%                     Will store params in fullfile(vcd_rootPath,'workspaces','info');
%                     Default: true
%
% OUTPUT:
%  condition_master : N rows x M table that lists every unique image presented
%                       in VCD-core, for every task crossing separated by
%                       individual trial.
%                    Each row is an individual image.
%                    For Gabors/RDKs/Dots/ComplexObj, two subsequent rows
%                    define a unique trial. For NS, there is only 1 row per
%                    unique trial, as there is only one central image.
%
%                    N rows: nr unique trials/2 * nr of repeats * nr of tasks.
%                    M columns:
%                       1: unique im nr (integer)
%                       2: trial nr (every trial occupies two rows)
%                       2: thickening_dir (1=left, 2=right, 3=both)
%                       3:
%                       4: stim loc (1=left, 2=right, 3=central)
%                       5: cue status (0=uncued, 1 =cued)
%                       ... rows 6 to M-3 are stimClass specific
%                       M-2: paired_stim (learned associate for LTM)
%                       M-1: lure (for LTM)
%                       M: repeat nr
%
%  There are two hierarchical structs, which are the building blocks of
%  the condition master table, listing the unique images and unique 
%  conditions for each stimulus class:
%  * all_unique_im  : matrix labeling the unique image features
%  * all_conds      : stimulus "condition master" matrix describing the
%                    pseudo-random shuffled order of trials for one
%                    stimulus class, across all sessions.
% [TODO: document what the abbreviatio is, what the range is, what the units are, if it is an index or not]
% [TODO: Consider using a table instead of column]
%    * gabor : N trials x 12 cols matrix with the following stim columns
%       6: angle
%       7: contrast
%       8: phase
%       9: ref_delta (query stim for WM)
%
%   * rdk: N trials x 11 cols matrix with the following stim columns
%       6: motion_dir
%       7: coherence
%       8: ref_delta (query stim for WM)
%
%   * dot: N trials x 10 cols matrix with the following stim columns
%       6: dot_angle
%       7: ref_delta (query stim for WM)
%
%   * cobj: N trials x 13 cols matrix with the following stim columns
%       6: super_cat
%       7: basic_cat
%       8: sub_cat
%       9: facing_dir
%       10: ref_delta (query stim for WM)
%
%   *  ns: N trials x 12 cols matrix with the following stim columns
%       6: super_cat
%       7: basic_cat
%       8: sub_cat
%       9: change_blindness scenes (query stim for WM)
%
%
% Written by Eline Kupers November 2024 (kupers [at] umn [dot] edu)


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('p'             , @isstruct);
p0.addParameter('load_params'  , true, @islogical);
p0.addParameter('store_params' , true, @islogical);

% Parse inputs
p0.parse(p,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0

%% Load params if requested and we can find the file
if load_params
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('trials_%s*.mat',p.disp.name)));
    fprintf('[%s]: Found %d trial .mat file(s)\n',mfilename,length(d));
    if ~isempty(d)
        if length(d) > 1
            warning('[%s]: Multiple .mat files! Will pick the most recent one', mfilename);
        end
        fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d.name);
        load(fullfile(d(end).folder,d(end).name));
    else
        error('[%s]: Can''t find trial params file!', mfilename)
    end
else % Recreate conditions and miniblocks and trials
    
    % Bookkeeping structs
    all_unique_im     = struct();
    all_cond          = struct();
    condition_master  = table();
    
    
    %% Define miniblock content for each stimulus class
    for stimClass_idx = 1:length(p.exp.stimClassLabels)
        
        stim_table = table();
        
        % Get stimclass name
        stimClass_name = p.exp.stimClassLabels{stimClass_idx};
        
        % Get corresponding task crossings
        bsc_idx = cell2mat(cellfun(@(x) strcmp(p.exp.stimClassLabels,x), {stimClass_name}, 'UniformOutput', false));
        task_crossings = find(p.exp.crossings(bsc_idx,:));
        
        % Define the unique images for Gabors
        [t_cond, n_unique_cases] = vcd_defineUniqueImageNr(p, stimClass_name);
        
        all_unique_im.(stimClass_name) = t_cond;
        
        % Add stimclass name and idx to temporary table
        t_cond.stim_class_name = repmat({stimClass_name},size(t_cond,1),1);
        t_cond.stim_class = repmat(stimClass_idx,size(t_cond,1),1);
        
        % Loop over each task crossing for this stim class
        for curr_task = 1:length(task_crossings)
            
            % Add task name to temp table
            taskClass_name = p.exp.taskClassLabels{task_crossings(curr_task)};
            t_cond.task_class_name = repmat({taskClass_name},size(t_cond,1),1);
            t_cond.task_class      = repmat(task_crossings(curr_task),size(t_cond,1),1);
            
            % Is this a fixation task? if so, we need to double nr
            % of unique trials
            use_fix_flag = strcmp(taskClass_name,'fix');
            
            % Get the number of trials and miniblocks, which depend on task
            if p.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                n_trials_per_block    = p.exp.miniblock.n_trials_single_epoch;
                n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
            else
                n_trials_per_block    = p.exp.miniblock.n_trials_double_epoch;
                n_miniblocks_per_task = ceil(n_unique_cases/n_trials_per_block);
            end
            
            % Create condition master table, where unique images are
            % shuffled and distributed across miniblocks according to the
            % stimulus feature of interest that receive priority. For
            % example, we want to make sure we present all 3 gabor
            % contrasts within a miniblock at least once.
            tbl = vcd_createConditionMaster(p, t_cond, n_trials_per_block);
            
            %% Allocate trials to miniblocks
            if unique(tbl.stim_class)<5
                miniblock_trials = reshape(1:2:size(tbl,1),n_trials_per_block,[]);
            
                % preallocate space
                tbl.miniblock_nr = NaN(size(tbl,1),1);
                tbl.miniblock_local_trial_nr = NaN(size(tbl,1),1);
            
                % assign miniblock nr
                for mm = 1:size(miniblock_trials,2)
                    selected_trials = miniblock_trials(:,mm);

                        tbl.miniblock_nr(selected_trials)   = mm;
                        tbl.miniblock_nr(selected_trials+1) = mm;
                        tbl.miniblock_local_trial_nr(selected_trials) = 1:length(selected_trials);
                        tbl.miniblock_local_trial_nr(selected_trials+1) = 1:length(selected_trials);
                end
            else
                 miniblock_trials = reshape(1:size(tbl,1),n_trials_per_block,[]);
                 
                 tbl.miniblock_nr = NaN(size(tbl,1),1);
                 tbl.miniblock_local_trial_nr = NaN(size(tbl,1),1);
                 
                 % assign miniblock nr
                 for mm = 1:size(miniblock_trials,2)
                     selected_trials = miniblock_trials(:,mm);
                     
                     tbl.miniblock_nr(selected_trials)   = mm;
                     tbl.miniblock_local_trial_nr(selected_trials) = 1:length(selected_trials);
                 end
            end
            
            % Concatete single stim-task crossing table to master table
            stim_table = cat(1,stim_table,tbl);
            
        end % class idx
        
        % Accumulate condition master and info
        all_cond.(stimClass_name) = stim_table;
        
        condition_master = cat(1,condition_master,stim_table);
    end % stim idx
    
    
    %% Shuffle stimuli for SCC task
    condition_master = vcd_shuffleStimForSCC(condition_master, p.exp.miniblock.n_trials_single_epoch);
    
    
    %% At last, do some checks:
    
    % Stim class nr should match with params
    assert(isequal(unique(condition_master.stim_class)',1:length(p.exp.stimClassLabels)))
    % Task class nr should match with params
    assert(isequal(unique(condition_master.task_class)',1:length(p.exp.taskClassLabels)))
    
    % Cued vs uncued stimuli should be matched
    assert(isequal(sum(condition_master.iscued==1),sum(condition_master.iscued==0)))
    
    % Thickening direction should following cuing status
    assert(isequal(sum(condition_master.thickening_dir==1),sum(condition_master.thickening_dir==2)))
    assert(isequal(sum(condition_master.thickening_dir==1),sum(condition_master.iscued==1)))
    assert(isequal(sum(condition_master.thickening_dir==2),sum(condition_master.iscued==1)))
    
    % Now dive into each stimulus class
    N_tbl = NaN(5,30);
    for ii = unique(condition_master.stim_class)'
        [N, N_edges] = histcounts(condition_master.unique_im_nr(condition_master.stim_class==ii));
        N_tbl(ii,1:length(N)) = N;
        assert(length(unique(N))==1);
        assert(isequal(size(N,2), max(all_unique_im.(p.exp.stimClassLabels{ii}).unique_im_nr)))
        
        % Check nr of unique stimuli per stimulus class
        if ii == max(unique(condition_master.stim_class))
            assert(isequal(sum(condition_master.stim_class==ii),...
                max(condition_master.unique_im_nr(condition_master.stim_class==ii))*p.exp.n_unique_trial_repeats*(sum(p.exp.crossings(ii,:)))))
            unique_im_from_table(ii) = max(condition_master.unique_im_nr(condition_master.stim_class==ii));
        else
            assert(isequal(sum(condition_master.stim_class==ii),...
                max(condition_master.unique_im_nr(condition_master.stim_class==ii))*p.exp.n_unique_trial_repeats*(sum(p.exp.crossings(ii,:))-0.5)*p.exp.trial.stim_LR_loc_cue(ii)*2))
            unique_im_from_table(ii) = max(condition_master.unique_im_nr(condition_master.stim_class==ii));
        end
    end
    
    % Now dive into each task class
    crossings_unique_im = max(p.exp.crossings.*unique_im_from_table',2);
    
    M_tbl = NaN(10,30);
    for ii = unique(condition_master.task_class)'
        [M, edges_M] = histcounts(condition_master.unique_im_nr(condition_master.task_class==ii));
        M_tbl(ii,1:length(M)) = M;
    end
    
    % Store structs if requested
    if p.store_params
        fprintf('[%s]:Storing trial data..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir,sprintf('trials_%s_%s.mat',p.disp.name,datestr(now,30))),'condition_master','all_unique_im','all_cond')
    end
    
    
    % Plot figures to check stimulus order
    if p.verbose
        
        %%
        vcd_visualizeMasterTable(condition_master, p.store_imgs);
        
        
        figure; set(gcf,'Position',[1,1,1200,300]);
        histogram(condition_master.unique_im_nr)
        xlabel('unique condition nr'); ylabel('trial count')
        set(gca,'XTick',[1:30]); title('Total unique image nr')
        box off;
        if p.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs');
            filename = sprintf('vcd_totaluniqueim.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        figure; set(gcf,'Position',[1,1,1200,300]); imagesc(N_tbl)
        xlabel('unique im nr'); ylabel('Stim class');
        C = parula(20); colormap(C); set(gca,'CLim',[0 max(N_tbl(~isnan(N_tbl)))])
        cb = colorbar; cb.Label.String = 'unique im count';
        cb.Ticks = [0,unique(N_tbl(~isnan(N_tbl)))'];
        set(gca,'YTick',[1:5],'YTickLabel', p.exp.stimClassLabels);  set(gca,'XTick',[1:30]); title('Sum of unique im per stimulus class')
        box off;
        if p.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs');
            filename = sprintf('vcd_uniqueim_per_stimclass.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        
        figure; set(gcf,'Position',[1,1,1200,300]); imagesc(M_tbl)
        xlabel('unique condition nr'); ylabel('Task class');
        C = parula(20); colormap(C); set(gca,'CLim',[0 max(M_tbl(~isnan(M_tbl)))])
        cb = colorbar; cb.Label.String = 'unique im count';
        cb.Ticks = [0,unique(M_tbl(~isnan(M_tbl)))'];
        set(gca,'YTick',[1:10],'YTickLabel', p.exp.taskClassLabels);  set(gca,'XTick',[1:30]); title('Sum of unique im per task class')
        box off;
        if p.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs');
            filename = sprintf('vcd_uniqueim_per_taskclass.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        foo = condition_master.stim_class((condition_master.task_class==3),:);
        figure; set(gcf,'Position',[1,1,1200,300])
        imagesc(reshape(foo,p.exp.miniblock.n_trials_single_epoch,[]))
        xlabel('SCC miniblock nr'); ylabel('Trial nr within miniblock');
        C = parula(4); colormap(C); set(gca,'CLim',[min(foo),max(foo)])
        cb = colorbar; cb.Label.String = 'stim class nr';
        cb.Ticks = [1:max(foo)];
        set(gca,'YTick',[1:p.exp.miniblock.n_trials_single_epoch]); title('SCC stimulus class distribution')
        box off;
        if p.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs');
            filename = sprintf('vcd_scc_stimulusclass_trialdistr.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
    end
end % load params

return





