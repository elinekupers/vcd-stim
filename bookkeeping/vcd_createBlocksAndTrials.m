function [p, condition_master, all_unique_im, all_cond] = vcd_createBlocksAndTrials(p,varargin)
% VCD function to define, shuffle, and organize trials into blocks:
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
%   5. Allocate trials to blocks
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
%   Columns contains the following information:
%   1: {'unique_im_nr'   } unique images (parafoveal stimulus patches: 1 trial occupies 2 rows)
%   2: {'stimloc'        } stimulus location relative to fixation dot: 1 = left, 2 = right, 3 = central
%   3: {'stimloc_name'   } same as column two but then in text 
%   4: {'orient_dir'     } orientation (gabor), motion direction (rdk), or angle (dot), or facing direction (obj) in deg
%   5: {'contrast'       } stimulus contrast (Michelson fraction)
%   6: {'gbr_phase'      } gabor stimulus phase  (NaN for non gabor stim)
%   7: {'rdk_coherence'  } rdk stimulus dot coherence (fraction of dots, NaN for non rdk stim)
%   8: {'super_cat'      } superordinate object category level (for obj and ns) 
%   9: {'basic_cat'      } basic object category level (for obj and ns)
%   10: {'sub_cat'        } subordinate object category level (for obj and ns)
%   11: {'super_cat_name' } same as 8, but then in text
%   12: {'basic_cat_name' } same as 9, but then in text
%   13: {'sub_cat_name'   } same as 10, but then in text
%   14: {'stim_class_name'} stimulus class name (gabor, rdk, dot, obj, ns)
%   15: {'stim_class'     } stim class nr (1:5)
%   16: {'task_class_name'} task class name (fix/cd/scc/pc/wm/ltm/img/what/where/how)
%   17: {'task_class'     } task class nr (1:10)
%   18: {'iscued'         } cue status: 0 = uncued, 1 = cued
%   19: {'unique_trial_nr'} unique trial nr
%   20: {'thickening_dir' } thickening direction of fixation dot rim (for spatial attention cue) 1 = left, 2 = right, 3 = both/neutral
%   21: {'stim2_delta'    } for double epoch tasks, what predefined delta between stim 1 and stim 2 did we choose from p.stim.(<stim_class_name>).delta_ref
%   22: {'stim2'          } same as 21, but then update stim feature of stim2 (e.g., 80 deg orientation for a stim1: 95 - 15 deg)
%   23: {'ltm_stim_pair'  } for LTM task: each unique image nr is associated with
%                           another stimulus
%   24: {'islure'         } for LTM task: whether we use lure stimulus or not
%   25: {'repeat_nr'      } Keep track how many times has this unique image has
%                           been repeated thusfar in the experiment
%
%  There are two hierarchical structs, which are the building blocks of
%  the condition master table, listing the unique images and unique 
%  conditions for each stimulus class:
%  * all_unique_im  : matrix labeling the unique image features
%  * all_conds      : stimulus "condition master" matrix describing the
%                    pseudo-random shuffled order of trials for one
%                    stimulus class, across all sessions.

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
else % Recreate conditions and blocks and trials
    
    % Bookkeeping structs
    all_unique_im     = struct();
    all_cond          = struct();
    condition_master  = table();
    
    
    %% Define block content for each stimulus class
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
            
            % Get the number of trials and blocks, which depend on task
            if p.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                n_trials_per_block    = p.exp.block.n_trials_single_epoch;
            else
                n_trials_per_block    = p.exp.block.n_trials_double_epoch;
            end
            
            %% ---- IMPORTANT FUNCTION: Create condition master table ---- %%
            % Create condition master table, where unique images are
            % shuffled and distributed across blocks according to the
            % stimulus feature of interest that receive priority. For
            % example, we want to make sure we present all 3 gabor
            % contrasts within a block at least once.
            tbl = vcd_createConditionMaster(p, t_cond, n_trials_per_block);
            
            %% --- Allocate trials to blocks ---
            block_trials = reshape(1:size(tbl,1),n_trials_per_block,[]);
            
            tbl.stim_class_unique_block_nr = NaN(size(tbl,1),1);
            tbl.block_local_trial_nr = NaN(size(tbl,1),1);
            
            % assign block nr
            for mm = 1:size(block_trials,2)
                selected_trials = block_trials(:,mm);
                
                tbl.stim_class_unique_block_nr(selected_trials)   = mm;
                tbl.block_local_trial_nr(selected_trials) = 1:length(selected_trials);
            end

            % Concatete single stim-task crossing table to master table
            stim_table = cat(1,stim_table,tbl);
            
        end % class idx
        
        % Accumulate condition master and info
        all_cond.(stimClass_name) = stim_table;
        
        condition_master = cat(1,condition_master,stim_table);
    end % stim idx
    
    
    %% ---- IMPORTANT FUNCTION: Shuffle stimuli for SCC task ----
    condition_master = vcd_shuffleStimForSCC(condition_master, p.exp.block.n_trials_single_epoch);
%     
%     condition_master = vcd_shuffleStimForLTM(condition_master, p.exp.block.n_trials_double_epoch);

    
    %% Store structs if requested
    if p.store_params
        fprintf('[%s]:Storing trial data..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir,sprintf('trials_%s_%s.mat',p.disp.name,datestr(now,30))),'condition_master','all_unique_im','all_cond')
    end
    
    
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
        imagesc(reshape(foo,p.exp.block.n_trials_single_epoch,[]))
        xlabel('SCC block nr'); ylabel('Trial nr within block');
        C = parula(4); colormap(C); set(gca,'CLim',[min(foo),max(foo)])
        cb = colorbar; cb.Label.String = 'stim class nr';
        cb.Ticks = [1:max(foo)];
        set(gca,'YTick',[1:p.exp.block.n_trials_single_epoch]); title('SCC stimulus class distribution')
        box off;
        if p.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs');
            filename = sprintf('vcd_scc_stimulusclass_trialdistr.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
    end
end % load params

return





