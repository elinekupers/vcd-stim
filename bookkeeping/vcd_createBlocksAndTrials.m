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
%  [session_type] : (optional) Are we dealing with MRI or BEHAVIOR sessions
%                     Default: 'MRI'
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
%   1: {'unique_im_nr'   } unique image nr"
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
%   18: {'is_cued'         } cue status: 0 = uncued, 1 = cued
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
p0.addParameter('session_type' , 'MRI', @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));

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
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('trials_%s_%s*.mat',p.disp.name,session_type)));
    fprintf('\n[%s]: Found %d trial .mat file(s)\n',mfilename,length(d));
    if ~isempty(d)
        if length(d) > 1
            warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
        end
        fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
        load(fullfile(d(end).folder,d(end).name));
    else
        error('[%s]: Can''t find trial params file!\n', mfilename)
    end
else % Recreate conditions and blocks and trials
    if ~isfield(p,'exp')
        warning('[%s]: params.exp doesn''t exist!)',mfilename)
        warning('[%s]: Will try to loading exp params from file.\n', mfilename);
        d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('exp_%s*.mat',p.disp.name)));
        if ~isempty(d)
            fprintf('\n[%s]: Found %d exp params .mat file(s)\n',mfilename,length(d));

            if length(d) > 1
                warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
            end
            
            fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
            tmp = load(fullfile(d(end).folder,d(end).name));
            p.exp = tmp.exp;
        else
            error('[%s]: Can''t find exp params file!\n', mfilename)
        end
    end
    
    % Fixation order and fixation
    fixsoafun = @() round(p.stim.fix.dotmeanchange);
    
    % Contrast decrement gaussian time window onset
    cdsoafun = @() round(p.stim.cd.meanchange + p.stim.cd.changeplusminus*(2*(rand-.5)));
    
    % Bookkeeping structs
    all_unique_im   = struct();
    all_cond          = struct();
    condition_master  = table();
    
    
    % Define block content for each stimulus class
    for stimClass_idx = 1:length(p.exp.stimclassnames)
        
        stim_table = table();
        
        % Get stimclass name
        stimClass_name = p.exp.stimclassnames{stimClass_idx};
        
        % Get corresponding task crossings
        bsc_idx = cell2mat(cellfun(@(x) strcmp(p.exp.stimclassnames,x), {stimClass_name}, 'UniformOutput', false));
        
        if strcmp(session_type,'MRI')
            task_crossings = find( p.exp.crossings(bsc_idx,:));
        elseif strcmp(session_type,'BEHAVIOR')
            task_crossings = find(p.exp.n_unique_trial_repeats_behavior(bsc_idx,:)>0);
        end
                
        % Loop over each task crossing for this stim class
        for curr_task = 1:length(task_crossings)

            % Define the unique images for Gabors
            [t_cond, n_unique_cases] = vcd_defineUniqueImages(p, stimClass_name);
            
            all_unique_im.(stimClass_name) = t_cond;
            
            % Add stimclass name and idx to temporary table
            t_cond.stim_class_name = repmat({stimClass_name}, size(t_cond,1),1);
            t_cond.stim_class      = repmat(stimClass_idx,    size(t_cond,1),1);
            
            % Add task name to temp table
            taskClass_name         = p.exp.taskclassnames{task_crossings(curr_task)};
            t_cond.task_class_name = repmat({taskClass_name},size(t_cond,1),1);
            t_cond.task_class      = repmat(task_crossings(curr_task),size(t_cond,1),1);
            
            % If IMG or LTM, we only want to use a subset of unique images
            if strcmp(taskClass_name,'ltm') || strcmp(taskClass_name,'img')
                t_cond = t_cond(t_cond.is_special_core,:);
            end
            
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
            tbl = vcd_createConditionMaster(p, t_cond, n_trials_per_block,session_type);
            
            % --- Allocate trials to blocks ---
            tbl.stim_class_unique_block_nr = NaN(size(tbl,1),1);
            tbl.trial_nr = NaN(size(tbl,1),1);
            
            % assign block nr
            block_start_idx = 1:n_trials_per_block:size(tbl,1);
            for mm = block_start_idx
                selected_trials = mm:(mm+(n_trials_per_block-1));
                if mm==block_start_idx(end)
                    selected_trials = selected_trials(selected_trials<=size(tbl,1));
                end
                block_nr = find(mm==block_start_idx);
                tbl.stim_class_unique_block_nr(selected_trials) = block_nr;
                tbl.trial_nr(selected_trials) = 1:length(selected_trials);
            end
            
            if strcmp(taskClass_name,'cd')
                if strcmp(stimClass_name, 'ns'), nsides = 1; else, nsides = 2; end
                % shuffle [1,2] such that there is a 50% chance on either side 
                % we will actually apply the contrast decrement change to
                % stimulus.  1=yes change, 2=no change
                cd_change = repmat((1:(1/p.stim.cd.prob)),2,(ceil(size(tbl,1)/2)*2)/(1/p.stim.cd.prob));
                cd_change = cd_change([randperm(length(cd_change),length(cd_change));randperm(length(cd_change),length(cd_change))]);
                cd_change = cd_change';
                % Get onset of contrast decrement within the
                % stimulus period and log in the table
                for stimrow = 1:size(tbl,1)
                    if ~tbl.is_catch(stimrow,:)
                        for side = 1:nsides
                            if cd_change(stimrow,side) == 1
                                tbl.cd_start(stimrow,side) = feval(cdsoafun);
                                
                            elseif cd_change(stimrow,side) == 2
                                tbl.cd_start(stimrow,side) = 0;
                            end
                        end
                    else
                        tbl.cd_start(stimrow,2) = NaN;
                    end
                end
                if nsides == 1, tbl.cd_start(:,2) = NaN; end

            else
                tbl.cd_start = NaN(size(tbl,1),2);
            end
            
            % Concatete single stim-task crossing table to master table
            stim_table = cat(1,stim_table,tbl);
            
        end % task class idx
        
        % Accumulate condition master and info
        all_cond.(stimClass_name) = stim_table;
        
        condition_master = cat(1,condition_master,stim_table);
    end % stim idx
    
    
    % ---- IMPORTANT FUNCTION: Shuffle stimuli for SCC and LTM task ----
    condition_master = vcd_shuffleStimForTaskClass(p,  'scc', condition_master, p.exp.block.n_trials_single_epoch,session_type);
    if strcmp(session_type,'MRI')
        condition_master = vcd_shuffleStimForTaskClass(p,  'ltm', condition_master, p.exp.block.n_trials_double_epoch,session_type);
    end
    
    % ---- IMPORTANT FUNCTION: Add correct button press ----
    button_response = NaN(size(condition_master,1),1);
    for ii = 1:size(condition_master,1)
        if  condition_master.is_catch(ii)
            button_response(ii) = NaN;
        else
            button_response(ii) = vcd_getCorrectButtonResponse(p, condition_master(ii,:));
        end
    end
    condition_master.correct_response = button_response;


    %% Store structs if requested
    if p.store_params
        fprintf('[%s]:Storing trial data..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        save(fullfile(saveDir,sprintf('trials_%s_%s_%s.mat',p.disp.name,session_type,datestr(now,30))),'condition_master','all_unique_im','all_cond')
    end
    
    
    %% At last, do some checks:
    
    % Stim class nr should match with params
    assert(isequal(unique(condition_master.stim_class)',[1:length(p.exp.stimclassnames),99]))
    % Task class nr should match with params
    if strcmp(session_type,'MRI')
        assert(isequal(unique(condition_master.task_class)',1:length(p.exp.taskclassnames)));
    elseif strcmp(session_type,'BEHAVIOR')
        assert(isequal(unique(condition_master.task_class)', find(any(p.exp.n_unique_trial_repeats_behavior,1))))
    end
    % Cued vs uncued stimuli should be matched
    assert(isequal(sum(condition_master.is_cued==1),sum(condition_master.is_cued==2)))
    
    % Now dive into each stimulus class
    unique_im_from_table = [];
    N_tbl = NaN(5,30); unique_stim_class_names_in_table = unique(condition_master.stim_class_name(:,1),'stable')';
    adjusted_crossings = p.exp.nr_unique_trials_per_crossing;
    for ii = 1:length(unique_stim_class_names_in_table)
        
        % Check nr of unique stimuli per stimulus class
        if ii == 5
            tmp = condition_master(~condition_master.is_catch & strcmp(condition_master.stim_class_name(:,1), unique_stim_class_names_in_table{ii}),:);
        
            [N, N_edges] = histcounts(tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), unique_stim_class_names_in_table{ii})));
        
            N_tbl(ii,1:length(N)) = N;
            assert(isequal(size(N,2), length(all_unique_im.(p.exp.stimclassnames{ii}).unique_im_nr)))

            unique_im_from_table(ii) = length(unique(tmp.stim_nr_left));
            
            if strcmp(session_type,'MRI')
                % lower IMG/LTM block contribution (we have less unique images)
                unique_trial_repeats = p.exp.n_unique_trial_repeats_mri;
            elseif strcmp(session_type,'BEHAVIOR') % NO LTM/IMG
                adjusted_crossings(ii,6) = 0;
                adjusted_crossings(ii,7) = 0;
                unique_trial_repeats = p.exp.n_unique_trial_repeats_behavior;
            end
            assert(isequal(size(tmp,1), sum(adjusted_crossings(ii,:) .* unique_trial_repeats(ii,:),'omitnan')))
            
        else
            tmp = condition_master(any(~condition_master.is_catch & strcmp(condition_master.stim_class_name, unique_stim_class_names_in_table{ii}),2),:);
            
            [N, N_edges] = histcounts(cat(1,tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), unique_stim_class_names_in_table{ii})),....
                                            tmp.stim_nr_right(strcmp(tmp.stim_class_name(:,2), unique_stim_class_names_in_table{ii}))));
            
            N_tbl(ii,1:length(N)) = N;
            assert(isequal(size(N,2), length(all_unique_im.(p.exp.stimclassnames{ii}).unique_im_nr)))
            colsLR = [tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), unique_stim_class_names_in_table{ii})), ...
                     tmp.stim_nr_right(strcmp(tmp.stim_class_name(:,2), unique_stim_class_names_in_table{ii}))];
                                            
            unique_im_from_table(ii) = length(unique(colsLR(:)));

            if strcmp(session_type,'MRI')
                adjusted_crossings(ii,6) = 2*adjusted_crossings(ii,6);
                adjusted_crossings(ii,7) = 2*adjusted_crossings(ii,7);
                unique_trial_repeats = p.exp.n_unique_trial_repeats_mri;
            elseif strcmp(session_type,'BEHAVIOR') % NO LTM/IMG
                adjusted_crossings(ii,6) = 0;
                adjusted_crossings(ii,7) = 0;
                unique_trial_repeats = p.exp.n_unique_trial_repeats_behavior;
            end
            assert(isequal(size(colsLR,1), sum(adjusted_crossings(ii,:) .* unique_trial_repeats(ii,:),'omitnan')));
        end
    end
    
    % Now dive into each task class
    crossings_unique_im = max(adjusted_crossings.*unique_im_from_table',[],2);
    
    M_tbl = NaN(10,112);
    for ii = unique(condition_master.task_class)'
        [M, edges_M] = histcounts([condition_master.stim_nr_left(~condition_master.is_catch & (condition_master.task_class==ii)); ...
                                    condition_master.stim_nr_left(~condition_master.is_catch & (condition_master.task_class==ii))],[0:1:111]);
        M_tbl(ii,1:length(M)) = M;
    end
    

    
    
    % Plot figures to check stimulus order
    if p.verbose
        
        %%
        vcd_visualizeMasterTable(condition_master, p.store_imgs,session_type);
        
        
        figure; set(gcf,'Position',[1,1,1200,300]);
        histogram([condition_master.stim_nr_left;condition_master.stim_nr_right],'numbins',111)
        xlabel('unique condition nr'); ylabel('trial count')
         title('Total unique image nr')
%          set(gca,'XTick',[1:30]);
        box off;
        if p.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_type));
            if ~exist(saveFigsFolder,'dir'); mkdir(saveFigsFolder); end
            filename = sprintf('vcd_totaluniqueim.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        figure; set(gcf,'Position',[1,1,1200,300]); imagesc(N_tbl)
        xlabel('unique im nr'); ylabel('Stim class');
        C = parula(20); colormap(C); set(gca,'CLim',[0 max(N_tbl(~isnan(N_tbl)))])
        cb = colorbar; cb.Label.String = 'unique im count';
        cb.Ticks = [0,unique(N_tbl(~isnan(N_tbl)))'];
        set(gca,'YTick',[1:5],'YTickLabel', p.exp.stimclassnames);  set(gca,'XTick',[1:30]); title('Sum of unique im per stimulus class')
        box off;
        if p.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_type));
            filename = sprintf('vcd_uniqueim_per_stimclass.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        
        figure; set(gcf,'Position',[1,1,1200,300]); imagesc(M_tbl)
        xlabel('unique condition nr'); ylabel('Task class');
        C = parula(20); colormap(C); set(gca,'CLim',[0 max(M_tbl(~isnan(M_tbl)))])
        cb = colorbar; cb.Label.String = 'unique im count';
        cb.Ticks = unique(M_tbl(~isnan(M_tbl)))';
        set(gca,'YTick',[1:10],'YTickLabel', p.exp.taskclassnames); 
        set(gca,'XTick',[0:25:110]); title('Sum of unique im per task class')
        box off;
        if p.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_type));
            filename = sprintf('vcd_uniqueim_per_taskclass.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        % Stim class distribution in SCC blocks
        foo = [condition_master.stim_nr_left((condition_master.task_class==3),:), ...
            condition_master.stim_nr_right((condition_master.task_class==3),:)]';
        figure; set(gcf,'Position',[1,1,1200,300])
        imagesc(foo)
        xlabel('Trial nr'); ylabel('left vs right stim');
        C = parula(max(foo(:))); colormap(C); set(gca,'CLim',[min(foo(:)),max(foo(:))])
        cb = colorbar; cb.Label.String = 'unique im nr';
        cb.Ticks = [min(foo(:)):10:max(foo(:))];
        set(gca,'YTick',[1:p.exp.block.n_trials_single_epoch]); title('SCC stimulus class distribution')
        box off;
        if p.store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_type));
            filename = sprintf('vcd_scc_stimulusclass_trialdistr.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        if strcmp(session_type,'MRI')
            % Stim class distribution in LTM blocks
            foo = [condition_master.stim_nr_left((condition_master.task_class==6),:), ...
                condition_master.stim_nr_right((condition_master.task_class==6),:)]';
            figure; set(gcf,'Position',[1,1,1200,300])
            imagesc(foo)
            xlabel('Trial nr'); ylabel('left vs right stim');
            C = parula(max(foo(:))); colormap(C); set(gca,'CLim',[min(foo(:)),max(foo(:))])
            cb = colorbar; cb.Label.String = 'unique im nr';
            cb.Ticks = [min(foo(:)):10:max(foo(:))];
            set(gca,'YTick',[1:p.exp.block.n_trials_single_epoch]); title('LTM stimulus class distribution')
            box off;
            if p.store_imgs
                saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',session_type));
                filename = sprintf('vcd_ltm_stimulusclass_trialdistr.png');
                print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
            end
        end
    end
end % load params

return





