function [params, condition_master, all_unique_im, all_cond] = vcd_createConditions(params,varargin)
% VCD function to define, shuffle, and organize trials into blocks:
%
%   [condition_master, all_unique_im, all_cond] = ...
%       vcd_createConditions(params,'load_params',<load_params>,'store_params',<store_params>)
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
%   5. Allocate trials to blocks.
%   6. Calculate with CD-trials a contrast decrement onset and when this
%      dip will happen during the trial.
%
% INPUTS:
%  params         : (required) Struct with stimulus and session params
%  [load_params]  : (optional) Load prior stored parameters or not. Will
%                     search for a file name starting with trials_<dispname>*.mat
%                     in fullfile(vcd_rootPath,'workspaces','info').
%                     Default: true
%  [store_params] : (optional) Store generated parameters or not.
%                     Will store params in fullfile(vcd_rootPath,'workspaces','info');
%                     Default: true
%  [env_type]     : (optional) Are we dealing with MRI or BEHAVIOR sessions
%                     Default: ''
%
% OUTPUT:
%  condition_master : table with N rows for every trial in the in VCD-core
%                      behavioralor MRI experiment, by M columns specifying
%                      information about the individual trial.
%                    Each row is an individual image.
%                    For Gabors/RDKs/Dots/ComplexObj, two subsequent rows
%                    define a unique trial. For NS, there is only 1 row per
%                    unique trial, as there is only one central image.
%
%   Columns description:
%  1:  'session_nr'       : (double, integral) Session number, behavioral
%                             experiment has 1  session. MRI experiment has
%                             27 sessions: 1 wide and 26 deep sessions.
%  2:  'session_type'     : (double, integral) Session type, two MRI
%                             sessions (the first and last) are split
%                             between subjects, such that half of the
%                             subjects see version A (session type = 1) and
%                             the other half of the subjects see version B
%                             (session type = 2).
%  3:  'run_nr'           : (double, integral) Run number. Resets every
%                             session behavioral experiment has 1-12 runs
%                             for itssingle session. MRI experiment has 10
%                             runs per session (except for the last
%                             session, which has 5 runs).
%  4:  'block_nr'         : (double, integral) Stimulus block number,
%                             resets every run. Behavioral experiment has
%                             1-12 runs for its single session. The MRI
%                             experiment has 10 runs per session, except
%                             for the last session, which has 5 runs.
%  5:  'trial_nr'         : (double, integral) Trial number within a
%                             stimulus block, resets every block. Single-
%                             stimulus presentation blocks have 8 trials
%                             per block (fix/cd/scc/pc/what/where/how
%                             tasks). Double-stimulus presentation blocks
%                             have 4 trials per block (wm/ltm/img tasks).
%  6:  'global_run_nr'    : (double, integral) Cumulative run number across
%                             all the sessions, i.e., this number does not
%                             reset for every session)
%  7:  'global_block_nr'  : (double, integral) Cumulative block number
%                             across all the sessions, i.e., this number
%                             does not reset for every run or session)
%  8:  'global_trial_nr'  : (double, integral) Cumulative trial number
%                             across all the sessions, i.e., this number
%                             does not reset for block, run or session)
%  9:  'condition_nr'     : (double, integral) Condition number for this
%                             particular trial. Every combination of
%                             a unique stimulus x task crossing x cueing
%                             status (cued/uncued/neutral cued) gets their
%                             own condition number. We use the same
%                             condition number regardless of whether there
%                             an additional test stimulus was shown (in
%                             double stimulus-presentation trials) or not.
%  10: 'condition_name'   : (cell, char) condition name, same as column
%                             9: 'condition_nr' but human-readible format.
%  11: 'stim_class'       : (double, integral) stimulus class number, where
%                            1: 'gabor' (Gabors: grating windowed by a 2D circular gaussian)
%                            2: 'rdk'   (Random dot motion kinematograms)
%                            3: 'dot'   (Single small dot on a 4-degree iso-eccentric ring)
%                            4: 'obj'   (Grayscale objects cropped from the background)
%                            5: 'ns'    (Large, centrally presented colorful naturalistic scenes)
%  12: 'stim_class_name'  : (cell, char) stimulus class name, same as
%                            column 11: 'stim_class' but human-readible
%                            format.
%  13: 'task_class'       : (double, integral) task class number. Ranges
%                            from 1-10, where 1: fix, 2: cd, 3: scc, 4: pc,
%                            5: wm, 6: ltm, 7: img, 8: what, 9: where,
%                            10: how.
%  14: 'task_class_name'  : (cell, char) task class name, same as column
%                            13: 'task_class' but human-readible format.
%  15: 'crossing_nr'      : (double, integral) stimulus-task class crossing
%                            number. Ranges from 1 through 32, where
%                            1: 'fix-gabor',  2: 'cd-gabor',  3: 'scc-all'
%                            4: 'pc-gabor' ,  5: 'wm-gabor',  6: 'ltm-all'
%                            7: 'img-gabor',  8: 'fix-rdk',   9: 'cd-rdk'
%                            10: 'pc-rdk' ,  11: 'wm-rdk' ,  12: 'img-rdk'
%                            13: 'fix-dot',  14: 'cd-dot' ,  15: 'pc-dot'
%                            16: 'wm-dot',   17: 'img-dot',  18: 'fix-obj'
%                            19: 'cd-obj' ,  20: 'pc-obj' ,  21: 'wm-obj'
%                            22: 'img-obj',  23: 'what-obj', 24: 'how-obj'
%                            25: 'fix-ns' ,  26: 'cd-ns',    27: 'pc-ns'
%                            28: 'wm-ns',    29: 'img-ns',   30: 'what-ns'
%                            31: 'where-ns', 32: 'how-ns'.
%                            'all' in 'scc-all' means that those stimulus
%                            blocks mix all stimulus classes crossed with
%                            the 'scc' task (i.e., gabors, rdks, dots,
%                            objects).
%                            'all' in 'ltm-all' means that those stimulus
%                            blocks mix all stimulus classes crossed with
%                            the 'ltm' task (i.e., gabors, rdks, dots,
%                            objects, scenes).
%  16: 'crossing_name'    : (cell, char) task class name, same as column
%                            15: 'crossing_nr' but human-readible.
%  17: 'stim_nr_left'     : (double, integral) unique stimulus number for
%                            the stimulus location on the left of the
%                            central fixation circle for gabors, rdks, dots,
%                            and objects, OR the unique stimulus number of
%                            the centrally presented scene.
%  18: 'stim_nr_right'    : (double, integral) unique stimulus number for
%                            the stimulus location on the right of the
%                            central fixation circle for gabors, rdks, dots,
%                            and objects. Scenes only have NaN in this column
%                            as we only show one centrally presented
%                            stimulus.
%  19: 'is_cued'          : (double, integral) what stimulus location is
%                            cued with a covert spatial attention cue by
%                            thickening the fixation circle rim where
%                            1 = left side, 2 = right side, 3 = both
%                            sides/neutral cue.
%  20: 'is_catch'         : (double, 0/1) whether the trial is
%                            a catch trial where no stimulus showed up.
%                            1 = yes, this is a catch trial. 0 = no,
%                            this is a regular, NON catch trial.
%  21: 'correct_response' : (double, integral) what is the correct response
%                            for the cued stimulus this given trial?
%                            Depending on the stimulus-task crossing, the
%                            correct response can be 2AFC (1 or 2),
%                            3AFC (1, 2, or 3), or 4AFC (1, 2, 3, or 4).
%                            See vcd_getCorrectButtonResponse.m for more
%                            details.
%  22: 'orient_dir'       : (double, decimal) Depending on the stimulus
%                            class, this column describes the Gabor tilt
%                            orientation (degrees), RDK motion direction
%                            (degrees), or Dot angular position (degrees),
%                            or object facing direction (degrees).
%  23: 'contrast'         : (double, decimal) stimulus contrast (Michelson fraction)
%  24: 'gbr_phase'        : (double, integral) gabor stimulus phase  (NaN for non gabor stim)
%  25: 'rdk_coherence'    : (double, decimal) rdk stimulus dot coherence,
%                            i.e. fraction of dots that move in the same
%                            direction (NaN for non rdk stim)
%  26: 'super_cat'        : (double, integral) OBJ and NS superordinate
%                            object category level: 1:'human', 2:'animal',
%                            3:'object', 4:'food', 5:'place'.
%                            (NaN for gabor, rdk and dot stimulus classes).
%  27: 'super_cat_name'   : (cell, char) same as colum 26: 'super_cat', but
%                            human-readible format.
%  28: 'basic_cat'        : (double, integral) OBJ and NS basic object
%                            category level for every superordinate
%                            category level:
%                            OBJECT BASIC CATEGORY LABELS:
%                            * human:  1:'facemale', 2:'facefemale'
%                            * animal: 1:'small',    2:'large'
%                            * object: 1:'tool',     2:'vehicle'
%                            * food:   1:'manmade',  2:'produce'
%                            * place:  1:'building', (no 2nd label)
%                            SCENE BASIC CATEGORY LABELS:
%                            * human:  1:'indoor',  2:'outdoor'
%                            * animal: 1:'indoor',  2:'outdoor'
%                            * object: 1:'indoor',  2:'outdoor'
%                            * food:   1:'indoor',  2:'outdoor'
%                            * place:  1:'indoor',  2:'outdoor'
%                            (NaN for gabor, rdk and dot stimulus classes).
%  29: 'basic_cat_name'   : (cell, char) same as colum 28: 'basic_cat', but
%                            human-readible format.
%  30: 'sub_cat'          : (double, integral) OBJ and NS subordinate
%                            object category level for every basic category
%                            level:
%                            OBJECT SUBORDINATE CATEGORY LABELS:
%                            * facemale:   1:'damon'
%                            * facefemale: 1:'lisa',   2:'sophia'
%                            * small:      1:'parrot', 2:'cat'
%                            * large:      1:'bear',   2:'giraffe'
%                            * tool:       1:'drill',  2:'brush'
%                            * vehicle:    1:'bus',    2:'suv'
%                            * manmade:    1:'pizza',  2:'suv'
%                            * produce:    1:'bus',    2:'banana'
%                            * building:
%                            SCENE SUBORDINATE CATEGORY LABELS:
%                            * indoor:  1:'left', 2:'center', 3:'right
%                            * outdoor: 1:'left', 2:'center', 3:'right
%                            (NaN for gabor, rdk and dot stimulus classes).
%  31: 'sub_cat_name'     : (cell, char) same as colum 30: 'sub_cat', but
%                            human-readible format.
%  32: 'affordance_cat'   : (double, integral) OBJ and NS affordance
%                            category labels
%                            OBJECT AFFORDANCE LABELS:
%                            1:'greet', 2:'grasp',
%                            3:'enter', 4: 'observe/do nothing'
%                            SCENE SUBORDINATE CATEGORY LABELS:
%                            1:'greet', 2:'grasp',
%                            3:'enter', 4: 'observe/do nothing'
%                            (NaN for gabor, rdk and dot stimulus classes).
%  33: 'affordance_name'  : (cell, char) same as 32: 'affordance_cat' but
%                            human-readible format.
%  34: 'cd_start'         : (double, integral) for cd task crossings only
%                            (other task crossings will have NaN), the time
%                            frame relative to stimulus onset of that trial
%                            when the contrast decrement started. Only 20%
%                            of cd trials will contain a contrast
%                            decrement. If no contrast decrement occurred
%                            in a cd trial, then 'cd_start' is set to NaN.
%  35: 'stim2_delta'      : (double, decimal) for double-stimulus
%                            presentation trials only, depending on the
%                            stimulus class, this column describes the
%                            difference in: Gabor tilt orientation
%                            (degrees), RDK motion direction (degrees)
%                            Dot angular position (degrees), object facing
%                            direction (degrees) between the first stimulus
%                            and the second stimulus shown after the delay.
%                            Values are defined by the stimulus parameter:
%                            params.stim.(<stim_class_name>).delta_ref
%  36: 'stim2_im_nr'      : (double, integral) for double-stimulus
%                            presentation trials only, what is the unique
%                            stimulus number that is associated with the
%                            test stimulus (the second stimulus show after
%                            the 8-s delay) in WM, LTM, IMG task crossings.
%  37: 'stim2_orient_dir' : (double, decimal) for double-stimulus
%                            presentation trials only, depending on the
%                            stimulus class, this column describes the
%                            total Gabor tilt orientation
%                            (degrees), RDK motion direction (degrees)
%                            Dot angular position (degrees), object facing
%                            direction (degrees) between the first stimulus
%                            and the second stimulus shown after the delay.
%  38: 'is_special_core'  : (double, 0/1) Is this stimulus part of
%                            the special core stimulus set (1) or not
%                            (0)? Special core stimuli are a subset of
%                            core stimuli that are only used for LTM and
%                            IMG task crossings.
%  39: 'is_lure'          : (double, 0/1) For LTM task crossings:
%                            did we use a novel lure stimulus (1)
%                            or not (0).
%  40: 'repeat_nr'        : (double, integral): How many times has this exact trial been
%                            repeated thusfar in the experiment, across all
%                            the sessions and runs.
%  41: 'trial_type'       : (double, integral): Is this a single-stimulus
%                            presentation trial (1) or a double-stimulus
%                            presentation trial (2)?
%
%  There are also two hierarchical structs, which are the building blocks
%  of the condition master table, listing the unique images and unique
%  conditions for each stimulus class:
%  * all_unique_im  : matrix labeling the unique image features
%  * all_conds      : stimulus "condition master" matrix describing the
%                     pseudo-random shuffled order of trials for one
%                     stimulus class, across all sessions.
%
%
% Written by Eline Kupers November 2024 (kupers [at] umn [dot] edu)


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p0 = inputParser;
p0.addRequired('params'        , @isstruct);
p0.addParameter('load_params'  , false, @islogical);
p0.addParameter('store_params' , true,  @islogical);
p0.addParameter('store_imgs'   , false, @islogical);
p0.addParameter('verbose'      , false, @islogical);
p0.addParameter('env_type'     , '', @(x) any(strcmp(x,{'BEHAVIOR','MRI'})));

% Parse inputs
p0.parse(params,varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p0.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p0.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p0


%% %%%%%%%%%%%%% BOOKKEEPING  %%%%%%%%%%%%%

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

% Load params if requested and we can find the file
if load_params
    if ~isfield(params,'is_demo'), params.is_demo = false; end
    if ~isfield(params,'is_wide'), params.is_wide = false; end
    if strcmp(env_type, 'MRI')
        fname = sprintf('condition_master_%s%s%s*.mat', choose(params.is_wide,'wide_','deep_'), choose(params.is_demo,'demo_',''), params.disp.name);
    else
        fname = sprintf('condition_master_%s%s*.mat', choose(params.is_demo,'demo_',''), params.disp.name);
    end
    d = dir(fullfile(vcd_rootPath,'workspaces','info',fname));
    fprintf('\n[%s]: Found %d condition_master .mat file(s)\n',mfilename,length(d));
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
    if ~isfield(params,'exp')
        warning('[%s]: params.exp doesn''t exist!)',mfilename)
        warning('[%s]: Will try to loading exp params from file.\n', mfilename);
        d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('exp_%s*.mat',params.disp.name)));
        if ~isempty(d)
            if verbose
                fprintf('\n[%s]: Found %d exp params .mat file(s)\n',mfilename,length(d));
                if length(d) > 1
                    warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
                end
                fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
            end
            tmp = load(fullfile(d(end).folder,d(end).name));
            params.exp = tmp.exp;
        else
            if verbose
                warning('[%s]: Can''t find exp session .mat files! Will run vcd_getSessionParams.m', mfilename);
            end
            params.exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
        end
    end
    if ~isfield(params,'stim')
        if verbose
            warning('[%s]: params.stim doesn''t exist!)',mfilename)
            warning('[%s]: Will try to loading stim params from file.\n', mfilename);
        end
        d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('stim_%s*.mat',params.disp.name)));
        if ~isempty(d)
            if verbose
                fprintf('\n[%s]: Found %d stim params .mat file(s)\n',mfilename,length(d));
                if length(d) > 1
                    warning('[%s]: Multiple .mat files! Will pick the most recent one\n', mfilename);
                end
                fprintf('[%s]: Loading stim params .mat file: %s\n', mfilename, d(end).name);
            end
            tmp = load(fullfile(d(end).folder,d(end).name));
            params.exp = tmp.exp;
        else
            if verbose
                warning('[%s]: Can''t find stim session .mat files! Will run vcd_getStimParams.m', mfilename);
            end
            params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);
        end
    end
    
    %% %%%%%%%%%%%%% CREATE CONDITIONS %%%%%%%%%%%%%
    % condition_master table with all conditions shown across all sessions.
    % Note that block_nr counts continuously because we will shuffle the
    % block order later for each subject. Same holds for run_nr, these
    % numbers do not represent the final nr of runs as they will change
    % depending on the allocation of blocks to each run within a session.
    
    % Preallocate space and set up tables/structs
    if verbose
        tic
        fprintf('\n[%s]: Start creating conditions for %s experiment.. \n',mfilename,env_type);
    end
    
    % assume we don't run a demo run if parameter isn't specified;
    if ~isfield(params,'is_demo'), params.is_demo = false; end
    if ~isfield(params,'is_wide'), params.is_wide = false; end
    all_unique_im     = struct();  % details about unique core images present in VCD-core experiment. This information is also present in condition_master.
    all_cond          = struct();  % details about unique conditions (unique images x cueing condition) present in VCD-core experiment. This information is also present in condition_master.
    
    
    % Define block content for each stimulus class
    for stimClass_idx = 1:length(params.exp.stimclassnames)
        
        clear stim_table
        
        % Get stimclass name
        stimClass_name = params.exp.stimclassnames{stimClass_idx};
        
        % Get corresponding task crossings
        bsc_idx = cell2mat(cellfun(@(x) strcmp(params.exp.stimclassnames,x), {stimClass_name}, 'UniformOutput', false));
        
        if strcmp(env_type,'MRI')
            if params.is_wide
                task_crossings = find( params.exp.n_unique_trial_repeats_wide(bsc_idx,:)>0);
            else
                task_crossings = find( params.exp.n_unique_trial_repeats_deep(bsc_idx,:)>0);
            end
        elseif strcmp(env_type,'BEHAVIOR')
            if params.is_demo
                task_crossings = find(params.exp.n_unique_trial_repeats_demo(bsc_idx,:)>0);
            else
                task_crossings = find(params.exp.n_unique_trial_repeats_behavior(bsc_idx,:)>0);
            end
        end
        
        % Loop over each task crossing for this stim class
        for curr_task = 1:length(task_crossings)
            
            % Define the unique images for given stimulus class
            [t_cond, n_unique_cases] = vcd_defineUniqueImages(params, stimClass_name);
            
            all_unique_im.(stimClass_name) = t_cond;
            
            % Add stimclass name and idx to temporary table
            t_cond.stim_class_name = repmat({stimClass_name}, size(t_cond,1),1);
            t_cond.stim_class      = repmat(stimClass_idx,    size(t_cond,1),1);
            
            % Add task name to temp table
            taskClass_name         = params.exp.taskclassnames{task_crossings(curr_task)};
            t_cond.task_class_name = repmat({taskClass_name},size(t_cond,1),1);
            t_cond.task_class      = repmat(task_crossings(curr_task),size(t_cond,1),1);
            
            % If IMG or LTM, we only want to use a subset of unique images
            if strcmp(taskClass_name,'ltm') || strcmp(taskClass_name,'img')
                t_cond = t_cond(t_cond.is_special_core==1,:);
            end
            
            % Get the number of trials and blocks, which depend on task
            if params.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                n_trials_per_block    = params.exp.block.n_trials_single_epoch;
            else
                n_trials_per_block    = params.exp.block.n_trials_double_epoch;
            end
            
            %% ---- IMPORTANT FUNCTION: Create condition master table ---- %%
            % Create condition master table, where unique images are
            % shuffled and distributed across blocks (subject to certain
            % constraints)
            tbl = vcd_createConditionMaster(params, t_cond, env_type);
            
            if ~any(ismember(tbl.Properties.VariableNames,'unique_trial_nr')) % check for the first column in the table to make sure we actually generated conditions
                error('[%s]: No condition table was made! Hint: Check params.exp.n_unique_trial_repeats and/or inputs/outputs of vcd_createConditionMaster', mfilename)
            end
            % --- Allocate space for block nr ---
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
            
            % set Nans for cd onset (we will deal with it later)
            tbl.cd_start = NaN(size(tbl,1),1);
            
            % Concatenate single stim-task crossing table to master table
            if ~exist('stim_table','var')
                stim_table = tbl([],:);
            end
            stim_table = cat(1,stim_table,tbl);
            
        end % task class idx
        
        % Accumulate condition master and info
        all_cond.(stimClass_name) = stim_table;
        
        if ~exist('condition_master','var')
            condition_master = stim_table([],:); % create empty table with the same columns
        end
        condition_master = cat(1,condition_master,stim_table);
    end % stim idx
    
    % ---- Add object catch trials ----
    % How many trials are object catch given the 20% probility?
    pcobj_trial_idx = find(condition_master.stim_class == 4 & condition_master.task_class == 4 & condition_master.is_catch==0);
    if ~isempty(pcobj_trial_idx)
        nr_objectcatch_trials_per_rep = round(size(condition_master(pcobj_trial_idx,:),1)*params.exp.trial.pc.prob_objcatch);
        
        % Select object catch trials randomly from the total list (without
        % replacement), using the following constraints:
        % 1. Ensure we distribute object catch trials evenly across left/right cued locations
        % 2. Ensure we distribute object catch trials evenly across unique objects (to the extent possible)
        % 3. Ensure we distribute object catch trials evenly across possible correct button presses.
        % 4. Ensure we distribute object catch trials evenly across blocks
        while 1
            objcatch_ok = false(1,4); % 4 constraints
            
            objcatch_idx = datasample(pcobj_trial_idx,nr_objectcatch_trials_per_rep,'Replace',false);
            
            % Check the cued stimulus in those selected trials
            cued_loc_objcatch = condition_master.is_cued(objcatch_idx);
            if (abs(diff(histcounts(cued_loc_objcatch, [1:3]))) <= 1)
                objcatch_ok(1) = true;
            end
            
            % Check the object stimulus in those selected trials
            obj_objcatch = condition_master.sub_cat_name(objcatch_idx);
            if length(unique(obj_objcatch))==length(obj_objcatch) || ...
                    length(unique(obj_objcatch))==length(unique(condition_master.sub_cat_name(pcobj_trial_idx)))
                objcatch_ok(2) = true;
            end
            
            % Check the button response in those selected trials
            for nn = 1:size(condition_master(objcatch_idx,:),1)
                button_objcatch(nn) = vcd_getCorrectButtonResponse(params, condition_master(objcatch_idx(nn),:));
            end
            if abs(diff(histcounts(button_objcatch,[1:3])))<=1
                objcatch_ok(3) = true;
            end
            
            % Check the block in those selected trials
            objcblock_nr = unique(condition_master.stim_class_unique_block_nr(objcatch_idx));
            possible_block_nrs = unique(condition_master.stim_class_unique_block_nr(pcobj_trial_idx));
            if length(objcatch_idx) >= length(possible_block_nrs) % if we have more obj catch trials than obj-pc blocks
                if length(objcblock_nr)/length(possible_block_nrs) > (2/3) % if 2/3 object catch trials are distributed across blocks, we are happy..
                    objcatch_ok(4) = true;
                end
            elseif length(objcatch_idx) < length(possible_block_nrs) % if we have fewer obj catch trials than obj-pc blocks
                if length(objcblock_nr)/length(possible_block_nrs) > (1/2) % if 1/2 object catch trials are distributed across blocks, we are happy..
                    objcatch_ok(4) = true;
                end
            else
                error('wtf')
            end
            
            if sum(objcatch_ok)==length(objcatch_ok)
                break
            end
        end
        
        
        % Get the unique image numbers for object catch stimuli
        % (reshape to 16 objects x 18 catch rotations)
        catch_im_nr = reshape(params.stim.obj.unique_im_nrs_objcatch',[],params.stim.obj.num_unique_objects)';
        
        % Copy catch rotations such that we can remove them when used in a trial.
        possible_objcatch_rotations = params.stim.obj.catch_rotation; % dims: 16 x 18
        
        % Loop over object catch trials
        for cc = 1:length(objcatch_idx)
            % is this a right or left cued stimulus location?
            if cued_loc_objcatch(cc)==1 % if we are dealing with a left cued stimulus
                old_stim_nr       = condition_master.stim_nr_left(objcatch_idx(cc));
            elseif cued_loc_objcatch(cc)==2 % if we are dealing with a right cued stimulus
                old_stim_nr       = condition_master.stim_nr_right(objcatch_idx(cc));
            else
                error('[%s]: Stimulus location must be 1 or 2 for defining objectcatch!',mfilename)
            end
            
            % What core object (1-16) are we dealing with?
            old_stim_obj_nr   = (old_stim_nr==params.stim.obj.unique_im_nrs_core);
            
            % What are the possible object catch rotations we can use?
            possible_objcatch_rotations0 = possible_objcatch_rotations(old_stim_obj_nr,:);
            possible_objcatch_rotations1 = possible_objcatch_rotations0(~isnan(possible_objcatch_rotations0));
            
            % randomly select object catch rotation
            objcatch_rot = datasample(possible_objcatch_rotations1,1);
            
            % update rotation of object
            condition_master.orient_dir(objcatch_idx(cc),cued_loc_objcatch(cc)) = objcatch_rot;
            
            % remove rotation from list
            remove_me = ismember(possible_objcatch_rotations0,objcatch_rot);
            assert(isequal(possible_objcatch_rotations(old_stim_obj_nr,remove_me),condition_master.orient_dir(objcatch_idx(cc),cued_loc_objcatch(cc))))
            possible_objcatch_rotations(old_stim_obj_nr,remove_me) = NaN;
            
            % Update special core column (objectcatch stimuli can never be
            % special core stimuli)
            condition_master.is_special_core(objcatch_idx(cc),cued_loc_objcatch(cc)) = 0;
            
            % Insert object catch image number into "is_objectcatch" for now..
            % (we need to hold on to the core object number for the condition
            % label).
            condition_master.is_objectcatch(objcatch_idx(cc)) = catch_im_nr(old_stim_obj_nr,remove_me);
        end
        condition_master.is_objectcatch(setdiff(pcobj_trial_idx,objcatch_idx)) = 0;
    end
    
    % ---- IMPORTANT STEP: Shuffle stimuli for SCC and LTM task ----
    condition_master = vcd_shuffleStimForTaskClass(params,  'scc', condition_master, params.exp.block.n_trials_single_epoch,env_type);
    if strcmp(env_type,'MRI') && ~params.is_wide
        condition_master = vcd_shuffleStimForTaskClass(params,  'ltm', condition_master, params.exp.block.n_trials_double_epoch,env_type);
    end
    
    % ---- IMPORTANT STEP: Add correct button press ----
    button_response = NaN(size(condition_master,1),1);
    for ii = 1:size(condition_master,1)
        if condition_master.is_catch(ii)==1
            button_response(ii) = NaN;
        else
            button_response(ii) = vcd_getCorrectButtonResponse(params, condition_master(ii,:)); % note: PC-OBJ button press relies on orient_dir, which has been updated for objectcatch trials above.
        end
    end
    condition_master.correct_response = button_response;
    
    % ---- bookkeeping: See if we can achieve roughly equally sample scc images within a block
    
    % Check how many SCC we will actually allocate across sessions
    [all_sessions,session_types] = vcd_getSessionEnvironmentParams(params, env_type);
    
    scc_trials = find(condition_master.stim_class==99 & condition_master.task_class==3);
    scc_cued0 = condition_master.is_cued(scc_trials);
    
    nr_scc_blocks = sum(sum(all_sessions(:,3,:,:),3),4);
    nr_scc_blocks = nr_scc_blocks(nr_scc_blocks>0);
    if all(nr_scc_blocks<1)
        nr_scc_blocks = length(ceil(nr_scc_blocks));
    else
        nr_scc_blocks = sum(nr_scc_blocks);
    end
    
    scc_catch_trials    = condition_master.is_catch(scc_trials);
    nr_scc_catch_trials = sum(scc_catch_trials);
    if strcmp(env_type,'BEHAVIOR')
        bb_bpress = reshape(condition_master.correct_response(scc_trials(~scc_catch_trials))',nr_scc_blocks,[]); % (nr blocks x trials )
        expected_scc_trials = (params.exp.n_unique_trial_repeats_behavior(:,3).* params.exp.nr_unique_trials_per_crossing(:,3));
    elseif strcmp(env_type,'MRI')
        bb_bpress     = reshape(condition_master.correct_response(scc_trials(~scc_catch_trials))',[],params.exp.block.n_trials_single_epoch); % (nr blocks x trials )
        if params.is_wide
            expected_scc_trials = (params.exp.n_unique_trial_repeats_wide(:,3).* params.exp.nr_unique_trials_per_crossing(:,3));
            expected_catch_trials = sum(params.exp.nr_catch_trials_wide(:,3),'omitnan');
        else
            expected_scc_trials   = params.exp.n_unique_trial_repeats_deep(:,3).*params.exp.nr_unique_trials_per_crossing(:,3);
            expected_catch_trials = sum(params.exp.nr_catch_trials_deep(:,3), 'omitnan');
        end
    end
    bb_cued  = reshape(scc_cued0(~scc_catch_trials)', params.exp.block.n_trials_single_epoch, []); % (trials x nr blocks)
    n0 = histcounts(bb_bpress,[1:5]);  % 1=gabor, 2=dot, 3=rdk, 4=obj
    m0 = histcounts(bb_cued,[1:3]); % 1=left, 2=right, 3=neutral
    assert(isequal(m0(1),m0(2))); % assume equal left/right spatial cues
    expected_scc_trials = expected_scc_trials([1,3,2,4])'; % remove nan for NS; swap single dot and RDK stim class to match response order for subject
    assert(isequal(n0, expected_scc_trials))
    
    % Check unique trial nr
    assert(isequal(condition_master.unique_trial_nr(scc_trials),[1:length(condition_master.unique_trial_nr(scc_trials))]'));
    
    new_trial_nr = repmat([1:params.exp.block.n_trials_single_epoch]',ceil(length(condition_master.unique_trial_nr(scc_trials))/params.exp.block.n_trials_single_epoch),1);
    new_trial_nr = new_trial_nr(1:length(scc_trials));
    condition_master.trial_nr(scc_trials) = new_trial_nr;
    condition_master.stim_class_unique_block_nr(scc_trials) = new_trial_nr;

    %% ---- bookkeeping: Check if button presses are balanced to the extent possible
    condition_master = vcd_balanceButtonCorrectPresses(params, condition_master, env_type);

    %% ---- IMPORTANT FUNCTION: Allocate trials to blocks all unique trials and repeats of trials.
    % !!WARNING!! There is a randomization component involved in creating the
    % conditions (i.e., order of trials allocated to a block). If you don't
    % want this, set params.load_params = true to load an existing
    % condition_master.
    condition_master = vcd_allocateBlocksToRuns(params,condition_master,env_type);

    %% ---- IMPORTANT FUNCTION: Add contrast decrement
    % !!WARNING!! There is a randomization component involved in
    % determining the contrast decrement component!!
    condition_master = vcd_determineContrastDecrementChangeTrials(params, condition_master);

    %% -- clean up condition master catch trials --
    % set stimulus class name to NaN
    condition_master.stim_class_name(condition_master.is_catch==1,:) = repmat({NaN,NaN}, sum(condition_master.is_catch==1),1);
    
    %% Store condition_master if requested
    if store_params
        fprintf('[%s]: Storing condition_master..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        if strcmp(env_type,'MRI')
            fname = sprintf('condition_master_%s%s%s_%s.mat',choose(params.is_wide,'wide_','deep2_'), choose(params.is_demo,'demo_',''),params.disp.name,datestr(now,30));
        else
            fname = sprintf('condition_master_%s%s_%s.mat',choose(params.is_demo,'demo_',''),params.disp.name,datestr(now,30));
        end
        save(fullfile(saveDir,fname),'condition_master','all_unique_im','all_cond')
    end
    
    
    %% At last, do some checks:
    
    % Stim class nr should match with params
    assert(isequal(unique(condition_master.stim_class)',[1:length(params.exp.stimclassnames),99]))
    % Task class nr should match with defined stimulus-task crossings in params.exp
    if strcmp(env_type,'MRI')
        if params.is_wide
            assert(isequal(unique(condition_master.task_class)',find(nansum(params.exp.n_unique_trial_repeats_wide,1)>1)));
        else
            assert(isequal(unique(condition_master.task_class)',1:length(params.exp.taskclassnames)));
            
            % load in wide session to check A/B combined sessions
            d = dir(fullfile(vcd_rootPath,'workspaces','info','condition_master_wide_7TAS_BOLDSCREEN32_*.mat'));
            a = load(fullfile(d(end).folder,d(end).name),'condition_master');
            cm_wide = a.condition_master; clear a;
        end
    elseif strcmp(env_type,'BEHAVIOR')
        if params.is_demo
            assert(isequal(unique(condition_master.task_class)', find(any(params.exp.n_unique_trial_repeats_demo,1))))
        else
            assert(isequal(unique(condition_master.task_class)', find(any(params.exp.n_unique_trial_repeats_behavior,1))))
        end
    end
    
    % Cued vs uncued stimuli should be matched for non-demo sessions
    if ~params.is_demo
        %         if mod(sum(condition_master.is_catch == 1 & ismember(condition_master.is_cued,[1,2])),2)==0
        %             assert(isequal(sum(condition_master.is_cued==1 & condition_master.is_catch==0),sum(condition_master.is_cued==2 & condition_master.is_catch==0)))
        %         else
        if params.is_wide
            assert(ismember(sum(condition_master.is_cued==1),sum(condition_master.is_cued==2)+[-2:2]))
            assert(ismember(sum(condition_master.is_cued==2),sum(condition_master.is_cued==1)+[-2:2]))
        else
            assert(ismember(sum(condition_master.is_cued==1),sum(condition_master.is_cued==2)+[-12:12])) % ltm slop
            assert(ismember(sum(condition_master.is_cued==2),sum(condition_master.is_cued==1)+[-12:12])) % ltm slop
        end
        %         end
    end
    
    % Now dive into each stimulus class
    unique_im_from_table = [];
    N_tbl = NaN(5,30);
    adjusted_crossings = params.exp.nr_unique_trials_per_crossing;
    
    all_sessions = vcd_getSessionEnvironmentParams(params, env_type);
    scc_probability = (params.exp.nr_unique_trials_per_crossing(1:4,3)./8)';
    ltm_probability = (params.exp.nr_unique_trials_per_crossing(:,6)./4./2)';
    % What difference between predicted and observed scc trials with a
    % given stimulus class do we allow?
    if strcmp(env_type,'BEHAVIOR')
        scc_tolerance = 1*2;  % max 1 per block // 2 blocks total
    else
        if params.is_wide
            scc_tolerance = 2*3; % max 2 per block // 3 blocks total
        end
    end
    
    for st = [1,2]
        
        for ii = 1:length(params.exp.stimclassnames)
            
            if strcmp(env_type,'MRI') && ~params.is_wide
                scc_tolerance = 2*squeeze(sum(all_sessions(ii,3,:,:))); % max 2 per block session type A / 0 per block session type B // A: 104 / B: 0 blocks total
                ltm_tolerance = 2*squeeze(sum(all_sessions(ii,6,:,:))); % max 1 per block sessoin type A / 2 per block session type B // A: 729 / B: 19 blocks total
            end

            nr_ltm_blocks = sum(sum(squeeze(all_sessions(:,6,:,st))));
            assert(isequal(length(unique(condition_master.global_block_nr(condition_master.session_type==st & condition_master.task_class==6,:))),nr_ltm_blocks))
            
            % Check nr of unique stimuli per stimulus class
            if ii == 5 % NS is dealt with separately because it only has a single central stim per trial
                tmp = condition_master(condition_master.is_catch==0 & condition_master.session_type==st & isnan(condition_master.is_objectcatch) & strcmp(condition_master.stim_class_name(:,1), params.exp.stimclassnames{ii}),:);
                tmp_catch = condition_master(condition_master.is_catch==1 & condition_master.session_type==st & isnan(condition_master.is_objectcatch) & strcmp(condition_master.stim_class_name(:,1), params.exp.stimclassnames{ii}),:);
                tmp_objectcatch = condition_master(condition_master.is_catch==0 & condition_master.session_type==st & ~isnan(condition_master.is_objectcatch) & strcmp(condition_master.stim_class_name(:,1), params.exp.stimclassnames{ii}),:);
                
                [N, ~] = histcounts(tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})));
                N_catch = sum(tmp_catch.is_catch==1);
                N_tbl(ii,1:length(N)) = N;
                if ~params.is_demo
                    if st == 1
                        assert(isequal(size(N,2), length(all_unique_im.(params.exp.stimclassnames{ii}).unique_im_nr)))
                    elseif strcmp(env_type,'MRI') && params.is_wide == 0 && st == 2 % Deep has NS-LTM
                        cm_wideB = cm_wide(cm_wide.session_type==st,:);
                        [N_wideB, ~] = histcounts(cm_wideB.stim_nr_left(strcmp(cm_wideB.stim_class_name(:,1), params.exp.stimclassnames{ii})),....
                            [params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_core, max(params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_core)+1]);
                        N_B = N + N_wideB;
                        assert(isequal(sum(N_B>0), length(all_unique_im.(params.exp.stimclassnames{ii}).unique_im_nr)));
                        
                        colsL_wideB      = cm_wideB.stim_nr_left(~ismember(cm_wideB.task_class,6) & isnan(cm_wideB.is_objectcatch) & strcmp(cm_wideB.stim_class_name(:,1), params.exp.stimclassnames{ii})); % exclude SCC, we deal with those trials separately
                    end
                end
                
                if strcmp(env_type,'MRI') && params.is_wide == 0 && st == 2
                    colsL = tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); 
                    
                    if st == 2
                        unique_im_from_table(ii) = length(unique([colsL;colsL_wideB])); % this should be 24 for GBR or RDK, 16 for DOT/OBJ
                        assert(isequal(unique_im_from_table(ii),params.exp.nr_unique_trials_per_crossing(ii,2)))
                    else
                        unique_im_from_table(ii) = length(unique(tmp.stim_nr_left)); % this should be 30 for NS
                        assert(isequal(unique_im_from_table(ii),params.exp.nr_unique_trials_per_crossing(ii,2)))
                    end
                end
                
                % Calculate the expected nr of trials from the params.exp.mri.ses_blocks table.
                expected_nr_of_trials = sum(squeeze(sum(all_sessions(ii,:,:,st).*params.exp.nr_trials_per_block(ii,:))));
                if params.is_demo
                    expected_nr_of_trials = 0.5*expected_nr_of_trials; % demo has only half the nr of trials per block.
                end
                empirical_nr_of_trials = size(tmp,1);
                if st == 2
                    N_catch = expected_nr_of_trials*0.2;
                    if ~(abs(floor(empirical_nr_of_trials)-expected_nr_of_trials)<= N_catch)
                        warning('[%s]: Mismatch between expected and empirical nr of NS blocks. YOU SHOULD DOUBLE CHECK!! that this is due to LTM mix of stimulus classes.', mfilename)
                    end
                else
                    assert(isequal(empirical_nr_of_trials, expected_nr_of_trials)) % nr of expected trials should match the empirical nr of trials
                end
                
            else % classic stimuli: GBR/RDK/DOT/OBJ
                if strcmp(env_type,'BEHAVIOR')
                    if ismember(ii,[1:3]) % if GBR/RDK/DOT
                        tmp = condition_master(condition_master.is_catch==0 & isnan(condition_master.is_objectcatch) & any(strcmp(condition_master.stim_class_name, params.exp.stimclassnames{ii}),2),:);
                    elseif ii == 4 % OBJ
                        tmp = condition_master(condition_master.is_catch==0 & condition_master.is_objectcatch==0 & any(strcmp(condition_master.stim_class_name, params.exp.stimclassnames{ii}),2),:);
                    end
                    % Get counts of stimulus numbers for left and right stimulus position.
                    [N, ~] = histcounts(cat(1,tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})),....
                        tmp.stim_nr_right(strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))));
                    N_tbl(ii,1:length(N)) = N;
                    if ~params.is_demo
                        assert(isequal(size(N,2), length(all_unique_im.(params.exp.stimclassnames{ii}).unique_im_nr)))
                    end
                    colsL = tmp.stim_nr_left(tmp.task_class ~=3 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); % exclude SCC, we deal with those trials separately
                    colsR = tmp.stim_nr_right(tmp.task_class ~=3 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii})); % exclude SCC, we deal with those trials separately
                    colsLR = [colsL(:),colsR(:)]; % if this fails, then we don't have equal left and right stimulus numbers for classic stimuli
                    unique_im_from_table(ii) = length(unique(colsLR)); % this should be 24 for GBR or RDK, 16 for DOT/OBJ
                    
                    colsLR_scc = [tmp.stim_nr_left(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); ... % deal with scc trials
                        tmp.stim_nr_right(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))];
                    
                    % Calculate the expected nr of trials from the params.exp.behavior.ses_blocks table.
                    expected_nr_of_trials = sum(sum(all_sessions(ii,:,:,:).*params.exp.nr_trials_per_block(ii,:)));
                    
                    if params.is_demo
                        expected_nr_of_trials = expected_nr_of_trials*0.5; % again: demo has half the trials per block
                    end
                    
                    if strcmp(params.exp.stimclassnames{ii},'obj')
                        empirical_nr_of_trials = ceil(sum([size(colsLR,1),size(unique(colsLR_scc),1)/2])+sum(condition_master.is_objectcatch==1)); % for OBJ: sum regular trials, scc trials and objectcatch trials
                    else
                        empirical_nr_of_trials = sum([size(colsLR,1),size(unique(colsLR_scc),1)/2]);  % for GBR/RDK/DOT: sum regular trials, scc trials
                    end
                    
                    assert((abs(floor(empirical_nr_of_trials)-expected_nr_of_trials))<=scc_tolerance); % we allow for 1-2 trials difference per stim class due to imbalanced nr of trials for SCC
                    
                    
                else % if we deal with MRI and version A/B
                    
                    for st = [1,2]
                        if ismember(ii,[1:3]) % if GBR/RDK/DOT
                            tmp = condition_master(condition_master.is_catch==0 & condition_master.session_type==st & isnan(condition_master.is_objectcatch) & any(strcmp(condition_master.stim_class_name, params.exp.stimclassnames{ii}),2),:);
                        elseif ii == 4 % OBJ
                            tmp = condition_master(condition_master.is_catch==0 & condition_master.session_type==st & any(strcmp(condition_master.stim_class_name, params.exp.stimclassnames{ii}),2),:);
                            tmp0 = cat(1,tmp(isnan(tmp.is_objectcatch),:),tmp(tmp.is_objectcatch==0,:)); tmp = tmp0; clear tmp0;
                            tmp_catch = condition_master(condition_master.is_catch==0 & condition_master.is_objectcatch==1 & condition_master.session_type==st & any(strcmp(condition_master.stim_class_name, params.exp.stimclassnames{ii}),2),:);
                        end
                        
                        % Get counts of stimulus numbers for left and right stimulus position.
                        [N, ~] = histcounts(cat(1,tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})),....
                            tmp.stim_nr_right(strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))), [params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_core, max(params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_core)+1]);
                        
                        N_tbl(ii,1:length(N)) = N;
                        if ~params.is_demo
                            if st==1
                                assert(isequal(sum(N>0), length(all_unique_im.(params.exp.stimclassnames{ii}).unique_im_nr)))
                            elseif strcmp(env_type,'MRI') && params.is_wide == 0 && st == 2
                                cm_wideB = cm_wide(cm_wide.session_type==st,:);
                                [N_wideB, ~] = histcounts(cat(1,cm_wideB.stim_nr_left(strcmp(cm_wideB.stim_class_name(:,1), params.exp.stimclassnames{ii})),....
                                    cm_wideB.stim_nr_right(strcmp(cm_wideB.stim_class_name(:,2), params.exp.stimclassnames{ii}))), [params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_core, max(params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_core)+1]);
                                N_B = N + N_wideB;
                                assert(isequal(sum(N_B>0), length(all_unique_im.(params.exp.stimclassnames{ii}).unique_im_nr)));
                                
                                colsL_wideB = cm_wideB.stim_nr_left(~ismember(cm_wideB.task_class,3) & isnan(cm_wideB.is_objectcatch) & strcmp(cm_wideB.stim_class_name(:,1), params.exp.stimclassnames{ii})); % exclude SCC, we deal with those trials separately
                                colsR_wideB  = cm_wideB.stim_nr_right(~ismember(cm_wideB.task_class,3) & isnan(cm_wideB.is_objectcatch) & strcmp(cm_wideB.stim_class_name(:,2), params.exp.stimclassnames{ii})); % exclude SCC, we deal with those trials separately
                                colsLR_wideB  = [colsL_wideB(:),colsR_wideB(:)]; % if this fails, then we don't have equal left and right stimulus numbers for classic stimuli
                                colsLR_scc_wideB  = [cm_wideB.stim_nr_left(cm_wideB.task_class ==3 & strcmp(cm_wideB.stim_class_name(:,1), params.exp.stimclassnames{ii})); ... % deal with scc trials
                                    cm_wideB.stim_nr_right(cm_wideB.task_class ==3 & strcmp(cm_wideB.stim_class_name(:,2), params.exp.stimclassnames{ii}))];
                            end
                        end
                        % WIDE only has SCC
                        if params.is_wide
                            colsL = tmp.stim_nr_left(tmp.task_class ~=3 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); % exclude SCC, we deal with those trials separately
                            colsR = tmp.stim_nr_right(tmp.task_class ~=3 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii})); % exclude SCC, we deal with those trials separately
                            colsLR = [colsL(:),colsR(:)]; % if this fails, then we don't have equal left and right stimulus numbers for classic stimuli
                            unique_im_from_table(ii) = length(unique(colsLR)); % this should be 24 for GBR or RDK, 16 for DOT/OBJ
                            assert(isequal(unique_im_from_table(ii),params.exp.nr_unique_trials_per_crossing(ii,2)))
                            colsLR_scc = [tmp.stim_nr_left(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); ... % deal with scc trials
                                tmp.stim_nr_right(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))];
                            colsLR_ltm = [];
                        else
                            % Deep has SCC an LTM
                            colsL = tmp.stim_nr_left(~ismember(tmp.task_class,[3,6,7]) & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); % exclude SCC/LTM/IMG, we deal with those trials separately
                            colsR = tmp.stim_nr_right(~ismember(tmp.task_class,[3,6,7]) & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii})); % exclude SCC/LTM/IMG, we deal with those trials separately
                            colsLR = [colsL(:),colsR(:)]; % if this fails, then we don't have equal left and right stimulus numbers for classic stimuli
                            if st == 2
                                unique_im_from_table(ii) = length(unique([colsLR;colsLR_wideB])); % this should be 24 for GBR or RDK, 16 for DOT/OBJ
                                assert(isequal(unique_im_from_table(ii),params.exp.nr_unique_trials_per_crossing(ii,2)))
                            else
                                unique_im_from_table(ii) = length(unique(colsLR)); % this should be 24 for GBR or RDK, 16 for DOT/OBJ
                                assert(isequal(unique_im_from_table(ii),params.exp.nr_unique_trials_per_crossing(ii,2)))
                            end
                            colsLR_scc = [tmp.stim_nr_left(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); ... % deal with scc trials
                                tmp.stim_nr_right(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))];
                            colsLR_ltm = [tmp.stim_nr_left(tmp.task_class ==6 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); ... % deal with ltm trials
                                tmp.stim_nr_right(tmp.task_class ==6 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))];
                            colsLR_img = [tmp.stim_nr_left(tmp.task_class ==7 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); ... % deal with ltm trials
                                tmp.stim_nr_right(tmp.task_class ==7 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))];
                        end
                        
                        % Calculate the expected nr of trials from the params.exp.mri.ses_blocks table.
                        expected_nr_of_trials = sum(squeeze(sum(all_sessions(ii,[1,2,4:10],:,st).*params.exp.nr_trials_per_block(ii,[1,2,4:10])))) + ...
                            sum(all_sessions(ii,3,:,st)).*scc_probability(ii) + sum(all_sessions(ii,6,:,st)).*ltm_probability(ii);
                        
                        if params.is_demo
                            expected_nr_of_trials = expected_nr_of_trials*0.5; % again: demo has half the trials per block
                        end
                        if params.is_wide
                            empirical_nr_of_trials = ceil(sum([size(colsLR,1),size(unique(colsLR_scc),1)/2, size(unique(colsLR_ltm),1)/2])); % for GBR/RDK/DOT/OBJ: sum regular trials, scc trials (OBJ: we treat objectcatch trials as regular trials)
                        else
                            empirical_nr_of_trials = ceil(sum([size(colsLR,1),size(colsLR_scc,1), size(colsLR_ltm,1), size(colsLR_img,1)])); % for GBR/RDK/DOT/OBJ: sum regular trials, scc trials (OBJ: we treat objectcatch trials as regular trials)

                            if ii == 4 % add object catch trials
                                empirical_nr_of_trials = empirical_nr_of_trials + size(tmp_catch,1);
                            end
                        end
                        if st==2
                            if ~(abs(floor(empirical_nr_of_trials)-expected_nr_of_trials)<=(scc_tolerance(st)+ltm_tolerance(st))) % we allow for x trials difference per stim class due to imbalanced nr of trials for SCC
                                warning('[%s]: %s: there is a difference of %d between the number empirical and expected trials for version B!', mfilename, params.exp.stimclassnames{ii}, abs(floor(empirical_nr_of_trials)-expected_nr_of_trials))
                            end
                        else
                            assert(abs(floor(empirical_nr_of_trials)-expected_nr_of_trials)<=(scc_tolerance(st)+ltm_tolerance(st))); % we allow for x trials difference per stim class due to imbalanced nr of trials for SCC
                        end
                    end
                end
            end
        end
    end
    % Now dive into each task class
    M_tbl = NaN(length(params.exp.taskclassnames),length(params.stim.all_core_im_nrs)+1);
    for ii = unique(condition_master.task_class)'
        [M, ~] = histcounts([condition_master.stim_nr_left(condition_master.is_catch==0 & (condition_master.task_class==ii)); ...
            condition_master.stim_nr_left(condition_master.is_catch==0 & (condition_master.task_class==ii))],[0:1:length(params.stim.all_core_im_nrs)+1]);
        M_tbl(ii,1:length(M)) = M;
    end
    
    
    
    
    %% Plot figures to check condition master content
    if verbose
        
        vcd_visualizeMasterTable(condition_master, store_imgs,env_type);
        
        figure; set(gcf,'Position',[1,1,1200,300]);
        histogram([condition_master.stim_nr_left;condition_master.stim_nr_right],'numbins',length(params.stim.all_core_im_nrs))
        xlabel('unique condition nr'); ylabel('trial count')
        title('Total unique image nr')
        box off;
        if store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',env_type));
            if ~exist(saveFigsFolder,'dir'); mkdir(saveFigsFolder); end
            filename = sprintf('vcd_totaluniqueim.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        figure; set(gcf,'Position',[1,1,1200,300]); imagesc(N_tbl)
        xlabel('unique im nr'); ylabel('Stim class');
        C = parula(20); colormap(C); set(gca,'CLim',[0 max(N_tbl(~isnan(N_tbl)))])
        cb = colorbar; cb.Label.String = 'unique im count';
        cb.Ticks = [0,unique(N_tbl(~isnan(N_tbl)))'];
        set(gca,'YTick',[1:5],'YTickLabel', params.exp.stimclassnames);  set(gca,'XTick',[1:30]); title('Sum of unique im per stimulus class')
        box off;
        if store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',env_type));
            filename = sprintf('vcd_uniqueim_per_stimclass.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        
        figure; set(gcf,'Position',[1,1,1200,300]); imagesc(M_tbl)
        xlabel('unique condition nr'); ylabel('Task class');
        C = parula(20); colormap(C); set(gca,'CLim',[0 max(M_tbl(~isnan(M_tbl)))])
        cb = colorbar; cb.Label.String = 'unique im count';
        cb.Ticks = unique(M_tbl(~isnan(M_tbl)))';
        set(gca,'YTick',[1:10],'YTickLabel', params.exp.taskclassnames);
        set(gca,'XTick',[0:25:110]); title('Sum of unique im per task class')
        box off;
        if store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',env_type));
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
        set(gca,'YTick',[1:params.exp.block.n_trials_single_epoch]); title('SCC stimulus class distribution')
        box off;
        if store_imgs
            saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',env_type));
            filename = sprintf('vcd_scc_stimulusclass_trialdistr.png');
            print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
        end
        
        if strcmp(env_type,'MRI')
            % Stim class distribution in LTM blocks
            foo = [condition_master.stim_nr_left((condition_master.task_class==6),:), ...
                condition_master.stim_nr_right((condition_master.task_class==6),:)]';
            figure; set(gcf,'Position',[1,1,1200,300])
            imagesc(foo)
            xlabel('Trial nr'); ylabel('left vs right stim');
            C = parula(max(foo(:))); colormap(C); set(gca,'CLim',[min(foo(:)),max(foo(:))])
            cb = colorbar; cb.Label.String = 'unique im nr';
            cb.Ticks = [min(foo(:)):10:max(foo(:))];
            set(gca,'YTick',[1:params.exp.block.n_trials_single_epoch]); title('LTM stimulus class distribution')
            box off;
            if store_imgs
                saveFigsFolder = fullfile(vcd_rootPath,'figs',sprintf('condition_master0_%s',env_type));
                filename = sprintf('vcd_ltm_stimulusclass_trialdistr.png');
                print(gcf,'-dpng','-r300',fullfile(saveFigsFolder,filename));
            end
        end
    end % if verbose
end % if load params

return





