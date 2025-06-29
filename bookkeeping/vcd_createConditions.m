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
%  20: 'is_catch'         : (logical, true/false) whether the trial is
%                            a catch trial where no stimulus showed up.
%                            true = yes, this is a catch trial. false = no,
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
%  38: 'is_special_core'  : (logical, true/false) Is this stimulus part of
%                            the special core stimulus set (true) or not 
%                            (false)? Special core stimuli are a subset of 
%                            core stimuli that are only used for LTM and
%                            IMG task crossings.
%  39: 'is_lure'          : (logical, true/false) For LTM task crossings: 
%                            did we use a novel lure stimulus (true)
%                            or not (false).
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
p0.addParameter('store_params' , true, @islogical);
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
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('condition_master_%s%s_%s*.mat',choose(params.is_demo,'demo_',''), params.disp.name,env_type)));
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
    % Preallocate space and set up tables/structs
    if verbose
        tic
        fprintf('\n[%s]: Start creating conditions for %s experiment.. \n',mfilename,env_type);
    end
    
    % assume we don't run a demo run if parameter isn't specified;
    if ~isfield(params,'is_demo'), params.is_demo = false; end
    
    all_unique_im     = struct();  % details about unique core images present in VCD-core experiment. This information is also present in condition_master.
    all_cond          = struct();  % details about unique conditions (unique images x cueing condition) present in VCD-core experiment. This information is also present in condition_master.
    
    % master table with all conditions shown across all sessions. 
    % Note that block_nr counts continuously because we will shuffle the 
    % block order later for each subject. Same holds for run_nr, these
    % numbers do not represent the final nr of runs as they will change
    % depending on the allocation of blocks to each run within a session.
    condition_master  = table();   
    
    % Define block content for each stimulus class
    for stimClass_idx = 1:length(params.exp.stimclassnames)
        
        stim_table = table();
        
        % Get stimclass name
        stimClass_name = params.exp.stimclassnames{stimClass_idx};
        
        % Get corresponding task crossings
        bsc_idx = cell2mat(cellfun(@(x) strcmp(params.exp.stimclassnames,x), {stimClass_name}, 'UniformOutput', false));
        
        if strcmp(env_type,'MRI')
            task_crossings = find( params.exp.crossings(bsc_idx,:));
        elseif strcmp(env_type,'BEHAVIOR')
            if params.is_demo
                task_crossings = find(params.exp.n_unique_trial_repeats_demo(bsc_idx,:)>0);
            else
                task_crossings = find(params.exp.n_unique_trial_repeats_behavior(bsc_idx,:)>0);
            end
        end
                
        % Loop over each task crossing for this stim class
        for curr_task = 1:length(task_crossings)

            % Define the unique images for Gabors
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
                t_cond = t_cond(t_cond.is_special_core,:);
            end
            
            % Get the number of trials and blocks, which depend on task
            if params.exp.trial.single_epoch_tasks(task_crossings(curr_task))
                n_trials_per_block    = params.exp.block.n_trials_single_epoch;
            else
                n_trials_per_block    = params.exp.block.n_trials_double_epoch;
            end

            %% ---- IMPORTANT FUNCTION: Create condition master table ---- %%
            % Create condition master table, where unique images are
            % shuffled and distributed across blocks according to the
            % stimulus feature of interest that receive priority. For
            % example, we want to make sure we present all 3 gabor
            % contrasts within a block at least once.
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
            stim_table = cat(1,stim_table,tbl);
            
        end % task class idx
        
        % Accumulate condition master and info
        all_cond.(stimClass_name) = stim_table;
        
        condition_master = cat(1,condition_master,stim_table);
    end % stim idx
    
    % ---- Add object catch trials ----
    % How many trials are object catch given the 20% probility?
    pcobj_trial_idx = find(condition_master.stim_class == 4 & condition_master.task_class == 4);
    nr_objectcatch_trials_per_rep = round(size(condition_master(pcobj_trial_idx,:),1)*params.exp.trial.pc.prob_objcatch);
    
    % Select object catch trials randomly from the total list (without
    % replacement). Ensure we distribute object catch trials across
    % left/right cued locations and unique objects (as much as possible)
    while 1
        objcatch_ok = false(1,3);
        
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
        
        if sum(objcatch_ok)==length(objcatch_ok)
            break
        end
    end
    
        
    % Get the unique image numbers for object catch stimuli
    % (reshape to 18 catch rotations x 16 objects)
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
        end
        
        % What core object (1-16) are we dealing with?
        old_stim_obj_nr   = (old_stim_nr==params.stim.obj.unique_im_nrs_core);
        
        % What are the possible object catch rotations we can use?
        possible_objcatch_rotations0 = possible_objcatch_rotations(old_stim_obj_nr,:);
        possible_objcatch_rotations1 = possible_objcatch_rotations0(~isnan(possible_objcatch_rotations0));
        
        % randomly select object catch rotation
        objcatch_rot = randi(length(possible_objcatch_rotations1),1);
        
        % update rotation of object
        condition_master.orient_dir(objcatch_idx(cc),cued_loc_objcatch(cc)) = possible_objcatch_rotations1(objcatch_rot);
        
        % remove rotation from list
        possible_objcatch_rotations(old_stim_obj_nr,objcatch_rot) = NaN;

        % Update special core column (objectcatch stimuli can never be 
        % special core stimuli)
        condition_master.is_special_core(objcatch_idx(cc),cued_loc_objcatch(cc)) = false;

        % Insert object catch image number into "is_objectcatch" for now..
        % (we need to hold on to the core object number for the condition
        % label).
        condition_master.is_objectcatch(objcatch_idx(cc)) = catch_im_nr(old_stim_obj_nr,objcatch_rot);

    end    
    
    
    % ---- IMPORTANT STEP: Shuffle stimuli for SCC and LTM task ----
    condition_master = vcd_shuffleStimForTaskClass(params,  'scc', condition_master, params.exp.block.n_trials_single_epoch,env_type);
    if strcmp(env_type,'MRI')
        condition_master = vcd_shuffleStimForTaskClass(params,  'ltm', condition_master, params.exp.block.n_trials_double_epoch,env_type);
    end
    
    % ---- IMPORTANT STEP: Add correct button press ----
    button_response = NaN(size(condition_master,1),1);
    for ii = 1:size(condition_master,1)
        if  condition_master.is_catch(ii)
            button_response(ii) = NaN;
        else
            button_response(ii) = vcd_getCorrectButtonResponse(params, condition_master(ii,:));
        end
    end
    condition_master.correct_response = button_response;

    % ---- bookkeeping: See if we can achieve roughly equally sample scc images within a block
    scc_trials0 = find(condition_master.stim_class==99);
    cued0 = condition_master.is_cued(scc_trials0);
    scc_trials = scc_trials0; % make a copy
    if ~isempty(scc_trials)
        while 1
            bb_bpress = reshape(condition_master.correct_response(scc_trials),params.exp.block.n_trials_single_epoch,[]);
            bb_cued   = reshape(cued0,params.exp.block.n_trials_single_epoch,[]);
            n0 = zeros(size(bb_bpress,2),4);
            m0 = zeros(size(bb_bpress,2),2);
            for bb = 1:size(bb_bpress,2)
                n0(bb,:) = histcounts(bb_bpress(:,bb),[1:5]);
                m0(bb,:) = histcounts(bb_cued(:,bb));
            end
            too_few = sum(n0==0,1);
            cued_ok = sum(m0==0,1);
            if all(too_few==false(1,4)) && all(cued_ok==false(1,2))
                break;
            else
                % shuffle trials while preserving 50:50 cue left/right
                cued00 = cued0;
                cueleft_idx0         = find(cued0==1);
                cueright_idx0        = find(cued0==2);
                cueleft_idx          = cueleft_idx0(randperm(length(cueleft_idx0),length(cueleft_idx0)));
                cueright_idx         = cueright_idx0(randperm(length(cueright_idx0),length(cueright_idx0)));
                cued0(cueleft_idx0)  = cued0(cueleft_idx); 
                cued0(cueright_idx0) = cued0(cueright_idx);
                assert(isequal(cued00,cued0)); % cued status should not change
                scc_trials(cueleft_idx0)  = scc_trials(cueleft_idx);
                scc_trials(cueright_idx0) = scc_trials(cueright_idx0);
            end
        end
        % reshuffle order
        condition_master(scc_trials0,:) = condition_master(scc_trials,:);
        % fix unique trial nr
        scc_trials = find(condition_master.stim_class==99);
        condition_master.unique_trial_nr(scc_trials) = [1:length(condition_master.unique_trial_nr(scc_trials))]';
        condition_master.trial_nr(scc_trials) = repmat([1:params.exp.block.n_trials_single_epoch]',length(condition_master.unique_trial_nr(scc_trials))/params.exp.block.n_trials_single_epoch,1);
        condition_master.stim_class_unique_block_nr(scc_trials) = repelem([1:length(condition_master.unique_trial_nr(scc_trials))/params.exp.block.n_trials_single_epoch],params.exp.block.n_trials_single_epoch)';
    end
    
    
    % ---- bookkeeping: Check if button presses are balanced to the extent possible
    curr_sc = unique(condition_master.stim_class);
    curr_tc = unique(condition_master.task_class);
    curr_tc = curr_tc(curr_tc~=1); % exclude fixation
    for jj = 1:length(curr_sc)
        for mm = 1:length(curr_tc)
            % check if we have such a crossing
            if sum(ismember(condition_master.task_class,curr_tc(mm)) & ismember(condition_master.stim_class,curr_sc(jj)))>0
                if ismember(curr_tc(mm),[2,4,5,6,7]) % cd, pc, wm, ltm, img have 2 response options
                    nr_responses = 2;
                elseif ismember(curr_tc(mm),[3,8,10]) % scc, what, how have 4 response options
                    nr_responses = 4;
                elseif ismember(curr_tc(mm),[9]) % where has 3 response options
                    nr_responses = 3;
                end
                resp = condition_master.correct_response( ...
                    condition_master.stim_class==curr_sc(jj) & ...
                    condition_master.task_class==curr_tc(mm));
                n = histcounts(resp,1:(nr_responses+1));
                % For WHAT/HOW tasks we go by category info, because we combine
                % object and foods into one button press..
                if ismember(curr_sc(jj),[4,5]) && ismember(curr_tc(mm),8)
                    supercat = condition_master.super_cat(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                    if curr_sc(jj) == 4 % objects
                        cued0  = condition_master.is_cued(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
                        n0     = histcounts([supercat(cued0==1,1);supercat(cued0==2,2)],1:6);
                    elseif curr_sc(jj) == 5 % scenes, no left/right, only center
                        n0     = histcounts(supercat,1:6);
                    end
                    n1 = [n0([1,2]), n0(3)+n0(4), n0(5)];
                    assert(all(n==n1))
                elseif ismember(curr_sc(jj),4) && ismember(curr_tc(mm),4) % PC-OBJ
                    resp2 = resp(condition_master.is_objectcatch( ...
                        condition_master.stim_class==curr_sc(jj) & ...
                        condition_master.task_class==curr_tc(mm))==0);
                    n2 = histcounts(resp2,1:(nr_responses+1));
                    if mod(length(resp2),2)==1 % if we have an uneven nr of trials after we removed objectcatch trials
                        if ~(diff(n2)<=1) % we allow for a difference of one, but not more
                            error('[%s]: Button response counts diverge more than 1 and are considered unbalanced across options! Please rerun vcd_createConditions.m',mfilename);
                        end
                    else % assume we have an even nr of trials after we removed objectcatch trials
                        if ~(all(n==n0)) % we don't allow for a difference of one
                            error('[%s]: Uneven nr of button responses! Please rerun vcd_createConditions.m',mfilename);
                        end
                    end
                elseif ismember(curr_sc(jj),5) && ismember(curr_tc(mm),9) % NS-WHERE
                    subcat = condition_master.sub_cat(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                    n0     = histcounts(subcat,1:(nr_responses+1));
                    assert(all(n==n0))
                elseif ismember(curr_sc(jj),[4,5]) && ismember(curr_tc(mm),10) % OBJ-HOW & NS-HOW
                    affordcat = condition_master.affordance_cat(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                    if curr_sc(jj) == 4
                        cued0 = condition_master.is_cued(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
                        n1    = histcounts([affordcat(cued0==1,1);affordcat(cued0==2,2)],1:(nr_responses+1));
                    elseif curr_sc(jj) == 5
                        n1     = histcounts(affordcat,1:(nr_responses+1));
                    end
                    assert(all(n==n1))
                elseif ismember(curr_sc(jj),99) || ismember(curr_tc(mm),3) % scc-all
                    % Check if we sample all core images across trials
                    stmclass0 = condition_master.stim_class_name(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm),:);
                    cued0     = condition_master.is_cued(condition_master.stim_class==curr_sc(jj) & condition_master.task_class==curr_tc(mm));
                    stmclass0_cued = [stmclass0(cued0==1,1);stmclass0(cued0==2,2)];
                    [~,stmclass_cued_i] = ismember(stmclass0_cued,params.exp.stimclassnames([1,3,2,4]));
                    n1 = histcounts(stmclass_cued_i,1:(nr_responses+1));
                    assert(all(n==n1))
                elseif all(n==0) && curr_tc(mm)~=2 && ~ismember(curr_tc(mm),[8,9,10])
                    error('[%s]: No correct responses found?!',mfilename);
                else
                    assert(all(diff(n)==0))
                end
            end
        end
    end
    
    %% ---- IMPORTANT FUNCTION: Allocate trials to blocks all unique trials and repeats of trials.
    % !!WARNING!! There is a randomization component involved in creating the
    % conditions (i.e., order of trials allocated to a block). If you don't
    % want this, set params.load_params = true to load an existing
    % condition_master.
    condition_master = vcd_allocateBlocksToRuns(params,condition_master,env_type);
    
    
    %% ---- IMPORTANT FUNCTION: Add contrast decrement
    condition_master = vcd_determineContrastDecrementChangeTrials(params, condition_master);
    
    % Convert is_objectcatch vector in logical.. 
    condition_master.is_objectcatch = logical(condition_master.is_objectcatch);
    
    %% Store condition_master if requested
    if store_params
        fprintf('[%s]:Storing condition_master..\n',mfilename)
        saveDir = fullfile(vcd_rootPath,'workspaces','info');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end
        if params.is_demo
            fname = sprintf('condition_master_demo_%s_%s.mat',params.disp.name,datestr(now,30));
        else
            fname = sprintf('condition_master_%s_%s.mat',params.disp.name,datestr(now,30));
        end
        save(fullfile(saveDir,fname),'condition_master','all_unique_im','all_cond')
    end
    
    
    %% At last, do some checks:
    
    % Stim class nr should match with params
    assert(isequal(unique(condition_master.stim_class)',[1:length(params.exp.stimclassnames),99]))
    % Task class nr should match with params
    if strcmp(env_type,'MRI')
        assert(isequal(unique(condition_master.task_class)',1:length(params.exp.taskclassnames)));
    elseif strcmp(env_type,'BEHAVIOR')
        if params.is_demo
            assert(isequal(unique(condition_master.task_class)', find(any(params.exp.n_unique_trial_repeats_demo,1))))
        else
            assert(isequal(unique(condition_master.task_class)', find(any(params.exp.n_unique_trial_repeats_behavior,1))))
        end
    end
    
    % Cued vs uncued stimuli should be matched for non-demo sessions
    if ~params.is_demo
        assert(isequal(sum(condition_master.is_cued==1),sum(condition_master.is_cued==2)))
    end
    
    % Now dive into each stimulus class
    unique_im_from_table = [];
    N_tbl = NaN(5,30); 
    adjusted_crossings = params.exp.nr_unique_trials_per_crossing;
                
    all_sessions = vcd_getSessionEnvironmentParams(params, env_type);
    for ii = 1:length(params.exp.stimclassnames)
        
        % Check nr of unique stimuli per stimulus class
        if ii == 5
            tmp = condition_master(~condition_master.is_catch & ~condition_master.is_objectcatch & strcmp(condition_master.stim_class_name(:,1), params.exp.stimclassnames{ii}),:);
        
            [N, ~] = histcounts(tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})));
        
            N_tbl(ii,1:length(N)) = N;
            if ~params.is_demo
                assert(isequal(size(N,2), length(all_unique_im.(params.exp.stimclassnames{ii}).unique_im_nr)))
            end
            unique_im_from_table(ii) = length(unique(tmp.stim_nr_left));

        
%             if strcmp(env_type,'MRI')
%                 % lower IMG/LTM block contribution (we have less unique images)
%                 unique_trial_repeats = params.exp.n_unique_trial_repeats_mri;
%             elseif strcmp(env_type,'BEHAVIOR') % NO LTM/IMG
%                 adjusted_crossings(ii,6) = 0;
%                 adjusted_crossings(ii,7) = 0;
%                 half_blocks = all_sessions>0 & all_sessions<1;
%                 if sum(half_blocks)>0
%                     unique_trial_repeats(half_blocks)=all_sessions(half_blocks);
%                 end
%             end

            expected_nr_of_trials = sum(sum(all_sessions(ii,:,:).*params.exp.nr_trials_per_block(ii,:,:)));
            if params.is_demo
                expected_nr_of_trials = 0.5*expected_nr_of_trials;
            end
            empirical_nr_of_trials = size(tmp,1);
            assert(isequal(empirical_nr_of_trials, expected_nr_of_trials))
            
        else
            tmp = condition_master(~condition_master.is_catch & ~condition_master.is_objectcatch & any(strcmp(condition_master.stim_class_name, params.exp.stimclassnames{ii}),2),:);
            
            [N, ~] = histcounts(cat(1,tmp.stim_nr_left(strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})),....
                                            tmp.stim_nr_right(strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))));
            
            N_tbl(ii,1:length(N)) = N;
            if ~params.is_demo
                assert(isequal(size(N,2), length(all_unique_im.(params.exp.stimclassnames{ii}).unique_im_nr)))
            end
            colsL = tmp.stim_nr_left(tmp.task_class ~=3 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii}));
            colsR = tmp.stim_nr_right(tmp.task_class ~=3 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}));
            colsLR = [colsL(:),colsR(:)];
            unique_im_from_table(ii) = length(unique(colsLR));
            colsLR_scc = [tmp.stim_nr_left(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,1), params.exp.stimclassnames{ii})); ...
                            tmp.stim_nr_right(tmp.task_class ==3 & strcmp(tmp.stim_class_name(:,2), params.exp.stimclassnames{ii}))];
            
            expected_nr_of_trials = sum(sum(all_sessions(ii,:,:).*params.exp.nr_trials_per_block(ii,:)));
            if params.is_demo
                expected_nr_of_trials = expected_nr_of_trials*0.5;
            end
            
            if strcmp(params.exp.stimclassnames{ii},'obj')
                empirical_nr_of_trials = sum([size(colsLR,1),size(unique(colsLR_scc),1)/2])+sum(condition_master.is_objectcatch);
            else
                empirical_nr_of_trials = sum([size(colsLR,1),size(unique(colsLR_scc),1)/2]);
            end
            assert((abs(floor(empirical_nr_of_trials)-expected_nr_of_trials))<=2); % we allow for 1-2 trials difference per stim class due to imbalanced nr of trials for SCC
        end
    end
    
    % Now dive into each task class
    M_tbl = NaN(length(params.exp.taskclassnames),length(params.stim.all_core_im_nrs)+1);
    for ii = unique(condition_master.task_class)'
        [M, ~] = histcounts([condition_master.stim_nr_left(~condition_master.is_catch & (condition_master.task_class==ii)); ...
                                    condition_master.stim_nr_left(~condition_master.is_catch & (condition_master.task_class==ii))],[0:1:length(params.stim.all_core_im_nrs)+1]);
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
    end
end % load params

return





