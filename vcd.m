function outputs = vcd(varargin)
% VCD function to call general information about experimental design, 
% including the number of stimulus and task class names, numbers, unique
% image number, full information for each unique stimulus.
%
%   outputs = vcd('varname1',var1,'varname2',var2, ...'varnameN',varN);
%
% Input argument names options, corresponding input type, and an example input type:
%  *  'stimulusclassnames',        [stimclass_number] or []                 -- example input: [1] -> output: 'gabor'
%  *  'stimulusclassnumbers',      {'stimclass_name'} or {}                 -- example input: {'gabor'} -> output: [1]
%  *  'taskclassnames',            [taskclass_number] or []                 -- example input: [4] -> output: 'pc'
%  *  'taskclassnumbers',          {'taskclass_name'} or {}                 -- example input: {'pc'} -> output: [4]
%  *  'crossingnames',             [crossing_number] or []                  -- example input: [32] -> output: 'gabor'
%  *  'crossingnumbers',           {'crossing_name'} or {}                  -- example input: {'how'} -> output: 'gabor'
%  *  'allstimulusnumbers',        {'stimclass_name'} or {}                 -- example input: {'rdk'} -> output: 'gabor'
%  *  'stimulusnumberstonames',    [stim_number]                            -- example input: [56] -> output: 'gabor'
%  *  'stimulusnamestonumbers',    {stim_name}                              -- example input: {'DOT-0056-L'} -> output: 'gabor'
%  *  'allcore',                   [stim_number], {stim_name}, [], or {}    -- example input: [4] or {'obj'} -> output: 'gabor'
%  *  'specialcore',               [stim_number], {stim_name}, [], or {}    -- example input: [3] or {'dot'} -> output: 'gabor'
%  *  'stimtostimclassname',       [stim_number]                            -- example input: [1400] -> output: 'gabor'
%  *  'stimtostimclassnumber',     [stim_number]                            -- example input: [800] -> output: 'gabor'
%  *  'stimtotaskclassname',       [stim_number]                            -- example input: [21] -> output: 'gabor'
%  *  'stimtotaskclassnumber',     [stim_number]                            -- example input: [480] -> output: 'gabor'
%  *  'fullinfo',                  [stim_number]                            -- example input: [3] -> output: 'gabor'
%  *  'conditionnumbertoname',     [stim_number]                            -- example input: [340] -> output: 'gabor'
%  *  'conditionnametonumber',     {'condition_name'}                       -- example input: {'OBJ-0065-L-UNCUED-WHAT'} -> output: 'gabor'
%  *  'allwmteststimulusnumbers',  [stimclass_number], {'stimclass_name'}, [], or {} -- example: [1] or {'gabor'} -> output: 'gabor'
%  *  'allltmteststimulusnumbers', [stimclass_number], {'stimclass_name'}, [], or {} -- example: [5] or {'ns'} -> output: 'gabor'
%  *  'allimgteststimulusnumbers', [stimclass_number], {'stimclass_name'}, [], or {} -- example: [4] or {'obj'} -> output: 'gabor'
%
%
% Note that 'fullinfo', 'stimulusnumberstonames', and
% 'stimulusnamestonumbers' requires the condition_master table. This
% variable is stored here: 
%   vcd_rootPath/workspaces/info/condition_master_<dispname>_YYYYMMDDTHHMMSS.mat.
%
% You can generate a new condition_master.mat file by running: 
%   [~, condition_master] = vcd_createConditions(params);
% See also s_createDesignMatrix.m
% 
% % --- STIMULUS CLASSES ---
% Abbreviations for stimulus classes are: 
%  1:'gabor'  - Gabor (grating subject to Gaussian window).
%  2:'rdk'    - Random dot motion kinetogram
%  3:'dot'    - Single dot
%  4:'obj'    - Objects
%  5:'ns'     - Natural scenes
%
% % --- TASK CLASSES ---
% Abbreviations for task classes are: 
%   1:'fix'   - Fixation brightness task
%   2:'cd'    - Contrast change task
%   3:'scc'   - Categorization task
%   4:'pc'    - Perceptual categorization tasks (tilt, motion direction, 
%                            dot position, object rotation, indoor/outdoor)
%   5:'wm'    - Working memory tasks (tilt, motion direction, position, 
%                            rotation, scene)
%   6:'ltm'   - Matching task
%   7:'img'   - Imagery task
%   8:'what'  - "What?" task
%   9:'where' - "Where?" task
%   10:'how'  - "How?" task
% 
% % --- STIMULUS NUMBERS ---
% Any stimulus used in VCD has its own unique integer number between 1-1421.
% # 1-110   are core VCD stimuli (used in all stimulus-task crossing)
% # 111-422 are WM test stimuli  (used in the second stimulus array after 
%                                the delay, for working memory task
%                                crossings only).
% # 1-110 are core VCD stimuli (used in all stimulus-task crossing)
%
% % --- CONDITION NAMES AND NUMBERS ---
% Stimulus numbers are different from Condition numbers.
% Condition numbers refer to the combination of:
% (1) a unique stimulus (e.g., "GBR-0001-L"). Note that each stimulus has
%  a fixed spatial location: either left or right (for gabors, RDKs, single
%  dots, and objects) or center (for scenes).
% (2) subject to a particular a cueing state ("CUED","UNCUED", or "NCUED" 
%                           the latter refers for neutral cue)
% (3) subject to a particular task crossing (e.g., "PC" or "FIX")
% All possible conditions in VCD are stored in vcd_getConditionNames().
%
% Possible inputs:
% 'stimulusclassnames'    : Provide all stimulus class names (use []),
%                           or translate stimulus class number(s) to
%                           name(s). Stimulus class number should be an
%                           integral number between 1 and 5, can be a
%                           single number or vector. Output will preserve
%                           the order specified by the user's input(s), or
%                           ascending if input is empty.
% 'stimulusclassnumbers'  : Provide all stimulus class numbers (use []),
%                           or translate provided stimulus class name(s) to
%                           number(s). Stimulus class name should be a
%                           string or character class, and should be part
%                           of the 5 stimulus classes, where 1:'gabor',
%                           2:'rdk', 3:'dot',4:'obj',5:'ns'. Name can be
%                           lower case or upper case, and can be a single
%                           string char 'rdk' or cell with one or multiple
%                           names {'rdk','gabor'}). Output will preserve
%                           the order specified by the user's input(s), or
%                           ascending if input is empty.
% 'taskclassnames'        : Provide all task class names (use [] or {})
%                           or translate task class number(s) to name(s).
%                           Task class number should be an integral number
%                           between 1 and 10, can be a single number or
%                           vector. Output will preserve the order
%                           specified by the user's input(s) or ascending
%                           if input is empty, where 1:'fix', 2:'cd',
%                           3:'scc', 4:'pc', 5:'wm', 6:'ltm', 7:'img',
%                           8:'what', 9:'where', 10:'how'.
% 'taskclassnumbers'      : Provide all task class numbers (use [] or {}),
%                           or translate task class name(s) to number(s). 
%                           Task class name should be a string or character
%                           class, and should be part of the 10 task
%                           classes: 'fix','cd','scc','pc','wm', 'ltm'
%                           'img','what','where','how'. Name can be lower
%                           case or upper case. Name can be a single name
%                           ('fix') or cell with multiple names
%                           {'scc','wm'}). Output will preserve the order
%                           specified by the user's input(s). or ascending
%                           if input is empty.
% 'crossingnames'         : Provide all task-stimulus class crossing names,
%                           (use [] or {}) or translate crossing
%                           numbers(s) to crossing names(s). Crossing
%                           numbers should be an integral number between 1
%                           and 32. Inputs can be empty vector, a single
%                           number or a list of numbers in a vector. Output
%                           will preserve the order specified by the user's
%                           input(s) or ascending if input is empty.
% 'crossingnumbers'       : Provide all crossings between task and stimulus
%                           classes used in VCD (use [] or []), 
%                           or translate, crossing name(s) to crossing 
%                           number(s). Crossing name should be a string or 
%                           character class in the format <'task-stim'>, 
%                           and should be part of the 32 crossings:
%                           'fix-gabor','cd-gabor','pc-gabor','wm-gabor',
%                           'img-gabor','fix-rdk','cd-rdk','pc-rdk',
%                           'wm-rdk','img-rdk','fix-dot','cd-dot','pc-dot',
%                           'wm-dot','img-dot''fix-obj','cd-obj','pc-obj',
%                           'wm-obj','img-obj','what-obj','how-obj',
%                           'fix-ns','cd-ns','pc-ns','wm-ns','img-ns',
%                           'what-ns','where-ns','how-ns'. Note that for
%                           'ltm' and 'scc' task classes, we mix classic
%                           stimulus classes and use crossing names:
%                           'ltm-all','scc-all'. Name can be lower case or
%                           upper case. Name can be a single name
%                           ('fix-gabor') or cell with multiple names
%                           {'scc-all','wm-obj'}). Output will preserve the
%                           order specified by the user's input(s), or
%                           ascending if input is empty,
% 'allstimulusnumbers'	  : Provide all stimulus numbers associated with
%                           given stimulus class name or number. To get all
%                           stimulus numbers of all stimulus classes, use
%                           [] or {}. For stimulus numbers
%                           associated with one stimulus class, the name
%                           should be a string or character class (e.g.,
%                           'rdk'), and should be part of the 5 stimulus
%                           classes: 'gbr','rdk','dot','obj','ns'. Name can
%                           be lower case or upper case. To get stimulus
%                           numbers associated with multiple stimulus
%                           classes, use a cell vector with multiple names
%                           , e.g.: {'rdk','gabor'}).  Output will be
%                           integral numbers, range between 1-1550, and
%                           preserve the order specified by the user's
%                           input (or ascending if input is {}).
% 'stimulusnumberstonames' : Provide stimulus name(s) associated with 
%                           given stimulus number(s). Stimulus number should
%                           be an integral number and range between 1-1550.
%                           Input can be a single number or a vector. Input
%                           cannot be an empty array.
%                           Output is a string or cell vector with strings
%                           in format:
%                            <stimclassname>-<stimnumber>-<stimlocation>.
%                           Note that for stimclassname: "gabor" is
%                           shortened to "GBR". Output will preserve the
%                           order specified by the user's input(s).
% 'stimulusnamestonumbers' : Provide unique stimulus numbers(s) associated 
%                           with stimulus name(s). Input stimulus name(s) 
%                           should be a string in the format:
%                           <stimclassname>-<stimnumber>-<stimlocation>,
%                           and can also be cell array of strings. Note
%                           that for stimclassname: "gabor" is shortened to
%                           "GBR". Output is a vector with image numbers. If
%                           the stimulus class name and/or stimulus
%                           location in the stimulus name does not
%                           correspond to the actual stimulus class and/or
%                           location, the function will throw a warning and
%                           return NaN. Output will preserve the order
%                           specified by the user's input(s). Note that
%                           stimulus names are different from condition
%                           names, as condition names also include the cued
%                           state and task-crossing (i.e.:
%                           <stimclassname>-<stimnumber>-<stimlocation>-<cuestate>-<taskclassname>)
% 'stimtostimclassname'   : Provide stimulus class name(s) for given
%                           stimulus number(s). Input stimulus number(s) should be
%                           be integral and range between 1-1550. Input can
%                           be a single number or a vector, and cannot be empty. Output is a
%                           string or cell vector with stimulus class
%                           names: 'gabor','rdk','dot','obj','ns'. Output
%                           will preserve the order specified by the user's
%                           input(s).
% 'stimtostimclassnumber' : Provide stimulus class number(s) for given
%                           stimulus number(s). Input stimulus number
%                           should be be a single integral number or a
%                           vector of integral numbers that range between
%                           1-1550. Input cannot be an empty vector or
%                           array. Output is a singleton or vector with
%                           stimulus class numbers ranging from 1-5 where
%                           1:'gabor', 2:'rdk',3:'dot',4:'obj',5:'ns'.
%                           Output will preserve the order specified by the
%                           user's input(s).
% 'stimtotaskclassname'   : Provide task class name(s) for given
%                           stimulus number(s). Stimulus number should
%                           be be integral and range between 1-1550. Input
%                           can be a single number or a vector. Output is a
%                           string or cell vector with task class names: 
%                           'fix','cd','scc','pc','wm', 'ltm','img','what',
%                           'where','how'. Output will preserve the order
%                           specified by the user's input(s).
% 'stimtotaskclassnumber' : Provide task class number(s) for given
%                           stimulus number(s). Stimulus number should be
%                           be integral and range between 1-1550. Input can
%                           be a single number or a vector with multiple
%                           numbers. Input cannot be an empty vector or
%                           array. Output is a singleton or vector with
%                           task class numbers. ranging from 1-10, where
%                           1:'fix',
%                           2:'cd',3:'scc',4:'pc',5:'wm',6:'ltm',7:'img',
%                           8:'what', 9:'where',10:'how'. Output will
%                           preserve the order specified by the user's
%                           input(s).
% 'conditionnametonumber' : Get the corresponding condition number for a
%                           given condition name. Condition names should be
%                           char (or a cell with list of char names) in the 
%                           format: 
%                           <stimclassname>-<stimnumber>-<stimlocation>-<cuestate>-<taskclassname>
%                           Output is a vector of integers, preserving the
%                           order specified by the user's input(s). The
%                           first three components use a similar format as
%                           described for "stimulus name", where the
%                           <stimclassname> is either: 'GBR, 'RDK', 'DOT',
%                           'OBJ', 'NS', or 'ALL' (latter is for 'SCC' and
%                           'LTM' task classes only). 
%                           <stimnumber> is an integer between 1-1550.  
%                           <stimlocation> is either: 'L' for left,'R' for
%                           right, 'C' for central (NS only), or 'X' (catch
%                           trial, no stimulus).
%                           <cuestate> is either: 'CUED','UNCUED','NCUED', 
%                           for trials with stimuli. If trial is a catch, 
%                           then there are blanks for stimuli and the 
%                           cuestate refers to the location that is cued: 
%                           'LCUED', 'RCUED', 'NCUED'   
%                           <taskclassname> can be 'FIX, 'CD','SCC','PC',
%                           'WM','LTM','IMG','WHAT','WHERE','HOW'. 
% 'conditionnumbertoname' : Get the corresponding condition name for a
%                           given condition number. Condition numbers should
%                           be integral ranging from 1-1200. Output is a 
%                           char (or a cell with list of char names) in the 
%                           format: 
%                           <stimclassname>-<stimnumber>-<stimlocation>-<cuestate>-<taskclassname>
%                           and preserves the order specified by the user's 
%                           input(s).
% 'allcore'               : Provide all core stimulus number(s) for one (or
%                           more) stimulus class(es). Input names should be
%                           lower case and any of the stimulus class names:
%                           'gabor','rdk','dot','obj','ns'. Output is a
%                           vector of integers, preserving the order
%                           specified by the user's input(s).
% 'specialcore'           : Provide subset of core stimulus number(s) used 
%                           in LTM-[stimclass] and IMG-[stimclass]
%                           crossings. Input names should be lower case and
%                           any of the stimulus class names: 'gabor',
%                           'rdk', 'dot', 'obj', 'ns'. Output is a vector
%                           of integers, preserving the order specified 
%                           by the user's input(s).
% 'allwmtestimages'       : Provide all (non-core) stimulus number(s) of
%                           test images used in WM-[stimclass] crossings.
%                           Input names should be lower case and any of the
%                           stimulus class names:'gabor','rdk','dot','obj',
%                           'ns'. Output is a vector of integers,
%                           preserving the order specified by the user's
%                           input(s).
% 'allimgtestimages'      : Provide all (non-core) stimulus number(s) of 
%                           test images used in LTM-[stimclass] crossings.
%                           Input names should be lower case and any of the
%                           stimulus class names:'gabor','rdk','dot','obj',
%                           'ns'. Output is a vector of integers,
%                           preserving the order specified by the user's
%                           input(s).
% 'allltmtestimages'      : Provide all (non-core) stimulus number(s) of 
%                           test images used in IMG-[stimclass] crossings.
%                           Input names should be lower case and any of the
%                           stimulus class names:'gabor','rdk','dot','obj',
%                           'ns'. Output is a vector of integers,
%                           preserving the order specified by the user's
%                           input(s).
% 'fullinfo'              : Provide all info (data dump) for a given 
%                           stimulus number. Input should be an integral
%                           number between 1-1550.
%
%
% More examples:
% vcd('stimulusclassnames',[3 2 2 2 1])  
% vcd('stimulusclassnames',[])  
% vcd('taskclassnames',[1 5 3 2 2 2 1]) 
% vcd('taskclassnames',{}) 
% vcd('crossingnames',[2 3 4 32])   
% vcd('crossingnumbers',{'pc-gabor' 'scc-all'})
% vcd('stimulusnumberstonames',[12 35 3 2 2 2 1])        
% vcd('taskclassnumbers','FIX')
% vcd('stimulusclassnumber',{'gabor','obj'})
% vcd('allstimulusnumbers',{})
% vcd('stimtostimclassname',12)
% vcd('stimulusnamestonumbers',{'GBR-023-L'})
% vcd('allcore',{'gabor' 'rdk'})                  
% vcd('specialcore','rdk')
% vcd('allwmteststimulusnumbers',{'rdk' 'ns'})
% vcd('allimgteststimulusnumbers','rdk')
% vcd('allltmteststimulusnumbers','ns')
% vcd('allstimulusnumbers','rdk')
% vcd('allstimulusnumbers',{'gabor' 'rdk' 'dot' 'obj' 'ns'})
% vcd('fullinfo',35) 
% vcd('stimulusclassname',[3 2 2 2 1],'taskclassnames',[1 5 3 2 2 2 1]) 
% vcd('stimulusnumberstonames',[10])
% vcd('conditionnumbertoname',[3, 400])
% vcd('conditionnametonumber',{'GBR-0023-L-UNCUED-PC', 'NS-0082-C-NCUED-WHERE'})
%
% Written by Eline Kupers @ UMN 2025/04


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%

% Get experimental and stimulus parameters.
exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);
all_condition_names = vcd_getConditionNames;

p0 = inputParser;

f_stimclassname = @(x) all(ismember(lower(x), exp.stimclassnames));
f_taskclassname = @(x) all(ismember(lower(x), exp.taskclassnames));
f_crossname     = @(x) all(ismember(lower(x), exp.crossingnames));
f_conditionname = @(x) all(vcd_conditionName2Number(x)~=0);

% General
p0.addParameter('verbose'               , true, @islogical);  

% Broad stimulus class & task class names and numbers
p0.addParameter('stimulusclassnames'    , [], @(x) choose_isempty(x) || isstimclassnr(x));    
p0.addParameter('stimulusclassnumbers'  , [], @(x) choose_isempty(x) || f_stimclassname(x)); 
p0.addParameter('taskclassnames'        , [], @(x) choose_isempty(x) || istaskclassnr(x));           
p0.addParameter('taskclassnumbers'      , [], @(x) choose_isempty(x) || f_taskclassname(x)); 
p0.addParameter('crossingnames'         , [], @(x) choose_isempty(x) || iscrossnr(x));           
p0.addParameter('crossingnumbers'       , [], @(x) choose_isempty(x) || f_crossname(x));  

% Stimulus specific
p0.addParameter('allstimulusnumbers'    , [], @(x) choose_isempty(x) || (isstimclassnr(x) || f_stimclassname(x)));
p0.addParameter('stimulusnumberstonames', [], @(x) isstimnr(x));   
p0.addParameter('stimulusnamestonumbers', [], @(x) all(iscell(x)) || ischar(x));
p0.addParameter('allcore'               , [], @(x) choose_isempty(x) || f_stimclassname(x));                 
p0.addParameter('specialcore'           , [], @(x) choose_isempty(x) || f_stimclassname(x));                 
p0.addParameter('stimtostimclassname'   , [], @(x) isstimnr(x));   
p0.addParameter('stimtostimclassnumber' , [], @(x) isstimnr(x));   
p0.addParameter('stimtotaskclassname'   , [], @(x) isstimnr(x));   
p0.addParameter('stimtotaskclassnumber' , [], @(x) isstimnr(x));   
p0.addParameter('fullinfo'              , [], @(x) isstimnr(x)); 

% Condition specific
p0.addParameter('conditionnumbertoname' , [], @(x) isconditionnr(x));
p0.addParameter('conditionnametonumber' , [], @(x) f_conditionname(x));

% Test image specific (only for WM, LTM, IMG task crossings)
p0.addParameter('allwmteststimulusnumbers' , [], @(x) choose_isempty(x) || f_stimclassname(x));   
p0.addParameter('allltmteststimulusnumbers', [], @(x) choose_isempty(x) || f_stimclassname(x));  
p0.addParameter('allimgteststimulusnumbers', [], @(x) choose_isempty(x) || f_stimclassname(x));  

% Parse inputs
p0.parse(varargin{:});
verbose     = p0.Results.verbose;
% displayname = p0.Results.displayname;

% Check which info was requested
requested_info = {};
for ivar = 1:length(p0.Parameters)
    if ~ismember(p0.Parameters(ivar),p0.UsingDefaults)
        requested_info{length(requested_info)+1} = p0.Parameters{ivar};  %#ok<AGROW>
        requested_info{length(requested_info)+1} = p0.Results.(p0.Parameters{ivar}); %#ok<AGROW>
    end
end
 
% Check if info table is cached, if not, we load it
global vcd_info;

if isempty(vcd_info) || ~exist('vcd_info','var')
    if ~exist('vcd_rootPath','file')
        error('[%s]: Please navigate to vcd-stim folder',mfilename)
    end
    
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('condition_master_*.mat')));
    dcell = struct2cell(d); filenames = dcell(1,:); clear dcell
    is_DEEP_MRI_file  = find(~cellfun(@isempty, (cellfun(@(x) regexp(x,'condition_master_deep_7TAS\w*','ONCE'), filenames, 'UniformOutput', false))));
    is_WIDE_MRI_file = find(~cellfun(@isempty, (cellfun(@(x) regexp(x,'condition_master_wide_7TAS\w*','ONCE'), filenames, 'UniformOutput', false))));
    is_BEH_file  = find(~cellfun(@isempty, (cellfun(@(x) regexp(x,'condition_master_PPROOM\w*','ONCE'), filenames, 'UniformOutput', false))));    
    is_DEMO_file = find(~cellfun(@isempty, (cellfun(@(x) regexp(x,'condition_master_demo\w*','ONCE'), filenames, 'UniformOutput', false))));    
    
    % remove demo file from list of condition_master files; we don't want to use that for vcd.m
    if ~isempty(is_DEMO_file), filenames(is_DEMO_file) = {NaN}; end

    % Check if we have a condition_master filename to load
    if isempty(filenames{1}) || all(cell2mat(cellfun(@(x) isequalwithequalnans(x,NaN), filenames, 'UniformOutput', false)))
        error('[%s]: Can''t find vcd_info! Looking for a file called: "../workspaces/info/condition_master*.mat"',mfilename)
    else
        % if we have multiple condition_master files, take the most recent one.
        if length(filenames) > 1
            if verbose
                fprintf('[%s]: *** WARNING *** \n', mfilename);
                fprintf('[%s]: Found %d vcd info files with the same name! Will pick the most recent one. (We ignore demo files.) \n', mfilename,length(filenames));
            end
            % if we have a MRI condition_master file, we pick MRI file, as
            % it is more accurate, otherwise we use PPROOM
            if ~isempty(is_DEEP_MRI_file)
                d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('condition_master_deep_7TAS*.mat')));
            elseif ~isempty(is_WIDE_MRI_file)
                d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('condition_master_wide_7TAS*.mat')));
            elseif ~isempty(is_BEH_file)
                d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('condition_master_PPROOM*.mat')));
            end
            % pick the most recent one
            fname = d(end).name;
        end
        
        % ensure we have a filename to load
        assert(~isempty(fname))

        if verbose
            fprintf('[%s]: Loading condition_master .mat file: %s\n', mfilename, fname);
        end
        vcd_info = load(fullfile(d(end).folder,fname),'condition_master');
    end
end

% If we haven't done so already, create a vector of stimulus class numbers or names for each unique image nr
% *** Number of stimuli and their stimulus numbers: 
% 	 1710 total (001-1710) 
% 	 110 core (001-110) 
% 	 312 WM test (111-422) 
% 	 940 IMG test (423-1362) 
% 	 60 LTM novel lures (1363-1422) 
% 	 288 OBJ catch images (1423-1710) 
% 	 47 special core
all_core_im_stimclassnumbers       = ...
    cat(2, ones(1, numel(stim.gabor.unique_im_nrs_core)), 2*ones(1, numel(stim.rdk.unique_im_nrs_core)), ...
    3*ones(1, numel(stim.dot.unique_im_nrs_core)), 4*ones(1, numel(stim.obj.unique_im_nrs_core)), ...
    5*ones(1, numel(stim.ns.unique_im_nrs_core)));
all_wm_test_im_stimclassnumbers    = ...
    cat(2, ones(1, numel(stim.gabor.unique_im_nrs_wm_test)), 2*ones(1, numel(stim.rdk.unique_im_nrs_wm_test)), ...
    3*ones(1, numel(stim.dot.unique_im_nrs_wm_test)), 4*ones(1, numel(stim.obj.unique_im_nrs_wm_test)), ...
    5*ones(1, numel(stim.ns.unique_im_nrs_wm_test)));
all_img_test_im_stimclassnumbers   = ...
    cat(2, ones(1, numel(stim.gabor.unique_im_nrs_img_test)), 2*ones(1, numel(stim.rdk.unique_im_nrs_img_test)), ...
    3*ones(1, numel(stim.dot.unique_im_nrs_img_test)), 4*ones(1, numel(stim.obj.unique_im_nrs_img_test)), ...
    5*ones(1, numel(stim.ns.unique_im_nrs_img_test)));
all_ltm_lure_im_stimclassnumbers    = 5*ones(1, numel(stim.ns.unique_im_nrs_novel_ltm_lures));
all_objectcatch_im_stimclassnumbers = 4*ones(1, numel(stim.obj.unique_im_nrs_objcatch));
all_im_stimclassnumbers            = cat(2,all_core_im_stimclassnumbers,all_wm_test_im_stimclassnumbers,all_img_test_im_stimclassnumbers,all_objectcatch_im_stimclassnumbers);

% Stimulus class names
all_core_im_stimclassnames         = exp.stimclassnames(all_core_im_stimclassnumbers);
all_wm_test_im_stimclassnames      = exp.stimclassnames(all_wm_test_im_stimclassnumbers);
all_img_test_im_stimclassnames     = exp.stimclassnames(all_img_test_im_stimclassnumbers);
all_ltm_lure_im_stimclassnames     = exp.stimclassnames(all_ltm_lure_im_stimclassnumbers);
all_objectcatch_im_stimclassnames  = exp.stimclassnames(all_objectcatch_im_stimclassnumbers);
all_im_stimclassnames              = exp.stimclassnames(all_im_stimclassnumbers);
 
assert(isequal(length(stim.all_im_nrs),length(all_im_stimclassnames)))

% Now check requested_info by user
for ii = 1:2:length(requested_info)
    
    info_name = requested_info{ii};
    info_val  = requested_info{ii+1};
    
    switch info_name
        
        case 'stimulusclassnames'
            if isempty(info_val) 
                out = exp.stimclassnames; % {'gabor','rdk','dot','obj','ns'};
            elseif iscell(info_val) && strcmp(info_val{:},'all')
                out = exp.stimclassnames; % {'gabor','rdk','dot','obj','ns'};
            elseif iscell(info_val) && all(cellfun(@isnumeric, info_val))
                info_val = cell2mat(info_val);
                out = exp.stimclassnames(info_val);
            elseif isnumeric(info_val)
                out = exp.stimclassnames(info_val);
            end
            
        case 'stimulusclassnumbers'
                        
            if choose_isempty(info_val)                                % provide all stimulus class numbers: 1-5
                out = 1:length(exp.stimclassnames); 
                
            elseif ischar(info_val)
                info_val = lower(info_val);                     % use lower case characters
                
                if strcmp(info_val,'all')                       % provide all stimulus class numbers: 1-5
                    out = 1:length(exp.stimclassnames); 
                elseif any(strcmp(info_val,exp.stimclassnames)) % provide specific stimulus class number(s)
                    [~,out] = ismember(info_val,exp.stimclassnames);
                end
 
            elseif iscell(info_val)
                if strcmp(info_val,'all')                    % provide all stimulus class numbers: 1-5
                    out = 1:length(exp.stimclassnames); 
                else                                            % provide specific stimulus class number(s)
                    [~,out] = ismember(lower(info_val(:)),exp.stimclassnames);
                end
            end
            
            if size(out,1)>size(out,2)
                out = out';
            end

            
        case 'taskclassnames'
            
            if choose_isempty(info_val)
                out = exp.taskclassnames;                    % provide all task class names: {'fix','cd','scc','pc','wm','ltm','img','what','where','how'};
            elseif iscell(info_val) && all(cellfun(@isnumeric, info_val))
                info_val = cell2mat(info_val);
                out = exp.taskclassnames(info_val);    
            elseif isnumeric(info_val)                          % provide specific task class name(s)
                out = exp.taskclassnames(info_val);
            end
            
        case 'taskclassnumbers'
            
            if choose_isempty(info_val)                               % provide all task class numbers: 1-10
                out = 1:length(exp.taskclassnames);         
                
            elseif ischar(info_val)
                info_val = lower(info_val); % use lower case characters
                
                if strcmp(info_val,'all')                      % provide all task class numbers: 1-10
                    out = 1:length(exp.taskclassnames);      
                elseif any(strcmp(info_val,exp.taskclassnames)) % provide specific task class number(s)
                    [~,out] = ismember(info_val,exp.taskclassnames);
                end  
            elseif iscell(info_val)
                if strcmp(info_val(:),'all')                    % provide all task class numbers: 1-10
                    out = 1:length(exp.taskclassnames);
                else                                            % provide specific task class number(s)
                    [~,out] = ismember(lower(info_val(:)),exp.taskclassnames);
                end
            end
            if size(out,1) > size(out,2)
                out = out';
            end
            
        case 'crossingnames'
            
            if choose_isempty(info_val)                % provide all task-stimulus class crossing names: {'fix-gbr',...,'how-ns'};
                out = exp.crossingnames; 
            elseif iscell(info_val) && all(cellfun(@isnumeric, info_val)) 
                info_val = cell2mat(info_val);
                out = exp.crossingnames(info_val);    
            elseif isnumeric(info_val)           % provide specific task-stimulus class crossing name(s)
                out = exp.crossingnames(info_val);
            end
            if size(out,1) > size(out,2)
                out = out';
            end
            
        case 'crossingnumbers'
            
            if choose_isempty(info_val) % provide all task-stimulus class crossing numbers (1-32), if user doesn't specify
                out = 1:length(exp.crossingnames);           
                
            elseif ischar(info_val)
                info_val = lower(info_val); % use lower case characters
                
                if strcmp(info_val,'all')
                    out = 1:length(exp.crossingnames);       % provide all task-stimulus class crossing numbers: 1-32
                    
                elseif any(strcmp(info_val,exp.crossingnames))  % provide specific task-stimulus class crossing number(s)
                    [~,out] = ismember(info_val,exp.crossingnames);
                end
                
            elseif iscell(info_val)
                if strcmp(info_val(:),'all')                                % provide all task-stimulus class crossing numbers: 1-32
                    out = 1:length(exp.crossingnames);
                else                                                        % provide specific task-stimulus class crossing number(s)
                    [~,out] = ismember(lower(info_val(:)), exp.crossingnames);
                end
            end
            
            if size(out,1) > size(out,2)
                out = out';
            end
            
        case 'allstimulusnumbers'
            % provide all stimulus numbers associated with given stimulus class name or number
            if choose_isempty(info_val) % provide all image nrs 1-1550, if user doesn't specifiy
                out = stim.all_im_nrs;
            else
                if iscell(info_val)
                    if all(cellfun(@isnumeric, info_val))
                        info_val = cell2mat(info_val);
                    elseif all(cellfun(@ischar, info_val))
                        % do nothing
                    end
                end
                
                if isnumeric(info_val) % if user provides stim class number
                    [~,idx] = ismember(info_val,1:length(exp.stimclassnames),'legacy');
                    stimclassname = exp.stimclassnames(idx);
                elseif ischar(info_val) % if user provides stim class name 
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names 
                    stimclassname = info_val(:);
                end
                
                out = [];
                for sc = 1:length(stimclassname)
                    if strcmp(stimclassname{sc},'ns')
                        out = cat(2, out, stim.ns.unique_im_nrs_core, stim.ns.unique_im_nrs_wm_test, stim.ns.unique_im_nrs_img_test,stim.ns.unique_im_nrs_ltm_lures);
                    else
                        out = cat(2, out, stim.(stimclassname{sc}).unique_im_nrs_core, stim.(stimclassname{sc}).unique_im_nrs_wm_test, stim.(stimclassname{sc}).unique_im_nrs_img_test);
                    end
                end
            end
            
        case 'stimulusnumberstonames'
            % provide stimulus name(s) associated with given stimulus number(s).
            
            % use cell vector if we deal with multiple stimulus numbers
            if length(info_val) > 1
                out = cell(1,length(info_val));
            else
                out = [];
            end
            
            if isnumeric(info_val) % if user provides singleton or vector
                stim_nr = info_val;

            elseif iscell(info_val)  % if user provides cell
                stim_nr = cell2mat(info_val);
            end
            
            for nn = 1:length(stim_nr)
                
                % if we deal with a core stimulus
                if ismember(stim_nr(nn), stim.all_core_im_nrs)  
                    
                    [~,idx] = ismember(stim_nr(nn),vcd_info.condition_master.stim_nr_left,'legacy');
                    
                    if idx == 0
                        [~,idx] = ismember(stim_nr(nn), vcd_info.condition_master.stim_nr_right,'legacy');
                        im_loc = 2; % right stim
                        stimnumber   = vcd_info.condition_master.stim_nr_right(idx);
                        stimlocation = 'R';
                    else
                        im_loc = 1; % left stim
                        stimnumber = vcd_info.condition_master.stim_nr_left(idx);
                        if vcd_info.condition_master.stim_class(idx) == 5
                            stimlocation  = 'C'; % central stim
                        else
                            stimlocation = 'L';
                        end
                    end
                    
                    

                elseif ismember(stim_nr(nn), stim.all_test_im_nrs)  
                    % find unique stimulus number for test image
                    [~,idx] = ismember(stim_nr(nn),vcd_info.condition_master.stim2(:,1),'legacy');
                    if idx == 0
                        [~,idx] = ismember(stim_nr(nn), vcd_info.condition_master.stim2(:,2),'legacy');
                        im_loc = 2; % right stim
                        stimnumber   = vcd_info.condition_master.stim_nr_right(idx);
                        stimlocation = 'R';
                    else
                        im_loc = 1; % left stim
                        stimnumber = vcd_info.condition_master.stim_nr_left(idx);
                        if vcd_info.condition_master.stim_class(idx) == 5
                            stimlocation  = 'C'; % central stim
                        else
                            stimlocation = 'L';
                        end
                    end
                end
                stimclassname = upper(vcd_info.condition_master.stim_class_name{idx,im_loc});
                if strcmp(stimclassname,'GABOR')
                    stimclassname = regexprep('GABOR','[AO]','');
                end
                stimnumber = num2str(stimnumber);
                
                if ischar(stimnumber)
                    stimnumber = sprintf('%03d',str2double(stimnumber));
                elseif isnumeric(stimnumber)
                    stimnumber = sprintf('%03d',stimnumber);
                end
                if iscell(out)
                    out{nn} = sprintf('%s-%s-%s',stimclassname,stimnumber,stimlocation);
                else
                    out = sprintf('%s-%s-%s',stimclassname,stimnumber,stimlocation);
                end
            end
            
        case 'stimulusnamestonumbers'
            % provide stimulus number(s) associated with given stimulus name(s).
            if ischar(info_val)
                info_val = {info_val};
            end 
            
            out = zeros(1,length(info_val));
            
            for nn = 1:length(info_val)    
                
                % extract stim nr from string
                input_name0 = strsplit(info_val{nn},'-');
                input_stim_class_name = lower(input_name0{1});
                input_stim_nr         = str2double(input_name0{2});
                input_stim_pos        = input_name0{3};
                
                % check if stim nr matches stim class and stim location
                [~,idx] = ismember(input_stim_nr,stim.all_im_nrs,'legacy');
                true_stim_class_name = all_im_stimclassnames(idx);
                
                if ismember(true_stim_class_name,'gabor')
                    true_stim_class_name = 'GBR';
                end
                
                if ~ismember(input_stim_class_name,true_stim_class_name)
                    warning('[%s]: Stimulus class name "%s" is not associated with this stimulus number %03d!',mfilename,input_stim_class_name,input_stim_nr)
                    out(nn) = NaN;
                else
                    
                    % Check position of stimulus
                    [~,idx_l] = ismember(input_stim_nr, vcd_info.condition_master.stim_nr_left,'legacy');
                    [~,idx_r] = ismember(input_stim_nr, vcd_info.condition_master.stim_nr_right,'legacy');

                    if vcd_info.condition_master.stim_class(idx_l,:)==5
                        true_stim_pos = 'C'; % central stim
                    else
                        if (idx_l == 0) && (idx_r > 0) % if right stimulus nr matches
                            true_stim_pos = 'R'; % right stim
                        elseif (idx_l > 0) && (idx_r == 0) % if left stimulus nr matches
                            true_stim_pos = 'L'; % left stim
                        elseif (idx_l == 0) && (idx_r == 0) % if none matches
                            true_stim_pos = NaN;
                        end
                    end

                    if ~strcmp(input_stim_pos,true_stim_pos) || isnan(true_stim_pos)
                        warning('[%s]: Stimulus position "%s" is not associated with stimulus number %03d!',mfilename,input_stim_pos,input_stim_nr)
                        out(nn) = NaN;
                    end
                end

                if out(nn) == 0
                    out(nn) = input_stim_nr;
                end
            end   
            
        case 'allcore' 
           % provide all core stimulus number(s) for one (or more) stimulus class(es)
           if choose_isempty(info_val) % provide all coreimage nrs (1-110), if user doesn't specifiy stimulus class
                out = stim.all_core_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    [~,idx] = ismember(info_val,1:length(exp.stimclassnames),'legacy');
                    stimclassname = exp.stimclassnames(idx);
                elseif ischar(info_val) % if user provides stim class name 
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names 
                    stimclassname = info_val;
                end
                
                if size(stimclassname,1) > size(stimclassname,2)
                    stimclassname = stimclassname';
                end
                
                out = [];
                for sc = 1:length(stimclassname)
                    out = cat(2, out, stim.(stimclassname{sc}).unique_im_nrs_core);
                end
           end 
            
        case 'specialcore'
            % provide all special core stimulus number(s) for one (or more) stimulus class(es)
            if isempty(info_val) % provide all special 47 core image nrs, if user doesn't specifiy stimulus class
                out = stim.all_specialcore_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    [~,idx] = ismember(info_val,1:length(exp.stimclassnames),'legacy');
                    stimclassname = exp.stimclassnames(idx);
                elseif ischar(info_val) % if user provides stim class name
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names
                    stimclassname = info_val;
                end
                
                if size(stimclassname,1) > size(stimclassname,2)
                    stimclassname = stimclassname';
                end
                
                out = [];
                for sc = 1:length(stimclassname)
                    out = cat(2, out, stim.(stimclassname{sc}).unique_im_nrs_specialcore);
                end
            end
            
        case 'stimtostimclassname'
            % translate stimulus number(s) to stimulus class name(s)
            
            if isnumeric(info_val) % if user provides singleton or vector
                stim_nr = info_val;

            elseif iscell(info_val)  % if user provides cell
                stim_nr = cell2mat(info_val);
            end
            if numel(stim_nr)==1
                out = all_im_stimclassnames{stim_nr};
            else
                out = all_im_stimclassnames(stim_nr);
            end
            % if needed, transpose to mantain size
            if ~isequal(size(info_val,1) > size(info_val,2),size(out,1) > size(out,2))
                out = out';
            end
        case 'stimtostimclassnumber'
            % translate stimulus number(s) to stimulus class number(s)
            if isnumeric(info_val) % if user provides singleton or vector
                stim_nr = info_val;
            elseif iscell(info_val)  % convert cell to vector if user provides cell
                stim_nr = cell2mat(info_val);
            end
            
            out = all_im_stimclassnumbers(stim_nr);
            
        case 'stimtotaskclassname'   
            % translate stimulus number(s) to task class name(s)
            if isnumeric(info_val) % if user provides single stim nr
                stim_nr = info_val;
            elseif iscell(info_val)  % convert cell to vector if user provides cell
                stim_nr = cell2mat(info_val);
            end
            
            % get stimulus class nr
            stim_class_nr = all_im_stimclassnumbers(stim_nr);
            
            if length(stim_class_nr) == 1
                out = exp.taskclassnames(exp.crossings(stim_class_nr,:));
            else
                out = cell(1,length(stim_class_nr));
                for nn = 1:length(stim_class_nr)
                    out{nn} = exp.taskclassnames(exp.crossings(stim_class_nr(nn),:));
                end
            end            
        case 'stimtotaskclassnumber' 
            % translate stimulus number(s) to task class number(s)
            if isnumeric(info_val) % if user provides single stim nr
                stim_nr = info_val;
            elseif iscell(info_val)  % convert cell to vector if user provides cell
                stim_nr = cell2mat(info_val);
            end
            
            % get stimulus class nr
            out = all_im_stimclassnumbers(stim_nr);
            
        case 'allwmteststimulusnumbers'       
            % provide all (non-core) stimulus number(s) of test images used in WM-[stimclass] crossings 
           if choose_isempty(info_val) % provide all wm test image nrs (1-110), if user doesn't specifiy stimulus class
                out = stim.all_wm_test_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    [~,idx] = ismember(info_val,1:length(exp.stimclassnames),'legacy');
                    stimclassname = exp.stimclassnames(idx);
                elseif ischar(info_val) % if user provides stim class name 
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names 
                    stimclassname = info_val;
                end
                
                out = [];
                for sc = 1:length(stimclassname)
                    if isfield(stim.(stimclassname{sc}),'unique_im_nrs_wm_test')
                        out = cat(2, out, stim.(stimclassname{sc}).unique_im_nrs_wm_test);
                    end
                end
           end 
            
            
        case 'allltmteststimulusnumbers'      
            % provide all (non-core) stimulus number(s) of test images used in LTM-[stimclass] crossings 
            if isempty(info_val) % provide all wm test image nrs (1-110), if user doesn't specifiy stimulus class
                out = stim.all_ltm_lure_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    [~,idx] = ismember(info_val,1:length(exp.stimclassnames),'legacy');
                    stimclassname = exp.stimclassnames(idx);
                elseif ischar(info_val) % if user provides stim class name
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names
                    stimclassname = info_val(:);
                end
                
                out = [];
                for sc = 1:length(stimclassname)
                    if isfield(stim.(stimclassname{sc}), 'unique_im_nrs_ltm_lures')
                        out = cat(2, out, stim.(stimclassname{sc}).unique_im_nrs_ltm_lures);
                    end
                end
            end
            
        case 'allimgteststimulusnumbers'      
            % provide all (non-core) stimulus number(s) of test images used in IMG-[stimclass] crossings 
            if isempty(info_val) % provide all wm test image nrs (1-110), if user doesn't specifiy stimulus class
                out = stim.all_img_test_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    [~,idx] = ismember(info_val,1:length(exp.stimclassnames),'legacy');
                    stimclassname = exp.stimclassnames(idx);
                elseif ischar(info_val) % if user provides stim class name
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names
                    stimclassname = info_val(:);
                end
                
                out = [];
                for sc = 1:length(stimclassname)
                    if isfield(stim.(stimclassname{sc}), 'unique_im_nrs_img_test')
                        out = cat(2, out, stim.(stimclassname{sc}).unique_im_nrs_img_test);
                    end
                end
            end

        case 'fullinfo'              
            % provide all info (data dump) for given stimulus number (1-1421)
            if isnumeric(info_val) % if user provides singleton or vector
                stim_nr = info_val;
            elseif iscell(info_val)  % if user provides cell, turn into vector
                stim_nr = info_val{:};
            end
            
            % prepare output
            if length(stim_nr) > 1
                out = [];
            else
                out = [];
            end
            
            for nn = 1:length(stim_nr)
                % preallocate stim_info struct for each stimulus nr
                stim_info = struct('condition_nr',[], 'condition_name',[],...
                    'stim_class', [], 'stim_class_name', [], 'task_class', [], ...
                    'task_class_name', [], 'crossing_nr', [], 'crossing_name', [], ...
                    'stim_nr', [], 'stim_loc', [], 'stim_loc_name', [],'is_objectcatch', [], ...
                    'orient_dir', [], 'contrast', [], 'gbr_phase', [], ...
                    'rdk_coherence', [], 'super_cat', [], 'super_cat_name', [], ...
                    'basic_cat', [], 'basic_cat_name', [], 'sub_cat', [], ...
                    'sub_cat_name', [], 'affordance_cat', [], 'affordance_name', [], ...
                    'is_special_core', [],'is_lure', []);
                
                % if we deal with a core stimulus
                if ismember(stim_nr(nn), stim.all_core_im_nrs,'legacy')
                    % find unique stim number for left core image 
                    [~,idx] = ismember(stim_nr(nn),vcd_info.condition_master.stim_nr_left,'legacy');
                    if isempty(idx) || idx == 0 % if we can't find it, we assume this is a stimulus on the right
                        [~,idx] = ismember(stim_nr(nn), vcd_info.condition_master.stim_nr_right,'legacy');
                        im_loc = 2; % right stim
                        stim_info.stim_nr = vcd_info.condition_master.stim_nr_right(idx);
                        
                        assert(isequal(stim_nr(nn),vcd_info.condition_master.stim_nr_right(idx)))
                        assert(~isempty(regexp(vcd_info.condition_master.condition_name{idx,2},'\w+-R-\w+')))
                    else
                        im_loc = 1; % left/central stim
                        stim_info.stim_nr = vcd_info.condition_master.stim_nr_left(idx);
                        
                        assert(isequal(stim_nr(nn),vcd_info.condition_master.stim_nr_left(idx)))
                        assert(~isempty(regexp(vcd_info.condition_master.condition_name{idx,1},'\w+-[LC]-\w+')))
                    end
                % if we deal with a test image (i.e., wm test, img test, or ltm lure)
                elseif ismember(stim_nr(nn), stim.all_test_im_nrs)
                    % find unique stim number for left test image 
                    [~,idx] = ismember(stim_nr(nn),vcd_info.condition_master.stim2_im_nr(:,1),'legacy');
                    if isempty(idx) || idx == 0 % if we can't find it, we assume this is a test stimulus on the right
                        [~,idx] = ismember(stim_nr(nn), vcd_info.condition_master.stim2_im_nr(:,2),'legacy');
                        im_loc = 2; % right stim
                        stim_info.stim_nr = vcd_info.condition_master.stim2_im_nr(idx,2);
                        
                        assert(isequal(stim_nr(nn),vcd_info.condition_master.stim2_im_nr(idx,2)))
                        assert(~isempty(regexp(vcd_info.condition_master.condition_name{idx,2},'\w+-R-\w+')))
                    else
                        im_loc = 1; % left stim
                        assert(isequal(stim_nr(nn),vcd_info.condition_master.stim2_im_nr(idx,1)))
                        stim_info.stim_nr = vcd_info.condition_master.stim2_im_nr(idx,1);
                        assert(~isempty(regexp(vcd_info.condition_master.condition_name{idx,1},'\w+-[LC]-\w+')))
                    end
                % if we deal with an object catch image    
                elseif ismember(stim_nr(nn),stim.all_objectcatch_im_nrs)
                    % find unique stim number for test image
                    [~,idx] = ismember(stim_nr(nn),vcd_info.condition_master.stim2_im_nr(:,1),'legacy');
                    if isempty(idx) || idx == 0 % if we can't find it, we assume this is an object catch stimulus on the right
                        [~,idx] = ismember(stim_nr(nn), vcd_info.condition_master.stim2_im_nr(:,2),'legacy');
                        im_loc = 2; % right stim
                        stim_info.stim_nr = vcd_info.condition_master.stim_nr_right(idx);
                        
                        assert(isequal(stim_nr(nn),vcd_info.condition_master.stim_nr_right(idx)))
                        assert(~isempty(regexp(vcd_info.condition_master.condition_name{idx,2},'\w+-R-\w+')))
                    else
                        im_loc = 1; % left stim
                        stim_info.stim_nr = vcd_info.condition_master.stim_nr_left(idx);
                        
                        assert(isequal(stim_nr(nn),vcd_info.condition_master.stim_nr_left(idx)))
                        assert(~isempty(regexp(vcd_info.condition_master.condition_name{idx,1},'\w+-[LC]-\w+')))
                    end
                else
                    error('[%s]: Not sure what type of image this is..',mfilename);
                end
                
                % If stimulus number exists in condition_master..
                if idx~=0
                    % get stim info
                    tmp_info = vcd_info.condition_master(idx,:);
                    
                    % add stim loc
                    
                    % Check if we deal with a NS..
                    if any(ismember(stim_info.stim_nr,[stim.ns.unique_im_nrs_core, ...
                                stim.ns.unique_im_nrs_wm_test,stim.ns.unique_im_nrs_novel_ltm_lures, stim.ns.unique_im_nrs_img_test])) ...
                                || (~isempty(stim_info.stim_class) && stim_info.stim_class == 5)
                        stim_info.stim_loc = 3;
                        stim_info.stim_loc_name = 'center';
                        assert(isnan(tmp_info.condition_nr(2)));
                    else
                        stim_info.stim_loc = im_loc;
                        stim_info.stim_loc_name = choose(im_loc==1,'left','right');
                    end
                    
                    % loop over fieldnames over stim info that we haven't
                    % filled yet..
                    fn = fieldnames(stim_info);
                    fn = fn(~ismember(fn,{'stim_nr','stim_loc','stim_loc_name',...
                        'condition_nr','task_class_name','crossing_name'}));
                    for ff = 1:length(fn)
                        if strcmp(fn{ff},'condition_name')
                            tmp_condname0 = tmp_info.(fn{ff});
                            tmp_condname1 = strsplit(tmp_condname0{1},'-');
                            tmp_condname2 = [tmp_condname1{1} '-'  tmp_condname1{2} '-' tmp_condname1{3}];
                            tmp_idx = ~cellfun(@isempty, regexp(all_condition_names,[tmp_condname2 '\w*']));
                            stim_info.(fn{ff}) = cat(2,all_condition_names(tmp_idx));
                            stim_info.condition_nr = vcd_conditionName2Number(stim_info.(fn{ff}));
                        elseif strcmp(fn{ff},'task_class')
                            % check for mixed stim class 
                            if isempty(stim_info.stim_class) % if we haven't defined this field yet
                                if tmp_info.stim_class == 99
                                    stim_idx = find(ismember(exp.stimclassnames,tmp_info.stim_class_name(stim_info.stim_loc)));   
                                else
                                    stim_idx = tmp_info.stim_class;
                                end
                            else % we use it and assume it is correct
                                if stim_info.stim_class == 99
                                    stim_idx = find(ismember(exp.stimclassnames,stim_info.stim_class_name));   
                                else
                                    stim_idx = stim_info.stim_class;
                                end
                            end
                            stim_info.task_class = find(exp.crossings(stim_idx,:));
                            stim_info.task_class_name = exp.taskclassnames(exp.crossings(stim_idx,:));
                            
                            % check if stim is not a special core, then we
                            % remove LTM/IMG crossings
                            if ~ismember(stim_info.stim_nr, stim.all_specialcore_im_nrs)
                                if any(ismember(stim_info.task_class,[6,7]))
                                    stim_info.task_class(ismember(stim_info.task_class,[6,7])) = [];
                                    stim_info.task_class_name(ismember(stim_info.task_class,[6,7])) = [];
                                end
                            end
                        elseif strcmp(fn{ff},'stim_class')
                            if tmp_info.stim_class == 99
                                stim_idx = find(ismember(exp.stimclassnames,tmp_info.stim_class_name(stim_info.stim_loc)));   
                            else
                                stim_idx = tmp_info.stim_class;
                            end
                            stim_info.stim_class = stim_idx;
                            stim_info.stim_class_name = exp.stimclassnames(stim_idx);
  
                        elseif strcmp(fn{ff},'crossing_nr')
                            stim_info.crossing_nr = find(~cellfun(@isempty, regexp(exp.crossingnames,[stim_info.stim_class_name{1} '\w*'])));
                            stim_info.crossing_name = exp.crossingnames(stim_info.crossing_nr);
                        else
                            if size(tmp_info.(fn{ff}),2) > 1
                                stim_info.(fn{ff}) = tmp_info.(fn{ff})(im_loc);
                            else
                                stim_info.(fn{ff}) = tmp_info.(fn{ff});
                            end
                        end
                    end
                    
                    if length(stim_nr) > 1
                        out = [out,stim_info];
                    else
                        out = stim_info;
                    end
                end
            end
            
        case 'conditionnumbertoname'
            % provide condition number(s) for given condition name(s) 
            if isnumeric(info_val) % if user provides single integer or vector of integers
                cond_nr = info_val;
            elseif iscell(info_val)  % if user provides cell, turn into vector
                if ischar(info_val{1})
                    cond_nr = str2double(char2mat(info_val))';
                else
                    cond_nr = info_val{:};
                end
            end
           
            out = vcd_conditionNumber2Name(cond_nr);
            
        case 'conditionnametonumber'
            if ischar(info_val)
                info_val = {info_val};
            end 
            
            out = vcd_conditionName2Number(info_val);
            
    end % end of switch/case statement
    
    % Check outputs
    if length(requested_info) > 2
        if ii == 1
            % create output struct if we have more than one info request
            outputs = struct();
        end
        % Add output to outputs struct
        outputs.(info_name) = out;
    else
        outputs = out;
    end
end

 
return


