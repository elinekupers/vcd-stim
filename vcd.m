function outputs = vcd(varargin)
% VCD function to call general information about experimental design, 
% including the number of stimulus and task class names, numbers, unique
% image number, full information for each unique stimulus.
%
%   outputs = vcd('varname1',var1,'varname2',var2, ...'varnameN',varN);
%
% Note that 'fullinfo', 'stimulusnumberstonames', and
% 'stimulusnamestonumbers' requires the condition_master table. This
% variable is stored here: 
%   vcd_rootPath/workspaces/info/trials_7TASBOLDSCREEN_YYYYMMDDTHHMMSS.mat.
% You can re-generate a new trials_*.mat file by running (see also
% s_createDesignMatrix.m.):
% [~, condition_master] = vcd_createBlocksAndTrials(params, ...
%   'load_params', false, 'store_params', true);
% 
% Abbreviations for stimulus classes are: 
%  1:'gabor'  - Gabor (grating subject to Gaussian window).
%  2:'rdk'    - Random dot motion kinetogram
%  3:'dot'    - Single dot
%  4:'obj'    - Objects
%  5:'ns'     - Natural scenes
%
% Abbreviations for task classes are: 
%   1:'fix'   - Fixation brightness task
%   2:'cd'    - Contrast change task
%   3:'scc'   - Categorization task
%   4:'pc'    - Perceptual categorization tasks (tilt, motion direction, dot position, object rotation, indoor/outdoor)
%   5:'wm'    - Working memory tasks (tilt, motion direction, position, rotation, scene)
%   6:'ltm'   - Matching task
%   7:'img'   - Imagery task
%   8:'what'  - "What?" task
%   9:'where' - "Where?" task
%   10:'how'  - "How?" task
%
% Possible inputs:
% 'stimulusclassnames'    : Provide all stimulus class names (use []), or
%                           translate stimulus class number(s) to name(s)
%                           Stimulus class number should be an integral 
%                           number between 1 and 5, can be a single number
%                           or vector. Output will be sorted by order of 
%                           input nr, or ascending if input is empty.
% 'stimulusclassnumbers'  : Provide all stimulus class numbers (use []), or
%                           translate stimulus class name(s) to number(s).
%                           Stimulus class name should be a string or
%                           character class, and should be part of the 5 
%                           stimulus classes, where 1:'gabor', 2:'rdk',
%                           3:'dot',4:'obj',5:'ns'.
%                           Name can be lower case or upper case, and can 
%                           be a single string ('rdk') or cell with multiple 
%                           names {'rdk','gabor'}). 
%                           Output will be sorted by order of input name, 
%                           or ascending if input is empty.                   
% 'taskclassnames'        : Provide all task class names (use []), or 
%                           translate task class number(s) to name(s). 
%                           Task class number should be an integral 
%                           number between 1 and 10, can be a single number
%                           or vector. Output will be sorted by order of 
%                           input nr or ascending number if input 
%                           is empty, where 1:'fix',2:'cd',:'scc',4:'pc',
%                           5:'wm',6:'ltm',7:'img',8:'what', 9:'where',
%                           10:'how'.
% 'taskclassnumbers'      : Provide all task class numbers, or translate 
%                           task class name(s) to number(s).
%                           Task class name should be a string or
%                           character class, and should be part of the 10 
%                           task classes: 'fix','cd','scc','pc','wm', 'ltm'
%                           'img','what','where','how'.
%                           Name can be lower case or upper case. Name can 
%                           be a single name ('fix') or cell with multiple 
%                           names {'scc','wm'}). Output will be sorted by 
%                           order of input nr or names in ascending order 
%                           if input is empty.
% 'crossingnames'         : Provide all task-stimulus class crossing names,
%                           (use []) or translate crossing numbers(s) to  
%                           crossing names(s). Crossing numbers should be an  
%                           integral number between 1 and 32. Inputs can be  
%                           empty vector, a single number or a list of 
%                           numbers in a vector. Output will be sorted by
%                           order of input nr, or names in ascending order 
%                           if input is empty.
% 'crossingnumbers'       : Provide all crossings between task and stimulus   
%                           classes used in VCD (use []), or translate 
%                           crossing name(s) to crossing number(s).
%                           Crossing name should be a string or
%                           character class in the format <'task-stim'>, 
%                           and should be part of the 32 crossings:
%                           'fix-gabor','cd-gabor','pc-gabor','wm-gabor', 
%                           'img-gabor','fix-rdk','cd-rdk','pc-rdk',
%                           'wm-rdk','img-rdk','fix-dot','cd-dot','pc-dot',
%                           'wm-dot','img-dot''fix-obj','cd-obj','pc-obj',
%                           'wm-obj','img-obj','what-obj','how-obj',
%                           'fix-ns','cd-ns','pc-ns','wm-ns','img-ns',
%                           'what-ns','where-ns','how-ns'. 
%                           Note that for 'ltm' and 'scc' task classes, 
%                           we mix classic stimulus classes and use
%                           crossing names: 'ltm-all','scc-all'.
%                           Name can be lower case or upper case. Name can 
%                           be a single name ('fix-gabor') or cell with 
%                           multiple names {'scc-all','wm-obj'}). Output 
%                           will be sorted by order of input name, of
%                           ascending crossing number if input is empty.
% 'allstimulusnumbers'	  : Provide all stimulus numbers associated with 
%                           given stimulus class name or number. Stimulus
%                           class name should be a string or character
%                           class, and should be part of the 5 stimulus
%                           classes: 'gbr','rdk','dot','obj','ns'. Name can
%                           be lower case or upper case, and can be a
%                           single string ('rdk') or cell with multiple
%                           names {'rdk','gabor'}). Stimulus class number
%                           should be an integral number between 1 and 5,
%                           can be a single number or vector. Output will
%                           be integral, sorted by ascending stimulus 
%                           number, and range between 1-1550.
% 'stimulusnumberstonames' : Provide stimulus name(s) associated with 
%                           given image number(s). Stimulus number should
%                           be be integral and range between 1-1550. Input
%                           can be a single number or a vector. Output is a
%                           string or cell vector with strings in format:
%                            <stimclassname>-<stimnumber>-<stimlocation>. 
%                           Note that for stimclassname: "gabor" is 
%                           shortened to "GBR".
% 'stimulusnamestonumbers' : Provide stimulus numbers(s) associated with 
%                           stimulus name(s). Input stimulus name(s) should
%                           be a string in the format:
%                           <stimclassname>-<stimnumber>-<stimlocation>,
%                           and can also be cell array of strings. Note
%                           that for stimclassname: "gabor" is shortened to
%                           "GBR". Output is a vector with image numbers. If
%                           the stimulus class name and/or stimulus
%                           location in the stimulus name does not
%                           correspond to the actual stimulus class and/or
%                           location, the function will throw a warning and
%                           return NaN.
% 'stimtostimclassname'   : Provide stimulus class name(s) for given 
%                           stimulus number(s). Stimulus number should
%                           be be integral and range between 1-1550. Input
%                           can be a single number or a vector. Output is a
%                           string or cell vector with stimulus class
%                           names: 'gabor','rdk','dot','obj','ns'.
% 'stimtostimclassnumber' : Provide stimulus class number(s) for given 
%                           stimulus number(s). Stimulus number should be
%                           be integral and range between 1-1550. Input can
%                           be a single number or a vector. Output is a
%                           singleton or vector with stimulus class numbers
%                           ranging from 1-5 corresponding to 1:'gabor',
%                           2:'rdk',3:'dot',4:'obj',5:'ns'.
% 'stimtotaskclassname'   : Provide task class name(s) for given 
%                           stimulus number(s). Stimulus number should
%                           be be integral and range between 1-1550. Input
%                           can be a single number or a vector. Output is a
%                           string or cell vector with task class names: 
%                           'fix','cd','scc','pc','wm', 'ltm','img','what',
%                           'where','how'.
% 'stimtotaskclassnumber' : Provide task class number(s) for given 
%                           stimulus number(s). Stimulus number should be
%                           be integral and range between 1-1550. Input can
%                           be a single number or a vector. Output is a
%                           singleton or vector with task class numbers.
%                           ranging from 1-10 corresponding to 1:'fix',
%                           2:'cd',3:'scc',4:'pc',5:'wm',6:'ltm',7:'img',
%                           8:'what', 9:'where',10:'how'.
% 'allcore'               : Provide all core stimulus number(s) for one (or
%                           more) stimulus class(es). Input names should be
%                           lower case and any of the stimulus class names:
%                           'gabor','rdk','dot','obj','ns'. Output is a
%                           vector of integers, sorted in ascending order.
% 'specialcore'           : Provide subset of core stimulus number(s) used 
%                           in LTM-[stimclass] and IMG-[stimclass] crossings.
%                           Input names should be lower case and any of the
%                           stimulus class names:'gabor','rdk','dot','obj',
%                           'ns'. Output is a vector of integers, sorted 
%                           in ascending order.
% 'allwmtestimages'       : Provide all (non-core) stimulus number(s) of
%                           test images used in WM-[stimclass] crossings.
%                           Input names should be lower case and any of the
%                           stimulus class names:'gabor','rdk','dot','obj',
%                           'ns'. Output is a vector of integers, sorted 
%                           in ascending order.
% 'allimgtestimages'      : Provide all (non-core) stimulus number(s) of 
%                           test images used in LTM-[stimclass] crossings.
%                           Input names should be lower case and any of the
%                           stimulus class names:'gabor','rdk','dot','obj',
%                           'ns'. Output is a vector of integers, sorted 
%                           in ascending order.
% 'allltmtestimages'      : Provide all (non-core) stimulus number(s) of 
%                           test images used in IMG-[stimclass] crossings. 
%                           Input names should be lower case and any of the
%                           stimulus class names:'gabor','rdk','dot','obj',
%                           'ns'. Output is a vector of integers, sorted 
%                           in ascending order.
% 'fullinfo'              : Provide all info (data dump) for a given 
%                           stimulus number. Input should be an integral
%                           number between 1-1550.
%
%
% Examples:
% vcd('stimulusclassnames',[3 2 2 2 1])  
% vcd('taskclassnames',[1 5 3 2 2 2 1])               
% vcd('crossingnames',[2 3 4 32])   
% vcd('crossingnumbers',{'pc-gabor' 'scc-all'})
% vcd('stimulusnumberstonames',[12 35 3 2 2 2 1])        
% vcd('taskclassnumbers','FIX')
% vcd('stimulusclassnumber',{'gabor','obj'})
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
%
% Written by Eline Kupers @ UMN 2025/04


%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%

% Get experimental and stimulus parameters.
exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);

p0 = inputParser;
% General
p0.addParameter('verbose'               , true, @islogical);  
p0.addParameter('displayname'           , '7TAS_BOLDSCREEN32', @(x) any(strcmp(x,{'7TAS_BOLDSCREEN32', 'KKOFFICE_AOCQ3277', 'PPROOM_EIZOFLEXSCAN', 'EKHOME_ASUSVE247'})));                   

% Broad stimulus class & task class names and numbers
p0.addParameter('stimulusclassnames'    , [], @(x) isempty(x) | all(isnumeric(x) & (x>=1 & x <=5)));                            
p0.addParameter('stimulusclassnumbers'  , [], @(x) isempty(x) | all(ismember(lower(x), exp.stimclassnames))); 
p0.addParameter('taskclassnames'        , [], @(x) isempty(x) | all(isnumeric(x) & (x>=1 & x <=10)));           
p0.addParameter('taskclassnumbers'      , [], @(x) isempty(x) | all(ismember(lower(x), exp.taskclassnames))); 
p0.addParameter('crossingnames'         , [], @(x) isempty(x) | all((isnumeric(x) & (x>=1 & x <=32))));           
p0.addParameter('crossingnumbers'       , [], @(x) isempty(x) | all(ismember(lower(x), exp.crossingnames)));  

% Stimulus specific
p0.addParameter('allstimulusnumbers'    , [], @(x) isempty(x) | all(ismember(lower(x), exp.stimclassnames))); 
p0.addParameter('stimulusnumberstonames', [], @(x) all(isnumeric(x) & (stim.all_im_nrs(1) >=1 & x <= stim.all_im_nrs(end))));   
p0.addParameter('stimulusnamestonumbers', [], @(x) all(iscell(x)) | ischar(x));
p0.addParameter('allcore'               , [], @(x) isempty(x) | all(ismember(lower(x), exp.stimclassnames)));                 
p0.addParameter('specialcore'           , [], @(x) isempty(x) | all(ismember(lower(x), exp.stimclassnames)));                 
p0.addParameter('stimtostimclassname'   , [], @(x) all(isnumeric(x) & (stim.all_im_nrs(1) >=1 & x <= stim.all_im_nrs(end))));   
p0.addParameter('stimtostimclassnumber' , [], @(x) all(isnumeric(x) & (stim.all_im_nrs(1) >=1 & x <= stim.all_im_nrs(end))));   
p0.addParameter('stimtotaskclassname'   , [], @(x) all(isnumeric(x) & (stim.all_im_nrs(1) >=1 & x <= stim.all_im_nrs(end))));   
p0.addParameter('stimtotaskclassnumber' , [], @(x) all(isnumeric(x) & (stim.all_im_nrs(1) >=1 & x <= stim.all_im_nrs(end))));   
p0.addParameter('fullinfo'              , [], @(x) all(isnumeric(x) & (stim.all_im_nrs(1) >=1 & x <= stim.all_im_nrs(end))));  

% Test image specific (only for WM, LTM, IMG task crossings)
p0.addParameter('allwmteststimulusnumbers' , [], @(x) isempty(x) | all(ismember(lower(x), exp.stimclassnames)));   
p0.addParameter('allltmteststimulusnumbers', [], @(x) isempty(x) | all(ismember(lower(x), exp.stimclassnames)));  
p0.addParameter('allimgteststimulusnumbers', [], @(x) isempty(x) | all(ismember(lower(x), exp.stimclassnames)));  

% Parse inputs
p0.parse(varargin{:});
verbose     = p0.Results.verbose;
displayname = p0.Results.displayname;

% Check which info was requested
requested_info = {};
for ivar = 1:length(p0.Parameters)
    if ~ismember(p0.Parameters(ivar),p0.UsingDefaults)
        requested_info{length(requested_info)+1} = p0.Parameters{ivar}; 
        requested_info{length(requested_info)+1} = p0.Results.(p0.Parameters{ivar});
    end
end
 
% Check if info table is cached, if not, we load it
global vcd_info;

if isempty(vcd_info) || ~exist('vcd_info','var')
    if ~exist('vcd_rootPath','file')
        error('[%s]: Please navigate to vcd-stim folder',mfilename)
    end
    
    d = dir(fullfile(vcd_rootPath,'workspaces','info',sprintf('trials*%s*.mat',displayname)));
    
    if isempty(d)
        error('[%s]: Can''t find vcd_info! Looking for a file called: "vcd_rootPath/workspaces/info/trials*%s.mat"',mfilename,displayname)
    else
        if verbose
        fprintf('[%s]: Found %d vcd info file(s)\n',mfilename,length(d));
            if length(d) > 1
                warning('[%s]: Multiple files with the same name exist! Will pick the most recent one', mfilename);
            end
            fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
        end
        vcd_info = load(fullfile(d(end).folder,d(end).name),'condition_master');
    end
end

% If we haven't done so already, create a vector of stimulus class numbers or names for each unique image nr
if ~exist('all_core_im_stimclassnumbers','var') || isempty(all_core_im_stimclassnumbers)
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
    
    all_ltm_lure_im_stimclassnumbers   = 5*ones(1, numel(stim.ns.unique_im_nrs_ltm_lures));
    
    all_core_im_stimclassnames         = exp.stimclassnames(all_core_im_stimclassnumbers);
    all_wm_test_im_stimclassnames      = exp.stimclassnames(all_wm_test_im_stimclassnumbers);
    all_img_test_im_stimclassnames     = exp.stimclassnames(all_img_test_im_stimclassnumbers);
    all_ltm_lure_im_stimclassnames     = exp.stimclassnames(all_ltm_lure_im_stimclassnumbers);
    all_im_stimclassnames              = cat(2,all_core_im_stimclassnames,all_wm_test_im_stimclassnames,all_img_test_im_stimclassnames,all_ltm_lure_im_stimclassnames);
    
    assert(isequal(length(stim.all_im_nrs),length(all_im_stimclassnames)))
end

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
                idx    = sort(info_val,'ascend');
                out = exp.stimclassnames(idx);
            elseif isnumeric(info_val)
                idx    = sort(info_val,'ascend');
                out = exp.stimclassnames(idx);
            end
            
        case 'stimulusclassnumbers'
                        
            if isempty(info_val)                                % provide all stimulus class numbers: 1-5
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
            
            if isempty(info_val)
                out = exp.taskclassnames;                    % provide all task class names: {'fix','cd','scc','pc','wm','ltm','img','what','where','how'};
            elseif iscell(info_val) && all(cellfun(@isnumeric, info_val))
                info_val = cell2mat(info_val);
                idx    = info_val;
                out = exp.taskclassnames(idx);    
            elseif isnumeric(info_val)                          % provide specific task class name(s)
                idx    = info_val;
                out = exp.taskclassnames(idx);
            end
            
        case 'taskclassnumbers'
            
            if isempty(info_val)                               % provide all task class numbers: 1-10
                out = 1:length(exp.taskclassnames);         
                
            elseif ischar(info_val)
                info_val = lower(info_val); % use lower case characters
                
                if strcmp(info_val,'all')                      % provide all task class numbers: 1-10
                    out = 1:length(exp.taskclassnames);      
                elseif any(strcmp(info_val,exp.taskclassnames)) % provide specific task class number(s)
                    [~,out] = ismember(info_val,exp.taskclassnames);
                end  
            elseif iscell(info_val)
                if strcmp(info_val{:},'all')                    % provide all task class numbers: 1-10
                    out = 1:length(exp.taskclassnames);
                else                                            % provide specific task class number(s)
                    [~,out] = ismember(lower(info_val(:)),exp.taskclassnames);
                end
            end
            
        case 'crossingnames'
            
            if isempty(info_val)                % provide all task-stimulus class crossing names: {'fix-gbr',...,'how-ns'};
                out = exp.crossingnames; 
            elseif iscell(info_val) && all(cellfun(@isnumeric, info_val)) 
                info_val = cell2mat(info_val);
                idx    = info_val;
                out = exp.crossingnames(idx);    
            elseif isnumeric(info_val)           % provide specific task-stimulus class crossing name(s)
                idx    = info_val;
                out = exp.crossingnames(idx);
            end
            if size(out,1) > size(out,2)
                out = out';
            end
            
        case 'crossingnumbers'
            
            if isempty(info_val) % provide all task-stimulus class crossing numbers (1-32), if user doesn't specify
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
            if isempty(info_val) % provide all image nrs 1-1550, if user doesn't specifiy
                out = stim.all_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
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
                out = sort(out,'ascend');
            end
            
        case 'stimulusnumberstonames'
            % provide stimulus name(s) associated with given stimulus number(s).
            
            % use cell vector if we deal with multiple stimulus numbers
            if length(info_val) > 1
                out = cell(1,length(info_val));
            else
                out = NaN(1,length(info_val));
            end
            
            if isnumeric(info_val) % if user provides singleton or vector
                stim_nr = info_val;

            elseif iscell(info_val)  % if user provides cell
                stim_nr = cell2mat(info_val);
            end
            
            for nn = 1:length(stim_nr)
                
                % if we deal with a core stimulus
                if ismember(stim_nr(nn), stim.all_core_im_nrs)  
                    
                    [~,idx] = ismember(stim_nr(nn),vcd_info.condition_master.stim_nr_left);
                    if idx == 0
                        [~,idx] = ismember(stim_nr(nn), vcd_info.condition_master.stim_nr_right);
                        im_loc = 2; % right stim
                        stimnumber = vcd_info.condition_master.stim_nr_right(idx);
                    else
                        im_loc = 1; % left stim
                        stimnumber = vcd_info.condition_master.stim_nr_left(idx);
                    end

                elseif ismember(stim_nr(nn), stim.all_test_im_nrs)  
                    % find unique stimulus number for test image
                    [~,idx] = ismember(stim_nr(nn),vcd_info.condition_master.stim2(:,1));
                    if idx == 0
                        [~,idx] = ismember(stim_nr(nn), vcd_info.condition_master.stim2(:,2));
                        im_loc = 2; % right stim
                        stimnumber = vcd_info.condition_master.stim_nr_right(idx);
                    else
                        im_loc = 1; % left stim
                        stimnumber = vcd_info.condition_master.stim_nr_left(idx);
                    end
                end
                stimclassname = upper(vcd_info.condition_master.stim_class_name{idx,im_loc});
                if strcmp(stimclassname,'GABOR')
                    stimclassname = regexprep('GABOR','[AO]','');
                end
                stimnumber = num2str(stimnumber);
                stimlocation  = choose((im_loc==1),'L','R');
                if ischar(stimnumber)
                    stimnumber = sprintf('%03d',str2double(stimnumber));
                elseif isnumeric(stimnumber)
                    stimnumber = sprintf('%03d',stimnumber);
                end
                if iscell(out)
                    out{nn} = sprintf('%s-%s-%s',stimclassname,stimnumber,stimlocation);
                else
                    out(nn) = sprintf('%s-%s-%s',stimclassname,stimnumber,stimlocation);
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
                input_stim_class_name = input_name0{1};
                input_stim_nr         = str2double(input_name0{2});
                input_stim_pos        = input_name0{3};
                
                % check if stim nr matches stim class and stim location
                true_stim_class_name = all_im_stimclassnames(input_stim_nr==stim.all_im_nrs);
                
                if strcmp(true_stim_class_name,'gabor')
                    true_stim_class_name = 'GBR';
                end
                
                if ~strcmp(input_stim_class_name,true_stim_class_name)
                    warning('[%s]: Stimulus class name "%s" is not associated with this stimulus number %03d!',mfilename,input_stim_class_name,input_stim_nr)
                    out(nn) = NaN;
                end
               
                % Check position of stimulus
                [~,idx_l] = ismember(input_stim_nr, vcd_info.condition_master.stim_nr_left);
                [~,idx_r] = ismember(input_stim_nr, vcd_info.condition_master.stim_nr_right);
                if (idx_l == 0) && (idx_r > 0) % if right stimulus nr matches
                    true_stim_pos = 'R'; % right stim
                elseif (idx_l > 0) && (idx_r == 0) % if left stimulus nr matches
                    true_stim_pos = 'L'; % left stim
                elseif (idx_l == 0) && (idx_r == 0) % if none matches
                    true_stim_pos = NaN;
                end

                if ~strcmp(input_stim_pos,true_stim_pos) || isnan(true_stim_pos)
                    warning('[%s]: Stimulus position "%s" is not associated with stimulus number %03d!',mfilename,input_stim_pos,input_stim_nr)
                    out(nn) = NaN;
                end
                
                if out(nn) == 0
                    out(nn) = input_stim_nr;
                end
            end   
            
        case 'allcore' 
           % provide all core stimulus number(s) for one (or more) stimulus class(es)
           if isempty(info_val) % provide all coreimage nrs (1-110), if user doesn't specifiy stimulus class
                out = stim.all_core_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
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
                out = sort(out,'ascend');
           end 
            
        case 'specialcore'
            % provide all special core stimulus number(s) for one (or more) stimulus class(es)
            if isempty(info_val) % provide all special 47 core image nrs, if user doesn't specifiy stimulus class
                out = stim.all_specialcore_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
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
                out = sort(out,'ascend');
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
           if isempty(info_val) % provide all wm test image nrs (1-110), if user doesn't specifiy stimulus class
                out = stim.all_wm_test_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
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
                out = sort(out,'ascend');
           end 
            
            
        case 'allltmteststimulusnumbers'      
            % provide all (non-core) stimulus number(s) of test images used in LTM-[stimclass] crossings 
            if isempty(info_val) % provide all wm test image nrs (1-110), if user doesn't specifiy stimulus class
                out = stim.all_ltm_lure_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
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
                out = sort(out,'ascend');
            end
            
        case 'allimgteststimulusnumbers'      
            % provide all (non-core) stimulus number(s) of test images used in IMG-[stimclass] crossings 
            if isempty(info_val) % provide all wm test image nrs (1-110), if user doesn't specifiy stimulus class
                out = stim.all_img_test_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
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
                out = sort(out,'ascend');
            end

        case 'fullinfo'              
            % provide all info (data dump) for given stimulus number (1-1550)
            if isnumeric(info_val) % if user provides singleton or vector
                stim_nr = info_val;
            elseif iscell(info_val)  % if user provides cell, turn into vector
                stim_nr = info_val{:};
            end
            
            out = cell(1,length(stim_nr));
            for nn = 1:length(stim_nr)
                
                % if we deal with a core stimulus
                if ismember(stim_nr(nn), stim.all_core_im_nrs)  
                    
                    [~,idx] = ismember(stim_nr(nn),vcd_info.condition_master.stim_nr_left);
                    if idx == 0
                        [~,idx] = ismember(stim_nr(nn), vcd_info.condition_master.stim_nr_right);
                        im_loc = 2; % right stim
                    else
                        im_loc = 1; % left stim
                    end

                elseif ismember(stim_nr(nn), stim.all_test_im_nrs)  
                    % find unique stim number for test image
                    [~,idx] = ismember(stim_nr(nn),vcd_info.condition_master.stim2(:,1));
                    if idx == 0
                        [~,idx] = ismember(stim_nr(nn), vcd_info.condition_master.stim2(:,2));
                        im_loc = 2; % right stim
                    else
                        im_loc = 1; % left stim
                    end
                end
                    
                % Get stim info (reference first)
                stim_info = table2struct(vcd_info.condition_master(idx,:));
                field_sz  = struct2array(structfun(@numel, stim_info, 'UniformOutput', false));
                fnames = fieldnames(stim_info);
                for ff = find(field_sz==2)
                    if iscell(stim_info.(fnames{ff})(im_loc))
                        stim_info.(fnames{ff}) = stim_info.(fnames{ff}){im_loc};
                    else
                        stim_info.(fnames{ff}) = stim_info.(fnames{ff})(im_loc);
                    end
                end

                % add stim loc           
                if stim_info.stim_class == 5
                    stim_info.stimloc = 3;
                    stim_info.stimloc_name = 'center';
                else
                    stim_info.stimloc = im_loc;
                    stim_info.stimloc_name = choose(im_loc==1,'left','right');
                end
                    
                % remove irrelevant stimulus location , fields with nans, and trial info
                info_to_delete = {'is_cued','repeat_nr','task_class_name','task_class',...
                    'unique_trial_nr','stim_class_unique_block_nr','block_local_trial_nr'};
                if im_loc==1
                    info_to_delete = cat(2,info_to_delete,'stim_nr_right');
                    if ismember(stim_nr(nn), stim.all_test_im_nrs)
                        stim_info.stim1 = stim_info.stim_nr_left;
                    end
                elseif im_loc==2
                    info_to_delete = cat(2,info_to_delete,'stim_nr_left');
                    if ismember(stim_nr(nn), stim.all_test_im_nrs)
                        stim_info.stim1 = stim_info.stim_nr_right;
                    end
                end
                    
                for jj = 1:length(info_to_delete)
                    stim_info = rmfield(stim_info,(info_to_delete{jj}));
                end
                
                % remove fields with NaNs
                fnames = fieldnames(stim_info); % get updated fieldnames
                field_nan  = cellfun(@(x) all(x==1), struct2cell(structfun(@isnan, stim_info, 'UniformOutput', false)));
                field_nan  = find(field_nan)';
                for ff = field_nan
                    stim_info = rmfield(stim_info,(fnames{ff}));
                end
                
                out{nn} = stim_info;
            end
    end
    
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


