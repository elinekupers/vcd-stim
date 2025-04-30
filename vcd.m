function outputs = vcd(varargin)
% VCD function to call general information about experimental design, 
% including the number of stimulus and task class names, numbers, unique
% image number, full information for each unique stimulus.
%
% There are 5 stimulus classes: 'gbr' 'rdk' 'dot' 'obj' 'ns'

% 'stimulusclassnames'               : provide all stimulus class names, or translate stimulus class number(s) to name(s)
% 'stimulusclassnumbers'	         : provide all stimulus class numbers, or translate stimulus class name(s) to number(s)
% 'taskclassnames'                   : provide all task class names, or translate task class number(s) to name(s)
% 'taskclassnumbers'                 : provide all task class numbers, or translate stimulus class name(s) to number(s)
% 'crossingnames'                    : provide all task-stimulus class crossing names, or translate crossing numbers(s) to crossing names(s)	
% 'crossingnumber'                   : provide all task-stimulus class crossing numbers, or translate crossing name(s) to crossing number(s)	
% 'allimages'	                     : provide all stimulus numbers associated with given stimulus class name or number  
% 'stimulusname'                     : []

% 'stimtostimclassname'              : translate integers to stimulus class names
% 'stimtostimclassnumber'            : translate integers to task names
% 'stimtotaskclassname'              : translate integers to stimulus class names
% 'stimtotaskclassnumber'            : translate integers to task names
% 'allcore'                          : provide all core stimulus number(s) for one (or more) stimulus class(es)
% 'specialcore'                      : provide subset of core stimulus number(s) used in LTM-[stimclass] and IMG-[stimclass] crossings 
% 'allwmtestimages'                  : provide all (non-core) stimulus number(s) of test images used in WM-[stimclass] crossings 
% 'allimgtestimages'                 : provide all (non-core) stimulus number(s) of test images used in LTM-[stimclass] crossings 
% 'allltmtestimages'                 : provide all (non-core) stimulus number(s) of test images used in IMG-[stimclass] crossings 
% 'fullinfo'                         : provide all info (data dump) for given stimulus number (1-1550)


% Examples:
% vcd('stimulusclassname',[1156 35 3 2 2 2 1])  
% vcd('taskname',[1 5 3 2 2 2 1])               
% vcd('crossingname',[2 3 4 32])                
% vcd('stimulusname',[12 35 3 2 2 2 1])        
% vcd('tasknumber','FIX')
% vcd('stimulusclassnumber',{'GBR','OBJ'})
% vcd('crossingnumber',{'PC-GBR' 'SCC-ALL'})
% vcd('stimulusnumbertostimulusname',{'12-GBR0-CON50-'})
% vcd('stimulusnametostimulusnumber',12)
% vcd('allcore',{'GBR' 'RDK'})                  
% vcd('specialcore','RDK')
% vcd('allwmtestimages',{'RDK' 'NS'})
% vcd('allimgtestimages','RDK')
% vcd('allltmtestimages',RDK')
% vcd('allimages','RDK',)
% vcd('allimages',{'GBR' 'RDK' 'DOT' 'OBJ' 'NS'})
% vcd('fullinfo',35) 



%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%

% Get experimental and stimulus parameters.
exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);

p0 = inputParser;
% General
p0.addParameter('stimulusclassnames'    , [], @(x) isempty(x) || validateattributes(x, 'numeric', {'real','finite','>=',1,'<=',5}));   % provide all stimulus class names, or translate stimulus class number(s) to name(s)
p0.addParameter('stimulusclassnumbers'  , [], @(x) isempty(x) || validatestring(lower(x), exp.stimclassnames));                        % provide all stimulus class numbers, or translate stimulus class name(s) to number(s)
p0.addParameter('taskclassnames'        , [], @(x) isempty(x) || validateattributes(x, 'numeric',{'real','finite','>=', 1, '<=', 10}));% provide all task class names, or translate task class number(s) to name(s)
p0.addParameter('taskclassnumbers'      , [], @(x) isempty(x) || validatestring(lower(x), exp.taskclassnames));                        % provide all task class numbers, or translate stimulus class name(s) to number(s)
p0.addParameter('crossingnames'         , [], @(x) isempty(x) || validateattributes(x, 'numeric',{'real','finite','>=', 1, '<=', 32}));% provide all task-stimulus class crossing names, or translate crossing number(s) to crossing name(s)
p0.addParameter('crossingnumbers'       , [], @(x) isempty(x) || validatestring(lower(x), exp.crossingnames));                         % provide all task-stimulus class crossing numbers, or translate crossing name(s) to crossing number(s)
p0.addParameter('allimages'             , [], @(x) isempty(x) || validatestring(lower(x), exp.stimclassnames));                        % provide all stimulus numbers associated with given stimulus class name or number

% Stimulus specific
p0.addParameter('allcore'               , [], @(x) isempty(x) || validatestring(lower(x), exp.stimclassnames));                                            % provide all core stimulus number(s) for one (or more) stimulus class(es)
p0.addParameter('specialcore'           , [], @(x) isempty(x) || validatestring(lower(x), exp.stimclassnames));                                            % provide subset of core stimulus number(s) used in LTM-[stimclass] and IMG-[stimclass] crossings 
p0.addParameter('stimtostimclassname'   , [], @(x) validateattributes(x, 'numeric', {'real','finite','>=', stim.all_im_nrs(1), '<=', stim.all_im_nrs(end)}));   % translate stimulus number(s) to stimulus class name(s)
p0.addParameter('stimtostimclassnumber' , [], @(x) validateattributes(x, 'numeric', {'real','finite','>=', stim.all_im_nrs(1), '<=', stim.all_im_nrs(end)}));   % translate stimulus number(s) to stimulus class number(s)
p0.addParameter('stimtotaskclassname'   , [], @(x) validateattributes(x, 'numeric', {'real','finite','>=', stim.all_im_nrs(1), '<=', stim.all_im_nrs(end)}));   % translate stimulus number(s) to task class name(s)
p0.addParameter('stimtotaskclassnumber' , [], @(x) validateattributes(x, 'numeric', {'real','finite','>=', stim.all_im_nrs(1), '<=', stim.all_im_nrs(end)}));   % translate stimulus number(s) to task class number(s)
p0.addParameter('fullinfo'              , [], @(x) validateattributes(x, 'numeric', {'real','finite','>=', stim.all_im_nrs(1), '<=', stim.all_im_nrs(end)}));   % provide all info (data dump) for given stimulus number (1-1550)

% Test image specific (only for WM, LTM, IMG task crossings)
p0.addParameter('allwmtestimages'       , [], @(x) isempty(x) || validatestring(lower(x), exp.stimclassnames));   % provide all (non-core) stimulus number(s) of test images used in WM-[stimclass] crossings 
p0.addParameter('allltmtestimages'      , [], @(x) isempty(x) || validatestring(lower(x), exp.stimclassnames));   % provide all (non-core) stimulus number(s) of test images used in LTM-[stimclass] crossings 
p0.addParameter('allimgtestimages'      , [], @(x) isempty(x) || validatestring(lower(x), exp.stimclassnames));   % provide all (non-core) stimulus number(s) of test images used in IMG-[stimclass] crossings 

% Parse inputs
p0.parse(varargin{:});

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
    if ~exist('vcd_rootPath',2)
        error('[%s]: Please navigate to vcd-stim folder',mfilename)
    end
    
    d = dir(fullfile(vcd_rootPath,'workspaces','info','trials*.mat'));
    
    if isempty(d)
        error('[%s]: Can''t find vcd_info! Looking for a file called: "vcd_rootPath/workspaces/info/trials*.mat"',mfilename)
    else
        if verbose
        fprintf('[%s]: Found %d vcd info file(s)\n',mfilename,length(d));
            if length(d) > 1
                warning('[%s]: Multiple files with the same name exist! Will pick the most recent one', mfilename);
            end
            fprintf('[%s]: Loading exp params .mat file: %s\n', mfilename, d(end).name);
        end
        tmp = load(fullfile(d(end).folder,d(end).name));
        vcd_info = tmp; clear tmp;
    end
end

% If we haven't done so already, create a vector of stimulus class numbers or names for each unique image nr
if ~isfield(vcd_info,'all_core_im_stimclassnumbers')
    vcd_info.all_core_im_stimclassnumbers       = ...
        cat(2, ones(1, numel(stim.gabor.unique_im_nrs_core)), 2*ones(1, numel(stim.rdk.unique_im_nrs_core)), ...
        3*ones(1, numel(stim.dot.unique_im_nrs_core)), 4*ones(1, numel(stim.obj.unique_im_nrs_core)), ...
        5*ones(1, numel(stim.ns.unique_im_nrs_core)));
    
    vcd_info.all_wm_test_im_stimclassnumbers    = ...
        cat(2, ones(1, numel(stim.gabor.unique_im_nrs_wm_test)), 2*ones(1, numel(stim.rdk.unique_im_nrs_wm_test)), ...
        3*ones(1, numel(stim.dot.unique_im_nrs_wm_test)), 4*ones(1, numel(stim.obj.unique_im_nrs_wm_test)), ...
        5*ones(1, numel(stim.ns.unique_im_nrs_wm_test)));
    
    vcd_info.all_img_test_im_stimclassnumbers   = ...
        cat(2, ones(1, numel(stim.gabor.unique_im_nrs_img_test)), 2*ones(1, numel(stim.rdk.unique_im_nrs_img_test)), ...
        3*ones(1, numel(stim.dot.unique_im_nrs_img_test)), 4*ones(1, numel(stim.obj.unique_im_nrs_img_test)), ...
        5*ones(1, numel(stim.ns.unique_im_nrs_img_test)));
    
    vcd_info.all_ltm_lure_im_stimclassnumbers   = 5*ones(1, numel(stim.ns.unique_im_nrs_ltm_lures));
    
    vcd_info.all_core_im_stimclassnames         = exp.stimclassnames(all_core_im_stimclassnumbers);
    vcd_info.all_wm_test_im_stimclassnames      = exp.stimclassnames(all_wm_test_im_stimclassnumbers);
    vcd_info.all_img_test_im_stimclassnames     = exp.stimclassnames(all_img_test_im_stimclassnumbers);
    vcd_info.all_ltm_lure_im_stimclassnames     = exp.stimclassnames(all_ltm_lure_im_stimclassnumbers);
    vcd_info.all_im_stimclassnames              = cat(2,all_core_im_stimclassnames,all_wm_test_im_stimclassnames,all_img_test_im_stimclassnames,all_ltm_lure_im_stimclassnames);
    
    assert(isequal(length(stim.all_im_nrs),length(all_im_stimclassnames)))
end

% Now check requested_info by user
for ii = 1:2:length(requested_info)
    
    info_name = requested_info{ii};
    info_val  = requested_info{ii+1};
    
    switch info_name
        
        case 'stimulusclassnames'
            if isempty(info_val) || strcmp(info_val{:},'all')
                output = exp.stimclassnames; % {'gabor','rdk','dot','obj','ns'};
                
            elseif isnumeric(info_val)
                idx    = sort(info_val,'ascend');
                output = exp.stimclassnames(idx);
            end
            
        case 'stimulusclassnumbers'
                        
            if isempty(info_val)                                % provide all stimulus class numbers: 1-5
                output = 1:length(exp.stimclassnames); 
                
            elseif ischar(info_val)
                info_val = lower(info_val);                     % use lower case characters
                
                if strcmp(info_val,'all')                       % provide all stimulus class numbers: 1-5
                    output = 1:length(exp.stimclassnames); 
                elseif any(strcmp(info_val,exp.stimclassnames)) % provide specific stimulus class number(s)
                    [~,output] = ismember(info_val,exp.stimclassnames);
                end
 
            elseif iscell(info_val)
                if strcmp(info_val{:},'all')                    % provide all stimulus class numbers: 1-5
                    output = 1:length(exp.stimclassnames); 
                else                                            % provide specific stimulus class number(s)
                    [~,output] = ismember(lower(info_val(:)),exp.stimclassnames);
                end
            end
            
            output = sort(output,'ascend');
            
        case 'taskclassnames'
            
            if isempty(info_val)
                output = exp.taskclassnames;                    % provide all task class names: {'fix','cd','scc','pc','wm','ltm','img','what','where','how'};
                
            elseif isnumeric(info_val)                              % provide specific task class name(s)
                idx    = sort(info_val,'ascend');
                output = exp.taskclassnames(idx);
            end
            
        case 'taskclassnumbers'
            
            if isempty(info_val)                               % provide all task class numbers: 1-10
                output = 1:length(exp.taskclassnames);         
                
            elseif ischar(info_val)
                info_val = lower(info_val); % use lower case characters
                
                if strcmp(info_val,'all')                      % provide all task class numbers: 1-10
                    output = 1:length(exp.taskclassnames);      
                elseif any(strcmp(info_val,exp.taskclassnames)) % provide specific task class number(s)
                    [~,output] = ismember(info_val,exp.taskclassnames);
                end  
            elseif iscell(info_val)
                if strcmp(info_val{:},'all')                    % provide all task class numbers: 1-10
                    output = 1:length(exp.taskclassnames);
                else                                            % provide specific task class number(s)
                    [~,output] = ismember(lower(info_val(:)),exp.taskclassnames);
                end
            end
            
            output = sort(output,'ascend');

            
        case 'crossingnames'                                    % provide all task-stimulus class crossing names: {'fix-gbr',...,'how-ns'};
            if isempty(info_val)
                output = exp.crossingnames; 
                
            elseif isnumeric(info_val)                              % provide specific task-stimulus class crossing name(s)
                idx    = sort(info_val,'ascend');
                output = exp.crossingnames(idx);
            end

        case 'crossingnumbers'
            
            if isempty(info_val) % provide all task-stimulus class crossing numbers (1-32), if user doesn't specify
                output = 1:length(exp.crossingnames);           
                
            elseif ischar(info_val)
                info_val = lower(info_val); % use lower case characters
                
                if strcmp(info_val,'all')
                    output = 1:length(exp.crossingnames);       % provide all task-stimulus class crossing numbers: 1-32
                    
                elseif any(strcmp(info_val,exp.crossingnames))  % provide specific task-stimulus class crossing number(s)
                    [~,output] = ismember(info_val,exp.crossingnames);
                end
                
            elseif iscell(info_val)
                if strcmp(info_val{:},'all')                                % provide all task-stimulus class crossing numbers: 1-32
                    output = 1:length(exp.crossingnames);
                else                                                        % provide specific task-stimulus class crossing number(s)
                    [~,output] = ismember(lower(info_val(:)), exp.crossingnames);
                end
            end
            
            output = sort(output,'ascend');

        case 'allimages'
            % provide all stimulus numbers associated with given stimulus class name or number
            if isempty(info_val) % provide all image nrs 1-1550, if user doesn't specifiy
                output = stim.all_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
                elseif ischar(info_val) % if user provides stim class name 
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names 
                    stimclassname = info_val(:);
                end
                
                output = [];
                for sc = 1:length(stimclassname)
                    if strcmp(stimclassname{sc},'ns')
                        output = cat(2, output, stim.ns.unique_im_nrs_core, stim.ns.unique_im_nrs_WM_test, stim.ns.unique_im_nrs_IMG_test,stim.ns.unique_im_nrs_LTM_lures);
                    else
                        output = cat(2, output, stim.(stimclassname{sc}).unique_im_nrs_core, stim.(stimclassname{sc}).unique_im_nrs_WM_test, stim.(stimclassname{sc}).unique_im_nrs_IMG_test);
                    end
                end
                output = sort(output,'ascend');
            end
            
        case 'allcore' 
           % provide all core stimulus number(s) for one (or more) stimulus class(es)
           if isempty(info_val) % provide all coreimage nrs (1-110), if user doesn't specifiy stimulus class
                output = stim.all_core_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
                elseif ischar(info_val) % if user provides stim class name 
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names 
                    stimclassname = info_val(:);
                end
                
                output = [];
                for sc = 1:length(stimulusclassname)
                    output = cat(2, output, stim.(stimclassname{sc}).unique_im_nrs_core);
                end
                output = sort(output,'ascend');
           end 
            
        case 'specialcore'
            % provide all special core stimulus number(s) for one (or more) stimulus class(es)
            if isempty(info_val) % provide all special 47 core image nrs, if user doesn't specifiy stimulus class
                output = stim.all_specialcore_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
                elseif ischar(info_val) % if user provides stim class name
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names
                    stimclassname = info_val(:);
                end
                
                output = [];
                for sc = 1:length(stimulusclassname)
                    output = cat(2, output, stim.(stimclassname{sc}).unique_im_nrs_specialcore);
                end
                output = sort(output,'ascend');
            end
            
        case 'stimtostimclassname'
            % translate stimulus number(s) to stimulus class name(s)
            
            if isnumeric(info_val) % if user provides singleton or vector
                
                stim_nr = info_val;

            elseif iscell(info_val)  % if user provides cell
                stim_nr = info_val{:};
            end
            
            output = all_im_stimclassnames(stim_nr);
            
        case 'stimtostimclassnumber'
            % translate stimulus number(s) to stimulus class number(s)
            if isnumeric(info_val) % if user provides singleton or vector
                
                stim_nr = info_val;
                
            elseif iscell(info_val)  % if user provides cell
                stim_nr = info_val{:};
            end
            
            output = all_im_stimclassnumbers(stim_nr);
            
        case 'stimtotaskclassname'   
            % translate stimulus number(s) to task class name(s)
            if isnumeric(info_val) % if user provides singleton or vector
                
                stim_nr = info_val;
                
            elseif iscell(info_val)  % if user provides cell
                stim_nr = info_val{:};
            end
            
            % get stimulus class nr
            stim_class_nr = all_im_stimclassnumbers(stim_nr);
            
            output = cell(1,length(stim_class_nr));
            for nn = 1:length(stim_class_nr)
                output{nn} = exp.taskclassnames(exp.crossings(stim_class_nr(nn),:));
            end
            
        case 'stimtotaskclassnumber' 
            % translate stimulus number(s) to task class number(s)
            if isnumeric(info_val) % if user provides singleton or vector
                
                stim_nr = info_val;
                
            elseif iscell(info_val)  % if user provides cell
                stim_nr = info_val{:};
            end
            
            % get stimulus class nr
            output = all_im_stimclassnumbers(stim_nr);
            
        case 'allwmtestimages'       
            % provide all (non-core) stimulus number(s) of test images used in WM-[stimclass] crossings 
           if isempty(info_val) % provide all wm test image nrs (1-110), if user doesn't specifiy stimulus class
                output = stim.all_wm_test_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
                elseif ischar(info_val) % if user provides stim class name 
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names 
                    stimclassname = info_val(:);
                end
                
                output = [];
                for sc = 1:length(stimulusclassname)
                    if isfield(stim.(stimclassname{sc}),'unique_im_nrs_wm_test')
                        output = cat(2, output, stim.(stimclassname{sc}).unique_im_nrs_wm_test);
                    end
                end
                output = sort(output,'ascend');
           end 
            
            
        case 'allltmtestimages'      
            % provide all (non-core) stimulus number(s) of test images used in LTM-[stimclass] crossings 
            if isempty(info_val) % provide all wm test image nrs (1-110), if user doesn't specifiy stimulus class
                output = stim.all_ltm_lure_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
                elseif ischar(info_val) % if user provides stim class name
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names
                    stimclassname = info_val(:);
                end
                
                output = [];
                for sc = 1:length(stimulusclassname)
                    if isfield(stim.(stimclassname{sc}), 'unique_im_nrs_ltm_lures')
                        output = cat(2, output, stim.(stimclassname{sc}).unique_im_nrs_ltm_lures);
                    end
                end
                output = sort(output,'ascend');
            end
            
        case 'allimgtestimages'      
            % provide all (non-core) stimulus number(s) of test images used in IMG-[stimclass] crossings 
            if isempty(info_val) % provide all wm test image nrs (1-110), if user doesn't specifiy stimulus class
                output = stim.all_img_test_im_nrs;
            else
                if isnumeric(info_val) % if user provides stim class number
                    stimclassname = exp.stimclassnames(ismember(info_val,1:length(exp.stimclassnames)));
                elseif ischar(info_val) % if user provides stim class name
                    stimclassname = {info_val};
                elseif iscell(info_val)  % if user provides multiple stim class names
                    stimclassname = info_val(:);
                end
                
                output = [];
                for sc = 1:length(stimulusclassname)
                    if isfield(stim.(stimclassname{sc}), 'unique_im_nrs_img_test')
                        output = cat(2, output, stim.(stimclassname{sc}).unique_im_nrs_img_test);
                    end
                end
                output = sort(output,'ascend');
            end

        case 'fullinfo'              
            % provide all info (data dump) for given stimulus number (1-1550)
            if isnumeric(info_val) % if user provides singleton or vector
                stim_nr = info_val;
            elseif iscell(info_val)  % if user provides cell, turn into vector
                stim_nr = info_val{:};
            end
            
            output = struct();
            for nn = 1:length(stim_nr)
                
                % if we deal with a core stimulus
                if ismember(stim_nr(nn), stim.all_core_im_nrs)  
                    
                    [~,idx] = ismember(stim_nr(nn),condition_master.stim_nr_left);
                    if idx == 0
                        [~,idx] = ismember(stim_nr(nn), condition_master.stim_nr_right);
                        im_loc = 2; % right stim
                    else
                        im_loc = 1; % left stim
                    end

                elseif ismember(stim_nr(nn), stim.all_test_im_nrs)  
                    % find unique stim number for test image
                    [~,idx] = ismember(stim_nr(nn),condition_master.stim2(:,1));
                    if idx == 0
                        [~,idx] = ismember(stim_nr(nn), condition_master.stim2(:,2));
                        im_loc = 2; % right stim
                    else
                        im_loc = 1; % left stim
                    end
                end
                    
                % Get stim info (reference first)
                stim_info = table2struct(condition_master(idx,:));
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
                info_to_delete = {'is_cued','repeat_nr','stim_class_unique_block_nr','block_local_trial_nr'};
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
                field_nan  = cellfun(@(x) all(x==1), struct2cell(structfun(@isnan, stim_info, 'UniformOutput', false)));
                
                for ff = find(field_nan)'
                    stim_info = rmfield(stim_info,(fnames{ff}));
                end
                
                    
               
                    
                
                    

                output(nn).stim_info;
            end
    end
    
    % Check outputs
    if length(requested_info) > 2
        if ii == 1
            % create output struct if we have more than one info request
            outputs = struct();
        end
        % Add output to outputs struct
        outputs.(info_name) = output;
    else
        outputs = output;
    end
end

 
return


