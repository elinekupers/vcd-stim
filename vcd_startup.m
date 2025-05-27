function [] = vcd_startup(exp_env)
% 
% Script to add paths and toolbox dependencies prior to running the
% experiment. Requires user to tell what environment we are in: 7TAS (1),
% psychophysics room (2) or other (3). if other, user needs to be change 
% the path in Matlab to the root of their local vcd-stim code folder

if ~exist('exp_env','var')
    exp_env = []; % can be 1:'7tas', 2:'pproom', 3:'other'
end

% Ask the user what environment we are in
if isempty(exp_env)
   exp_env = input('Where are we running the core VCD experiment?  1: 7TAS   2: psychophysics room   3: other \n')
elseif ~isempty(exp_env) 
    % Check if inputs are either 1, 2, or 3
    assert(isequal(length(exp_env),1))
    assert(ismember(exp_env,[1,2,3]))
end
    

% Define paths
switch exp_env
    
    case 1 % 7TAS
        vcdcode_dir  = '/Users/7tasuser/Desktop/cvnlab/VCD/vcd-stim';             % where does vcd stim code live?
        knkutils_dir = '/Users/7tasuser/Desktop/cvnlab/kendrick/knkutils-master'; % where does knk utils code live?
        ptb_dir      = '/Applications/Psychtoolbox';                              % where does ptb code live?
        edf2asc_func = '/Applications/Eyelink/EDF_Access_API/Example';
    case 2 % psychophysics room
        vcdcode_dir  = '/Users/psphuser/Desktop/cvnlab/VCD/vcd-stim';
        knkutils_dir = '/Users/psphuser/Desktop/cvnlab/VCD/knkutils';
        ptb_dir      = '/Applications/Psychtoolbox';
        edf2asc_func = '/Applications/Eyelink/EDF_Access_API/Example';
        
    case 3 % other
        if exist('vcd_rootPath','file')
            vcdcode_dir = vcd_rootPath;
        else
            error('Please navigate to vcd-stim folder on your local machine or define it in vcd_startup.m !\n')
        end
        % assuming knkutils and ptb code lives in the same folder tree as 
        % vcd-stim code
        knkutils_dir = fullfile(fileparts(vcd_rootPath),'knkutils');
        ptb_dir      = fullfile(fileparts(vcd_rootPath),'Psychtoolbox-3');
        edf2asc_func = []; % assume we don't do eyetracking
end

fprintf('The following folders have been added to the path: \n')

% Add knkutils toolbox
if ~exist('pton','file')
    addpath(genpath(knkutils_dir))
    fprintf('%s\n',knkutils_dir)
else
    fprintf('%s already is added to paths\n',knkutils_dir)
end

% Add vcd-stim toolbox
if ~exist('vcd_rootPath','file')
    addpath(genpath(vcdcode_dir))
    fprintf('%s\n',vcdcode_dir)
else
    fprintf('%s already is added to paths\n',vcdcode_dir)
end

% Add PTB toolbox
if ~exist('Screen','file')
    addpath(genpath(ptb_dir))
    fprintf('%s\n',ptb_dir)
else
    fprintf('%s already is added to paths\n',ptb_dir)
end

% Add location of edf2asc converter 
if ~exist('edf2asc','file') && ~isempty(edf2asc_func)
    pth = strcat([edf2asc_func ':/usr/bin/:/bin:/usr/sbin:/sbin']);
    if ~isequal(getenv('PATH'),pth)
        setenv('PATH',pth);
        fprintf('%s\n',edf2asc_func)
    else
        fprintf('%s already is added to paths\n',pth)
    end
else
    fprintf('%s already is added to paths\n',edf2asc_func)
end
    
    
% Go to root of vcd-stim folder
cd(vcdcode_dir)








