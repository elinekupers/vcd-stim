% vcd_startup.m
% 
% Script to add paths and toolbox dependencies prior to running the
% experiment. Requires user to tell what environment we are in: 7TAS (1),
% psychophysics room (2) or other (3). if other, user needs to be change 
% the path in Matlab to their local
% vcd-stim code folder or
exp_env = []; % can be 'pproom', '7tas', 'other'

% Ask the user what environment we are in
if ~exist('exp_env','var') || isempty(exp_env)
  exp_env = input('Where are we running the core VCD experiment?  1: 7TAS   2: psychophysics room   3: other \n')
end

switch exp_env
    
    case 1 % 7TAS
        vcdcode_dir  = '/Users/7tasuser/Desktop/cvnlab/VCD/vcd-stim';             % where does vcd stim code live?
        knkutils_dir = '/Users/7tasuser/Desktop/cvnlab/kendrick/knkutils-master'; % where does knk utils code live?
        ptb_dir      = '/Applications/Psychtoolbox';                              % where does ptb code live?
        
    case 2 % psychophysics room
        vcdcode_dir  = '/Users/psphuser/Desktop/cvnlab/VCD/vcd-stim';
        knkutils_dir = '/Users/psphuser/Desktop/cvnlab/VCD/knkutils';
        ptb_dir      = '/Applications/Psychtoolbox';
        
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
end

% Add vcd-stim subfolders
cd(vcdcode_dir)
addpath(genpath(vcdcode_dir))
addpath(genpath(knkutils_dir))

fprintf('The following folders have been added to the path: \n')
fprintf('%s\n',vcdcode_dir)
fprintf('%s\n',knkutils_dir)

% Check PTB toolbox
if ~exist('Screen','file')
    addpath(genpath(ptb_dir))
    fprintf('%s:\n',ptb_dir)
end






