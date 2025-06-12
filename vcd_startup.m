function [] = vcd_startup(exp_env)
% 
% Script to add paths and toolbox dependencies prior to running the
% experiment. Requires user to tell what environment we are in: 7TAS (1),
% CMRR psychophysics room (2),  NYU psychophysics room (3), or other environments like KK office (4). if other, user needs to be change 
% the path in Matlab to the root of their local vcd-stim code folder

if ~exist('exp_env','var')
    exp_env = []; % can be 1:'7tas', 2:'cmrr pproom', 2:'nyu pproom', 4:'other'
end

% Ask the user what environment we are in
if isempty(exp_env)
   exp_env = input('Where are we running the core VCD experiment?  1: 7TAS   2: CMRR psychophysics room   3: NYU psychophysics room   4: other \n');
elseif ~isempty(exp_env) 
    % Check if inputs are either 1, 2, or 3
    assert(isequal(length(exp_env),1))
    assert(ismember(exp_env,[1,2,3,4]))
end
    

% Define paths
switch exp_env
    
    case 1 % 7TAS
        vcdcode_dir  = '/Users/7tasuser/Desktop/cvnlab/VCD/vcd-stim';             % where does vcd stim code live?
        knkutils_dir = '/Users/7tasuser/Desktop/cvnlab/kendrick/knkutils-master'; % where does knk utils code live?
        ptb_dir      = '/Applications/Psychtoolbox';                              % where does ptb code live?
        edf2asc_func = '/Applications/Eyelink/EDF_Access_API/Example';            % where does edf2asc converter function live?
    case 2 % CMRR psychophysics room
        vcdcode_dir  = '/Users/psphuser/Desktop/cvnlab/VCD/vcd-stim';
        knkutils_dir = '/Users/psphuser/Desktop/cvnlab/VCD/knkutils';
        ptb_dir      = '/Applications/Psychtoolbox';
        edf2asc_func = '/Applications/Eyelink/EDF_Access_API/Example';
    case 3 % NYU psychophysics room
        vcdcode_dir  = '/Users/zlu/Documents/MATLAB/vcd-stim-0.21';
        knkutils_dir = '/Users/zlu/Documents/MATLAB/knkutils-master';
        ptb_dir      = '/Users/zlu/Documents/MATLAB/Psychtoolbox';
        edf2asc_func = '/Applications/Eyelink/EDF_Access_API/Example';
    case 4 % other
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


added_paths = {};

% Add knkutils toolbox
if ~exist('pton','file')
    addpath(genpath(knkutils_dir))
    added_paths = cat(1,added_paths,{knkutils_dir});
end

% Add vcd-stim toolbox
if ~exist('vcd_singleRun','file')
    addpath(genpath(vcdcode_dir))
    added_paths = cat(1,added_paths,{vcdcode_dir});
end

% Add PTB toolbox
if ~exist('Screen','file')
    addpath(genpath(ptb_dir))
    added_paths = cat(1,added_paths,{ptb_dir});
end

% Add location of edf2asc converter 
if ~exist('edf2asc','file') && ~isempty(edf2asc_func)
    pth = strcat([edf2asc_func ':/usr/bin/:/bin:/usr/sbin:/sbin']);
    if ~isequal(getenv('PATH'),pth)
        setenv('PATH',pth);
        added_paths = cat(1,added_paths,{pth});
    end    
end

% deal with setting Bits Sharp settings
if exp_env == 2 % CMRR psychophysics room
    smatch = matchfiles('/dev/tty.usbmodem*');
    assert(length(smatch)==1);
    s1 = serial(smatch{1}); fopen(s1);
    fprintf(s1,['$BitsPlusPlus' 13]);
    fprintf(s1,['$enableGammaCorrection=[greyLums.txt]' 13]);
    % MONO PLUS MODE FOR REFERENCE
    %fprintf(s1,['$monoPlusPlus' 13]);
    %fprintf(s1,['$enableGammaCorrection=[13bitLinearLUT.txt]' 13]);
    %fprintf(s1,['$Help' 13]);
    %output=fscanf(s1)
    fclose(s1); delete(s1); clear smatch s1;
end

% Go to root of vcd-stim folder
cd(vcdcode_dir)

% Tell user what we did. If we didn't add any paths, we don't say anything..
if ~isempty(added_paths)
    fprintf('The following folders have been added to the path: \n')
    fprintf('%s\n',added_paths{:})
end




