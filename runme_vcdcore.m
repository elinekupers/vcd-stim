function runme_vcdcore(subj_nr,ses_nr,ses_type,run_nr, varargin)
% __FUNCTION__
%  runme_vcdcore(subj_nr,sesID,run_nr, varargin)
%
% __PURPOSE__
%  Display a single 5.5-min fMRI run of the core Visual Cognition Dataset
%  (VCD) experiment.
% 
% __INPUTS__
%  Input parser requires first three inputs (subj_nr, sesID, run_nr).
%  The other input arguments are optional and will be set to default if 
%  undefined. To define optional input arguments, use: 'vararg', <val>.
%
%   subj_nr           : subject number 
%   ses_nr            : session number (numeric integer between 1-27)
%   ses_type          : session type number (numeric integer: 1=A and 2=B)
%   run_nr            : run number     (numeric integer between 1-10)
%   [dispName]        : name of display (string), this affects stimulus size
%                       Choose: '7TAS_BOLDSCREEN32' - BOLD screen at the 7TAS MRI scanner
%                               'KKOFFICE_AOSQ3277' - external monitor in kendrick's CMRR office
%                               'EKHOME_ASUSVE247'  - external monitor at Eline's home
%                               'PPROOM_EIZOFLEXSCAN' - CMRR's Psychophysics room monitor 
%                        Default: '7TAS_BOLDSCREEN32'                      
%   [debugmode]       : if true, there is no eyetracking, no waiting for external
%                       trigger from scanner, no monitor synctest. 
%                        Default: false
%   [loadparams]      : if true, load stored parameter values. 
%                        Default: true
%   [storeparams]     : if true, store created parameter values. 
%                        Default: true
%   [loadstimtiming]  : if true, load temporary stimulus file that was store right before ptb flipping (to save time and rerun the same run) 
%                        Default: true
%   [savestimtiming]  : if true, store temporary stimulus file prior to ptb flipping (to save time and rerun the same run) 
%                        Default: true
%   [savestim]        : if true, store stimuli and timing in temporary file prior to ptb flipping. File is ~1.5 GB!
%                        Default: false
%   [offsetpix]       : offset of center [x,y]-coordinate in pixels. 
%                        Default: No offset [0,0]
%   [movieflip]       : flip presented stimulus display left-right (first argument) or up-down (second argument). 
%                        Default: no flipping: [0,0].
%   [stimdir]         : folder where pre-made stimuli are stored.
%                        Default: fullfile(vcd_rootPath,'workspaces','info')
%   [savedatadir]     : folder where subject's button presses and stim timing are stored. 
%                        Default: [], which will turn into: sprintf('vcd_subj%03d_ses%02d',subj_nr,sesID);
%   [subjfilename]    : pre-fix for behavioral mat-file. 
%                        Default: [], which will turn into: sprintf('%s_vcd_subj%03d_ses%02d_run%02d.mat',ts0,subj_nr,sesID,run_nr);
%   [wanteyetracking] : if true, we will initialize eyetracking. 
%                        Default: false;
%
% __OUTPUTS__
%  None
%
% __EXAMPLE__
%  runme_vcdcore(1, 1, 1,1, 'debugmode', true)
%  runme_vcdcore(1, 1, 1,1, 'debugmode', true, 'dispName','KKOFFICE_AOCQ3277')
%  runme_vcdcore(1, 1, 1,1, 'debugmode', true, 'dispName','KKOFFICE_AOCQ3277', 'savetempstimuli', true)
%  runme_vcdcore(1, 1, 1,1, 'debugmode', true, 'dispName','EKHOME_ASUSVE247', 'savetempstimuli', true)
%
% __DEPENDENCIES__
%  * Psychtoolbox-3 
%  7TAS: v. 3.0.16? or lower?.
%      commit ef093cbf296115badddb995fa06452e34c8c7d02 (HEAD -> master, origin/master, origin/HEAD)
%      Date:   Tue Nov 17 19:46:08 2020 +0100
%
%  PP room: v. 3.0.14 — build date 0ct 3, 2017
% 
% __Author__
%  Written by Eline Kupers @ UMN (kupers@umn.edu)
%
% __HISTORY__
%  2024/12/01 - first version committed to github 
%  YYYY/MM/DD - 
%  

%% %%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%
p = inputParser;
p.addRequired ('subj_nr'        , @isnumeric); % subject number 
p.addRequired ('ses_nr'         , @isnumeric); % session number 
p.addRequired ('ses_type'       , @isnumeric); % session type (1=A or 2=B) 
p.addRequired ('run_nr'         , @isnumeric); % nun number
p.addParameter('dispName'       , '7TAS_BOLDSCREEN32' , @(x) any(strcmp(x, {'7TAS_BOLDSCREEN32','KKOFFICE_AOCQ3277','PPROOM_EIZOFLEXSCAN','EKHOME_ASUSVE247'})))
p.addParameter('debugmode'      , false, @islogical);
p.addParameter('loadparams'     , true, @islogical);
p.addParameter('storeparams'    , true, @islogical);
p.addParameter('savestim'       , false, @islogical);
p.addParameter('loadstimfromrunfile', false, @islogical);
p.addParameter('offsetpix'      , [0 0], @isnumerical); % [x,y]
p.addParameter('movieflip'      , [0 0], @isnumerical); % up/down, left/right
p.addParameter('stimDir'        , fullfile(vcd_rootPath,'workspaces','info'), @ischar);
p.addParameter('savedatadir'    , [], @ischar);
p.addParameter('subjfilename'   , [], @ischar);
p.addParameter('wanteyetracking', false, @islogical);

% Parse inputs
p.parse(subj_nr, ses_nr, ses_type, run_nr, varargin{:});

% Rename variables into general params struct
rename_me = fieldnames(p.Results);
for ff = 1:length(rename_me)
    eval([sprintf('%s = p.Results.%s;', rename_me{ff},rename_me{ff})]);
end
clear rename_me ff p

cd(vcd_rootPath);
addpath(genpath(pwd));

%% Setup display resolution:
% Conservative estimate of what subject can see is h: 800 px, w: 1024 pix.
% This results in atan(800/1024*37.9/2/(5.5+96.1))/pi*180*2 = 16.58 deg 
% total (jointly across two eyes). In the vertical dimension, 
% atan(500/768*28.5/2/(5.5+96.1))/pi*180*2 = 10.43 deg total.
% Note: nsdheightdeg = 12.70 in degrees (BOLDscreen vertical, Nova1x32)

%%%%%%%%%%%%%%%% anon functions
% pix2deg = @(res_pix, res_cm, dist_cm) (atan( ((res_pix/res_cm) / 2) / dist_cm) / pi*180*2);

%%%%%%%%%%%%%%%% path stuff

% addpath(genpath(strrep(which('runvcdcore'),'runvcdcore.m','knkutils')));
% stimulusdir = strrep(which('runvcdcore'),'runvcdcore.m','stimulusfiles');

% how do we map subject ID to subject numbers 1-50?
subjfun = @(x) find(ismember([1:50],x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DON'T EDIT BELOW

% history:
% - 2018/12/25 - initial runnsd version 
% - 2024/11/06 - initial vcd version with some stuff adopted from runnsd

close all; clc;

%%%%%%%%%%%%%%%% need to get info from user
% ask the user what to run
if ~exist('subj_nr','var') || isempty(subj_nr)
  subj_nr = input('What is the subj number? ')
end

if ~exist('ses_nr','var') || isempty(ses_nr)
    % ask the user more stuff
    ses_nr  = input('What session nr (1-27)? ');
end

if ismember(ses_nr,[1,27])
    if ~exist('ses_type','var') || isempty(ses_type)
        % ask the user more stuff
        ses_type  = input('What session type (1 = A and 2 = B)? ');
    end
end

if ~exist('run_nr','var') || isempty(run_nr)
    run_nr = input('What run number (1-12)? ')
end

% give some hints to the operator
if isempty(savedatadir)
    savedatadir = fullfile(vcd_rootPath,'data',sprintf('vcd_subj%03d_ses%02d',subj_nr, ses_nr));
end
if ~exist(savedatadir, 'dir'); mkdir(savedatadir); end


tempfiles = matchfiles(fullfile(savedatadir,sprintf('*_vcd_subj%03d_ses%02d_%s_run%02d.mat',subj_nr, ses_nr, choose(ses_type==1,'A','B'), run_nr)));
fprintf('\nIt appears you have completed the following runs already:\n');
tempfiles


if ~exist('wanteyetracking','var')
    wanteyetracking  = input('Do you want to use eyetracking? (press 1 for yes, 0 for no) ')
end

if ~exist('dispName','var')
    dispName = input('What monitor are you using? (press 1 for 7TAS, 2 for Psychophysics room, 3 for Office) ')
end

fprintf('[%s]: Running subj_nr %03d - session %02d %s - run %02d \n', mfilename, subj_nr,ses_nr,choose(ses_type==1,'A','B'),run_nr)
fprintf('[%s]: %s eyetracking \n', mfilename, choose(wanteyetracking,'YES','NO'))
fprintf('[%s]: Running experiment with images optimized for %s\n', mfilename, dispName)


% Deal with folders and filenames
instructionsDir = fullfile(vcd_rootPath, 'workspaces','instructions');

ts0 = gettimestring;
if isempty(subjfilename)
    behavioralfilename    = sprintf('behavior_%s_vcd_subj%03d_ses%02d_%s_run%02d.mat',ts0,subj_nr,ses_nr,choose(ses_type==1,'A','B'),run_nr);
    
    if wanteyetracking
        eyelinkfilename = fullfile(vcd_rootPath,sprintf('eye_%s_vcd_subj%03d_ses%02d_%s_run%02d.edf',ts0,subj_nr,ses_nr,choose(ses_type==1,'A','B'),run_nr));
    else
        eyelinkfilename = '';
    end
        
else
    behavioralfilename    = fullfile(savedatadir,sprintf('%s_%s.mat',ts0,subjfilename));
    
    if wanteyetracking
        eyelinkfilename = fullfile(savedatadir,sprintf('%s_%s.edf',ts0,subjfilename));
    else
        eyelinkfilename = '';
    end
end
    

% Experiment-specific stuff
switch ses_nr
    
    case num2cell(1:40)
        
        assert(~isempty(subjfun(subj_nr)));
        assert(subjfun(subj_nr)>=1 && subjfun(subj_nr)<=50);
        assert(ses_type>=1 && ses_type<=2);
        assert(ses_nr>=1 && ses_nr<=27);
        assert(run_nr>=1 && run_nr<=10);
        
        setnum = [131 subjfun(subj_nr) ses_nr ses_type run_nr];
        
    otherwise
        error(sprintf('%d is a bad experiment number!!!',expnum));
        
end



% Deal with debug mode
if debugmode
    wanteyetracking = false;
end

% We will load stimuli from scratch 
scan = struct();

%% run experiment
[data, params, getoutearly] = vcd_singleRun(subj_nr, ses_nr, ses_type, run_nr, ... % mandatory inputs
    'debugmode',debugmode, ...
    'behaviorfile',behavioralfilename, ... % optional inputs use: 'var',<val>
    'eyelinkfile',eyelinkfilename, ...
    'savedatadir', savedatadir, ...
    'wanteyetracking', wanteyetracking, ...
    'loadparams', loadparams, ...
    'storeparams', storeparams, ...
    'dispName', dispName, ... 
    'offsetpix', offsetpix, ...
    'movieflip', movieflip, ...
    'instrtextdir',instructionsDir, ...
    'scan', scan, ...
    'savestim', savestim, ...
    'loadstimfromrunfile', loadstimfromrunfile); 

% analyze dat
% runvcdbehavioralanalysis(behaviorfile,3000,[299.8 300.2],[stimulusdir filesep 'vcd_expdesign.mat']);
