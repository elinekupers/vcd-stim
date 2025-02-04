%% runvcdcore.m

% Setup display resolution:
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




% how do we map subject ID to subject numbers 1-8?
subjfun = @(x) find(ismember([1 2 3 4 5 6 7 8],x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DON'T EDIT BELOW

% history:
% - 2018/12/25 - initial version from runnsd
% - 20248/11/06 - initial version from runnsd

%%%%%%%%%%%%%%%% historical purposes only:

%%%%%%%%%%%%%%%% some Kendrick MATLAB prep stuff

close all;
kendrickstartup;

%%%%%%%%%%%%%%%% need to get info from user

% ask the user what to run
if ~exist('subjnum','var') 
  subjnum = input('What is the subj id? ')
end

% give some hints to the operator
tempfiles = matchfiles(fullfile(vcd_rootPath,sprintf('*_vcd_subj%d_run*_exp*.mat',subjnum)));
fprintf('\nIt appears you have completed the following runs already:\n');
tempfiles

% ask the user more stuff
expnumNEW = input(['What experiment (session number 1-40)? '])
runnum    = input('What run number (1-12)? ')
trialnum  = input('What trial number to start at (1-75)? (press ENTER for full run) ')

% massage
if isempty(trialnum)
  trialnum = 0;
end

%%%%%%%%%%%%%%%% deal with caching stuff

% NOPE. no caching.
expnum = expnumNEW;

%%%%%%%%%%%%%%%% continue

% there is no pre-loading, so just set these to []
images = [];
maskimages = [];

% calc
ts0 = gettimestring;
filename    = fullfile(vcd_rootPath,sprintf('%s_vcd_subj%d_run%02d_exp%02d.mat',ts0,subjnum,runnum,expnum));
eyefilename = fullfile(vcd_rootPath,sprintf('eye_%s_vcd_subj%d_run%02d_exp%02d.edf',ts0,subjnum,runnum,expnum));

% experiment-specific stuff
switch expnum
    
    
case num2cell(1:40)

  assert(~isempty(subjfun(subjnum)));
  assert(subjfun(subjnum)>=1 && subjfun(subjnum)<=8);
  assert(expnum>=1 && expnum<=40);
  assert(runnum>=1 && runnum<=12);
  assert(trialnum>=0 && trialnum<=75);

  setnum = [131 subjfun(subjnum) expnum runnum trialnum];

otherwise
  error(sprintf('%d is a bad experiment number!!!',expnum));

end

% deal with eyetracking
if wanteyetracking
    tfun = tfunEYE;
    eyefilename = eyefilename;
else
    tfun = tfunNOEYE;
    eyefilename = [];
end

% run experiment
vcd_singleRun(subjnum, runnum)

% single_vcd_run (filename,offset,movieflip,frameduration,fixdotcol,dotdiam_pix,tfun, ...
%                ptonparams,soafun,0,images,setnum,[],grayval,iscolor,[],[], ...
%                [],dres,triggerkey,[],trialparams,eyefilename,maskimages,specialoverlay,stimulusdir, ...
%                [],[],setupscript);

% analyze data
runvcdbehavioralanalysis(filename,3000,[299.8 300.2],[stimulusdir filesep 'vcd_expdesign.mat']);
