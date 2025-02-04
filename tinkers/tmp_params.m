% 
% if iscolor
%     [timeframes,timekeys,digitrecord,trialoffsets] = ptviewmovie(...
%         images, ... .mat  uint8 file A x B x 1/3 x N matrix
%         frameorder,... vector of positive integers refer to specific images in <images>. 0=no image (just the background and fixation).
%         framecolor,... default is 255*ones(size(<frameorder>,2),3)
%         frameduration,... how many monitor refreshes you want a single movie frame to last.  default: 15.
%         fixationorder,... check options. default: ones(1,1+size(<frameorder>,2)+1).
%         fixationcolor,...
%         fixationsize,...
%         grayval,...
%         [],...
%         [],...
%         offset,...
%         choose(con==100,[],1-con/100),...
%         movieflip,...
%         scfactor,...
%         [],....
%         triggerfun,...
%         framefiles,...
%         [], ...
%         triggerkey,...
%         specialcon,...
%         trialtask,...
%         maskimages,...
%         specialoverlay);
% else
%     % OLD AND WASTEFUL: reshape(cat(3,images{:}),size(images{1},1),size(images{1},2),1,[])
%     [timeframes,timekeys,digitrecord,trialoffsets] = ptviewmovie(images, ...
%         frameorder,framecolor,frameduration,fixationorder,fixationcolor,fixationsize,grayval,[],[], ...
%         offset,choose(con==100,[],1-con/100),movieflip,scfactor,[],triggerfun,framefiles,[], ...
%         triggerkey,specialcon,trialtask,maskimages,specialoverlay);
% end



% function [images,maskimages] = showmulticlass(outfile,offset,movieflip,frameduration,fixationinfo,fixationsize, ...
%   triggerfun,ptonparams,soafun,skiptrials,images,setnum,isseq,grayval,iscolor, ...
%   numrep,con,existingfile,dres,triggerkey,framefiles,trialparams,eyelinkfile,maskimages,specialoverlay)
%
% <outfile> is the .mat file to save results to
% <offset> is horizontal and vertical offset for display purposes (see ptviewmovie.m)
% <movieflip> is flip to apply to movie (see ptviewmovie.m)
% <frameduration> is number of monitor refreshes for one movie frame (see ptviewmovie.m)
% <fixationinfo> is
%   {A B C} where A is the base color (uint8 1x3) for the fixation dot (see fixationcolor in ptviewmovie.m)
%                 B is alpha for the fixation dot in the default case
%                 C is alpha for the fixation dot when flips happen
%   {D E} where D is a set of colors (uint8 Nx3) (see the negative-integers case for fixationcolor in ptviewmovie.m)
%               E is an alpha value in [0,1]
%   {F} where F is the {A B C D E F G} case of <fixationorder> in ptviewmovie.m
%     note that F should really only be {A B C D E F} since we will add G if <existingfile> is supplied.
% <fixationsize> is size in pixels for fixation dot or an entire alpha image for the fixation (see ptviewmovie.m)
% <triggerfun> is the trigger function to call when starting the movie (see ptviewmovie.m)
% <ptonparams> is a cell vector with parameters for pton.m
% <soafun> is a function that returns a (presumably stochastic) stimulus-onset asynchrony for the
%   fixation dot (in number of movie frames).  the output must be a positive integer.
%   <soafun> is ignored when <fixationinfo> is {F}.
% <skiptrials> is number of trials to skip at the beginning
% <images> (optional) is a speed-up (no dependencies).  <images> can be reused only if <setnum>
%   stays within [1 2 3 7 8]
% <setnum> (optional) is
%   1 means the original 31 stimulus classes [15 frames, 3s / 3s]
%   default: 1.
% <isseq> (optional) is whether to do the special sequential showing case.  should be either 0 which
%   means do nothing special, or a positive integer indicating which frame to use.  if a positive
%   integer, <fixationinfo> and <soafun> are ignored.  default: 0.
% <grayval> (optional) is the background color as uint8 1x1 or 1x3 (see ptviewmovie.m).  default: uint8(127).
% <iscolor> (optional) is whether to expect that the images are color.  default: 0.
% <numrep> (optional) is number of times to repeat the movie.  note that fixation stuff
%   is not repeated, but stochastically generated.  default: 1.
% <con> (optional) is the contrast in [0,100] to use (achieved via <moviemask> in ptviewmovie.m).
%   default: 100.  NOTE: this currently has a slow implementation (the initial setup time is long).
% <existingfile> (optional) is an old <outfile>.  if supplied, we pull the 'framedesign',
%   'classorder', 'fixationorder', 'trialoffsets', and 'digitrecord' variables from this old file
%   instead of computing it fresh.  prior to May 13 2013, the trialtask and the digit-stream
%   were not preserved.  now, the 'trialoffsets' and 'digitrecord' means that they are preserved!
% <dres> (optional) is
%   [A B] where this is the desired resolution to imresize the images to (using bicubic interpolation).
%     if supplied, the imresize takes place immediately after loading the images in, and this imresized
%     version is what is cached in the output of this function.
%  -C where C is the <scfactor> input in ptviewmovie.m
%   default is [] which means don't do anything special.
% <triggerkey> (optional) is the input to ptviewmovie.m
% <framefiles> (optional) is the input to ptviewmovie.m
% <trialparams> (optional) is the {B E F G H} of the <trialtask> inputs
%   to ptviewmovie.m.  specify when <setnum> is 51 or 52 or 53 or 54 or 56,57,58, 59,60,61,
%   62,63,64,65, 66,  78,79
% <eyelinkfile> (optional) is the .edf file to save eyetracker data to.
%   default is [] which means to do not attempt to use the Eyelink.
% <maskimages> (optional) is a speed-up (no dependencies).  <maskimages> is applicable and can be reused only
%   if <setnum> stays within [73 74 75 76 77] or [78 79] or [89 90 91 92 93].
% <specialoverlay> (optional) is the input to ptviewmovie.m
%
% show the stimulus and then save workspace (except the variable 'images') to <outfile>.

% history:

outfile,offset,movieflip,frameduration,fixationinfo,fixationsize, ...
    triggerfun,ptonparams,soafun,skiptrials,images,setnum,isseq,grayval,iscolor, ...
    numrep,con,existingfile,dres,triggerkey,framefiles,trialparams,eyelinkfile,maskimages,specialoverlay)

%%%%%%%%%%%%% input

% if ~exist('setnum','var') || isempty(setnum)
%     setnum = 1;
% end
% if ~exist('isseq','var') || isempty(isseq)
%     isseq = 0;
% end
% if ~exist('grayval','var') || isempty(grayval)
%     grayval = uint8(127);
% end
% if ~exist('iscolor','var') || isempty(iscolor)
%     iscolor = 0;
% end
% if ~exist('numrep','var') || isempty(numrep)
%     numrep = 1;
% end
% if ~exist('con','var') || isempty(con)
%     con = 100;
% end
% if ~exist('existingfile','var') || isempty(existingfile)
%     existingfile = [];
% end
% if ~exist('dres','var') || isempty(dres)
%     dres = [];
% end
% if ~exist('triggerkey','var') || isempty(triggerkey)
%     triggerkey = [];
% end
% if ~exist('framefiles','var') || isempty(framefiles)
%     framefiles = [];
% end
% if ~exist('trialparams','var') || isempty(trialparams)
%     trialparams = [];
% end
% if ~exist('eyelinkfile','var') || isempty(eyelinkfile)
%     eyelinkfile = [];
% end
% if ~exist('maskimages','var') || isempty(maskimages)
%     maskimages = [];
% end
% if ~exist('specialoverlay','var') || isempty(specialoverlay)
%     specialoverlay = [];
% end
% if ~isempty(existingfile)
%     efile = load(existingfile,'framedesign','classorder','fixationorder','trialoffsets','digitrecord');
% end

% numinclass = cellfun(@(x) size(x,choose(iscolor,4,3)),images);  % a vector with number of images in each class

%%%%%%%%%%%%% perform run-specific randomizations (NOTE THAT THERE ARE HARD-CODED CONSTANTS IN HERE)

% load in some aux info
% switch setnum(1)
% case {6}
%   load(infofile_version2,'lettersix','numbersix','polygonsix');
% case {1 2 3 4 5 7 8}
%   load(infofile,'lettersix','numbersix','polygonsix');
% end

% figure out frame assignment for each class
% if ~isempty(existingfile)
%     framedesign = efile.framedesign;
% else
%     switch setnum(1)
%         case 1  %
%             if isseq
%                 framedesign = {};
%                 for p=1:73
%                     framedesign{p} = isseq;
%                 end
%             end
%     end
% end

%%%%%%%%%%%%% more preparations

% load in some aux info
% if isseq
%     switch setnum(1)
%         case 1
%             trialpattern = eye(31);
%             onpattern = [1];
%     end
% end

% skip trials?
% if exist('trialpattern','var')
%     trialpattern = trialpattern(skiptrials+1:end,:);
% end

% decide assignment of classes to the events in the experimental design
% if ~isempty(existingfile)
%     classorder = efile.classorder;
% else
%     switch setnum(1)
%         case 1
%             classorder = [1 2 3 4 5 6 7 9 11 8 10 12 13:31];  % e.g., event 1 will be classorder(1)
%     end
% end