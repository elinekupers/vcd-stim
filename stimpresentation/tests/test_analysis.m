[file0,pathname0] = uigetfile('*.mat','Select .mat file');
bfile = fullfile(pathname0,file0);

[file0,pathname0] = uigetfile('*.edf','Select .edf file');
efile = fullfile(pathname0,file0);

b1 = load(bfile);
ptviewmoviecheck(b1.data.timing.timeframes,b1.data.timeKeys,[],{'5' 't'});
behresults = vcdbehavioralanalysis(bfile);
eyeresults = vcdeyetrackingpreprocessing(efile,bfile,behresults);
