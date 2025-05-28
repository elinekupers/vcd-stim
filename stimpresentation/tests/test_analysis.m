bfile = 'behavior_20250528104429_vcd_subj001_ses01_A_run01.mat';
efile = 'eye_20250528104429_vcd_subj001_ses01_A_run01.edf';
behresults = vcdbehavioralanalysis(bfile);
eyeresults = vcdeyetrackingpreprocessing(efile,bfile,behresults);
