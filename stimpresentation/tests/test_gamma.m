%% test_gamma.m

% Script for a quick and dirty check to see if monitor acts linearly or
% not.

curr_path = pwd; % store current path

% go to root of cvnlab folder with ./calibrationimages
cd('~/Desktop/cvnlab/')

% add paths
vcd_startup

% run calibration script
calibrate

cd(curr_path); % go back to initial path