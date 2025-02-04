function rootPath = vcd_rootPath()
%
% Return the path to root vcd project folder
%
% This function must reside in the directory base of this code repository.
% It is used to determine the location of various subdirectories
%
% Example:
%   fullfile(vcdRootPath, 'stimulus')
%
% By Eline Kupers 2024 @ UMN


rootPath = which('vcd_rootPath');
rootPath = fileparts(rootPath);

return

