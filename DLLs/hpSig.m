function [] = hpSig(string1)

% HELPFILE for hpSig(string1)
% this is a mexfile dll
%
% Sends the argin string to the 
% HP35665A signal analyzer.
% A trailing "\n" (line feed) is automatically added.
% Maximum string length is 256 characters.
% String must be in a row vector - 
% multiple rows are read columnwise and
% made into a single long row.

