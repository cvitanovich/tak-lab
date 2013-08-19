function [h,header]=readsi1(filename);
% READSR1 reads SI1 file
%         [h,header]=READSI1(filename);

[fid,msg]=fopen(filename,'r');
if fid<0 error(msg); end;
header=fread(fid,256,'char');
h=fread(fid,inf,'short','ieee-le');
fclose(fid);
