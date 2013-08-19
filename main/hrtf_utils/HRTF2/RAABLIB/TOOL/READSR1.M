function [h,header]=readsr1(filename);
% READSR1 reads SR1 file
%         [h]=READSR1(filename);

[fid,msg]=fopen(filename,'r','ieee-le');
if fid<0 error(msg); end;
header=fread(fid,256,'char');
h=fread(fid,inf,'float32');
fclose(fid);
