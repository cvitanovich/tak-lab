function [hl,hr]=readsr2(filename);
% READSR2 reads SR2 file
%         [hl,hr]=READSR2(filename);

[fid,msg]=fopen(filename,'r','ieee-le');
if fid<0 error(msg); end;
fseek(fid,256,'bof');
h=fread(fid,inf,'float32');
hl=h(1:1024); hr=h(1025:2048);
fclose(fid);
