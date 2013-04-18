function [hl,hr]=readsig2(filename);
% READSIG2 reads SIG2 file
%          [hl,hr]=READSIG2(filename);

[fid,msg]=fopen(filename,'r','ieee-be');
if fid<0 error(msg); end;
h=fread(fid,inf,'float32');
hl=h(1:4096); hr=h(4097:8192);
fclose(fid);