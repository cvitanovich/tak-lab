function []=writesi1(data,filename,header);
% WRITESI1  Write data to file using SI1 format.
%           []=WRITESI1(data,filename,header);

[fid,msg]=fopen(filename,'w');
if fid<0 error(msg); end;
l=fwrite(fid,header,'char');
if l~=length(header) error('header: write fault'); end;
l=fwrite(fid,data,'short');
if l~=length(data) error('data: write fault'); end;
fclose(fid);
