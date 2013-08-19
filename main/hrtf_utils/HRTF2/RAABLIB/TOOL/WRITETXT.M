function writetxt(fname,str)
% WRITETXT writes a text string to a file.
%          WRITETXT(fname,str)

[fid,msg]=fopen(fname,'wt');
if fid<0 error(msg); end;
fprintf(fid,'%s',str);
fclose(fid);
