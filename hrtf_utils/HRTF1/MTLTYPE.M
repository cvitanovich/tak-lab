function filetype = mtlfiletype(filename)

%  filetype = mtlfiletype(filename)


disp(['file: ',filename]);
[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

fseek(fid,79,'bof');
eof = fread(fid,1,'char');
eof_int = abs(eof);
if (eof_int ~= 26)
  fclose(fid);
  error('choosen file isn''t readable');
end;

fseek(fid,80,'bof');
filetype = fread(fid,1,'short');
if ( (filetype<2)|(filetype>11) )
 fclose(fid);
 error('choosen file includes unknown filefiletype');
end;

return;
