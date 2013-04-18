function [filetype,info_blocks,...
          n_channels,n_lines, sampling_rate,...
          first_line,last_line,n_directions,...
          comment1, comment2] = readh(filename)



%function [filetype,info_blocks,...
%          n_channels,n_lines, sampling_rate,...
%          first_line,last_line,n_directions,...
%          comment1, comment2] = readh(filename)
%
% READ_H Reads the header of a SIGMAP and MTLSIG file.


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
 error('choosen file includes unknown filetype');
end;


fseek(fid,82,'bof');
info_blocks = fread(fid,1,'short');


fseek(fid,84,'bof');
n_channels = fread(fid,1,'short');


fseek(fid,128,'bof');
n_lines = fread(fid,1,'long');



if (filetype<7)
  fseek(fid,132,'bof');
  sampling_rate = fread(fid,1,'float32');
  fseek(fid,136,'bof');
  first_line = fread(fid,1,'float32');
  fseek(fid,140,'bof');
  last_line = fread(fid,1,'float32');
end;


fseek(fid,220,'bof');
n_directions = fread(fid,1,'short');

if (filetype<=3)
  fseek(fid,0,'bof');
  comment1 = setstr(fread(fid,80,'uchar'));    
end;

if (filetype >3)
  fseek(fid,222,'bof');
  comment2 = setstr(fread(fid,200,'uchar'));
end;

%added new lines Dec 1, 1995:
fseek(fid,422,'bof');
birdnum=fread(fid,1,'short');
disp(birdnum);
fseek(fid,512,'bof');
elev=fread(fid,1,'short');
disp(elev);
fseek(fid,516,'bof');
az=fread(fid,1,'short');
disp(az);
%end of added lines

fclose(fid);
return;

