function  writeh(filename,filetype,info_blocks,...
          n_channels,n_lines, sampling_rate,...
          first_line,last_line,n_directions,...
          comment1, comment2)



%function  writeh(filename,filetype,info_blocks,...
%          n_channels,n_lines, sampling_rate,...
%          first_line,last_line,n_directions,...
%          comment1, comment2)
% WRITEH writes the header of a SIGMAP and MTLSIG file.


disp(['file: ',filename]);
[fid,message] = fopen(filename,'w+','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

d=zeros(1,512);
fwrite(fid,d,'uchar');

message = 'Fehler beim Schreiben';
disp(message);

if (filetype<=3)
  status=fseek(fid,0,'bof');
  
   n = length (comment1);
   comment1=[comment1 zeros(1,79-n)];
   fwrite(fid,comment1,'uchar');
end;

status=fseek(fid,79,'bof');
eof=26;
count = fwrite(fid,eof,'char');



status=fseek(fid,80,'bof');
count = fwrite(fid,filetype,'short');



status=fseek(fid,82,'bof');
count = fwrite(fid,info_blocks,'short');


status=fseek(fid,84,'bof');
count = fwrite(fid,n_channels,'short');


status=fseek(fid,128,'bof');
count = fwrite(fid,n_lines,'long');




if (filetype<7)
  status=fseek(fid,132,'bof');
  count = fwrite(fid,sampling_rate,'float32');
  end;
  status=fseek(fid,136,'bof');
  count = fwrite(fid,first_line,'float32');
  
  status=fseek(fid,140,'bof');
  count = fwrite(fid,last_line,'float32');
  
end;


status=fseek(fid,220,'bof');
count = fwrite(fid,n_directions,'short');




if (filetype >3)
  status=fseek(fid,222,'bof');
  comment2 = setstr(fread(fid,200,'uchar'));
  fwrite(fid,comment2,'uchar');
end;

status= fseek(fid,511,'bof');
fwrite(fid,'a','char');
n=fclose(fid);
return;
