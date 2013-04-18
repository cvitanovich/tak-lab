function direction_matrix=wrt_d(filename,direction_matrix)


%function direction_matrix=wrt_d(filename,direction_matrix)

[fid,message] = fopen(filename,'a','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

sn=size(direction_matrix);
n_directions=sn(2);
n_byte=n_directions*2*4;
n_blocks=ceil(n_byte/256);

fseek(fid,512,'bof');
d=zeros(1,256*n_blocks);
fwrite(fid,d,'uchar');

fseek(fid,512,'bof');
fwrite(fid,direction_matrix,'float32');

fclose(fid);

return;
