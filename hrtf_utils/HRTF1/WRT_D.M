function direction_matrix=read_dir(filename,n_directions)

[fid,message] = fopen(filename,'w');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

fseek(fid,512,'bof');
fwrite(fid,direction_matrix,'float32');

fclose(fid);

return;

