function direction_matrix=read_dir(filename,n_directions)

% function direction_matrix=read_dir(filename,n_directions)


[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

fseek(fid,512,'bof');
direction_matrix=fread(fid,2*n_directions,'float32');
reshape(direction_matrix,2,n_directions);

fclose(fid);

return;
