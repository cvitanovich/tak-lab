function status = mtlwdel(filename,delay_l,delay_r)

%function status = mtlwdel(filename,delay_l,delay_r)

binaural=1;
if (nargin < 3) binaural=0; end;

%[fid,message] = fopen(filename,'a','ieee-le');
if exist(filename)==2		% file exists
   [fid,message] = fopen(filename,'r+','l');
else				% file does not exist
   [fid,message] = fopen(filename,'w+','l');
end;


if fid==-1          % file can't be opened
  disp(message);
  return;
end;

n_directions=length(delay_r);
n_byte=n_directions*4*(binaural+1);
n_blocks=ceil(n_byte/256);

fseek(fid,512+n_directions*2*4,'bof');
fwrite(fid,[delay_l;delay_r],'float32');


fclose(fid);

return;
