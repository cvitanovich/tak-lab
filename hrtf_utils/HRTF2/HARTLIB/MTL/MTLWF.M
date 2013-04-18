function status = mtlwf(filename,frq)

%function status = mtlwf(filename,frq)

binaural=1;

if (nargin < 3) binaural=0; end;

           [filetype,info_blocks,...
           n_channels,n_lines, sampling_rate,...
           first_line,last_line,n_directions,...
         comment1, comment2] = mtlrh(filename);

[fid,message] = fopen(filename,'a','ieee-le')
if fid==-1          % file can't be opened
  disp(message);
  return;
end;


status=fseek(fid,512+n_directions*2*4,'bof');

fwrite(fid, frq, 'float32');

fclose(fid);

return;
