function status = mtlwdel2(filename,delay)

% status = mtlwdel2(filename,delay)
%
% MTLWDEL writes delay to catalogue
%
%         filename: string with filename
%         delay: array with delay in lines for each channel
% 

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 21.4.1995


[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 start,stop,n_directions,comment1, comment2] = mtlrh(filename);

if (file_type ~=2) error('MTLWDEL2: wrong file type.'); status=-1; return; end;

n_chan=length(delay);
if (n_chan ~= n_channels) error('MTLWDEL2: wrong number of channels in delay.'); status=-1;
return; end;

n_byte=n_chan*4;
n_blocks=ceil(n_byte/256);


%[fid,message] = fopen(filename,'a','ieee-le');
if exist(filename)==2		% file exists
   [fid,message] = fopen(filename,'r+','l');
else				% file does not exist
   [fid,message] = fopen(filename,'w+','l');
end;


if fid==-1          % file can't be opened
  disp(message);
  status=-1;
  return;
end;


status=fseek(fid,512+n_directions*2*4,'bof');
if (status<0) error('MTLWDEL2: fseek failed.'); status=-1; return; end;

count=fwrite(fid,delay,'float32');
if (count ~= n_chan) error('MTLWDEL2: fwrite failed.'); status=-1; return; end;

status=fclose(fid);
if (status<0) error('MTLWDEL2: fclose failed.'); status=-1; return; end;

return;

end;
