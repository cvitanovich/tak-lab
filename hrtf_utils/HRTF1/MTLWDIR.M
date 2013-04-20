function status=mtlwdir(filename,direction_matrix)

% status=mtlwdir(filename,direction_matrix)
%
% MTLWDIR writes direction matrix of a mtlsig file
%
% status: negative values indicate error
%
% filename: name of the mtlsig file
% direction_matrix: matrix with n-directions columns
%                   first line: elevation
%                   second line: azimuth
%

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 21.4.1995

status=0;

%[fid,message] = fopen(filename,'a+','ieee-le');

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

[m,n]=size(direction_matrix);
num_val=m*n;
status=fseek(fid,512,'bof');
if (status<0) error('ERROR: mtlwdir, fseek to end of info_block');status=-1;return, end;
count=fwrite(fid,direction_matrix,'float32');
if (count<num_val) error('ERROR: mtlwdir, fwrite');status=-1;return, end;

status=fclose(fid);
if (count<num_val) error('ERROR: mtlwdir, fclose');status=-1;return, end;

return;
end;