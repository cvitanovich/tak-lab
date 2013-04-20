function direction_matrix=mtlrdir(filename,n_directions)

% direction_matrix=mtlrdir(filename,n_directions)
%
% MTLRDIR reads direction matrix of a mtlsig file
%
% direction_matrix: matrix with n-directions columns
%                   first line: elevation
%                   second line: azimuth
%
% filename: name of the mtlsig file
%
% optional arguments:
%
% n_directions: number of directions
%
%
% See also: MTLWH, MTLRCH, MTLWCH

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995
% Date: 21.4.1995


[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  direction_matrix=-1;
  return
end

if (nargin < 2)
   status=fseek(fid,220,'bof');
   if (status<0) error('ERROR: mtlrdir, fseek n_directions'); 
       direction_matrix=-1;
       return, 
   end
   n_directions = fread(fid,1,'short');
   if (n_directions<=0) direction_matrix=-1;
       return; 
   end
end


status=fseek(fid,512,'bof');
if (status<0) 
    error('ERROR: mtlrdir, fseek n_directions'); 
    direction_matrix=-1;
    return; 
end

direction_matrix=fread(fid,2*n_directions,'float32');

status=fclose(fid);
if (status<0) error('ERROR: mtlrdir, fseek n_directions'); 
    direction_matrix=-1;
    return; 
end

direction_matrix=reshape(direction_matrix,2,n_directions);