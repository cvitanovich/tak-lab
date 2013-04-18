function map=mtlrdir(filename)

% map=mtlrdir(filename)
%
% MTLRDIR reads direction matrix of a mtlsig file
%
% direction_matrix: matrix with n-directions columns
%                   first line: elevation
%                   second line: azimuth
%
% filename: name of the mtlsig file
%
%
%
%
% See also: MTLWH, MTLRCH, MTLWCH

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 19.4.1995
% Date: 21.4.1995
% Date: 19.1.1996
% Date: 12.5.1996
% status mask

BIFI_MEASURED     = 1;
BIFI_X_TESTSIGNAL = 2;
BIFI_SPEAKER_EQ   = 4;
BIFI_TOTAL_EQ     = 8;

BIFI_RES_NUMBER = 5;

BIFI_RES_AS_MEASURED         = 1;
BIFI_RES_CONSTANT_AZIMUTH    = 2;
BIFI_RES_CONSTANT_DENSITY    = 3;
BIFI_RES_OLD_FULL_SPHERE     = 4;
BIFI_RES_WITH_DIRECTION_INFO = 5;



[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  direction_matrix=-1;
  return;
end;

% read number of channels

stat=fseek(fid,84,'bof');
if (stat<0) error('ERROR: mtlrh, fseek n_channels'); filetype=-1;return, end;
n_channels = fread(fid,1,'ushort');

% read bifi res_type
stat=fseek(fid,234,'bof');
if (stat<0) error('ERROR: mtlrh, fseek res type'); filetype=-1;return, end;
res_type = fread(fid,1,'ushort');

if (res_type == BIFI_RES_WITH_DIRECTION_INFO)

 % read channels per direction
 stat=fseek(fid,238,'bof');
 if (stat<0) error('ERROR: mtlrh, fseek channels per direction');     filetype=-1;return, end;
 channels_per_dir = fread(fid,1,'ushort');

n_dir = n_channels/channels_per_dir

if (n_dir == 1)
   map = zeros(2,1);
   stat=fseek(fid,240,'bof');
   if (stat<0) error('ERROR: mtlrh, fseek map for single direction'); filetype=-1;return, end;
   map = fread(fid,2,'float');
   map = reshape(map,2,1);
end

if (n_dir > 1)
   map =  zeros(2,n_dir);
   stat=fseek(fid,512,'bof');
   if (stat<0) error('ERROR: mtlrh, fseek map for multiple directions'); filetype=-1;return, end;
   map = fread(fid,n_dir*2,'float');
   map = reshape(map,2,n_dir);
 end
 
else
 error(' mtlrdir: File is not BIFI_RES_WITH_DIRECTION_INFO'); return;
end % of BIFI_RES_TYPE


status=fclose(fid);

if (status<0) error('ERROR: mtlrdir, close '); map=-1;return, end;

return;