function status=mtlwdir(filename,map)

% status=mtlwdir(filename,map)
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
% Date: 19.1.1996

SIGM_BLOCK_LEN = 256;
DOS_EOF = 26;

SO_EOF = 79;
SO_TYPE = 80;
SO_N = 128;
SO_CHA = 84;

INT_TYPE = 1;
FLOAT_TYPE = 2;
COMPLEX_TYPE = 3;
ADDA_TYPE = 4;


BIFI_MAGIC_NUMBER = hex2dec('7caffee1');

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


status=0;

[fid,message] = fopen(filename,'a+','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  status=-1;
  return;
end;

% read number of info_blocks

stat=fseek(fid,82,'bof');
if (stat < 0) error('ERROR: mtlrh, fseek info_blocks'); filetype=-1;return, end;

info_blocks = fread(fid,1,'ushort')

info_blocks
[m,n]=size(map);

n_dir = n;

if (n_dir == 1) 
 stat=fseek(fid,240,'bof');
 if (stat < 0) error('ERROR: mtlwh, map for single direction, fseek'); return;  end;
 count = fwrite(fid,map(:,1),'float');
 if (count < length(map(:,1))) error('ERROR: mtlwh, map for single direction , fwrite'); return; end;
else

 n_dir 
 info_blocks
 if ((n_dir > 1) & (info_blocks == 2 + ceil(n_dir*8/SIGM_BLOCK_LEN)))
    stat=fseek(fid,512,'bof');
    if (stat < 0) error('ERROR: mtlwh, map for multiple directions, fseek'); return;  end;
   count = fwrite(fid,map(:,1:n_dir),'float');
   if (count < length(map(:,1:n_dir))) error('ERROR: mtlwh, map for multiple directions , fwrite'); return; end;
 else
   error('ERROR: mtlwh, map for multiple directions, not enough info_blocks'); return;  end;
 end
end


status=fseek(fid,512,'bof');
if (status<0) error('ERROR: mtlwdir, fseek to end of info_block');status=-1;return, end;

num_val=m*n;
count=fwrite(fid,map,'float32');
if (count<num_val) error('ERROR: mtlwdir, fwrite');status=-1;return, end;

status=fclose(fid);
if (status<0) error('ERROR: mtlwdir, fclose');status=-1;return, end;

return;
end;
