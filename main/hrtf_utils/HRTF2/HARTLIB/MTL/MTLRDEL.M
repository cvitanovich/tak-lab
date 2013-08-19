function [delay]=mtlrdel(filename)

% function [delay_l]=mtlrdel(filename)
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

BIFI_RES_NUMBER = 5;

BIFI_RES_AS_MEASURED         = 1;
BIFI_RES_CONSTANT_AZIMUTH    = 2;
BIFI_RES_CONSTANT_DENSITY    = 3;
BIFI_RES_OLD_FULL_SPHERE     = 4;
BIFI_RES_WITH_DIRECTION_INFO = 5;

 


[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(filename);


if (res_type == BIFI_RES_WITH_DIRECTION_INFO)
  n_dir = n_channels/channels_per_dir;
  pos = n_dir * 4 * 2 + 512;
else
  pos = 512;
end;

[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

fseek(fid,pos,'bof');

delay=fread(fid,n_channels,'float32');
fclose(fid);


return;
