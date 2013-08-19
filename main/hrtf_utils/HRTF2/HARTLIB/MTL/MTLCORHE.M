function mtlcorhe(filename)


%function mtlcorhe(filename)

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 21.4.1995


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




sampling_rate=24000;
n_lines=2048;
subject_no=802;
lo_cutoff=100;
hi_cutoff=12000;
no_of_aver=100;
res_type=BIFI_RES_WITH_DIRECTION_INFO;
channels_per_dir = 2;
n_dir = n_channels/channels_per_dir;

[fid,message] = fopen(filename,'a+','ieee-le');
if (fid==-1)         % file can't be opened
  disp(message);
  filetype=-1;
  return;
end;

map =  zeros(2,n_dir);
stat=fseek(fid,512,'bof');
if (stat<0) error('ERROR: mtlrh, fseek map for multiple directions'); filetype=-1;return, end;
map = fread(fid,n_dir*2,'float');

map = reshape(map,2,n_dir);

% stat=fseek(fid,240,'bof');
% if (stat < 0) error('ERROR: mtlwh, map for single direction, fseek'); return;  %end;
% count = fwrite(fid,map(:,1),'float');
% if (count < length(map(:,1))) error('ERROR: mtlwh, map for single direction , %fwrite'); return; end;

mode='a';

mtlwh(filename,...
comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2,mode);


return;


