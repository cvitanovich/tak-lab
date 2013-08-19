function  channel = mtlrt4(filename, index)

%MTLRT4 reads multiplexed array (LRLRL) (type 4)

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 22.1.1996


if (nargin < 2) index=1; end;


type = 'short';   n_byte = 2;

[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(filename);



[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  error(message);
  channel=-1;
  return;
end;

pos=info_blocks*256+n_lines*n_byte*(index-1);


status=fseek(fid,pos,'bof');
if (status<0) error('ERROR: mtlrch,  fseek to channel not succesful' );channel=-1; return; end;


channel=fread(fid,n_lines*2,type);
%if (channel<0) error('ERROR: mtlrch, reading of channel not succesful' );channel=-1; end;

status = fclose(fid);

% demultiplexing array
l = length(channel);

channel = reshape(channel,2,l);

