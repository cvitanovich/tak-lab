function mtlwt4 (filename, STIM, index)
% mtlwt4 (filename, STIM, index)
%
% writes STIM in filename 
% PC format, 16 bit, LRLR

n_byte = 2;

[m,n] = size(STIM);

[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(filename);


if (nargin  < 3) pos = info_blocks*256; index = 1:m; end;
if (nargin == 3) pos = info_blocks*256+(index-1)*n_lines*2*n_byte; end;

[fid,message] = fopen(filename,'a+','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  channel=-1;
  return;
end;

status=fseek(fid,pos,'bof');
if (status<0) error('ERROR: mtlwch,  fseek to channel position not succesful' ); status=-1;return; end;



multiplex = zeros(2*m,1);

for i=1:n/2
 multiplex(1:2:2*m,1) = STIM(:,2*i-1);
 multiplex(2:2:2*m,1) = STIM(:,2*i);
 count = fwrite(fid, multiplex(:,1), 'short');
 if (count<n_lines*2) error('ERROR: mtlwch, writing of channel not succesful');  status=-1; return; end;
end;

end

fclose(fid);

return
