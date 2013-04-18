function mtlifft(infile,outfile)

% function mtlifft(infile,outfile)
% 
% Performs IFFT on a BIFI-File
%
% Klaus Hartung
% Lehrtstuhl f. allg. Elektrotechnkik und Akustik
% Ruhr-Universitaet Bochum
% 1994
% Modified
% 12.5.1996 Klaus Hartung


if (nargin<2) 
outfile=strrep(infile,'.sc','.sr')
outfile=strrep(infile,'.SC','.sr')
end;


[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(infile);





if (filetype ~= 3) 
  error(['File: ' infile ' has wrong type']);
  return;
end;





filetype=2;

mtlwh(outfile,...
comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2);


for index=1:n_channels
clc
disp(['MTLIFFT: Calculating channel ', num2str(index), ' of ', num2str(n_channels)] );
channel = mtlrch(infile,index);
channel=channel(:);
t_channel = rifft(channel);
status = mtlwch(outfile,t_channel,index);

end;
return;
                                                                                                                                                                                                                                                          
