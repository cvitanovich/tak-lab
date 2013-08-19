function mtlmask(infile,outfile, lb, ub, rv)

% function mtlmask(infile,outfile)

% MTLMASK relplaces all values which are below lb (in dB) and ub (in dB) 
% by rv (in dB)
% 
% infile: name of input file
% outfile: name of output file
% lb: lower bound
% ub: upper bound
% rv: value to place instead of old values

% Klaus Hartung
% Lehrstuhl f. allg. Elektrotechnik und Akustik
% Ruhr-Universtitaet Bochum
% 12.5.1996


if (nargin<2) 
disp(['outfile not specified']);
  return
end;


[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(infile,0);

if (filetype ~= 3) 
  error(['File: ' infile ' has wrong type']);
  return;
end;



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

epsilon=db2abs(rv);

index=0;

for index=1:n_channels
	clc
	disp(['MTLFADE: Calculating channel ', num2str(index), ' of ', num2str(n_channels)] );
	channel = mtlrch(infile,index);
	channel=channel(:);
        logmag=abs2db(channel);
	p=find((logmag<=lv) | (logmag>=up)); channel(p)=ones(size(p))*epsilon;
	status = mtlwch(outfile,channel,index);
	end;
return;


                                                                                                                                                                                                                                                          