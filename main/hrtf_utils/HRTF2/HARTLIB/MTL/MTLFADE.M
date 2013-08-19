function mtlfade(infile,outfile)

% function mtlfade(infile,outfile)


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

epsilon=1E-8;

index=0;

for index=1:n_channels
	clc
	disp(['MTLFADE: Calculating channel ', num2str(index), ' of ', num2str(n_channels)] );
	channel = mtlrch(infile,index);
	channel=channel(:);
%	channel=fade(channel,171,271,839,939);
        logmag=abs2db(channel);
	p=find((logmag<=-150) | (logmag>=0)); channel(p)=ones(size(p))*epsilon;

status = mtlwch(outfile,channel,index);
	end;
return;


                                                                                                                                                                                                                                                          