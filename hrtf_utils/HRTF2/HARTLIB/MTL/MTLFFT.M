function mtlfft(infile,outfile,len)

% function mtlfft(infile,outfile,len)



[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(infile);


if (filetype ~= 2) 
  error(['File: ' infile ' has wrong type']);
  return;
end;


if(nargin<2) 
 filename=infile;
 strrep(filename,'.sr','.sc')
 strrep(filename,'.SR','.sc')
 outfile=filename;
end;


if(nargin<3) 
len=n_lines;
end;

n_lines=len;

filetype=3;

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
  disp(['MTLFFT: Calculating channel ', num2str(index), ' of ', num2str(n_channels)] );

   channel=mtlrch(infile,index);
   f_channel=rfft(channel,len);
   if (length(f_channel)>n_lines) disp('mtlfft: problem with channel_length'); end;
   status = mtlwch(outfile,f_channel,index);
end;
return;


																																																												


                                                                                                                                                                                                                                                          
