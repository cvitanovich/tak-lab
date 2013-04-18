function mtlsf(infile, outfile, new_sampling_rate)


% function mtlsf(infile, outfile, new_sampling_rate)

[comment1, file_type, info_blocks,...
           n_channels, day_of_date, month_of_date, year_of_date,...
           no_bl_per_chan, subject_no, ...
           subject_first_name, subject_last_name, ...
           n_lines, sampling_rate,...
           lo_cutoff, hi_cutoff, eff_s_to_n, no_of_aver, ... 
           no_of_qu_bits, chan_code, ampl_gain, ...
           max_ele, min_ele, status, res_type, res, ...
           channels_per_dir, map, comment2] = mtlrh(infile);

n_directions = n_channels/channels_per_dir;
if (file_type ~= 2) error('Wrong filetype'); return; end;

factor = new_sampling_rate/sampling_rate
n_lines
n_lines_out=n_lines*factor

end;

last_line_out=n_lines_out/sampling_rate;

index=0;

%if direction_matrix
if ((info_blocks-2) > 0)
   direction_matrix=mtlrdir(infile);
   if (index ~= 0) direction_matrix(:) = direction_matrix(:,index); end;
   % if time_delay exists read time_delay
   if (( (info_blocks-2)*256 - n_directions*8) >= (n_channels*4))
      delay=mtlrdel(infile);
     if (index ~= 0) 
        delay= delay(index); 
     end;
   end;
end;


 




% write output_header


% write directions;
% write delays;
if ((n_channels/n_directions) == 2) binaural=1; end;
n_dir=n_directions;
if (n_dir > 0)
  n_byte=n_dir*2*4;
  sn_del=length(delay);
  if (sn_del > 0) n_byte=n_byte+n_dir*4; end; 
  n_blocks=ceil(n_byte/256);
end;

info_blocks=2+n_blocks;


 mtlwh(outfile,...
           comment1, file_type, info_blocks,...
           n_channels, day_of_date, month_of_date, year_of_date,...
           no_bl_per_chan, subject_no, ...
           subject_first_name, subject_last_name, ...
           n_lines_out, new_sampling_rate,...
           lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
           no_of_qu_bits, chan_code, ampl_gain, ...
           max_ele, min_ele, status, res_type, res, ...
           channels_per_dir, map, comment2);

 


[fid,message] = fopen(outfile,'a','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

 if(length(delay) > 0) 
   mtlwdel(outfile, delay*factor);   
end;
end



for n=1:n_channels
%for n=1:4
close all
clc
disp(['MTLSF: : channel ', num2str(n), ' of ', num2str(n_channels)] );
 
channel=mtlrch(infile,n);
%irplot(channel,sampling_rate);
n_channel=resample(channel,new_sampling_rate,sampling_rate)/factor;
 mtlwch(outfile,n_channel,n);
%pause
end; 
