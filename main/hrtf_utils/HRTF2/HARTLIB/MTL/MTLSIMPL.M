function mtlsimplesf(infile, outfile, new_sampling_rate)


% function mtlsf(infile, outfile, new_sampling_rate)

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,comment1, comment2] = mtlrh(infile)


if (file_type ~= 2) error('Wrong filetype'); return; end;

factor = new_sampling_rate/sampling_rate
n_lines
n_lines_out=n_lines*factor
n_directions=n_channels/2;
end;
pause
last_line_out=n_lines_out/sampling_rate;

index=0;

% write output_header

mtlwh(outfile,file_type,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line_out,n_directions,...
       comment1, comment2);



% write directions;
% write delays;
if ((n_channels/n_directions) == 2) binaural=1; end;
n_dir=n_directions;


mtlwh(outfile,file_type,info_blocks,n_channels,n_lines_out, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);


[fid,message] = fopen(outfile,'a','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;



for n=1:n_channels
clc
disp(['MTLSF: : channel ', num2str(n), ' of ', num2str(n_channels)] );
 
channel=mtlrch(infile,n);
n_channel=resample(channel,new_sampling_rate,sampling_rate)/factor;
mtlwch(outfile,n_channel,n);
end; 
