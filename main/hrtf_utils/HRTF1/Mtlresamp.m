function mtlresamp(infile, outfile, factor)

%function mtlresamp(infile, outfile, factor)
% factor is a two element vector: [P Q] as in resample.m
% resamples at P/Q the original rate

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,comment1, comment2] = mtlrh(infile)

if (file_type ~= 2) error('Wrong filetype'); return; end;

if (nargin < 3) factor = [2 1]; end;
if length(factor) ~= 2 error('Factor incorrect.  Requires [P Q]');  end;
index=0;

if n_channels == n_directions
   binaural = 1;
else
   binaural =2;
end

%if direction_matrix
if ((info_blocks-2) > 0)
   direction_matrix=mtlrdir(infile,n_directions);
   if (index ~= 0) 
      direction_matrix(:) = direction_matrix(:,index);
   end;
   % if time_delay exists read time_delay
   %if (( (info_blocks-2)*256 - n_directions*8) == (n_channels*4))
   %   [delay_l,delay_r]=mtlrdel(infile,n_directions,binaural);
   %  if (index ~= 0) 
   %     delay_l = delay_l(index); 
   %     if (binaural == 1) delay_r = delay_r(index); end;
   %  end;
   end;
end;

% write output_header
sampling_rate = sampling_rate * factor(1)/factor(2);
comment2 = 'after resampling';

mtlwh(outfile,file_type,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);

% write directions;
% write delays;
if ((n_channels/n_directions) == 2) binaural=1; end;
n_dir=n_directions;
if (n_dir > 0)
  n_byte=n_dir*2*4;
  %if (exist(delay_l) ~= 0) n_byte=n_byte+n_dir*4*(binaural+1); end; 
  n_blocks=ceil(n_byte/256);
end;

info_blocks=2+n_blocks;

channel=mtlrch(infile,1);
temp = resample(channel,factor(1),factor(2));
n_lines = length(temp);
comment1 = 'after resampling';

mtlwh(outfile,file_type,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);


[fid,message] = fopen(outfile,'a','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;

if ((n_blocks) > 0)
 fseek(fid,512,'bof');
 fwrite(fid,zeros(n_blocks*256,1),'float32');
 fclose(fid);
end

%if direction_matrix
if ((n_blocks) > 0)
 status=mtlwdir(outfile,direction_matrix)
 %if(exist(delay_l) ~= 0) mtlwdel(outfile, delay_l,delay_r); end;
end

for n=1:n_channels
  clc
  disp(['MTLRESAMP: channel ', num2str(n), ' of ', num2str(n_channels)] );
  channel=mtlrch(infile,n);
  mtlwch(outfile,resample(channel,factor(1),factor(2)),n);
end; 
