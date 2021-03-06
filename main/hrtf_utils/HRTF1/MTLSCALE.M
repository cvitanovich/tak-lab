function mtlscale(infile, outfile, factor)

%function mtlscale(infile, outfile, factor)

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,comment1, comment2] = mtlrh(infile)

if (file_type ~= 2) error('Wrong filetype'); return; end;

if (nargin < 3) factor=1; end;

index=0;

%if direction_matrix
if ((info_blocks-2) > 0)
   direction_matrix=mtlrdir(infile,n_directions);
   if (index ~= 0) direction_matrix(:) = direction_matrix(:,index); end;
   % if time_delay exists read time_delay
   if (( (info_blocks-2)*256 - n_directions*8) == (n_channels*4))
      [delay_l,delay_r]=mtlrdel(infile,n_directions,binaural);
     if (index ~= 0) 
        delay_l = delay_l(index); 
        if (binaural == 1) delay_r = delay_r(index); end;
     end;
   end;
end;

% write output_header

mtlwh(outfile,file_type,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);

% write directions;
% write delays;
if ((n_channels/n_directions) == 2) binaural=1; end;
n_dir=n_directions;
if (n_dir > 0)
  n_byte=n_dir*2*4;
  sn_del=length(delay_l);
  if (sn_del > 0) n_byte=n_byte+n_dir*4*(binaural+1); end; 
  n_blocks=ceil(n_byte/256);
end;

info_blocks=2+n_blocks;
comment1 = [comment1 '; after scaling'];

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
 if(length(delay_l) > 0) mtlwdel(outfile, delay_l,delay_r); end;
end

for n=1:n_channels
  clc
  disp(['MTLSCALE: channel ', num2str(n), ' of ', num2str(n_channels)] );
  channel=mtlrch(infile,n);
  status = mtlwch(outfile,channel*factor,n)
end; 
fclose(fid);