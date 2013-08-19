function mtlscale(infile, outfile, factor)


% function mtlcut(infile, outfile, factor)

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,comment1, comment2] = mtlrh(infile)


if (file_type ~= 2) error('Wrong filetype'); return; end;



if (nargin < 4) n_lines_out=duration*sampling_rate; end;
last_line_out=n_lines_out/sampling_rate;

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


 
nl=duration;
%nl=duration*sampling_rate;

% get minimum delay time
for n=1:n_channels
%for n=1:4

clc
disp(['MTLCUT: First Run: channel  ', num2str(n), ' of ', num2str(n_channels)] );
 channel=mtlrch(infile,n);
 [dummy,t(n)]=h_max_en(channel,nl);

end; 



t_min=min(t)
t_start=t-t_min+1
[r,c]=size(t_start);
delay_l=t_start(1:2:c);
delay_r=t_start(2:2:c);
%pause
% write output_header

mtlwh(outfile,file_type,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line_out,n_directions,...
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


mtlwh(outfile,file_type,info_blocks,n_channels,n_lines_out, sampling_rate,...
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


%pause

for n=1:n_channels
%for n=1:4
close all
clc
disp(['MTLCUT: Second Run: channel ', num2str(n), ' of ', num2str(n_channels)] );
 
channel=mtlrch(infile,n);
%irplot(channel,sampling_rate);

n_channel=zeros(1,n_lines_out);
t(n)
 n_channel(t_start(n):t_start(n)+nl)=channel(t(n):t(n)+nl);
% irplot(n_channel',sampling_rate);
 mtlwch(outfile,n_channel,n);
%pause
end; 

