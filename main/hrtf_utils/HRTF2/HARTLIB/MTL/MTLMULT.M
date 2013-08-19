function mtlmult(file1, file2, outfile)


% function mtlmult(file1, file2, outfile)

% read head file1
[file_type1,info_blocks1,n_channels1,n_lines1, sampling_rate1,...
 first_line1,last_line1,n_directions1,comment1_1, comment2_1] = mtlrh(file1);

% read head file2
[file_type2,info_blocks2,n_channels2,n_lines2, sampling_rate2,...
 first_line2,last_line2,n_directions2,comment1_2, comment2_2] = mtlrh(file2);

% compare header
% if they do not fit error

if (file_type1 ~= file_type2) error('Filetypes are different'); return; end;
if (n_lines1 ~= n_lines2) error('N_LINES are different'); return; end;
if (sampling_rate1 ~= sampling_rate2) error('Sampling_rates are different'); return; end;
if (first_line1 ~= first_line2) error('First_line is different'); return; end;
if (last_line1 ~= last_line2) error('First_line is different'); return; end;

filename=file1
n_directions=n_directions1
filetype=file_type1
index=0
info_blocks=info_blocks1
n_channels=n_channels1
comment2=comment2_1
%pause

%if direction_matrix
if ((info_blocks1-2) > 0)
   direction_matrix=mtlrdir(filename,n_directions);
   if (index ~= 0) direction_matrix(:) = direction_matrix(:,index); end;
   % if time_delay exists read time_delay
   if (( (info_blocks-2)*256 - n_directions*8) == (n_channels*4))
      [delay_l,delay_r]=mtlrdel(filename,n_directions,binaural);
     if (index ~= 0) 
        delay_l = delay_l(index); 
        if (binaural == 1) delay_r = delay_r(index); end;
     end;
   end;
end;


% calculate output_header

file_type=file_type1;
info_blocks=max(info_blocks1,info_blocks2)
n_channels=max(n_channels1,n_channels2);
n_c2=min(n_channels1,n_channels2);
n_lines=n_lines1;
sampling_rate=sampling_rate1;
first_line=first_line1;
last_line=last_line2;
n_directions=max(n_directions1,n_directions2);
comment1=comment1_1;
comment2=comment1_2;

% if necessary swap file1 file2
if (n_channels1 < n_channels2)
 string=file2;
 file2=file1;
 file1=string;
end;
 
sn_dir=size(direction_matrix);
n_dir=sn_dir(2);

n_blocks=0;
n_dir=sn_dir(2);
if (n_dir > 0)
  n_byte=n_dir*2*4;
  sn_del=length(delay_l);
  if (sn_del > 0) n_byte=n_byte+n_dir*4*(binaural+1); end; 
  n_blocks=ceil(n_byte/256);
end;

info_blocks=2+n_blocks;


% write output_header


mtlwh(outfile,filetype,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);

if ((n_channels/n_directions) == 2) binaural=1; end;






%if direction_matrix
[fid,message] = fopen(outfile,'a','ieee-le');
if (fid==-1)          % file can't be opened
  disp(message);
  return;
end;

if ((n_blocks) > 0)
 fseek(fid,512,'bof');
  fwrite(fid,zeros(n_blocks*256,1),'float32');
 fclose(fid);

 status=mtlwdir(outfile,direction_matrix)
 if(length(delay_l) > 0) mtlwdel(outfile, delay_l,delay_r); 
  end;
end


chan_num=1;
for n=1:n_channels
clc
  disp(['MTLDIV: Calculating channel ', num2str(n), ' of ', num2str(n_channels)] );

 channel1=mtlrch(file1,n);
remainder = rem(n,n_c2);
 if (remainder == 0) n_2 = n_c2; 
    else n_2=remainder;
 end; 
  

  channel2=mtlrch(file2,n_2);
  channel_out = channel1.*channel2;

  mtlwch(outfile,channel_out,n);

end;
return;

