function mtlmean(infile,outfile)

% function mtlmean(infile, outfile)

% read head infile
[file_type1,info_blocks1,n_channels1,n_lines1, sampling_rate1,...
 first_line1,last_line1,n_directions1,comment1_1, comment2_1] = mtlrh(infile);


% compare header
% if they do not fit error

index=0;
%if direction_matrix
if ((info_blocks1-2) > 0)
   direction_matrix=mtlrdir(infile,n_directions1);
   if (index ~= 0) direction_matrix(:) = direction_matrix(:,index); end;
   % if time_delay exists read time_delay
   if (( (info_blocks1-2)*256 - n_directions1*8) == (n_channels1*4))
      [delay_l,delay_r]=mtlrdel(infile,n_directions1,binaural);
     if (index ~= 0) 
        delay_l = delay_l(index); 
        if (binaural == 1) delay_r = delay_r(index); end;
     end;
   end;
end;


% calculate output_header

filetype=file_type1;
info_blocks=info_blocks1;
n_channels=n_channels1/n_directions1;
n_lines=n_lines1;
sampling_rate=sampling_rate1;
first_line=first_line1;
last_line=last_line1;
n_directions=1;
comment1=comment1_1;
comment2=comment2_1;
filename=outfile;

 
% write output_header

mtlwh(filename,filetype,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);


% write directions;
% write delays;
if ((n_channels/n_directions) == 2) binaural=1; end;

sn_dir=size(direction_matrix);

n_dir=sn_dir(2);
if (n_dir > 0)
  n_byte=n_dir*2*4;
  sn_del=length(delay_l);
  if (sn_del > 0) n_byte=n_byte+n_dir*4*(binaural+1); end; 
  n_blocks=ceil(n_byte/256);
end;

info_blocks=2+n_blocks;


mtlwh(filename,filetype,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);


%if direction_matrix
[fid,message] = fopen(filename,'a','ieee-le');
if (fid==-1)          % file can't be opened
  disp(message);
  return;
end;

if ((n_blocks) > 0)
 fseek(fid,512,'bof');
 fwrite(fid,zeros(n_blocks*256,1),'float32');
 fclose(fid);
 
status=mtlwdir(filename,direction_matrix);
 if(length(delay_l) > 0) mtlwdel(filename, delay_l,delay_r); end;
end


m_ch1=zeros(n_lines,1);
m_ch2=zeros(n_lines,1);

for n=1:n_directions1
  clc
  disp(['MTLMEAN: Calculating direction ', num2str(n), ' of ', num2str(n_directions1)] );

 channel1=mtlrch(infile,(n-1)*(binaural+1)+1);
m_ch1=m_ch1+channel1;

 if (binaural==1)  
    channel2=mtlrch(infile,n*(1+binaural));
    m_ch2=m_ch2+channel2;
  end
end

m_ch1=m_ch1/n_directions1;
mtlwch(outfile,m_ch1,1);
if (binaural==1)
  m_ch2=m_ch2/n_directions1;
  mtlwch(outfile,m_ch2,2);
end;

return;
