function mtlifft(infile,outfile)

% function mtlifft(infile,outfile)


if (nargin<2) 
outfile=strrep(infile,'.sc','.sr');
outfile=strrep(infile,'.SC','.sr');
end;

[filetype,info_blocks,...
 n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,...
 comment1, comment2] = mtlrh(infile);


if (filetype ~= 3) 
  disp(['File: ' infile ' has wrong type']);
  return;
end;

mtltype=1;
if (last_line==(sampling_rate/2)) mtltype=0; end;

binaural=(n_channels/n_directions)-1;
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



n_dir=n_directions;
if (n_dir > 0)
  n_byte=n_dir*2*4;

  sn_del=length(delay_l);
  if (sn_del > 0) n_byte=n_byte+n_dir*4*(binaural+1); end; 
  n_blocks=ceil(n_byte/256);
  info_blocks=2+n_blocks;
end;

if (mtltype==0) 
   n_lines=(n_lines-1)*2;
else n_lines=n_lines*2;
end;
end;


first_line=0;
filetype=2;

last_line=(n_lines-1)*1/sampling_rate;
comment2 = 'passed through mtlifft';
mtlwh(outfile,filetype,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);



%[fid,message] = fopen(outfile,'a','ieee-le');
if exist(outfile)==2		% file exists
   [fid,message] = fopen(outfile,'r+','l');
else				% file does not exist
   [fid,message] = fopen(outfile,'w+','l');
end;


if fid==-1          % file can't be opened
  disp(message);
  return;
end;

 if ((n_blocks) > 0)
     status=fseek(fid,512,'bof');
     if (status<0) error('mtlifft not succesful -1' ); end;

      fwrite(fid,zeros(n_blocks*256,1),'float32');
    fclose(fid);
 end
 

%if direction_matrix
if ((n_blocks) > 0)
 status=mtlwdir(outfile,direction_matrix);
 if(length(delay_l) > 0) mtlwdel(outfile, delay_l,delay_r); end;
end


for index=1:n_channels
clc
disp(['MTLIFFT: Calculating channel ', num2str(index), ' of ', num2str(n_channels)] );


channel = mtlrch(infile,index);
channel=channel(1:n_lines/2);


t_channel = rifft(channel);

status = mtlwch(outfile,t_channel,index);
end;
return;
