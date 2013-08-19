function mtlinv(infile,outfile,f1, f2 ,f3 ,f4)

% function mtlinv(infile,outfile,f1, f2 ,f3 ,f4)

if(nargin < 2) error("Args are missing") return; end;

[filetype,info_blocks,...
 n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,...
 comment1, comment2] = mtlrh(infile)

if (filetype ~= 3) 
  disp(['File: ' infile ' has wrong type']);
  return;
end;


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
if (n_dir < 0)
  n_byte=n_dir*2*4;
  sn_del=length(delay_l);
  if (sn_del > 0) n_byte=n_byte+n_dir*4*(binaural+1); end; 
  n_blocks=ceil(n_byte/256);
  info_blocks=2+n_blocks;

end;


mtlwh(outfile,filetype,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);
%if direction_matrix
if ((n_blocks) > 0)
 status=mtlwdir(outfile,direction_matrix)
 if(length(delay_l) > 0) mtlwdel(outfile, delay_l,delay_r); end;
end

if_nargin==6
resolution=sampling_rate/2/n_lines;
l1=(f1-first_line)/resolution+1;
l2=(f2-first_line)/resolution+1;
l3=(f3-first_line)/resolution+1;
l4=(f4-first_line)/resolution+1;
end

for index=1:n_channels
channel = mtlrch(infile,index);
  channel = 1./channel;
  if (nargin==6) channel=fade(channel,l1,l2,l3,l4);  end;  
 status = mtlwch(outfile,channel, index)
end;
return;
                                                                                                                                                                                                                                                          
