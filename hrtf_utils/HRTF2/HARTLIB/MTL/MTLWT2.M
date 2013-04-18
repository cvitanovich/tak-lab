function status =mtlwt2(filename, file_format, file_type, h_l, h_r,...
                           direction_matrix, delay_r, delay_l, n_directions,...
                           n_channels, n_lines, start, stop);

% MTLWT2 reads file of type 2
% status = mtlwt2(filename, file_format, file_type, h_l,h_r,...
%                           direction_matrix, delay_r, delay_l, n_directions,...
%                           n_channels, n_lines, start, stop);
%
%          filename: string with filename
%                 
%


binaural=0;
n_blocks=0;

 
   fileformat='';
   filtetype=3;
   [n_directions,n_lines]=size(h_l); 
   if ([n_directions,n_lines]==size(hl))  binaural=1; end;

if ((n_channels/n_directions) == 2) binaural=1; end;

sn_dir=size(direction_matrix);

n_dir=sn_dir(2);
if (n_dir < 0)
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
if ((n_blocks) > 0)
 status=mtlwdir(filename,direction_matrix)
 if(length(delay_l) > 0) mtlwdel(filename, delay_l,delay_r); end;
end



[fid,msg]=fopen(filename,'a','ieee-le');
if fid<0 error(msg); end;


if (index >n_directions) error('Index not in catalogue'); end;

if (index == 0)
  first=1;
  last=n_directions;
 else
 first=index;
 last=index;
end;


for  n= first:last
  fseek(fid,info_blocks*256+n_lines*4*(binaural+1)*(n-1),'bof');
 
  h=fwrite(fid,[h_l(n,1:n_lines);h_r(n,1:n_lines)],'float32');

end;


fclose(fid);
return;
