function [file_format, file_type, H_l,H_r, direction_matrix,...
          delay_l, delay_r, n_directions, n_channels, n_lines,...
          first_frq, last_frq]=mtlrt3(filename, index);

% MTLRT3  reads file of type 3
%         [file_format, file_type, H_l,H_r, direction_matrix,...
%         delay_l, delay_r, n_directions, m_channels, n_lines,...
%         first_frq, last_frq]=mtlrt3(filename, index);
%
%         filename: string with filename
%         index: number or transfer function or direction to select        
%

if (nargin < 2) index=0; end;
binaural=0;

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 first_frq,last_frq,n_directions,comment1, comment2] = mtlrh(filename);

if ((n_channels/n_directions) == 2) binaural=1; end;

%if direction_matrix
if ((info_blocks-2) > 0)
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

if (index ~= 0)
 direction_matrix(:) = direction_matrix(:,index);
end;


[fid,msg]=fopen(filename,'r','ieee-le');
if fid<0 error(msg); end;


N=(binaural+1)*n_lines*2;
if (index >n_directions) error('Index not in catalogue'); end;

if (index == 0)
  first=1;
  last=n_directions;
 else
 first=index;
 last=index;
end;


% aloc memory
H_l=zeros(last-first+1,n_lines/2)+j*ones(last-first+1,n_lines/2);
if (binaural == 1)
   H_r=zeros(last-first+1,n_lines/2)+j*ones(last-first+1,n_lines/2);
end;

N=(binaural+1)*n_lines*2;

for  n= first:last
  disp('test')
  fseek(fid,info_blocks*256+n_lines*4*2*(binaural+1)*(n-1),'bof');
 
  H=fread(fid,N,'float32');
  H=H';

  Hl=H(1:n_lines*2); 
  if (binaural == 1) Hr=H(n_lines*2+1:4*n_lines); end;
  H_l(n,1:n_lines)=Hl(1:2:n_lines*2)+i*Hl(2:2:n_lines*2);
  if (binaural == 1) H_r(n,1:n_lines)=Hr(1:2:n_lines*2)+i*Hr(2:2:n_lines*2); end;

end;

fclose(fid);
return;
