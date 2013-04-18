function [file_format, file_type, H_i, direction_matrix,...
          frq_axis, n_directions, n_channels, n_lines,...
          first_frq, last_frq]=mtlrt11(filename, index);

% MTLRT11 reads file of type 11
%         [file_format, file_type, H_i, direction_matrix,...
%         frq_axis n_directions, m_channels, n_lines,...
%         first_frq, last_frq]=mtlrt3(filename, index);
%
%         filename: string with filename
%         index: number or transfer function or direction to select        
%

if (nargin < 2) index=0; end;
binaural=0;

[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 first_frq,last_frq,n_directions,comment1, comment2] = mtlrh(filename);

if (filetype ~=11) disp('Wrong Filetype'); return; end;


%if direction_matrix
if ((info_blocks-2) > 0)
   direction_matrix=mtlrdir(filename,n_directions);
   if (index ~= 0) direction_matrix(:) = direction_matrix(:,index); end;
   % if frq_axis exists read frq_axis
   if (( (info_blocks-2)*256 - n_directions*8) == (n_channels*4))
      frq_axis=mtlrfrq(filename,n_directions);
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
H_i=zeros(last-first+1,n_lines/2)+j*ones(last-first+1,n_lines/2);

N=(binaural+1)*n_lines*2;

for  n= first:last
  disp('test')
  fseek(fid,info_blocks*256+n_lines*4*2*(binaural+1)*(n-1),'bof');
 
  H=fread(fid,N,'float32');
  H=H';

  Hl=H(1:n_lines*2); 
  H_i(n,1:n_lines)=Hl(1:2:n_lines*2)+i*Hl(2:2:n_lines*2);
 end;

end;

fclose(fid);
return;
