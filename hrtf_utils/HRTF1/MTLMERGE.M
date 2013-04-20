function mtlmerge(file_list, outputfile)

% function mtlmerge(file_list, outputfile)

[n_files,str_len] = size(file_list)

[file_type_ref,info_blocks,n_channels,n_lines_ref, sampling_rate_ref,...
 first_line_ref,last_line_ref,n_directions,comment1, comment2] = mtlrh(file_list(1,:));

% set binaural
binaural_ref=0;
if ((n_channels/n_directions) == 2) binaural_ref=1; end;

total_channels=0;
total_directions=0;

for n=1:n_files
  clc
  disp(['MTLMERGE: Scanning file ', num2str(n), ' of ', num2str(n_files)] );

% read header (file_list(n,:=))

  [file_type,info_blocks,n_channels,n_lines, sampling_rate,...
   first_line,last_line,n_directions,comment1, comment2] = mtlrh(file_list(n,:));

% compare header
% if different error-message and end

  if (file_type ~= file_type_ref) error('different file_type'); return; end;
  if (n_lines ~= n_lines_ref) error('different number_of_lines'); return; end;

% this is not corrected until using mtlcorhe.m CHK 2/96:
% if (sampling_rate ~= sampling_rate_ref) ; error('different sampling_rate'); disp(sampling_rate);disp(sampling_rate_ref); pause; return; end;
  if (first_line ~= first_line_ref) error('different first_line'); return; end;

% this is not corrected until using mtlcorhe.m CHK 2/96:
% if (last_line ~= last_line_ref) ;disp(last_line);disp(last_line_ref);pause; return; end;

% compare binaural
  binaural=0;
  if ((n_channels/n_directions) == 2) binaural=1; end;
  if (binaural ~= binaural_ref) error('can t mix mono and stereo files'); return; end;

  total_channels=total_channels+n_channels;
  total_directions=total_directions+n_directions;

% if direction_matrix
  if ((info_blocks-2) > 0)
     direction_matrix=mtlrdir(file_list(n,:),n_directions);

  % if (n ~= 0) direction_matrix = direction_matrix(:,1); end;
   % if time_delay exists read time_delay
     if (( (info_blocks-2)*256 - n_directions*8) == (n_channels*4))
        [delay_l,delay_r]=mtlrdel(filename,n_directions,binaural);
        if (n ~= 0) 
          delay_l = delay_l(index); 
          if (binaural == 1) delay_r = delay_r(index); end;
        end;
     end;

% set_new_direction_matrix
  new_direction_matrix(1:2,n)=direction_matrix;
%moved end here from two lines up CHK 2/96
end;

% set new_delay
% new_delay_l(n) = delay_l;
% if (binaural==1)new_delay_r(n) = delay_r;end;

% added Feb 9, 1996 CHK:
if (any(isnan(mtlrch(file_list(n,:),1)))) 
  disp(file_list(n,:));
  error('MTLMERGE: BAD FILE: NAN');
  return;
end;

if (binaural == 1)
if (any(isnan(mtlrch(file_list(n,:),2)))) 
  disp(file_list(n,:));
  error('MTLMERGE: BAD FILE: NAN');
  return;
end;
end;

end;

% calculate header param

sn_dir=size(new_direction_matrix);

n_dir=sn_dir(2);
if (n_dir > 0)
  n_byte=n_dir*2*4;
  sn_del=length(delay_l);
  if (sn_del > 0) n_byte=n_byte+n_dir*4*(binaural+1); end; 
  n_blocks=ceil(n_byte/256);
  new_info_blocks=2+n_blocks;
else				%new 2/96 CHK
  new_info_blocks=2;
end;

% write header(outputfile)
disp('MTLMERGE WRITING HEADER OF OUTPUTFILE');

mtlwh(outputfile,file_type_ref,new_info_blocks,total_channels,n_lines_ref, sampling_rate_ref,...
       first_line_ref,last_line_ref,total_directions,...
       comment1, comment2);


% if direction, write dir(outputfile)
if ((n_blocks) > 0)
 status=mtlwdir(outputfile,new_direction_matrix);
% if type 2, write delay(outputfile)
% if(length(delay_l) > 0) mtlwdel(filename, delay_l,delay_r); end;

% new_direction_matrix is correct at this point

end;

out_chan_n=1;

for n=1:n_files
[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 start,stop,n_directions,comment1, comment2] = mtlrh(file_list(n,:));
%faster if you only read n_channels

clc
disp(['MTLMERGE: Appending file ', num2str(n), ' of ', num2str(n_files)] );

   for c=1:n_channels
       channel=mtlrch(file_list(n,:),c);
       mtlwch(outputfile,channel,out_chan_n);
       out_chan_n = out_chan_n+1;
   end;
end;


return;