function mtlcut(infile, outfile, lwindow, n_lines_out)

% mtlcut(infile, outfile, lwindow, n_lines_out)
%
%
% MTLCUT cuts the most important part of an impulse response.
% The preceding delay is evaluated and stored in the header of the new file.
%    
% infile: name of input file with long impulse response
% outfile: file with cut impulse response
% l_window: length of the window in samples
% n_lines_out: number of lines in new impulse response ( >= lwindow)

% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 21.4.1995


[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,comment1, comment2] = mtlrh(infile);


if (file_type ~= 2) error('MTLCUT: Wrong filetype'); return; end;
if ((n_channels/n_directions) == 2) binaural=1; end;
if (nargin < 4) n_lines_out=lwindow; end;
last_line_out=n_lines_out/sampling_rate;
if (lwindow > n_lines_out) error('MTLCUT: lwindow > n_lines_out'); return; end;

%if direction_matrix
if ((info_blocks-2) > 0)
   direction_matrix=mtlrdir(infile);

   % if time_delay exists read time_delay
   if (( (info_blocks-2)*256 - n_directions*8) == (n_channels*4))
       delay=mtlrdel2(infile);
   end;

end;

n_dir=n_directions;
if (n_dir > 0)
  n_byte=n_dir*2*4+n_channels*4;
  n_blocks=ceil(n_byte/256);
end;

info_blocks=2+n_blocks;


mtlwh(outfile,file_type,info_blocks,n_channels,n_lines_out, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);

if ((n_blocks) > 0) status=mtlwdir(outfile,direction_matrix); end;


nl=lwindow;

% get minimum delay time

for n=1:n_channels

      clc
      disp(['MTLCUT: First Run: channel  ', num2str(n), ' of ', num2str(n_channels)] );
      channel=mtlrch(infile,n);


% added Feb 9, 1996 CHK:
if (any(isnan(channel))) 
  error('MTLCUT: BAD FILE: NAN');
  disp(num2str(n));
  return;
end;
      [dummy,t(n)]=h_max_en(channel,nl);

end; 


t_min=min(t);
t_start=t-t_min+1;
mtlwdel2(outfile, t_start);

for n=1:n_channels
      %clc
      disp(['MTLCUT: Second Run: channel ', num2str(n), ' of ', num2str(n_channels)] );
      channel=mtlrch(infile,n);
      n_channel=zeros(1,n_lines_out);
      n_channel(t_start(n):t_start(n)+nl)=channel(t(n):t(n)+nl);
      mtlwch(outfile,n_channel,n);
end; 

return;

%end;