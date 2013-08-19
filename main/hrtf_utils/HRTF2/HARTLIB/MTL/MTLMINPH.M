function mtlminphase(infile, outfile)

% mtlminphase(infile, outfile)
%
%
% MTLMINPHASE calculates the minimum phase impulse response
% The delay of the envelope is preserved 
%    
% infile: name of input file with long impulse response
% outfile: file with cut impulse response
% Author: Klaus Hartung [Lehrstuhl AEA, Ruhr-Universitaet Bochum]
% Date: 28.7.1995


[file_type,info_blocks,n_channels,n_lines, sampling_rate,...
 first_line,last_line,n_directions,comment1, comment2] = mtlrh(infile);


if (file_type ~= 2) error('MTLCUT: Wrong filetype'); return; end;
if ((n_channels/n_directions) == 2) binaural=1; end;

last_line_out=n_lines/sampling_rate;


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


mtlwh(outfile,file_type,info_blocks,n_channels,n_lines, sampling_rate,...
       first_line,last_line,n_directions,...
       comment1, comment2);

if ((n_blocks) > 0) status=mtlwdir(outfile,direction_matrix); end;




%mtlwdel2(outfile, t_start);
if(length(delay)>0) mtlwdel2(outfile, delay); end;


for n=1:n_channels
	clc
	disp(['MTLMINPHASE channel ', num2str(n), ' of ', num2str(n_channels)] );
 
	channel=mtlrch(infile,n);
        h_min=minphase(channel,1,1);
 	mtlwch(outfile,h_min,n);
end; 


return;

end;
