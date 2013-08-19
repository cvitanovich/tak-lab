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


[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(infile);
filetype


if (filetype ~= 2) error('MTLCUT: Wrong filetype'); return; end;
n_dir = n_channels/channels_per_dir;

if (nargin < 4) n_lines_out=lwindow; end;

if (lwindow > n_lines_out) error('MTLCUT: lwindow > n_lines_out'); return; end;

info_blocks = 2 + ceil((n_dir*2*4+n_channels*4)/256)


mtlwh(outfile,...
comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines_out, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2);






nl=lwindow;


% get minimum delay time


glob_max=0;

for n=1:n_channels

	clc
	disp(['MTLCUT: First Run: channel  ', num2str(n), ' of ', num2str(n_channels)] );
 	channel=mtlrch(infile,n);
        max_amp=max(max(abs(channel)));
	if (max_amp > glob_max) glob_max=max_amp; end;
	
[dummy,t(n)]=h_max_en(channel,nl);
end; 


t_min=min(t);
t_start=t-t_min+1
mtlwdel(outfile, t_start);


for n=1:n_channels
	clc
	disp(['MTLCUT: Second Run: channel ', num2str(n), ' of ', num2str(n_channels)] );
 
	channel=mtlrch(infile,n);

	n_channel=zeros(n_lines_out,1);
 	n_channel(t_start(n):t_start(n)+nl,1)=channel(t(n):t(n)+nl,1);
        n_channel=n_channel/glob_max; 
	mtlwch(outfile,n_channel,n);
end; 


return;

