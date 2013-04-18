function status = mtlwch(filename, channel, index, line)

% MTLWCH write one channel to a file (type 1,2,3 or 4)
%
% status = mtlwch(filename,channel, index, line)
%         filename: string with filename
%         channel: vector with data
%         index: number of channel (index > 0)       
%         line: number of line
% See also: MTLRCH, MTLRH, MTLWH


% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%
% Date: 19.4.1995
% Mod:   25.4.1995
% Mod:   19.1.96 
% Mod:	 12.5.1996
% Modified: 8.7.96 Klaus Hartung
%     improved help
%     nargin check

% Calls 	Toolbox
% 
% mtlrh		hartung/mtl
% sizeoffile    hartung/tools

if (nargin < 1) error('No filename defined.'); end;
if (nargin < 2) error('No channel number defined.'); end;

if (nargin < 3) index=1; end;

if (nargin < 4)  line = -1; end;

sz=size(channel);

if (sz(2) ~=length(channel)) channel=channel.'; end;

iscomplex=0;

[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(filename,0);

if (filetype == 1) iscomplex=0; type = 'short';   n_byte = 2; end;
if (filetype == 2) iscomplex=0; type = 'float32'; n_byte = 4; end;
if (filetype == 3) iscomplex=1; type = 'float32'; n_byte = 4; end;
if (filetype == 4) iscomplex=0; type = 'short';   n_byte = 2; end;

end_pos=sizeoffile(filename);

[fid,msg]=fopen(filename,'r+','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  status=-1;
  return;
end;

%status=fseek(fid,0,'eof');
%if (status<0) error('ERROR: mtlwch,  fseek to end not succesful' ); status=-1;return; end;
%end_pos=ftell(fid);
%if (end_pos < 0) error('ERROR: mtlwch,  ftell at end not succesful' ); status=-1;return; end;
%status=fseek(fid,0,'bof');

pos=info_blocks*256+n_lines*n_byte*(index-1);

if (line > 0)
   pos =  pos + (line-1)*n_byte; 
   n_lines=1+iscomplex;
end;

if (pos > end_pos) error('ERROR: pos > end_pos'); status=-1; return; end;


status=fseek(fid,pos,'bof');
if (status<0) error('ERROR: mtlwch,  fseek to channel position not succesful' ); status=-1;return; end;

if ((filetype == 1) | (filetype == 2))
 ochannel=channel';
 o_lines=n_lines;
end;


if (iscomplex == 1)
 ochannel(1,:)= real(channel(1:n_lines/2));
 ochannel(2,2:n_lines/2)= imag(channel(2:n_lines/2));
 ochannel(2,1)= real(channel(n_lines/2+1));
 ochannel = reshape(ochannel,1,n_lines);
 o_lines = n_lines;
end


%if ((filetype == 4)
%ochannel(1:2:2*n_lines)= channel(1:n_lines);
%ochannel(1:2:2*n_lines)= channel2(1:n_lines);
%o_lines=n_lines*2;
%end

count=fwrite(fid,ochannel(1:n_lines),type);

if (count<o_lines) error('ERROR: mtlwch, writing of channel not succesful'); status=-1; return; end;

%end;

status=fclose(fid);
if (status < 0) error('ERROR: mtlwch, close file');status=-1; return; end;

%return;

%end;


% end of file mtlwch.m