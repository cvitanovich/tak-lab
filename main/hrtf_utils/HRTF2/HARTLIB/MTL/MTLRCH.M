function channel=mtlrch(filename, index, line, nyquist_option);

% MTLRCH reads one channel from a mtl-, sigmaps- or bifi-file (type 2 or 3)
%
%         channel=mtlrch(filename, index, line, nyquist_option, old_format)
%         channel=mtlrch(filename)      
%         channel=mtlrch(filename, index)
%         channel=mtlrch(filename, index, line);
%         channel=mtlrch(filename, index, line, nyquist_option);
%         filename: string with filename
%         index: number of channel to select (index > 0 )     
%         line: read only one line
%	  nyquist_option : 1 switches on nyquist moving from 0 to N/2
%
%         mtlrch works with real and complex data.
%
% See also: MTLWCH, MTLWH, MTLRH

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 

% Date: 19.4.1995
% Modified: 12.5.1996 Klaus Hartung
% Modified: 7.7.1996 Klaus Hartung
%           help improved

% if no index available read first channel

if (nargin < 1) error('Insufficient number of arguments. Filename is missing.'); end;
if (nargin < 2) index=1; end;
if (nargin < 3) line=0;  end;
if (nargin < 4) nyquist_option=1;  end;

iscomplex=0;


[comment1, filetype, info_blocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2] = mtlrh(filename);


n_byte=0;

if (filetype == 1) iscomplex=0; type = 'short';   n_byte = 2; end;
if (filetype == 2) iscomplex=0; type = 'float32'; n_byte = 4;end;
if (filetype == 3) iscomplex=1; type = 'float32'; n_byte = 4;end;
if (filetype == 4) iscomplex=0; type = 'short';   n_byte = 2; end;

filesize=sizeoffile(filename);
n_lines_v = (filesize-256*info_blocks)/(n_channels*4);


%if (n_lines ~= n_lines_v) error('mtlrch: something wrong with header or file length'); return; end;

n_lines=n_lines_v;

[fid,message] = fopen(filename,'r','ieee-le');
if fid==-1          % file can't be opened
  error(message);
  channel=-1;
  return;
end;

pos=info_blocks*256+n_lines*n_byte*(index-1);
 

r_lines=n_lines;

if (line > 0)
   if (iscomplex)
       pos= pos + (line*2-1)*n_byte;    
       r_lines=2; 
   end;
   if (iscomplex==0) 
       pos = pos + (line-1)*n_byte;
       r_lines=1; 
   end;
end;

% fseek to channel

status=fseek(fid,pos,'bof');
if (status<0) error('ERROR: mtlrch,  fseek to channel not succesful' );channel=-1; return; end;


%fread of channel
ichannel=fread(fid,r_lines,type);
%if (ichannel<0) error('ERROR: mtlrch, reading of channel not succesful' );channel=-1; end;


status = fclose(fid);
if (status < 0) error('ERROR: mtlrch, close file');channel=-1; return; end;

% creating a complex vector if necessary

if ((iscomplex==1) & (nyquist_option == 1))
   channel=zeros(r_lines/2+1,1)+ i*zeros(r_lines/2+1,1);
   channel(1:r_lines/2,1)=ichannel(1:2:r_lines)+i*ichannel(2:2:r_lines);
   channel(r_lines/2+1,1)=imag(channel(1));
   channel(1)=real(channel(1))+i*0;
end


if ((iscomplex==1) & (nyquist_option == 0))
   channel=zeros(r_lines/2,1) + i*zeros(r_lines/2,1);
   channel(1:r_lines/2,1)=ichannel(1:2:r_lines)+i*ichannel(2:2:r_lines);
end

if (iscomplex==0)
%   disp('not complex');
   channel=ichannel;
end;

channel=channel(:);


% end of file mtlrch.m 