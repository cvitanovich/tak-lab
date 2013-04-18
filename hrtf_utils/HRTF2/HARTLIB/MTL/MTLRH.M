function [comment1, filetype, info_blocks,...
          n_channels, day_of_date, month_of_date, year_of_date,...
          no_bl_per_chan, subject_no, ...
          subject_first_name, subject_last_name, ...
          n_lines, sampling_rate,...
          lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
          no_of_qu_bits, chan_code, ampl_gain, ...
          max_ele, min_ele, status, res_type, res, ...
          channels_per_dir, map, comment2] = mtlrh(filename, dispinfo)

% MTLRH  reads header of bifi-, sigmaps- and mtl-files.
%
%          [comment1, filetype, info_blocks,...
%          n_channels, day_of_date, month_of_date, year_of_date,...
%          no_bl_per_chan, subject_no, ...
%          subject_first_name, subject_last_name, ...
%          n_lines, sampling_rate,...
%          lo_cutoff, hi_cutoff, eff_s_to_n, no_of_aver, ... 
%          no_of_qu_bits, chan_code, ampl_gain, ...
%          max_ele, min_ele, status, res_type, res, ...
%          channels_per_dir, map, comment2] = MTLRH(filename,dispinfo)
%
% filetype: 
%           1 -> int, time-domain
%           2 -> float time-domain
%           3 -> complex frq-domain
%           4 -> int, AD-DA type
% info_blocks: number of info blocks
% n_channels: number of channels in file
% n_lines: number of lines per channel
% sampling_rate: sampling rate in kHz
% comment1: string with  79 characters
% comment2: string with 262 characters
% dispinfo: display information 0:no (default) 1:yes
% 
% See also: MTLWH, MTLRCH, MTLWCH

%
% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%
% Date: 19.4.1995
% Modified: 19.01.1996 Klaus Hartung 
%           Using bifi and sigmap structure now  
% Modified: 11.5.1996 Klaus Hartug
%	    dispinfo implemented
% Modified: 7.7.1996 Klaus Hartung
%	    help info improved

% Calls		toolbox
%
% no calls to other routines


if (nargin < 1) error('filename is missing'); end;
if (nargin < 2) dispinfo=0; end;
if (nargout < 1) dispfinfo=1; end;

comment1 = 'this is a test';
eof = -1;
filetype = -1;
info_blocks = -1;
n_channels = -1;
day_of_date = -1;
month_of_date = -1;
year_of_date = -1;
no_bl_per_chan = -1;
subject_no = -1;
subject_first_name = 'NON';
subject_last_name = 'NOMINATUR';
n_lines = -1;
sampling_rate = -1;
lo_cutoff = -1;
hi_cutoff = -1;
eff_s_to_n = 9999;
no_of_aver = -1;
no_of_qu_bits = -1;
chan_code = -1;
ampl_gain = ones(16,1) * (-1);
magic_nr = -1;
max_ele = 9999;
min_ele = 9999;
status = -1;
res_type = -1;
res = -1;
channels_per_dir = -1;
map = [-999; -999];
comment2 ='comment2';


SIGM_BLOCK_LEN = 256;
DOS_EOF = 26;

SO_EOF = 79;
SO_TYPE = 80;
SO_N = 128;
SO_CHA = 84;

INT_TYPE = 1;
FLOAT_TYPE = 2;
COMPLEX_TYPE = 3;
ADDA_TYPE = 4;


BIFI_MAGIC_NUMBER = hex2dec('7caffee1');

% status mask

BIFI_MEASURED     = 1;
BIFI_X_TESTSIGNAL = 2;
BIFI_SPEAKER_EQ   = 4;
BIFI_TOTAL_EQ     = 8;

BIFI_RES_NUMBER = 5;

BIFI_RES_AS_MEASURED         = 1;
BIFI_RES_CONSTANT_AZIMUTH    = 2;
BIFI_RES_CONSTANT_DENSITY    = 3;
BIFI_RES_OLD_FULL_SPHERE     = 4;
BIFI_RES_WITH_DIRECTION_INFO = 5;

%open file for reading

[fid,message] = fopen(filename,'r','ieee-le');
if (fid==-1)         % file can't be opened
  error(message);
  filetype=-1;
  return;
end;


% read comment1

stat=fseek(fid,0,'bof');
if (stat<0) error('ERROR: mtlrh, fseek comment1'); filetype=-1;return, end;
comment1 = setstr(fread(fid,79,'uchar')) ;
comment1=comment1';




% read eof

%test if data are stored in PC-format (EOF=26) of Workstation format(EOF=XX)

fseek(fid,79,'bof');
eof = fread(fid,1,'char');
eof_int = abs(eof);
if (eof_int ~= 26);
  fclose(fid);
  error('choosen file isn''t readable');
  filetype=-1;
  return;
end;


%check filetype

fseek(fid,80,'bof');
filetype = fread(fid,1,'short');
if ( (filetype<1)|(filetype>4) )
 fclose(fid);
 error('choosen file includes unknown filetype');
end;



% read number of info_blocks

stat=fseek(fid,82,'bof');
if (stat<0) error('ERROR: mtlrh, fseek info_blocks'); filetype=-1;return, end;
info_blocks = fread(fid,1,'ushort');

% read number of channels

stat=fseek(fid,84,'bof');
if (stat<0) error('ERROR: mtlrh, fseek n_channels'); filetype=-1;return, end;
n_channels = fread(fid,1,'ushort');


% read day_of_date

stat=fseek(fid,86,'bof');
if (stat<0) error('ERROR: mtlrh, fseek day_of_date'); filetype=-1;return, end;
day_of_date = fread(fid,1,'ushort');


% read month_of_date

stat=fseek(fid,88,'bof');
if (stat<0) error('ERROR: mtlrh, fseek month_of_date'); filetype=-1;return, end;
month_of_date = fread(fid,1,'ushort');

% read year_of_date

stat=fseek(fid,90,'bof');
if (stat<0) error('ERROR: mtlrh, fseek year_of_date'); filetype=-1;return, end;
year_of_date = fread(fid,1,'ushort');



% read number of blocks per channel

stat=fseek(fid,92,'bof');
if (stat<0) error('ERROR: mtlrh, fseek number of blocks per channel'); filetype=-1;return, end;
no_bl_per_chan = fread(fid,1,'ushort');


% the following data are optional


% read number of subject/animal

stat=fseek(fid,96,'bof');
if (stat<0) error('ERROR: mtlrh, fseek subject/animal number'); filetype=-1;return, end;
subject_no = fread(fid,1,'long');


% read name of subject/animal

stat=fseek(fid,100,'bof');
if (stat<0) error('ERROR: mtlrh, fseek subject first name'); filetype=-1;return, end;
subject_first_name = setstr(fread(fid,10,'uchar'))';

stat=fseek(fid,110,'bof');
if (stat<0) error('ERROR: mtlrh, fseek subject last_name name'); filetype=-1;return, end;
subject_last_name = setstr(fread(fid,10,'uchar'))';




% read number of points per channel
stat=fseek(fid,128,'bof');
if (stat<0) error('ERROR: mtlrh, fseek n_lines'); filetype=-1;return, end;
n_lines = fread(fid,1,'ulong');



% read sampling rate
stat=fseek(fid,132,'bof');
if (stat<0) error('ERROR: mtlrh, fseek sampling_rate'); filetype=-1;return, end;
sampling_rate = fread(fid,1,'float32');


% read lower cutoff frequency

stat=fseek(fid,136,'bof');
if (stat<0) error('ERROR: mtlrh, fseek lower cutoff_frequency'); filetype=-1;return, end;
lo_cutoff = fread(fid,1,'float32');


% read lower cutoff frequency

stat=fseek(fid,140,'bof');
if (stat<0) error('ERROR: mtlrh, fseek higher cutoff_frequency'); filtetype=-1; return, end;
hi_cutoff = fread(fid,1,'float32');


% read effective S/N ratio

stat=fseek(fid,144,'bof');
if (stat<0) error('ERROR: mtlrh, fseek effective S/N ratio'); filetype=-1;return, end;
eff_s_to_n = fread(fid,1,'float32');


% read number of averaging

stat=fseek(fid,148,'bof');
if (stat<0) error('ERROR: mtlrh, fseek number of averagings'); filetype=-1;return, end;
no_of_aver = fread(fid,1,'ulong');


% read number of quantization bits

stat=fseek(fid,152,'bof');
if (stat<0) error('ERROR: mtlrh, fseek number of quantization bits'); filetype=-1;return, end;
no_of_qu_bits = fread(fid,1,'ushort');


% read number of channel code y/n, give 0 or 1

stat=fseek(fid,154,'bof');
if (stat<0) error('ERROR: mtlrh, fseek channel code'); filetype=-1;return, end;
chan_code = fread(fid,1,'ushort');

% read amplification by input VCA

stat=fseek(fid,156,'bof');
if (stat<0) error('ERROR: mtlrh, fseek amplification by input VCA'); filetype=-1;return, end;
ampl_gain = fread(fid,16,'float32');



% now BIFI information starts

stat=fseek(fid,220,'bof');
if (stat<0) error('ERROR: mtlrh, fseek magic number'); filetype=-1;return, end;
magic_nr = fread(fid,1,'ulong');

if (magic_nr == BIFI_MAGIC_NUMBER)

% read max_ele
stat=fseek(fid,224,'bof');
if (stat<0) error('ERROR: mtlrh, fseek max ele'); filetype=-1;return, end;
max_ele = fread(fid,1,'float');

% read min_ele
stat=fseek(fid,228,'bof');
if (stat<0) error('ERROR: mtlrh, fseek min ele'); filetype=-1;return, end;
min_ele = fread(fid,1,'float');

% read bifi stat
stat=fseek(fid,232,'bof');
if (stat<0) error('ERROR: mtlrh, fseek status'); filetype=-1;return, end;
status = fread(fid,1,'ushort');


% read bifi res_type
stat=fseek(fid,234,'bof');
if (stat<0) error('ERROR: mtlrh, fseek res type'); filetype=-1;return, end;
res_type = fread(fid,1,'ushort');

% read bifi resolution
stat=fseek(fid,236,'bof');
if (stat<0) error('ERROR: mtlrh, fseek resolution'); filetype=-1;return, end;
res = fread(fid,1,'ushort');

% read channels per direction
stat=fseek(fid,238,'bof');
if (stat<0) error('ERROR: mtlrh, fseek channels per direction'); filetype=-1;return, end;
channels_per_dir = fread(fid,1,'ushort');


if ((n_channels > 0) & (channels_per_dir > 0)) 
   n_dir = n_channels/channels_per_dir;
else 
  n_dir = 1;
end;

if (res_type == BIFI_RES_WITH_DIRECTION_INFO)

	if (n_dir == 1)
   		map = zeros(2,1);
  		stat=fseek(fid,240,'bof');
   		if (stat<0) error('ERROR: mtlrh, fseek map for single direction'); filetype=-1;return, end;
   		map = fread(fid,2,'float');
   		map = reshape(map,2,1);
	end

	if (n_dir ~= 1) & (info_blocks > 2)
		map =  zeros(2,n_dir);
   		stat=fseek(fid,512,'bof');
   		if (stat<0) error('ERROR: mtlrh, fseek map for multiple directions'); filetype=-1;return, end;
   		map = fread(fid,n_dir*2,'float');
   		map = reshape(map,2,n_dir);
	end 


end % of BIFI_RES_TYPE

% end  of BIFI_MAGIC_NUMBER

else 
	stat = fseek(fid,220,'bof');
	if (stat<0) error('ERROR: mtlrh, fseek n_dir for non BIFI-file'); filetype=-1; return; end
	n_dir = fread(fid,1,'short');
	channels_per_dir=n_channels/n_dir;

	if (info_blocks>2)
		map=zeros(2,n_dir);
		stat=fseek(fid,512,'bof');
		if(stat<0) error('Error: mtlrh, fseek map for uuo Option'); filetype=-1; return; end
		map = fread(fid,n_dir*2,'float');
		map = reshape(map,2,n_dir);
	end
end


stat=fseek(fid,248,'bof');
if (stat<0) error('ERROR: mtlrh, fseek comment2'); filetype=-1;return, end;
comment2 = setstr(fread(fid,262,'uchar'));
comment2=comment2'; 


stat = fclose(fid);

if (filetype == 1) iscomplex=0; type = 'short';   n_byte = 2; end;
if (filetype == 2) iscomplex=0; type = 'float32'; n_byte = 4;end;
if (filetype == 3) iscomplex=1; type = 'float32'; n_byte = 4;end;
if (filetype == 4) iscomplex=0; type = 'short';   n_byte = 2; end;

filesize=sizeoffile(filename);
n_lines_v = (filesize-256*info_blocks)/(n_channels*4);


if ((n_lines ~= n_lines_v) & (n_lines_v ~=0) & magic_nr~=BIFI_MAGIC_NUMBER)
	disp ('mtlrch: something wrong with header or file length');
	disp (['Value in header: ' int2str(n_lines) ' Number of lines derived from file length: ' num2str(n_lines_v)]);
	disp ('Ignore value in header ');
	n_lines=n_lines_v;
end;


if (dispinfo==1)
disp(['Filename: ' filename]);
disp(['Comment 1: ' comment1]);
disp(['Filetype: ' int2str(filetype)]);
disp(['Number of Info Blocks: ' int2str(info_blocks)]);
disp(['Number of channels: ' int2str(n_channels)]);
disp(['Date: ' int2str(day_of_date) '.' int2str(month_of_date) '.' int2str(year_of_date)]);
disp(['Number of blocks per channel: ' int2str(no_bl_per_chan)]);
disp(['Subject number: ' int2str(subject_no) ' Subject first name: ' subject_first_name ' Subject last name: ' subject_last_name]);
disp(['Number of lines per channel: ' int2str(n_lines)]);
disp(['Sampling rate: ' int2str(sampling_rate) ' Hz']);
disp(['Cutoff Frequencies: Lo: ' int2str(lo_cutoff) ' Hi: ' int2str(hi_cutoff) ' Hz']);
disp(['effective SNR: ' int2str(eff_s_to_n) ' dB']);
disp(['Number of averages: ' int2str(no_of_aver)]);
disp(['Number of quantization bits: ' int2str(no_of_qu_bits)]);
disp(['Channel code: ' int2str(chan_code)]);
disp(['Amplifier gain: ' num2str(ampl_gain) ' dB']);
disp(['Maximum elevation: ' int2str(max_ele) ' Minimum elevation: ' int2str(min_ele)]);
disp(['Status: ' int2str(status)]);
disp(['Resolution Type: ' int2str(res_type) ]);
disp(['Resolution: ' int2str(res)]);
disp(['Channels per direction: ' int2str(channels_per_dir)]);
disp(['Comment 2: ' comment2]);
if (exist('map') & (n_dir==1) & ( (res_type == BIFI_RES_WITH_DIRECTION_INFO) | (info_blocks>2) ) ) disp(['Map: Ele ' num2str(map(1,1)) ', Azi: ' num2str(map(2,1))]); end;
end;

if (nargout <1) clear all; end; 

return;

% end of mtlrh.m


