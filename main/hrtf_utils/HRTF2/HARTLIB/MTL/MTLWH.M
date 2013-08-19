function  mtlwh(filename,...
          comment1, filetype, info_blocks,...
          n_channels, day_of_date, month_of_date, year_of_date,...
          no_bl_per_chan, subject_no, ...
          subject_first_name, subject_last_name, ...
          n_lines, sampling_rate,...
          lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
          no_of_qu_bits, chan_code, ampl_gain, ...
          max_ele, min_ele, status, res_type, res, ...
          channels_per_dir, map, comment2)

% MTLWH  write header of SIGMAP-, BIFI- or MTLSIG file.
%
%          mtlwh(filename,...
%          comment1, filetype, info_blocks,...
%          n_channels, day_of_date, month_of_date, year_of_date,...
%          no_bl_per_chan, subject_no, ...
%          subject_first_name, subject_last_name, ...
%          n_lines, sampling_rate,...
%          lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
%          no_of_qu_bits, chan_code, ampl_gain, ...
%          max_ele, min_ele, status, res_type, res, ...
%          channels_per_dir, map, comment2)
%
% 	   filetype: 3 -> frq-domain
%    	             2 -> time-domain
%          info_blocks: number of info blocks
% 	   n_channels: number of channels in file
% 	   n_lines: number of lines per channel
% 	   sampling_rate: sampling rate in kHz
%          comment1: string with  79 characters
%          comment2: string with 262 characters
% 
% See also: MTLRH, MTLRCH, MTLWCH

%
% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%
% Date: 19.4.1995
% Modified: 19.01.1996 Klaus Hartung 
%           Using bifi and sigmap structure now  
% Modified: 8.7.96 Klaus Hartung
%           improved help
%           nargin check


% Calls		toolbox
% 
% no calls


if (nargin < 1) error('Filename is missing'); end;

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

if (nargin >27) res_type=BIFI_RES_WITH_DIRECTION_INFO; end;

%open file for reading and appending new data

[fid,message] = fopen(filename,'w+','ieee-le');
if fid==-1          % file can't be opened
  disp(message);
  return;
end;


%erase old info blocks or create new info blocks
stat=fseek(fid,0,'bof');
if (stat < 0) error('ERROR: mtlwh, erase old info blocks, fseek'); return; end;
d=zeros(1,256*info_blocks);
count = fwrite(fid,d,'uchar');
if (count < (256*info_blocks)) error('ERROR: mtlwh, erase old info blocks, fwrite'); return; end;



% write comment

if (filetype<=3)
   stat=fseek(fid,0,'bof');
   if (stat < 0) error('ERROR: mtlwh, write comment, fseek'); return; end;
   n = length (comment1);
   comment1=[comment1 zeros(1,79-n)];
   count=fwrite(fid,comment1,'uchar');
   if (count < length(comment1)) error('ERROR: mtlwh, write comment, fwrite'); return; end;
end;


%write EOF

stat=fseek(fid,79,'bof');
if (stat < 0) error('ERROR: mtlwh, write EOF, fseek'); return; end;
eof=26;
count = fwrite(fid,eof,'char');
if (count < length(eof)) error('ERROR: mtlwh, write EOF, fwrite'); return; end;


%write filetype

stat=fseek(fid,80,'bof');
if (stat < 0) error('ERROR: mtlwh, write filetype, fseek'); return; end;
count = fwrite(fid,filetype,'short');
if (count < length(filetype)) error('ERROR: mtlwh, write filetype, fwrite'); return; end;


%write number of info_blocks

stat=fseek(fid,82,'bof');
if (stat < 0) error('ERROR: mtlwh, write info_blocks, fseek'); return; end;
count = fwrite(fid,info_blocks,'ushort');
if (count < length(info_blocks)) error('ERROR: mtlwh, write info_blocks, fwrite'); return; end;


%write number of channels 

stat=fseek(fid,84,'bof');
if (stat < 0) error('ERROR: mtlwh, write n_channels, fseek'); return; end;
count = fwrite(fid,n_channels,'short');
if (count < length(n_channels)) error('ERROR: mtlwh, write n_channels, fwrite'); return; end;


% write day_of_date

stat=fseek(fid,86,'bof');
if (stat < 0) error('ERROR: mtlwh, fseek day_of_date') ;return, end;
count = fwrite(fid,day_of_date,'ushort');
if (count < length(day_of_date)) error('ERROR: mtlwh, write day of date, fwrite'); return; end;


% write month_of_date

stat=fseek(fid,88,'bof');
if (stat < 0) error('ERROR: mtlwh, fseek month_of_date') ;return, end;
count = fwrite(fid,month_of_date,'ushort');
if (count < length(month_of_date)) error('ERROR: mtlwh, write month of date, fwrite'); return; end;

% write year_of_date

stat=fseek(fid,90,'bof');
if (stat < 0) error('ERROR: mtlwh, fseek year_of_date') ;return, end;
count = fwrite(fid,year_of_date,'ushort');
if (count < length(year_of_date)) error('ERROR: mtlwh, write year of date, fwrite'); return; end;


% write number of blocks per channel

stat=fseek(fid,92,'bof');
if (stat < 0) error('ERROR: mtlwh, fseek blocks per channel');return, end;
count = fwrite(fid,no_bl_per_chan,'ushort');
if (count < length(no_bl_per_chan)) error('ERROR: mtlwh, write number_of_blocks per channels, fwrite'); return; end;



% the following data are optional
if (nargin > 9)

stat=fseek(fid,96,'bof');
if (stat < 0) error('ERROR: mtlwh, fseek subject animal number');return, end;
count = fwrite(fid,subject_no,'long');
if (count < length(subject_no)) error('ERROR: mtlwh, write subject number, fwrite'); return; end;

end

if (nargin > 10)

stat=fseek(fid,100,'bof');
if (stat < 0) error('ERROR: mtlwh, fseek subject name');return, end;
count = fwrite(fid,subject_first_name,'char');
if (count < length(subject_first_name)) error('ERROR: mtlwh, write subject first_name, fwrite'); return; end;

end



if (nargin > 11)

stat=fseek(fid,100,'bof');
if (stat < 0) error('ERROR: mtlwh, fseek subject last  name');return, end;
count = fwrite(fid,subject_last_name,'char');
if (count < length(subject_last_name)) error('ERROR: mtlwh, write subject latst name, fwrite'); return; end;

end


if (nargin > 12)

stat=fseek(fid,128,'bof');
if (stat < 0) error('ERROR: mtlwh, write n_lines, fseek'); return; end;
count = fwrite(fid,n_lines,'long');
if (count < length(n_lines)) error('ERROR: mtlwh, write n_lines, fwrite'); return; end;

end



if (nargin > 13)

%write sampling rate 
  
stat=fseek(fid,132,'bof');
if (stat < 0) error('ERROR: mtlwh, write sampling_rate, fseek'); return; end;
count = fwrite(fid,sampling_rate,'float32'); 
if (count < length(sampling_rate)) error('ERROR: mtlwh, write sampling_rate, fwrite'); return; end;


%write first_line
if (nargin > 14)

stat=fseek(fid,136,'bof');
if (stat < 0) error('ERROR: mtlwh, write lower cutoff frequency, fseek'); return; end;
count = fwrite(fid,lo_cutoff,'float32');
if (count < length(lo_cutoff) )error('ERROR: mtlwh, write lower cutoff frequency, fwrite'); return; end;

end

%write last_line

if (nargin > 15)

stat=fseek(fid,140,'bof');
if (stat < 0) error('ERROR: mtlwh, write higher cutoff frequnency, fseek'); return; end;
count = fwrite(fid,hi_cutoff,'float32');
if (count < length(hi_cutoff)) error('ERROR: mtlwh, write higher cutoff frequency, fwrite'); return; end;

end;


if (nargin > 16)

stat=fseek(fid,144,'bof');
if (stat < 0) error('ERROR: mtlwh, write effective S/N ratio, fseek'); return; end;
count = fwrite(fid,eff_s_to_n,'float32');
if (count < length(eff_s_to_n)) error('ERROR: mtlwh, write effective S/N ratio, fwrite'); return; end;

end;

if (nargin > 17)

stat=fseek(fid,148,'bof');
if (stat < 0) error('ERROR: mtlwh, write number of averagings, fseek'); return; end;
count = fwrite(fid,no_of_aver,'ulong');
if (count < length(no_of_aver)) error('ERROR: mtlwh, number of averagings, fwrite'); return; end;

end;


if (nargin > 18)

stat=fseek(fid,152,'bof');
if (stat < 0) error('ERROR: mtlwh, number of quantization bits, fseek'); return; end;
count = fwrite(fid,no_of_qu_bits,'ushort');
if (count < length(no_of_qu_bits)) error('ERROR: mtlwh,  number of quantization bits, fwrite'); return; end;

end;


if (nargin > 19)

stat=fseek(fid,154,'bof');
if (stat < 0) error('ERROR: mtlwh, channel code, fseek'); return; end;
count = fwrite(fid,chan_code,'ushort');
if (count < length(chan_code)) error('ERROR: mtlwh,  channel code, fwrite'); return; end;

end;

if (nargin > 20)

stat=fseek(fid,154,'bof');
if (stat < 0) error('ERROR: mtlwh, amplification by input VCA, fseek'); return; end;
count = fwrite(fid,ampl_gain,'float32');
if (count < length(ampl_gain)) error('ERROR: mtlwh, amplification by input VCA, fwrite'); return; end;

end;


if (map(1,1) ~= -999) 
if (nargin > 21)

stat=fseek(fid,220,'bof');
if (stat < 0) error('ERROR: mtlwh, magic number, fseek'); return; end;
count = fwrite(fid,BIFI_MAGIC_NUMBER,'ulong');
if (count < length(BIFI_MAGIC_NUMBER)) error('ERROR: mtlwh, magic number, fwrite'); return; end;

end;

if (nargin > 22)

stat=fseek(fid,224,'bof');
if (stat < 0) error('ERROR: mtlwh, max ele, fseek'); return; end;
count = fwrite(fid,max_ele,'float');
if (count < length(max_ele)) error('ERROR: mtlwh, max ele, fwrite'); return; end;

end;



if (nargin > 23)

stat=fseek(fid,228,'bof');
if (stat < 0) error('ERROR: mtlwh, min ele, fseek'); return; end;
count = fwrite(fid,min_ele,'float');
if (count < length(min_ele)) error('ERROR: mtlwh, min ele, fwrite'); return; end;

end;


if (nargin > 24)

stat=fseek(fid,232,'bof');
if (stat < 0) error('ERROR: mtlwh, status, fseek'); return; end;
count = fwrite(fid,status,'ushort');
if (count < length(status)) error('ERROR: mtlwh, status , fwrite'); return; end;

end;


if (nargin > 25)

stat=fseek(fid,234,'bof');
if (stat < 0) error('ERROR: mtlwh, res_type, fseek'); return; end;
count = fwrite(fid,res_type,'ushort');
if (count < length(res_type)) error('ERROR: mtlwh, res_type , fwrite'); return; end;

end;


if (nargin > 26)

stat=fseek(fid,236,'bof');
if (stat < 0) error('ERROR: mtlwh, resolution, fseek'); return; end;
count = fwrite(fid,res,'ushort');
if (count < length(res)) error('ERROR: mtlwh, resolution , fwrite'); return; end;

end;


if (nargin > 27)

stat=fseek(fid,238,'bof');
if (stat < 0) error('ERROR: mtlwh, channels per direction, fseek'); return; end;
count = fwrite(fid,channels_per_dir,'ushort');
if (count < length(channels_per_dir)) error('ERROR: mtlwh, channels per direction , fwrite'); return; end;

end;


if (nargin > 28)
n_dir = n_channels/channels_per_dir;
if (n_dir == 1) 
 stat=fseek(fid,240,'bof');
 if (stat < 0) error('ERROR: mtlwh, map for single direction, fseek'); return;  end;
 count = fwrite(fid,map(:,1),'float');
 if (count < length(map(:,1))) error('ERROR: mtlwh, map for single direction , fwrite'); return; end;
else
% n_dir
%ceil(n_dir*8/SIGM_BLOCK_LEN)
 if ((n_dir > 1) & (info_blocks >= 2 + ceil(n_dir*8/SIGM_BLOCK_LEN)))
    stat=fseek(fid,512,'bof');
    if (stat < 0) error('ERROR: mtlwh, map for multiple directions, fseek'); return;  end;
   count = fwrite(fid,map(:,1:n_dir),'float');
   if (count < length(map(:,1:n_dir))) error('ERROR: mtlwh, map for multiple directions , fwrite'); return; end;
 else
   error('ERROR: mtlwh, map for multiple directions, not enough info_blocks'); return;  end;
 end
end

end;

end;


if (nargin > 29)


%write comment2

stat=fseek(fid,248,'bof');
if (stat < 0) error('ERROR: mtlwh, write comment2, fseek'); return; end;
if (length(comment2) > 262) comment2=comment2(1:262); end;
count= fwrite(fid,comment2,'uchar');
if (count < length(comment2)) error('ERROR: mtlwh, write comment2, fwrite'); return; end;

end;


%write something to fill end of info block

stat= fseek(fid,511,'bof');
if (stat < 0) error('ERROR: mtlwh, write something, fseek'); return; end;
fwrite(fid,'a','char');
if (count < 1) error('ERROR: mtlwh, write something, fwrite'); return; end;


%close file

stat=fclose(fid);
if (stat < 0) error('ERROR: mtlwh, close file'); return; end;

return;



% end of file mtlwh.m