function mat2mtl(file)

% MAT2MTL convert matlab file to mtl-, bifi- or sigmaps-file
%
%         mat2mtl(a) will read the matlab file a and create a
%         bifi-, sigmaps-, mtl-file a.srx.

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date 26.09.1995
% Modified: 7.7.1996
%       help modified
%       debug information removed

% Calls		Toolbox
% 
% mtlwch	hartung/mtl
% mtlwh		hartung/mtl


if (nargin<1) error('Filenname is not defined.'); end;

eval(['load ' file]);

maxamp=max(max(abs(HRIR)));

HRIR=HRIR/maxamp;

%N_LINES=N_LINES-1

filename = [file '.srx'];

filetype=2;
first_line=0;
last_line=Fs/2;

comment1=[];
comment2=[];

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

info_blocks = 2  + ceil(N_DIR*2*4/256);

[x,y]=size(HRIR);

comment1 = 'this is a test';
eof = -1;
filetype = 2;
n_channels = y;
day_of_date = -1;
month_of_date = -1;
year_of_date = -1;
no_bl_per_chan = -1;
subject_no = 0;
subject_first_name = 'NON';
subject_last_name = 'NOMINATUR';
n_lines = x;
sampling_rate = -1;
lo_cutoff = -1;
hi_cutoff = -1;
eff_s_to_n = 9999;
no_of_aver = -1;
no_of_qu_bits = -1;
chan_code = -1;
ampl_gain = ones(16,1) * (-1);
magic_nr = -1;
max_ele = -9999;
min_ele = -9999;
status = BIFI_TOTAL_EQ;
res_type = BIFI_RES_WITH_DIRECTION_INFO;
res = -1;
channels_per_dir = 2;
map = [-999; -999];
comment2 ='comment2';

sampling_rate =Fs;
channels_per_dir=N_CHANNELS/N_DIR;

min_ele=min(MAP(1,:));
max_ele=max(MAP(1,:));

mtlwh(filename,...
comment1, filetype, info_blocks,...
N_CHANNELS, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
N_LINES, Fs,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, MAP, comment2);
 



%status=mtlwdir(filename,MAP);

for n=1:N_CHANNELS
  mtlwch(filename,HRIR(:,n),n);
end;


% end of file mat2mtl.m