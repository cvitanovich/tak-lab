function wrstim(filename, STIM, subject_no, sampling_rate, map)

% writes STIM in filename 
%
%wrstim(filename, STIM, subject_no, sampling_rate, map)
%
% PC format, 16 bit



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


if (nargin < 3) subject_no = -1; end;
if (nargin < 4) sampling_rate = 44100; end;
if (nargin < 5) res_type = BIFI_RES_AS_MEASURED;
                res = -1;
 else    
                res_type = BIFI_RES_WITH_DIRECTION_INFO;
                res = -1;
end;

[m,n] = size(STIM);

comment1='ADDA_TYPE';
filetype=4;
infoblocks = 2;
if (res_type == 5) infoblocks= infoblocks+ceil((n/2*2*4)/256); end;
n_channels= n;
day_of_date = dofdate;
month_of_date = mofdate;
year_of_date = yofdate;
no_bl_per_chan = ceil(m/256);
% subject_no = 
subject_first_name ='';
subject_last_name ='';
n_lines = m;
%sampling_rate = -1;
lo_cutoff = -1;
hi_cutoff = -1;
eff_s_to_n = 9999;
no_of_aver = -1;
no_of_qu_bits = 16;
chan_code = -1;
ampl_gain = ones(16,1) * (-1);
magic_nr = BIFI_MAGIC_NUMBER;
max_ele = 9999;
min_ele = 9999;
status=-1;
res_type = BIFI_RES_WITH_DIRECTION_INFO;
res = -1;
channels_per_dir = 2;
%map = [-999; -999];
comment2 ='comment2';


mtlwh(filename,...
comment1, filetype, infoblocks,...
n_channels, day_of_date, month_of_date, year_of_date,...
no_bl_per_chan, subject_no, ...
subject_first_name, subject_last_name, ...
n_lines, sampling_rate,...
lo_cutoff,hi_cutoff, eff_s_to_n, no_of_aver, ... 
no_of_qu_bits, chan_code, ampl_gain, ...
max_ele, min_ele, status, res_type, res, ...
channels_per_dir, map, comment2);



multiplex = zeros(2*m,1);
index = 1:m;

for i=1:n/2
 mtlwt4(filename, STIM(:,2*i-1:2*i),i);
end

return
