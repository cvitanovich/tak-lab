function raw2clea(filename)

% RAW2CLEA  remove reflections and noise from impulse response 
%            (HRTF_RAW -> HRTF)
%
%            raw2clean(filename) will load the mat-file with HRTF_RAW 
%            shorten the impulse response to 3 ms and store the result 
%            in the matrix HRTF in the file filenamec.mat
%

% removal of dc part has to be included

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% 
% Date: ?
% Modified: 7.7.96
%      improved help
%      removal of debug information

% Calls 	Toolbox
% rifft		hartung/loafo2
% optwin	hartung/tools
% rfft		hartung/loafo2
%

if (nargin < 1) error('Insufficient numer of arguments. Specify filename.'); end;

disp('Load File')
eval(['load ' filename]);

disp('Calculation of HRIR_RAW')
HRIR_RAW = rifft(HRTF_RAW);

disp('Scaling of HRIR_RAW')
HRIR_RAW = HRIR_RAW/max(max(HRIR_RAW))*32767;


clear HRTF_RAW
% 3 ms
n=3e-3*Fs;

disp('Calculation of HRIR')
HRIR     = optwin(HRIR_RAW,n,512);

disp('Calculation of HRTF')
HRTF     = rfft(HRIR);

[m,n] = size(HRTF);

N_LINES = m;
N_CHANNELS = n;
N_DIR = N_CHANNELS/2;


filename=[filename 'c'];
disp('Saving data')


eval(['save ' filename ' Fs N_LINES N_DIR N_CHANNELS MAP HRTF HRIR ']);




% end of file raw2clean.m 
