function rawhead2(filename)

% RAWHEAD2 remove reflections and noise from headphone impulse response
%            (HPTF_RAW -> HPTF)
%
%            rawhead3(filename) will load the mat-file with HPTF_RAW 
%            shorten the impulse response to 128 samples and stores the result 
%            in the matrix HPTF in the file filenamec.mat
%


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
%
% rifft		hartung/loafo2
% optwin	hartung/tools
% rfft		hartung/loafo2
%

if (nargin < 1) error('Insufficient numer of arguments. Specify filename.'); end;

disp('Load File')
eval(['load ' filename]);

disp('Calculation of HRIR_RAW')
HPIR_RAW = rifft(HPTF_RAW);

disp('Scaling of HRIR_RAW')
HPIR_RAW = HPIR_RAW/max(max(HPIR_RAW))*32767;

%HPIR_RAW_MEAN=mean(HPIR_RAW);
%HPIR_RAW=HPIR_RAW-HPIR_RAW_MEAN;

disp('Calculation of HPIR')
HPIR     = optwin(HPIR_RAW,89,128);

disp('Calculation of HPTF')
HPTF     = rfft(HPIR);

[m,n] = size(HPTF);

N_LINES = m;
N_DIR=0;
N_CHANNELS=n;

filename=[filename 'c'];

eval(['save ' filename '    Fs N_LINES N_DIR N_CHANNELS HPTF HPIR']);

% end of file rawhead2.m
