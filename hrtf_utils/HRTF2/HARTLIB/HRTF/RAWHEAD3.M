function rawhead3(filename)

% RAWHEAD3 remove reflections and noise from headphone impulse response
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

HPIR_RAW_MEAN=mean(HPIR_RAW);

[k,j]=size(HPIR_RAW);
%for jc=1:j
%  HPIR_RAW(:,jc)=HPIR_RAW(:,jc)-HPIR_RAW_MEAN(jc)*ones(k,1);
%end

[m,p]=max(abs(HPIR_RAW));
min_p=min(p);
first_line=min_p-0.005*Fs;
last_line=first_line+1024-1;

HPIR_RAW2=zeros(1024,11);
HPIR_RAW2=HPIR_RAW(first_line:last_line,:);
HPIR_RAW=HPIR_RAW2;
clear HPIR_RAW2;


disp('Calculation of HPIR')
HPIR     = optwin2(HPIR_RAW,89,128);

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