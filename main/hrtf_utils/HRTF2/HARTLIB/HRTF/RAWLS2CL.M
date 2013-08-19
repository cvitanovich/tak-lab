function rawls2cl(filename)

% RAWLS2CL removes reflections and noise from impulse response (HRTF_RAW -> HRTF)

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum


% removal of dc part has to be included
disp('Load File')
eval(['load ' filename]);

disp('Calculation of LSIR_RAW')
LSIR_RAW = rifft(LSTF_RAW);


[m,p]=max(abs(LSIR_RAW));
min_p=min(p);
first_line=min_p-0.005*Fs;
last_line=first_line+1024-1;

LSIR_RAW2=zeros(1024,11);
LSIR_RAW2=LSIR_RAW(first_line:last_line,:);
LSIR_RAW=LSIR_RAW2;
clear LSIR_RAW2;


disp('Scaling of LSIR_RAW')
LSIR_RAW = LSIR_RAW/max(max(LSIR_RAW))*32767;

LSIR_RAW_MEAN=mean(LSIR_RAW);


clear LSTF_RAW
% 3 ms
n=3e-3*Fs;


disp('Calculation of LSIR')
LSIR     = optwin2(LSIR_RAW,n,512);


clear LSIR_RAW


disp('Calculation of LSTF')
LSTF     = rfft(LSIR);

[m,n] = size(LSTF);

N_LINES = m;
N_CHANNELS = n;
N_DIR = N_CHANNELS;


filename=[filename 'c'];
disp('Saving data')


eval(['save ' filename ' Fs N_LINES N_DIR N_CHANNELS LSTF LSIR ']);

