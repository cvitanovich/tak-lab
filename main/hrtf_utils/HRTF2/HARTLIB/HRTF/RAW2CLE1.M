function raw2cle1(filename)

% RAW2CLE1  remove reflections and noise from impulse response 
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

%disp('Limited version')
%HRTF_RAW=HRTF_RAW(:,92);

disp('Calculation of HRIR_RAW')
HRIR_RAW = rifft(HRTF_RAW);

disp('Scaling of HRIR_RAW')
HRIR_RAW = HRIR_RAW/max(max(HRIR_RAW))*32767;

[m,p]=max(abs(HRIR_RAW));
[min_p,pp]=min(p)
first_line=min_p-0.005*Fs;
last_line=first_line+1024-1;

HRIR_RAW_MEAN=mean(HRIR_RAW(first_line-80:first_line+80,:))
size(HRIR_RAW_MEAN)

HRIR_RAW2=zeros(1024,11);
HRIR_RAW2=HRIR_RAW(first_line:last_line,:);
HRIR_RAW =HRIR_RAW2;
clear HRIR_RAW2;

[n,n_chan]=size(HRIR_RAW);
for jc=1:n_chan
   HRIR_RAW(:,jc)=HRIR_RAW(:,jc)-HRIR_RAW_MEAN(jc)*ones(n,1);
end

clear HRTF_RAW
% 3 ms
n=3e-3*Fs;

disp('Calculation of HRIR')
HRIR     = optwin2(HRIR_RAW,n,256);

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