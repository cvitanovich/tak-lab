% OPTIHRTF ?

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%27.3.96


n_opt=31;

%HRIR=HRIR/(max(max(HRIR)));

HRTF_128_H=rfft(HRIR_128,4096);
HRIR_OPT=zeros(size(HRIR_128));
HRIR_GUESS=zeros(n_opt+1,244);
tofs=zeros(1,244);

for n=1:244
clc
disp(['Calculating channel ' int2str(n) '.']);
	[HRIR_GUESS(1:n_opt+1,n),tofs(n)]=firguess(HRIR_128(:,n),n_opt);
	b=firopt(HRTF_128_H(:,n),getsampl('m'),HRIR_GUESS(:,n),2,20);
        HRIR_OPT(tofs(n):tofs(n)+n_opt,n)=b';
end

HRTF_GUESS = rfft(HRIR_GUESS,4096);
HRTF_OPT = rfft(HRIR_OPT,4096);

