function eq_mic(file)

%filename=[file 'c']

%eval(['load ' filename])
load ../MAP
clear HRIR

filename=[file 'eq']

eval(['load ' filename])

load ../ls
load ../mic

HRIR_MIC=equ2(HRIR_LS,mic_ir_clean,[14500,16500]);

[m,n]=size(HRIR_MIC);

for k=1:n
 H(:,k)=conv(HRIR_MIC(:,k),LSIR(:,2));
 %irplot(HRIR_MIC(:,k),Fs);
 %irplot(H(:,k),Fs);
 %pause
end;

clear HRIR_MIC HRIR_LS mic_ir

HRIR=optwin(H,512,1024);
HRTF=rfft(HRIR);
%HRTF=HRTF/HRTF(2,1);
HRTF=HRTF/HRTF(9,1);
HRIR=rifft(HRTF);
outfile=[file 'mic'];

eval(['save ' outfile ' Fs  HRIR HRTF MAP' ]);