
Fs=30000; min_freq=2000; max_freq=12000; dur=1;
% make a white noise base signal
% [X,Xtime] = whnoise(min_freq,max_freq,Fs,dur);
nTAPS=256;
% snd=Xtime;
% make a gammatone
% [coefs]=makeGammaFIR(Fs,6000,'owl');
% snd=rand(1,Fs*dur);
% snd=filtfilt(coefs,1,snd);
% snd=snd./max(abs(snd));
%t=0:(1/Fs):dur;
%snd=sin(2*pi*6000.*t);

% load hrtf coeficients
cd /Users/cvitanovich/Documents/MATLAB/tak-lab/HRTFs
uiopen('*.mat');

azvals=-90:5:90; elvals=90:-5:-90;
Xes=repmat(azvals,length(elvals),1);
Yes=repmat(elvals',1,length(azvals));
ITDs=NaN*ones(length(azvals),length(elvals));
ILDs=NaN*ones(size(ITDs));
%hPlot=figure;
%hITD=subplot(2,1,1); hold on; axis square;
%hILD=subplot(2,1,2); hold on; axis square;
cnt=0;
hWait=waitbar(0,'crunching...');
for j=1:length(direc)
    cnt=cnt+1;
    waitbar((cnt/length(direc)),hWait);
    EL=direc(1,j);
    AZ=direc(2,j);
    LT=TF1(j,:);
    RT=TF2(j,:);
    idx1=find(Xes==AZ);
    idx2=find(Yes==EL);
    idx=intersect(idx1,idx2);
    %left_snd=filter(LT,1,snd);
    %right_snd=filter(RT,1,snd);
    tmp=xcorr(LT,RT);
    mx=find(tmp==max(tmp));
    ms=mx-(length(mx)/2);%(mx-Fs*dur)*(10^6)*(1/Fs);
    ITDs(idx)=ms;
    % calc ILD
    left_rms=sqrt(mean(LT.^2));
    right_rms=sqrt(mean(RT.^2));
    ILDs(idx)=right_rms-left_rms;   
end
close(hWait)
figure;
colormap hot
subplot(2,1,1); pcolor(Xes,Yes,ITDs); colorbar
subplot(2,1,2); pcolor(Xes,Yes,ILDs); colorbar
