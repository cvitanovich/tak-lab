% make a bunch of gammatones normalized to a +/- 1 range
% then run some stats
clear all
cd /Users/cvitanovich/Documents/MATLAB/tak-lab/HRTFs/Matlab_V6
uiopen('*.mat');

Fs=30000; cF=6000; species = 'owl';
dur=10; % 10 second sounds
target_rms=.25;
hWait=waitbar(0);
s=[];
nStates=100;
good_states=[];
j=1; cnt=0;
ramplen=5;
[coefs]=makeGammaFIR(Fs,cF,species);
while j < nStates
    cnt=cnt+1;
    rand('state',cnt);
    snd=rand(1,dur*Fs);
    snd=filtfilt(coefs,1,snd);
    snd=snd./max(abs(snd));
    mx=[]; r=[];
    for j=1:length(direc)
        EL=direc(1,j);
        AZ=direc(2,j);
        LT=TF1(j,:);
        RT=TF2(j,:);
        lt_snd=filter(LT,1,snd);
        rt_snd=filter(RT,1,snd);
        if(std(snd)>=target_rms)
            j=j+1;
            rms=std(snd);
            snd=(target_rms/rms).*snd;
            good_states=[good_states j];
        end
    end
    waitbar(j/nStates,hWait);
end
close(hWait)