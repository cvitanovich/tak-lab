% make a bunch of gammatones normalized to a +/- 1 range
% then run some stats
Fs=30000; cF=6000; species = 'owl';
dur=4; % 4 second sound
target_rms=.25;
hWait=waitbar(0);
s=[];
nStates=50;
good_states=[];
j=0; cnt=0;
ramplen=5;
[coefs]=makeGammaFIR(Fs,cF,species);
while j < nStates
    cnt=cnt+1;
    tic
    rand('state',cnt);
    snd=rand(1,dur*Fs);
    snd=filtfilt(coefs,1,snd);
    snd=snd./max(abs(snd));
    if(std(snd)<=target_rms)
        j=j+1;
        rms=std(snd);
        snd=(target_rms/rms).*snd;
        good_states=[good_states j];
    end
    ramped_snd = rampMySound(snd,ramplen,Fs);
    toc
    s(j)=std(snd);
    waitbar(j/nStates,hWait);
end
close(hWait)