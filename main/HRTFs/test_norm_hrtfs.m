% test normalized hrtfs (to determine largest possible rms value)
clear all
cd /Users/cvitanovich/Documents/MATLAB/tak-lab/HRTFs/Matlab_V6
fnames={'924AD_eq_ABLequal_normed.mat','925AD_eq_ABLequal_normed.mat', ...
    '1073AC_eq_ABLequal_normed.mat'};
mx=[]; R=[]; t=[];
nLocs=613;
nStates=16; tot=nStates*length(fnames); cnt=0;
figure; axis off;
for f=1:length(fnames)
    if(exist('txt0')>0)
        delete(txt0)
    end
    txt0=text(.5,.8,['BIRD ' num2str(f) ' of ' num2str(length(fnames))]);
    drawnow;
    load(fnames{f});
    SND=[];
    for k=1:nStates
        rand('state',k);
        
        % make a white noise base signal
        %[X,Xtime] = whnoise(2000,11000,30000,1);
        % normalize it
        % snd=Xtime./max(abs(Xtime));
        
        % % make a gammatone
        Fs=30000; cF=6000; species = 'owl'; dur=6;
        [coefs]=makeGammaFIR(Fs,cF,species);
        snd=rand(1,dur*Fs);
        snd=filtfilt(coefs,1,snd);
        snd=snd./max(abs(snd));
        SND=[SND snd'];
        % make octave band noise
        % [X,Xtime] = whnoise(4000,8000,30000,1);
        % %normalize it
        % snd=Xtime./max(abs(Xtime));
    end
    t=[];
    for j=1:length(direc)
        tic
        EL=direc(1,j);
        AZ=direc(2,j);
        LT=TF1(j,:);
        RT=TF2(j,:);
        L_SNDS=filter(LT,1,SND);
        R_SNDS=filter(RT,1,SND);
        M=max(abs([L_SNDS; R_SNDS]));
        for p=1:length(M)
            L_SNDS(:,p)=L_SNDS(:,p)./M(p);
            R_SNDS(:,p)=R_SNDS(:,p)./M(p);
        end
        R=[R (sqrt(mean(L_SNDS.^2))+sqrt(mean(R_SNDS).^2))./2];
        t(j)=toc;
        rem_time=mean(t)*(length(direc)-j);
        if(exist('txt')>0)
            delete(txt)
        end
        txt=text(.5,.5,['EST REM TIME (THIS BIRD): ' num2str(rem_time) ' secs']);
        drawnow;
    end
end

close(hWait)
max(r)