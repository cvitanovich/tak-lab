% script to produce a "wrap-around" (a.k.a. circularly filtered) g-tone sound
% params
sound_dir='C:/andrew/sounds/';
ADAPTOR.Fs=30000;
ADAPTOR.Dur=5; % seconds
ADAPTOR.cF=6000;
ADAPTOR.species='owl';
ADAPTOR.target_rms=0.1; % target rms for adaptor

% first produce a randomized sequence
X=rand(1,ADAPTOR.Dur*ADAPTOR.Fs);
len=length(X);

% second get gtone coeffs
B = makeGammaFIR(ADAPTOR.Fs,ADAPTOR.cF,ADAPTOR.species);

% third filter random sequence with coeffs (wrap around sequence at end)
N=length(B); % nTAPs
Y=NaN.*ones(1,len); % output vector
range=0:(N-1); % range for filter input
cnt=0; % loop counter
j=N-1; % idx of output val
% filter loop
h = waitbar(0,'Please wait...');
while(cnt < len)
    % calculate finite output value
    Y(j+1) = B * X(range+1)';
    % increment finite output value
    j=j+1; j=mod(j,len);
    % increment input range
    range=range+1; range=mod(range,len);
    % increment loop count
    cnt=cnt+1;
    waitbar(cnt/len,h)
end
close(h);

% normalize and then scale to desired rms
Y=0.9999.*(Y./max(abs(Y)));
rms=sqrt(mean(Y.^2));
S=(ADAPTOR.target_rms/rms).*Y;

% check if final sound is in proper range (+/- 0.9999)
if(max(abs(S))>0.9999)
    warndlg('Oops! Final sound outside desired range! (+/- 0.9999)');
else
    Zval=double('Z'); Aval=double('A');
    fname='gtone_A';
    CharVals=Aval:Zval;
    cnt=1; k=length(fname);
    while(exist([sound_dir fname '.mat']))
        fname(k)=char(CharVals(cnt));
        cnt=cnt+1;
        if(cnt>length(CharVals))
            cnt=1; k=k+1;
        end
    end
    disp(['Saving sound/params in file: ' sound_dir fname '.mat']);
    save([sound_dir fname],'S','ADAPTOR');
end