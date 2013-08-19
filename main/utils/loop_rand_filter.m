% initialize random no. generator

rand('state',0);
% make coefs
Fs=30000; cF=6000; species='owl';
[coefs]=makeGammaFIR(Fs,cF,species);
dur=Fs; % 1 sec buffers
rand('state',0); s=rand('state');
t=1;
nTAPS=length(coefs);
x=zeros(1,dur+nTAPS); % input buffer
y=zeros(1,dur); % out
tmp=zeros(1,dur+nTAPS);

tm=[];
for i=1:1000
    tic
    rand('state',s);
    data=rand(1,dur);
    s=rand('state');
    x=[x((end-nTAPS+1):end) data]; % circular buffer for data (shift by duration of each data buffer)
    tmp=conv(coefs,x);
    y=tmp(nTAPS+1:dur+nTAPS);
    t=t+dur;
    tm(i)=toc;
end
mean(tm)

