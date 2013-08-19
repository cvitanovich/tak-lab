function PINK_STIM = genpink(HRIR,duration,Fs)

% GENPINK generation of pink noise stimuli convolved with HRIR 
% 
% function STIM = genpink(HRIR, duration, Fs)
%
% convolution of pink noise with all impulse responses in HRIR.
% the stimulus has the lenght of duration (in ms)

ramp = 5; %  in ms 
Fl = 100 ; % Hz
Fu = 16000; % Hz
max_amp = 16000;

nramp = Fs * (5 /1000);
nsamp = Fs * (duration / 1000);
nbp = 128;

[m,n] = size(HRIR); 
nchannel = n;

b = fir1(nbp,[Fl Fu]*2/Fs);

%window1 = (sin( (0:nramp-1)* (pi/(2*nramp)) )).^2;
%window2 = ones((nsamp-2*nramp+1),1);
%window3 = (cos( (0:nramp-1)* (pi/(2*nramp)) )).^2;
%window = [ window1' ;window2;  window3' ];


bbnoise = randn(nsamp*1.2,1);
f=linspace(0,1,256);
m(2:256)=1./sqrt(f(2:256));
m(1)=1;

pink=fir2(256,f,m);
pnoise = conv(bbnoise,pink);

nbnoise = conv(pnoise,b);
%nbnoise = nbnoise(nbp:nbp+nsamp-1);
nbnoise = fade(nbnoise,nbp,nbp+nramp,nbp+nsamp-nramp,nbp+nsamp);

STIM = zeros(length(nbnoise)+length(HRIR(:,1))-1,nchannel);

for i=1:nchannel
 STIM(:,i) = conv(HRIR(:,i), nbnoise);
end

maximum = max(max(STIM))

PINK_STIM = STIM*(max_amp/maximum);

return


 
