function STIM = genstim(HRIR,duration,Fs)

% GENSTIM generation of white noise stimuli convolved with HRIR 
%
% function STIM = genstim(HRIR, duration, Fs)
%
% generation of stimuli 
% convolution of broadband noise with all impulse responses in HRIR.
% the stimulus has the lenght of duration (in ms)

ramp = 5; %  in ms 
Fl = 200 ; % Hz
Fu = 15000; % Hz
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
nbnoise = conv(bbnoise,b);
%nbnoise = nbnoise(nbp:nbp+nsamp-1);
nbnoise = fade(nbnoise,nbp,nbp+nramp,nbp+nsamp-nramp,nbp+nsamp);

STIM = zeros(length(nbnoise)+length(HRIR(:,1))-1,nchannel);

for i=1:nchannel
 STIM(:,i) = conv(HRIR(:,i), nbnoise);
end

maximum = max(max(STIM))

STIM = STIM*(max_amp/maximum);

return


 
