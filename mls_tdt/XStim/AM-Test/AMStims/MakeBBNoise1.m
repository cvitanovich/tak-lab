function [bbnoise] = MakeBBNoise1(Fs,dur)
%MakeBBNoise:	Create a BroadBand Noise (2-11 kHz)
%Fs:		Sampling rate (Hz)
%dur:		Stimulus duration (ms)

dur = dur/1000;
len = dur*Fs;
bandwidth = 200; %bandwidth (Hz) of chunk for ILD adjustment

minfreq = 2000; maxfreq = 11000;
minfreq = round(((minfreq + 1)/Fs) * len);
maxfreq = round(((maxfreq + 1)/Fs) * len);
range = maxfreq-minfreq+1;
len = round(len); 
   
% mag spectrum = 1 between set frequencies:
mag = zeros(len,1);
mag(minfreq:maxfreq) = ones(range,1);
   

% random phase spectrum between set frequencies:
rand('state',sum(100*clock));
phi = (rand(len,1) - 0.5) * (2*pi);
% combine phase and magnitude:
X = mag .* ( (cos(phi)) + (i .* sin(phi)) );


% convert to time domain:
bbnoise = real(ifft(X));
bbnoise = (bbnoise/max(abs(bbnoise)))';

return;

