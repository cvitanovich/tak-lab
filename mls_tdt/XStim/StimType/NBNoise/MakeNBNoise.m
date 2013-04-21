function [nbnoise] = MakeNBNoise(Fs,dur,bandlimits)
%MakeBBNoise:	Create a BroadBand Noise (2-11 kHz)
%Fs:		Sampling rate (Hz)
%dur:		Stimulus duration (ms)

dur = dur/1000;

len = dur*Fs;

minfreq = bandlimits(1); maxfreq = bandlimits(2);
minfreq = round(((minfreq + 1)/Fs) * len);
maxfreq = round(((maxfreq + 1)/Fs) * len);
range = maxfreq-minfreq+1;

% mag spectrum = 1 between set frequencies:
mag = zeros(len,1);
mag(minfreq:maxfreq) = ones(range,1);
mag = 50 * mag;

% random phase spectrum between set frequencies:
rand('state',sum(100*clock));
phi = (rand(len,1) - 0.5) * (2*pi);
% combine phase and magnitude:
X = mag .* ( (cos(phi)) + (i .* sin(phi)) );


% convert to time domain:
nbnoise = real(ifft(X));
nbnoise = (nbnoise/max(abs(nbnoise)))';

return;