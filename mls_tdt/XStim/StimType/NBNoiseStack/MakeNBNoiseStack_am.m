function [nbnoisestack] = MakeNBNoiseStack(Fs,dur,bandlimits)
%MakeNBNoiseStack:	Create a NarrowBand Noise Stack
%Fs:			Sampling rate (Hz)
%dur:			Stimulus duration (ms)
%bandlimits:	N X 2 array of bandlimits, LO and HI

dur = dur/1000;

len = dur*Fs;
mag = zeros(len,1);

for bandnum = 1:size(bandlimits,1)
   minfreq = bandlimits(bandnum,1); maxfreq = bandlimits(bandnum,2);
   minfreq = round(((minfreq + 1)/Fs) * len);
   maxfreq = round(((maxfreq + 1)/Fs) * len);
   range = maxfreq-minfreq+1;
   % mag spectrum = 1 between set frequencies:
   mag(minfreq:maxfreq) = ones(range,1);
end

mag = 50 * mag;

% random phase spectrum between set frequencies:
rand('state',sum(100*clock));
phi = (rand(len,1) - 0.5) * (2*pi);
% combine phase and magnitude:
X = mag .* ( (cos(phi)) + (i .* sin(phi)) );


% convert to time domain:
nbnoisestack = real(ifft(X));

modsig = tone(55,Fs);
mod_depth = input('Modulation depth (0 to 1): ')
modenv = (1 - mod_depth/2) + (mod_depth/2)*modsig;
nbnoisestack = nbnoisestack .* modsig;
nbnoisestack = (nbnoisestack/max(abs(nbnoisestack)))';

return;