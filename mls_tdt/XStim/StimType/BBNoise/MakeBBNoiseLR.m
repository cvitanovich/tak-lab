function [bbnoiseL,bbnoiseR] = MakeBBNoiseLR(Fs,dur,targ_freq,ild)
%MakeBBNoiseLR:	Create a BroadBand Noise (2-11 kHz) with an ILD-adjusted chunk
%Fs:		Sampling rate (Hz)
%dur:		Stimulus duration (ms)

dur = dur/1000;
len = dur*Fs;
bandwidth = 500; %bandwidth (Hz) of chunk for ILD adjustment

minfreq1 = targ_freq - (bandwidth/2);
maxfreq1 = targ_freq + (bandwidth/2);
minfreq1 = round(((minfreq1 + 1)/Fs) * len);
maxfreq1 = round(((maxfreq1 + 1)/Fs) * len);
range1 = maxfreq1-minfreq1+1;

minfreq2 = 2000; maxfreq2 = 11000;
minfreq2 = round(((minfreq2 + 1)/Fs) * len);
maxfreq2 = round(((maxfreq2 + 1)/Fs) * len);

% mag spectrum = 1 between set frequencies:
magL = zeros(len,1);

magL(minfreq2:minfreq1-1) = 1;
magL(maxfreq1+1:maxfreq2) = 1;
magR = magL;
magL(minfreq1:maxfreq1) = 1/sqrt(10^(ild/20));
magR(minfreq1:maxfreq1) = sqrt(10^(ild/20));

% random phase spectrum between set frequencies:
rand('state',sum(100*clock));
phi = (rand(len,1) - 0.5) * (2*pi);
% combine phase and magnitude:
XL = magL .* ( (cos(phi)) + (i .* sin(phi)) );
XR = magR .* ( (cos(phi)) + (i .* sin(phi)) );

% convert to time domain:
bbnoiseL = real(ifft(XL'));
bbnoiseR = real(ifft(XR'));

return;