function [X,Xtime] = seeded_whnoise(minfreq,maxfreq,Fs,duration,seedval,phi)
%A script to produce noise with a flat spectrum between specified frequencies.
%Usage: [freq domain, time domain] = whnoise(min frequency, max frequency,Fs,durationor)
%durationor is used to determine length; len = duration * Fs;
%Borrowed from Kip Keller w_rnd.m

if (nargin < 3) Fs = 30000; end
if (nargin < 4) duration = .1; end
if (nargin < 5) randn('seed',47);
else
    randn('seed',seedval);
end
if (nargin < 6) phase = (rand(duration*Fs,1) - 0.5) * (2*pi);
else
phase = phi;
end

len = duration*Fs;

minfreq = round(((minfreq+1)/Fs) * len);
maxfreq = round(((maxfreq+1)/Fs) * len);
range = maxfreq-minfreq+1;

% mag spectrum = 1 between set frequencies:
mag = zeros(len,1);
mag(minfreq:maxfreq) = ones(range,1);

% random phase spectrum between set frequencies:
			%-0.5 moves phases from 0to1 range to -0.5to0.5 range, then have to multiply by 2pi instead of pi
% combine phase and magnitude:
X = mag .* ( (cos(phase)) + (i .* sin(phase)) );

% convert to time domain:
Xtime = real(ifft(X));
Xtime=Xtime/max(abs(Xtime));						%have to normalize so that max and mins are between -1 and 1 otherwise will saturate = bad


while 0
valnoise = Xtime;
valnoise = abs(hilbert(valnoise));
valnoise = valnoise*10000;
valnoise = valnoise-mean(valnoise);
figure
% Use next highest power of 2 greater than or equal to length(x) to calculate FFT.
nfft= 2^(nextpow2(length(valnoise))); 
% Take fft, padding with zeros so that length(fftx) is equal to nfft 
fftx = fft(valnoise,nfft); 
% Calculate the numberof unique points
NumUniquePts = ceil((nfft+1)/2); 
% FFT is symmetric, throw away second half 
fftx = fftx(1:NumUniquePts); 
% Take the magnitude of fft of x and scale the fft so that it is not a function of the length of x
mx = abs(fftx)/length(valnoise); 
% Take the square of the magnitude of fft of x. 
mx = mx.^2; 
% Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
% The DC component and Nyquist component, if it exists, are unique and should not be multiplied by 2.
if rem(nfft, 2) % odd nfft excludes Nyquist point
  mx(2:end) = mx(2:end)*2;
else
  mx(2:end -1) = mx(2:end -1)*2;
end
% This is an evenly spaced frequency vector with NumUniquePts points. 
f = (0:NumUniquePts-1)*Fs/nfft; 
% normalize mx
mx = mx/max(mx);

mx = 20*log10(mx);
% Generate the plot, title and labels. 
plot(f,mx);
axis([0 20000 -100 100])
title('Power Spectrum'); 
xlabel('Frequency (Hz)');
ylabel('Power');
end

return;