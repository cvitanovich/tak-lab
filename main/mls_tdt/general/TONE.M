function [X,Xtime] = tone(freq,Fs,fact)
%A script to produce noise with a flat spectrum between specified frequencies.
%Usage: [freq domain, time domain] = whnoise(freq,sampling rate,factor)
%factor is used to determine length of sound; len = factor * Fs
%Borrowed from Kip Keller w_rnd.m

if (nargin < 2) Fs = 45000; end
if (nargin < 3) fact = 10; end

Nyq = Fs;
len = fact*Fs;

freq = round(((freq + 1)/Nyq) * len);

% mag spectrum = 1 between set frequencies:
mag = zeros(len,1);
mag(freq) = 1;
mag = 50 * mag;

% random phase spectrum between set frequencies:
phase = (rand(len,1) - 0.5) * (2*pi);
% combine phase and magnitude:
X = mag .* ( (cos(phase)) + (i .* sin(phase)) );


% convert to time domain:
Xtime = ifft(X);
Xtime = real(Xtime/max(abs(Xtime)));

return;