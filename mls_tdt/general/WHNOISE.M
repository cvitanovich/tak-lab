function [X,Xtime] = whnoise(minfreq,maxfreq,Fs,t,J)
%A script to produce noise with a flat spectrum between specified frequencies.
%Usage: [freq domain, time domain] = whnoise(min frequency, max frequency,Fs,t,seed)
%t is specified time in seconds;
%Borrowed from Kip Keller w_rnd.m

if nargin < 3 | isempty(Fs) Fs = 45000; end
if nargin < 4 | isempty(t) fact = 10; end
if nargin < 5 | isempty(J)
    J = sum(100*clock);
end

rand('state',J)

len = t*Fs;

minfreq = round(((minfreq + 1)/Fs) * len);
maxfreq = round(((maxfreq + 1)/Fs) * len);
range = maxfreq-minfreq+1;

% mag spectrum = 1 between set frequencies:
mag = zeros(len,1);
mag(minfreq:maxfreq) = ones(range,1);
mag = 50 * mag;

% random phase spectrum between set frequencies:
phase = (rand(len,1) - 0.5) * (2*pi);
% combine phase and magnitude:
X = mag .* ( (cos(phase)) + (i .* sin(phase)) );

% convert to time domain:
Xtime = real(ifft(X));

return;