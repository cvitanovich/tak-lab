function [forward, feedback]=make1_ERBfilt(cf,Fs,p)

% [forward, feedback]= make1_ERBfilt(cf,Fs,p)
% computes the filter coefficients for a single Gammatone filter.
% These filters were defined by Patterson and Holdworth for simulating 
% the cochlea.  The result is returned as a set of filter
% coefficients (forward and feedback), which
% can be passed to the MatLab "filter" function.
% March 8, 1999
% argins:
%		cf: center frequency for filtert
%		Fs: sampling frequency
%		p: 0) plot or 1) not

if nargin <2   	Fs = 30000;	end
if nargin <3 		p=0; 			end
fftpts=1024*8;

T=1/Fs;
% EarQ = 2.5 * ...    [reflects change from Q10dB]
% to something near critical band (ca. Q3dB ??)

EarQ = 1.9 * 0.074 * (cf .^ .504);
%EarQ = 2.5 * 0.074 * (cf .^ .504);

ERB = cf ./ EarQ;

B=1.019*2*pi*ERB;
gain = abs((-2*exp(4*i*cf*pi*T)*T + ...
				2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
				(cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
				sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
               sin(2*cf*pi*T))).* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) - ...
               sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
           (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
          (-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +  ...
           2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);
forward(1) = T^4 ./ gain;
forward(2) = -4*T^4*cos(2*cf*pi*T)./exp(B*T)./gain;
forward(3) = 6*T^4*cos(4*cf*pi*T)./exp(2*B*T)./gain;
forward(4) = -4*T^4*cos(6*cf*pi*T)./exp(3*B*T)./gain;
forward(5) = T^4*cos(8*cf*pi*T)./exp(4*B*T)./gain;
feedback(1) = ones(length(cf),1);
feedback(2) = -8*cos(2*cf*pi*T)./exp(B*T);
feedback(3) = 4*(4 + 3*cos(4*cf*pi*T))./exp(2*B*T);
feedback(4) = -8*(6*cos(2*cf*pi*T) + cos(6*cf*pi*T))./exp(3*B*T);
feedback(5) = 2*(18 + 16*cos(4*cf*pi*T) + cos(8*cf*pi*T))./exp(4*B*T);
feedback(6) = -8*(6*cos(2*cf*pi*T) + cos(6*cf*pi*T))./exp(5*B*T);
feedback(7) = 4*(4 + 3*cos(4*cf*pi*T))./exp(6*B*T);
feedback(8) = -8*cos(2*cf*pi*T)./exp(7*B*T);
feedback(9) = exp(-8*B*T);

% plot the response to a square pulse
if p > 0
impulse = [1 zeros(1,1023)];
y = filter(forward,feedback,impulse);
figure
response = 20*log10(abs(fft(y',fftpts)));
freqScale = (0:fftpts-1)/fftpts*Fs;
axis([2 4 -70 10]);
semilogx(freqScale(1:fftpts/2),response(1:fftpts/2,:));
ylabel ('dB magnitude')
xlabel ('Frequency Hz')
title ('ERB gammatone filter bank')
end