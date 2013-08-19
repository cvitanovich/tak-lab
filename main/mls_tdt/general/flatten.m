function [output] = flatten(input, p);

%function [output] = flatten(input, p);
%
% returns input, now with a flattened
% envelope.  Such that the original carrier
% spectrum is little changed
% output is normalized to the Vrms of input
% Mar 19, 1999
% argin p) 0) don't plot    1) plot env spectra 2) print moments

if nargin <2
   p=0;
end

% make sure input is centered at mean = zero
input = input - mean(input);

% get envelope of input
env = abs(hilbert(input));

% make & apply inverse envelope
Ienv = max(env)./ env;
output = input .* Ienv;

% normalize to RMS
output = output - mean(output);   
M = mom(output);
output = output / (1.414 * M(2));

% plot spectra of envelopes'
if p ==1
   freq = 0:15000/(4095):15000;
   FTin = fft(env,1024*8);
   Outenv = abs(hilbert(output));
   FTout = fft(Outenv,1024*8);
   figure
   plot(freq(1:500),20*log10(abs(FTout(1:500))),'r');
   hold on
   plot(freq(1:500),20*log10(abs(FTin(1:500))));
   xlabel('frequency of modulation (Hz)')
   ylabel('dB')
end

% display moments
if p ==2
	disp('Moments.  Before : ')
   M = mom(input)
   disp('After : ');
   M = mom(output)
end

