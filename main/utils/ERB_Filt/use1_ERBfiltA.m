function [output] = use1_ERBfiltA(input, Fs, cF, pflag)

%function [output] = use1_ERBfiltA(input, Fs, cF);
%
% uses the June 1998 implementation of gammatone filtering by Slaney
% to filter "input" through a series of four filters approximating
% the owl's critical bands
%
% input: vector to be filtered
% Fs: sampling frequency
% cF: center frequency for filter
% pflag: >0) plot impulse-response
%
% output has been filtered
% by a gammatone, whose impulse response has the form:
%time = 0:.0001:.02;
%bw = 125;			%1.019 * critical band;
%cF = 1000;			%centerFrequency;
%H = (time .^3) .* (exp(1) .^(-2*pi*bw*time)) .* cos(2 * pi * cF .* time);
%plot(H);
%April 2, 1999
%
% updated to use three-stage filtering process

if nargin < 4
   pflag =0;
end
if nargin < 3
   error('insufficient input arguments');
   return
end

% make filter
fcoefs=make_ERBFiltA(Fs,cF);

% apply each of four sequential filters
A0  = fcoefs(1);
A11 = fcoefs(2);
A12 = fcoefs(3);
A13 = fcoefs(4);
A14 = fcoefs(5);
A2  = fcoefs(6);
B0  = fcoefs(7);
B1  = fcoefs(8);
B2  = fcoefs(9);
gain = fcoefs(10);	

if ~isempty(input)
	y1 = filter([A0/gain A11/gain A2/gain], [B0 B1 B2], input);
	y2 = filter([A0 A12 A2], [B0 B1 B2], y1);
	y3 = filter([A0 A13 A2], [B0 B1 B2], y2);
	output = filter([A0 A14 A2], [B0 B1 B2], y3);
else
   disp('input vector is empty');
   output = [];
   return
end

if pflag > 0
	impulse = [1 zeros(1, 511)];
	y1 = filter([A0/gain A11/gain A2/gain], [B0 B1 B2], impulse);
	y2 = filter([A0 A12 A2], [B0 B1 B2], y1);
	y3 = filter([A0 A13 A2], [B0 B1 B2], y2);
   y4 = filter([A0 A14 A2], [B0 B1 B2], y3);
   figure
   pts = length(impulse)
   temp = fft(y4,pts*2);
   semilogx(0:(Fs/2)/(pts-1):(Fs/2),20*log10(abs(temp(1:pts))));
   title(['Impulse response of Gammatone filter with cF = ' num2str(cF)]);
   xlabel('Frequency')
   ylabel('dB magnitude')
end
