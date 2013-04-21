function IRout = altmagIR (IRin, minF, maxF, Fact)

%function IRout = altmagIR (IRin, minF, maxF, Fact)
% adds noise to magnitude spectrum of IR within frequency band
% argins:
%		IRin
%		minF
%		maxF
%		Fact (+/- dB noise to add)

if nargin < 4		Fact = 5;	end
if nargin < 3		maxF = 6000;	end
if nargin < 2		minF = 2000;	end
pflag = 1;
pflag1 = 0;
momflag = 0;            % equalizes by 2nd moment (leaving flat part of spectrum lowered)

fftpts = 255;
rand ('state',sum(100*clock));
freq = 0:15000/(fftpts/2 -1):15000;
minF = nearest_index(freq,minF);
maxF = nearest_index(freq,maxF);

if nargin < 1 | isempty('IRin')
	% read an HRTF in the time domain for testing purposes only
	FN = 'e:\hrtfdata\901\out8a_ad';
	IRin = mtlrch(FN,200*2-1);
end

% pts to change
Npts = maxF-minF+1;
FT_IR = fft(IRin,fftpts);
FT_new = FT_IR;
Xpart = FT_IR(minF:maxF);

temp = (rand(Npts,1) *2-1);
temp = temp .* tukeywin(length(temp),.1) *Fact;		% dB of change desired
temp = 10 .^ (temp /20);							% linear change desired

FT_new(minF:maxF) = Xpart .* temp;

if pflag
% plot magnitude
figure; hold on;
plot(freq,20*log10(abs(FT_new(1:fftpts/2))),'r');
plot(freq,20*log10(abs(FT_IR(1:fftpts/2))),'b');
title('HRTF before/after alteration')
xlabel('freq (Hz)')
ylabel('magnitude (dB)')
end
if pflag1
% plot phase
figure; hold on;
plot(freq,angle(FT_new(1:fftpts/2)),'r');
plot(freq,angle(FT_IR(1:fftpts/2)),'b');
title('HRTF before/after modification')
xlabel('freq (Hz)')
ylabel('phase angle (radians)')
end

% calc new IR
IRtemp = real(ifft(FT_new));
IRout = IRtemp(1:255);
IRout([1 255]) = [0 0];
if momflag
    IRout = IRout * (mom(IRin,2)/mom(IRout,2));
end
FT_IRout = fft(IRout,fftpts);

if pflag
% plot magnitude
figure; hold on;
plot(freq,20*log10(abs(FT_IRout(1:fftpts/2))),'r');
plot(freq,20*log10(abs(FT_IR(1:fftpts/2))),'b');
title('IR before (b)/after (r) alteration')
xlabel('freq (Hz)')
ylabel('magnitude (dB)')
end
if pflag1
% plot phase
figure; hold on;
plot(freq,angle(FT_IRout(1:fftpts/2)),'r');
plot(freq,angle(FT_IR(1:fftpts/2)),'b');
title('IR before/after modification')
xlabel('freq (Hz)')
ylabel('phase angle (radians)')

% plot time domain signal
figure; hold on;
plot(real(IRin))
plot(IRout,'r')
end