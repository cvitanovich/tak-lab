function IRout = altphaseIR2 (IRin, minF, maxF, Fact)

%function IRout = altphaseIR2 (IRin, minF, maxF, Fact)
% adds noise to phase spectrum of IR within frequency band
% argins:
%		IRin
%		minF
%		maxF
%		Fact (+/- radians to change, ranges from 0 to pi)
%

if nargin < 4		Fact = pi/2;	    end
if nargin < 3		maxF = 6000;	end
if nargin < 2		minF = 4000;	end

pflag = 0;
momflag = 1;      % equalizes by 2nd moment

if nargin < 1 | isempty('IRin')
	% read an HRTF in the time domain for testing purposes only
    FN = 'e:\kip\hrtfdata\839ad.eq';
	IRin = mtlrch(FN,200*2-1);
end

rand ('state',sum(100*clock));

% find the part of IRin to change
ind = find(IRin ~=0);
ind = 1:255;                        % for now, change them all
IRinShort = IRin(ind);
Npts = length(ind);
freq = 0:15000/(Npts/2 -1):15000;

% mag spectrum between broadly defined freqs
minF1 = nearest_index(freq,1000);
maxF1 = nearest_index(freq,12000);
mag = zeros(Npts,1);
mag(minF1:maxF1) = ones(maxF1-minF1+1,1);                       % positive freqs
mag(Npts+2-maxF1:Npts+2-minF1) = ones(maxF1-minF1+1,1);         % negative freqs

% phase is randomized between narrowly defined freqs
minF = nearest_index(freq,minF);
maxF = nearest_index(freq,maxF);
phase = zeros(Npts,1);
temp = (rand(maxF-minF+1,1) -.5)*2*Fact;		                % part of pi-change desired
temp = temp .* tukeywin(length(temp),.1);                       % taper the ends
phase(minF:maxF) = temp;                                        % positive freqs
phase(Npts+2-maxF:Npts+2-minF) = -flipud(temp);                 % negative freqs

% combine phase and magnitude spectra:
FT_X = (mag .*cos(phase)) + ((mag .*sin(phase))*sqrt(-1));

% multiply the two in freq domain
FT_IRshort = fft(IRinShort,Npts);
FT_X2 = FT_X .* FT_IRshort;

% convert to time domain:
IRout = zeros(255,1);
IRout(ind) = real(ifft(FT_X2,Npts));
IRout([1 255]) = [0 0];

if momflag
    IRout = IRout * (mom(IRin(ind),2)/mom(IRout,2));
end

if pflag
% plot freq domain
fftpts = 1024;
FT_IRout = fft(IRout,Npts);
FT_IRin = fft(IRin,Npts);
figure; hold on
plot(20*log10(abs(FT_IRin)))
plot(20*log10(abs(FT_IRout)),'r')
warning('off')
plot(20*log10(abs(FT_X2)),'k')
warning ('on')

figure; hold on
plot(angle(FT_IRin))
plot(angle(FT_IRout),'r')
plot(angle(FT_X2),'k')

% plot time domain
figure; hold on;
plot(IRin)
plot(IRout,'r');
end
