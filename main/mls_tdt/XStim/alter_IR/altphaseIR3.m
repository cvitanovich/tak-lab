function IRout = altphaseIR3 (IRin, minF, maxF, Fact)

%function IRout = altphaseIR3 (IRin, minF, maxF, Fact)
% adds noise to phase spectrum of IR within frequency band
% argins:
%		IRin
%		minF
%		maxF
%		Fact (+/- radians to change, ranges from 0 to pi)
%
% in progress to allow multiple frequency bands

if nargin < 4		Fact = pi/2;	end
if nargin < 3		maxF = 6000;	else maxF = sort(maxF(:));    end
if nargin < 2		minF = 4000;	else minF = sort(minF(:));    end

Nbands = min1([size(maxF,1) size(minF,1)]);
pflag = 0;
momflag = 1;      % equalizes by 2nd moment

if nargin < 1 | isempty('IRin')
	% read an HRTF in the time domain for testing purposes only
    FN = 'e:\kip\hrtfdata\904ad.eq';
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
phase = zeros(Npts,1);

for iband = 1:Nbands
    minFind = nearest_index(freq,minF(iband));
    maxFind = nearest_index(freq,maxF(iband));
    temp = (rand(maxFind-minFind+1,1) -.5)*2*Fact;		                % part of pi-change desired
    temp = temp .* tukeywin(length(temp),.1);                       % taper the ends
    phase(minFind:maxFind) = temp;                                        % positive freqs
    phase(Npts+2-maxFind:Npts+2-minFind) = -flipud(temp);                 % negative freqs
end

% combine phase and magnitude spectra:
FT_X = (mag .*cos(phase)) + ((mag .*sin(phase))*sqrt(-1));

% multiply the two in freq domain
FT_IRshort = fft(IRinShort,Npts);
FT_X2 = FT_X .* FT_IRshort;

% convert to time domain:
IRout = zeros(255,1);
IRout(ind) = real(ifft(FT_X2,Npts));
%IRout([1 255]) = [0 0];

if momflag
    IRout = IRout * (mom(IRin(ind),2)/mom(IRout,2));
end

if pflag
% plot freq domain
fftpts = Npts;
freq = 0:15000/(fftpts/2-1):15000;
FT_IRout = fft(IRout,fftpts);
FT_IRin = fft(IRin,fftpts);
figure; hold on
plot(freq,20*log10(abs(FT_IRin(1:round(fftpts/2)-1))))
plot(freq,20*log10(abs(FT_IRout(1:round(fftpts/2)-1))),'r')
warning('off')
plot(freq,20*log10(abs(FT_X2(1:round(fftpts/2)-1))),'k')
warning ('on')
title('blue: original, red: new, black: intermediate before ifft')

figure; hold on
plot(freq,angle(FT_IRin(1:round(fftpts/2)-1)))
plot(freq,angle(FT_IRout(1:round(fftpts/2)-1)),'r')
plot(freq,angle(FT_X2(1:round(fftpts/2)-1)),'k')
title('blue: original, red: new, black: intermediate before ifft')

% plot time domain
figure; hold on;
plot(IRin)
plot(IRout,'r');
title('blue: original, red: new')

end
