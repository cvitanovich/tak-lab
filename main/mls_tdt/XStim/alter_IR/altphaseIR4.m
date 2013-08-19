function IRout = altphaseIR4 (IRin, octlist, Fact)

%function IRout = altphaseIR4 (IRin, minF, maxF, Fact)
% adds noise to phase spectrum of IR within frequency band
% argins:
%		IRin
%		octlist = list of 1/3 octaves to change
%		Fact (+/- radians to change, ranges from 0 to pi)
%

if nargin < 3		Fact = pi/2;	end
if nargin < 2       octlist = [1:6];   end 

octlist = octlist * 2;

pflag = 0;

if nargin < 1 | isempty('IRin')
	% read an HRTF in the time domain for testing purposes only
    FN = 'e:\kip\hrtfdata\914AD.std';
	IRin = mtlrch(FN,200*2-1);
end

rand ('state',sum(100*clock));

% find the part of IRin to change
IRinShort = IRin;
Npts = length(IRinShort);
freq = 0:15000/(Npts/2 -1):15000;

% convert to freq domain
%FT_IRshort = fft(IRinShort,Npts);
%FT_new = FT_IRshort;

% mag spectrum between broadly defined freqs
minF1 = nearest_index(freq,1000);
maxF1 = nearest_index(freq,12000);
mag = zeros(Npts,1);
mag(minF1:maxF1) = ones(maxF1-minF1+1,1);                       % positive freqs
mag(Npts+2-maxF1:Npts+2-minF1) = ones(maxF1-minF1+1,1);         % negative freqs

% frequencies to change
cF = round(1000*exp(log(2)*([5:21]/6)));            % finds the 1/6 octaves
for i = 1:length(cF)
    ind(i) = nearest_index(freq, cF(i));            % creates an index of the x-values of freq that are close to the values of cF
end

noctlist = size(octlist,2);

% first frequencies to change (1/6 octave on either side of first and last 1/3 octave)
startpt = ind(octlist(1)-1);
endpt = ind(octlist(noctlist)+1);

%Xpart = FT_IRshort(startpt:endpt);
%Xpart2 = FT_IRshort(Npts+2-endpt:Npts+2-startpt);


% to change phase for the defined 1/3octaves
phase = zeros(Npts,1);
phase(ind(octlist)) = (rand(size(octlist))*2-1) * Fact;

for ioct = 4:2:size(ind,2)      % test every 1/3 octave point
    if (phase(ind(ioct - 2)) ~= 0) & (phase(ind(ioct)) ~=0)
        phase(ind(ioct-1)) = (phase(ind(ioct-2)) + phase(ind(ioct)))/2;
    end
end

% smoothly ramp between all ind-values
dphase = diff(phase(ind));    % difference between temps (y)
dphase2 = diff(ind);         % number of pts between indices (x)
for ioct = 1:size(ind,2)-1
    phase(ind(ioct):ind(ioct+1)) = phase(ind(ioct)) + dphase(ioct) .* ([0:dphase2(ioct)]/dphase2(ioct));
end

phase((((Npts+1)/2)+1):Npts) = -flipud(phase(2:((Npts+1)/2)));

% combine phase and magnitude spectra:
FT_X = (mag .*cos(phase)) + ((mag .*sin(phase))*sqrt(-1));

% multiply the two in freq domain
FT_IRshort = fft(IRinShort,Npts);
FT_X2 = FT_X .* FT_IRshort;

% convert to time domain:
IRout = real(ifft(FT_X2,Npts));
%IRout([1 255]) = [0 0];

if pflag
% plot freq domain
fftpts = Npts;
freq = 0:15000/(fftpts/2-1):15000;
FT_IRout = fft(IRout,fftpts);
FT_IRin = fft(IRin,fftpts);
figure; hold on
plot(freq,20*log10(abs(FT_IRin(1:round(fftpts/2)-1))),'linewidth',2)
plot(freq,20*log10(abs(FT_IRout(1:round(fftpts/2)-1)))+5,'r.','linewidth',1)
warning('off')
%plot(freq,20*log10(abs(FT_X2(1:round(fftpts/2)-1))),'k')
warning ('on')
%title('blue: original, red: new, black: intermediate before ifft')

figure; hold on
plot(freq,angle(FT_IRin(1:round(fftpts/2)-1)),'linewidth',2)
plot(freq,angle(FT_IRout(1:round(fftpts/2)-1)),'r--','linewidth',1)
%plot(freq,angle(FT_X2(1:round(fftpts/2)-1)),'k')
%title('blue: original, red: new, black: intermediate before ifft')

figure; hold on
plot(freq,unwrap(angle(FT_IRin(1:round(fftpts/2)-1))))
plot(freq,unwrap(angle(FT_IRout(1:round(fftpts/2)-1)),'r'),'linewidth',3)
plot(freq,unwrap(angle(FT_X2(1:round(fftpts/2)-1)),'k'),'linewidth',3)
title('blue: original, red: new, black: intermediate before ifft')

% plot time domain
figure; hold on;
plot(IRin,'linewidth',2)
plot(IRout,'r','linewidth',2);
%title('blue: original, red: new')

end
