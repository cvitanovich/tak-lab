function [IRout, IRin, ALTout] = altphaseIR5 (IRin, octlist, Fact, ALTin)

%function [IRout, IRin, ALTout] = altphaseIR5 (IRin, octlist, Fact, ALTin)
% adds noise to phase spectrum of IR within frequency band
% argins:
%		IRin
%		octlist = list of 1/12 octaves to change
%		Fact (+/- radians to change, ranges from 0 to pi)
%
%       ALTout is the altered mag spectrum between startpt and endpt
%           useful for changing future IRs in same way (supplied as ALTin)
%           ALTout/in MUST BE the length of freq (Npts/2)

if nargin < 3		Fact = pi/2;	end
if nargin < 2       octlist = [13:25];   end 

octlist = octlist * 2;

pflag = 0;

if nargin < 1 | isempty('IRin')
    % read an HRTF in the time domain for testing purposes only
    FN = 'e:\kip\hrtfdata\938AD.std.mat';
    dir = 0;
    eval(['load ' FN]);
    IRin = TF1(100,:)';
end
IRin = IRin(:);

rand ('state',sum(100*clock));

% find the part of IRin to change
IRinShort = IRin;
Npts = length(IRinShort);
freq = 0:15000/(Npts/2):15000;
cF = round(1000*exp(([23:84]/24)*log(2)))'; 
n_cF = length(cF);

% mag spectrum between broadly defined freqs
minF1 = nearest_index(freq,1000);
maxF1 = nearest_index(freq,12000);
%mag = zeros(Npts,1);
%mag(minF1:maxF1) = ones(maxF1-minF1+1,1);                       % positive freqs
%mag(Npts+2-maxF1:Npts+2-minF1) = ones(maxF1-minF1+1,1);         % negative freqs
mag = ones(Npts,1);

% frequencies to change
for i = 1:n_cF
    ind(i) = nearest_index(freq, cF(i));            % creates an index of the x-values of freq that are close to the values of cF
end

% first frequencies to change (1/24 octave on either side of first and last 1/3 octave)
% octlist has been doubled, so is now indexed to 24th octave cF
startpt = ind(octlist(1)-1);
endpt = ind(octlist(end)+1);

%Xpart = FT_IRshort(startpt:endpt);
%Xpart2 = FT_IRshort(Npts+2-endpt:Npts+2-startpt);

% make list of zeros for all freqs and change to "random" those 1/12 octave freqs requested
temp = zeros(size(freq));
temp(ind(octlist)) = rand(size(octlist))*2-1;

% change values of zero-list at 1/24 octave to each side of the requested 1/12 octave freqs
% first and last
pt0 = ind(octlist(1)); pt1 = ind(octlist(end));
temp(pt0-1) = temp(pt0)/2;
temp(pt1+1) = temp(pt1)/2;
% all the others
ind0 = find(temp(pt0:pt1)==0)+pt0-1;
ind0 = [pt0 ind0(find(diff(ind0)>1))+1 pt1];
for i = 1:length(ind0)-1
    temp(ind0(i)+1:ind0(i+1)-1) = mean(temp(ind0(i:i+1)));
end
temp((((Npts+1)/2)+1):Npts) = -flipud(temp(2:((Npts+1)/2)));
temp = temp' * Fact;
ALTout = temp;                              % save alteration for argout
if nargin >= 4                              % use alteration from argin
    temp = ALTin;
end

% combine phase and magnitude spectra:
FT_X = (mag .*cos(temp)) + ((mag .*sin(temp))*sqrt(-1));

% multiply the two in freq domain
FT_IRshort = fft(IRinShort,Npts);
FT_X2 = FT_X .* FT_IRshort;

% convert to time domain:
IRout = real(ifft(FT_X2,Npts));

% relics of trying longer HRIRs
IRout = IRout(1:255)';
IRin = IRin(1:255);
Npts = 255;
freq = 0:15000/(Npts/2):15000;

if pflag
    FT_IRout = fft(IRout,Npts);
    FT_IRin = fft(IRin,Npts);
    figure; hold on
    plot(fftshift([freq -fliplr(freq(2:end))]),fftshift(20*log10(abs(FT_IRin))),'b')
    plot(fftshift([freq -fliplr(freq(2:end))]),fftshift(20*log10(abs(FT_IRout))),'r')
    set(gca,'xlim',[0 15000])
    xlabel('Frequency (Hz)')
    ylabel('dB magnitude')
    title('blue: IRin, red: IRout')
    
    figure
    hold on
    FT_IRout = fft(IRout,Npts);
    FT_IRin = fft(IRin,Npts);
    plot(fftshift([freq -fliplr(freq(2:end))]), fftshift(unwrap(angle(FT_IRout))),'r')
    plot(fftshift([freq -fliplr(freq(2:end))]), fftshift(unwrap(angle(FT_IRin))),'b')
    set(gca,'xlim',[0 15000])
    xlabel('Frequency (Hz)')
    ylabel('unwrapped phase')
    title('blue: IRin    red: IRout ')
    
    figure
    hold on
    plot(fftshift([freq -fliplr(freq(2:end))]),fftshift(angle(FT_IRout)),'r')
    plot(fftshift([freq -fliplr(freq(2:end))]),fftshift(angle(FT_IRin)),'b')
    set(gca,'xlim',[0 15000])
    xlabel('Frequency (Hz)')
    ylabel('wrapped phase')
    title('blue: IRin    red: IRout ')
    
    figure
    hold on
    plot(IRout,'r')
    plot(IRin)
    title('blue: IRin    red: IRout ')
    
end
