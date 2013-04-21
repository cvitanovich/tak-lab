function IRout = altmagIR2 (IRin, minF, maxF, Fact)

%function IRout = altmagIR2 (IRin, minF, maxF, Fact)
% adds noise to magnitude spectrum of IR within frequency band
% argins:
%		IRin
%		minF
%		maxF
%		Fact (+/- dB noise to add)
%
% modified version to only change non-zero parts of IR

if nargin < 4		Fact = 5;	    end
if nargin < 3		maxF = 6000;	end
if nargin < 2		minF = 2000;	end

pflag = 0;
momflag = 1;      % equalizes by 2nd moment


if nargin < 1 | isempty('IRin')
	% read an HRTF in the time domain for testing purposes only
    FN = 'e:\kip\hrtfdata\839ad.eq';
	IRin = mtlrch(FN,200*2-1);
end

rand ('state',sum(100*clock));

ind = find(IRin ~=0);
%ind = 1:255;
IRinShort = IRin(ind);


Npts = length(IRinShort);

FT_IRshort = fft(IRinShort,Npts);
FT_new = FT_IRshort;

% calc freqs to change
freq = 0:15000/(Npts/2 -1):15000;
minF = nearest_index(freq,minF);
maxF = nearest_index(freq,maxF);

Xpart = FT_IRshort(minF:maxF);
Xpart2 = FT_IRshort(Npts+2-maxF:Npts+2-minF);
temp = (rand(maxF-minF+1,1) *2 - 1) ;
temp = temp .* tukeywin(length(temp),.1) *Fact;		% dB of change desired
temp = (10 .^ (temp /20));				        % linear change desired

FT_new(minF:maxF) = Xpart .* temp;                          % positive freqs
FT_new(Npts+2-maxF:Npts+2-minF) = Xpart2 .* flipud(temp);    % negative freqs

% calc new IR
IRnew = real(ifft(FT_new,Npts));
FT_new2 = fft(IRnew,Npts);

if momflag
    IRnew = IRnew * (mom(IRin(ind),2)/mom(IRnew,2));
end

% calc IRout
IRout = IRin;
IRout(ind) = IRnew;

if pflag
FT_IRout = fft(IRout,1024*2);
FT_IRin = fft(IRin,1024*2);
figure
hold on
plot(20*log10(abs(FT_IRout)),'k')
plot(20*log10(abs(FT_IRin)),'r')

figure
hold on
plot(20*log10(abs(FT_new)),'k')
plot(20*log10(abs(FT_new2)),'r')
plot(20*log10(abs(FT_IRshort)))

figure
hold on
plot(IRout,'r')
plot(IRin)
end