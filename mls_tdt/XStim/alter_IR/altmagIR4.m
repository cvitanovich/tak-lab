function [IRout, IRin] = altmagIR4 (IRin, octlist, Fact)

%function IRout = altmagIR4 (IRin, octlist, Fact)
% adds noise to magnitude spectrum of IR within frequency band
% argins:
%		IRin
%       octlist = list of 1/3
%       octaves to change
%		Fact (+/- dB noise to add)
%
% modified version to only change non-zero parts of IR
% allows multiple frequency bands

if nargin < 3		Fact = 20;	    end
if nargin < 2		octlist = [1 3:6] ; end

octlist = octlist*2;

pflag = 0;

if nargin < 1 | isempty('IRin')
	% read an HRTF in the time domain for testing purposes only
    FN = 'e:\kip\hrtfdata\904ad.eq';
	IRin = mtlrch(FN,100*2-1);
end

rand ('state',sum(100*clock));

% find non-zero parts of IRin
IRinShort = IRin;
Npts = length(IRinShort);

% convert to freq domain
FT_IRshort = fft(IRinShort,Npts);
FT_new = FT_IRshort;

% calc freqs to change
freq = 0:15000/(Npts/2 -1):15000;
cF = round(1000*exp(log(2)*([5:21]/6)));
% ind = [14   18    22    28    35    44    55    68    86   108];
for i = 1:length(cF)
    ind(i) = nearest_index(freq,cF(i));             % index to all 1/6th octave freqs
end

noctlist = size(octlist,2);

% get first and last frequencies to change (1/6 octave to each side of requested freqs)
startpt = ind(octlist(1)-1);
endpt = ind(octlist(noctlist)+1);

% get y values that correspond to these frequencies to change
Xpart = FT_IRshort(startpt:endpt);                      % positive freqs
Xpart2 = FT_IRshort(Npts+2-endpt:Npts+2-startpt);       % negative freqs

% make list of zeros for all freqs and change to "random" those 1/3 octave freqs requested
temp = zeros(size(freq));
temp(ind(octlist)) = rand(size(octlist))*2-1;

% change values of zero-list at 1/6 octave to each side of the requested 1/3 octave freqs
for ioct = 4:2:size(ind,2)      % test every 1/3 octave point
    if (temp(ind(ioct - 2)) ~= 0) & (temp(ind(ioct)) ~=0)
        temp(ind(ioct-1)) = (temp(ind(ioct-2)) + temp(ind(ioct)))/2;
    end
end

% smoothly ramp between all ind-values
dtemp = diff(temp(ind));    % difference between temps (y)
dtemp2 = diff(ind);         % number of pts between indices (x)
for ioct = 1:size(ind,2)-1
    temp(ind(ioct):ind(ioct+1)) = temp(ind(ioct)) + dtemp(ioct) .* ([0:dtemp2(ioct)]/dtemp2(ioct));
end

% just use temp between startpt and endpt
temp = temp(startpt:endpt);

temp = temp * Fact;
A = (10 .^(temp/20))';				        % linear change desired
FT_new(startpt:endpt) = Xpart .* A;                          % positive freqs
FT_new(Npts+2-endpt:Npts+2-startpt) = Xpart2 .* flipud(A);        % negative freqs


% calc new IR
IRnew = real(ifft(FT_new,Npts));
FT_new2 = fft(IRnew,Npts);


% calc IRout
IRout = IRnew;
%IRout([1 255]) = [0 0];

if pflag
    figure
    hold on
    FT_IRout = fft(IRout,Npts);
    FT_IRin = fft(IRin,Npts);
    plot(20*log10(abs(FT_IRout)),'r')
    plot(20*log10(abs(FT_IRin)),'b')
    title('the ENTIRE: blue: IRin    red: IRout ')

    figure
    hold on
    FT_IRout = fft(IRout,Npts);
    FT_IRin = fft(IRin,Npts);
    plot(unwrap(angle(FT_IRout)),'r')
    plot(unwrap(angle(FT_IRin)),'b')
    title('the ENTIRE: blue: IRin    red: IRout ')

    figure
    hold on
    plot(angle(FT_IRout),'r')
    plot(angle(FT_IRin),'b')
    title('the ENTIRE: blue: IRin    red: IRout ')
    
   
    figure
    hold on
    plot(IRout,'r')
    plot(IRin)
end