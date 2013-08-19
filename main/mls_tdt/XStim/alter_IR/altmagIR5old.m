function [IRout, IRin, ALTout] = altmagIR5 (IRin, octlist, Fact, ALTin)

%function [IRout, IRin, ALTout]  = altmagIR5 (IRin, octlist, Fact, ALTin)
% adds noise to magnitude spectrum of IR within frequency band
% argins:
%		IRin
%       octlist = list of 1/12 octaves to change
%		Fact (+/- dB noise to add)
%       ALTout is the altered mag spectrum between startpt and endpt
%           useful for changing future IRs in same way (supplied as ALTin)
%           ALTout/in MUST BE the length from startpt:endpt
%
% modified version to only change non-zero parts of IR
% allows multiple frequency bands

if nargin < 3		Fact = 10;	    end
if nargin < 2		octlist = [13:25] ; end

octlist = octlist*2;

pflag = 0;

if nargin < 1 | isempty('IRin')         % read an HRTF in the time domain for testing purposes only
    FN = 'e:\kip\hrtfdata\938AD.std.mat';
    dir = 0;
	eval(['load ' FN]);
    IRin = TF1(100,:);
end

rand ('state',sum(100*clock));

% find non-zero parts of IRin
IRinShort = IRin;
Npts = length(IRinShort);

% convert to freq domain
FT_IRshort = fft(IRinShort,Npts);
FT_new = FT_IRshort;

% calc freqs to change
freq = 0:15000/(Npts/2):15000;
cF = round(1000*exp(([23:84]/24)*log(2)))'; n_cF = length(cF);

% get index into freq vector for each cF
for i = 1:n_cF
    ind(i) = nearest_index(freq,cF(i));             % index to all 1/24th octave freqs
end

% get first and last frequencies to change (1/24 octave to each side of requested freqs)
% octlist has been doubled, so is now indexed to 24th octave cF
startpt = ind(octlist(1)-1);
endpt = ind(octlist(end)+1);

% get y values that correspond to these frequencies to change
Xpart = FT_IRshort(startpt:endpt);                      % positive freqs
Xpart2 = FT_IRshort(Npts+2-endpt:Npts+2-startpt);       % negative freqs

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

% just use temp between startpt and endpt
temp = temp(startpt:endpt);
temp = temp * Fact;
ALTout = temp;                              % save alteration for argout

if nargin >= 4                              % use alteration from argin
    temp = ALTin;
end

A = (10 .^(temp/20))';				        % linear change desired
FT_new(startpt:endpt) = Xpart .* A;                          % positive freqs
FT_new(Npts+2-endpt:Npts+2-startpt) = Xpart2 .* flipud(A);        % negative freqs

% calc new IR
IRnew = real(ifft(FT_new,Npts));
FT_new2 = fft(IRnew,Npts);

% calc IRout
IRout = IRnew';
%IRout([1 255]) = [0 0];

if pflag
    figure
    hold on
    FT_IRout = fft(IRout,Npts);
    FT_IRin = fft(IRin,Npts);
    plot([freq fliplr(freq(2:end))],20*log10(abs(FT_IRout)),'r')
    plot([freq fliplr(freq(2:end))],20*log10(abs(FT_IRin)),'b')
    title('blue: IRin    red: IRout ')

    figure
    hold on
    FT_IRout = fft(IRout,Npts);
    FT_IRin = fft(IRin,Npts);
    plot([freq fliplr(freq(2:end))],unwrap(angle(FT_IRout)),'r')
    plot([freq fliplr(freq(2:end))],unwrap(angle(FT_IRin)),'b')
    title('blue: IRin    red: IRout ')

    figure
    hold on
    plot([freq fliplr(freq(2:end))],angle(FT_IRout),'r')
    plot([freq fliplr(freq(2:end))],angle(FT_IRin),'b')
    title('blue: IRin    red: IRout ')
    
   
    figure
    hold on
    plot(IRout,'r')
    plot(IRin)
end