function IRout = altmagIR3 (IRin, minF, maxF, Fact)

%function IRout = altmagIR3 (IRin, minF, maxF, Fact)
% adds noise to magnitude spectrum of IR within frequency band
% argins:
%		IRin
%		minF
%		maxF
%		Fact (+/- dB noise to add)
%
% modified version to only change non-zero parts of IR
% allows multiple frequency bands

if nargin < 4		Fact = 0;	    end
if nargin < 3		maxF = 8000;	else maxF = sort(maxF(:));    end
if nargin < 2		minF = 2000;	else minF = sort(minF(:));    end

Nbands = min1([size(maxF,1) size(minF,1)]);
pflag = 0;
momflag = 0;      % equalizes by 2nd moment

if nargin < 1 | isempty('IRin')
	% read an HRTF in the time domain for testing purposes only
    FN = 'e:\kip\hrtfdata\904ad.eq';
	IRin = mtlrch(FN,200*2-1);
end

rand ('state',sum(100*clock));

% find non-zero parts of IRin
ind = find(IRin ~=0);
ind = 1:255;
IRinShort = IRin(ind);
Npts = length(IRinShort);

% convert to freq domain
FT_IRshort = fft(IRinShort,Npts);
FT_new = FT_IRshort;

% calc freqs to change
freq = 0:15000/(Npts/2 -1):15000;

for iband = 1:Nbands
    minFind = nearest_index(freq,minF(iband));
    maxFind = nearest_index(freq,maxF(iband));

    Xpart = FT_IRshort(minFind:maxFind);
    Xpart2 = FT_IRshort(Npts+2-maxFind:Npts+2-minFind);
    temp = (rand(maxFind-minFind+1,1) *2 - 1) ;
    temp = temp .* tukeywin(length(temp),.1) *Fact;		% dB of change desired
    A = (10 .^(temp/20));				        % linear change desired
    FT_new(minFind:maxFind) = Xpart .* A;                          % positive freqs
    FT_new(Npts+2-maxFind:Npts+2-minFind) = Xpart2 .* flipud(A);        % negative freqs
end


% calc new IR
IRnew = real(ifft(FT_new,Npts));
FT_new2 = fft(IRnew,Npts);

if momflag
    IRnew = IRnew * (mom(IRin(ind),2)/mom(IRnew,2));
end

% calc IRout
IRout = IRin;
IRout(ind) = IRnew;
IRout([1 255]) = [0 0];

if pflag
    figure
    hold on
    FT_IRout = fft(IRout,Npts);
    FT_IRin = fft(IRin,Npts);
    plot(20*log10(abs(FT_IRout)),'r')
    plot(20*log10(abs(FT_IRin)),'b')
    title('the ENTIRE: blue: IRin    red: IRout vs non-zero portion (cyan and magenta)')
    FT_IRout = fft(IRout(ind),Npts);
    FT_IRin = fft(IRin(ind),Npts);
    plot(20*log10(abs(FT_IRout)),'m')
    plot(20*log10(abs(FT_IRin)),'c')
    plot(minFind:maxFind,20*log10(abs(FT_IRin(minFind:maxFind))) + temp,'m*')
    plot(Npts+2-maxFind:Npts+2-minFind,20*log10(abs(FT_IRin(Npts+2-maxFind:Npts+2-minFind))) + flipud(temp),'m*')

    figure
    hold on
    FT_IRout = fft(IRout,Npts);
    FT_IRin = fft(IRin,Npts);
    plot(unwrap(angle(FT_IRout)),'r')
    plot(unwrap(angle(FT_IRin)),'b')
    title('the ENTIRE: blue: IRin    red: IRout ')

    
    figure
    hold on
    FT_IRout = fft(IRout(ind),Npts);
    FT_IRin = fft(IRin(ind),Npts);
    plot(20*log10(abs(FT_IRout)),'r')
    plot(20*log10(abs(FT_IRin)),'b')
    plot(minFind:maxFind,20*log10(abs(FT_IRin(minFind:maxFind))) + temp,'r*')
    plot(Npts+2-maxFind:Npts+2-minFind,20*log10(abs(FT_IRin(Npts+2-maxFind:Npts+2-minFind))) + flipud(temp),'r*')
    title('the non-zero portions of: blue: IRin         red: IRout')
    
    figure
    hold on
    plot(20*log10(abs(FT_new)),'k')
    plot(20*log10(abs(FT_new2)),'r')
    plot(20*log10(abs(FT_IRshort)))
    plot(minFind:maxFind,20*log10(abs(FT_IRshort(minFind:maxFind))) + temp,'r*')
    plot(Npts+2-maxFind:Npts+2-minFind,20*log10(abs(FT_IRshort(Npts+2-maxFind:Npts+2-minFind))) + flipud(temp),'r*')
    title('black: FT new       blue: FT IRshort           red: FT new2 (after ifft & fft of FT new)') 

    figure
    hold on
    plot(IRout,'r')
    plot(IRin)
end