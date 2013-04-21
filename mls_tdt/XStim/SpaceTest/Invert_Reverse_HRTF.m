function [Lout, Rout] = Invert_Reverse_HRTF (Lin, Rin,pflag)

%function [Lout, Rout] = Invert_Reverse_HRTF (Lin, Rin, pflag)
% reverses and inverts L and R magnitude spectra
% with only minimal effects on phase spectra
% based on Roian's test

if nargin < 3 | isempty('pflag')
    pflag = 0;
end

if nargin < 2 | isempty('Lin') | isempty('Rin')
    % read an HRTF in the time domain for testing purposes only
    locnum = 200;
    FN = 'e:\kip\hrtfdata\914ad.eq';
    Lin = mtlrch(FN,locnum*2-1);
    Rin = mtlrch(FN,locnum*2);  
end

Npts = length(Lin);

% convert to freq domain
FT_L = fft(Lin,Npts);
FT_R = fft(Rin,Npts);
FT_Lnew = FT_L;
FT_Rnew = FT_R;

% calc freqs to change (the range for mtlfade is [1700 2000 11000 11300])
freq = 0:15000/(Npts/2 -1):15000;
startpt = nearest_index(freq, 2000);
endpt = nearest_index(freq, 11000);

% get parts of FT to change
Lpart = FT_L(startpt:endpt);                      % positive freqs
Lpart2 = FT_L(Npts+1-endpt:Npts+1-startpt);       % negative freqs
Rpart = FT_R(startpt:endpt);                      % positive freqs
Rpart2 = FT_R(Npts+1-endpt:Npts+1-startpt);       % negative freqs

Lpart_abs = 20*log10(abs(Lpart));
Lpart2_abs = 20*log10(abs(Lpart2));
Rpart_abs = 20*log10(abs(Rpart));
Rpart2_abs = 20*log10(abs(Rpart2));

M =  max([max1(Lpart_abs) max1(Rpart_abs)]);
m = min([min1(Lpart_abs) min1(Rpart_abs)]);
midpt = M + ((m-M)/2);

% left IR phase & magnitude
Lpart_angle = angle(Lpart);
Lpart_abs = 10 .^(-Lpart_abs /20);                      % inverted magnitude

Lpart2_angle = angle(Lpart2);
Lpart2_abs = 10 .^(-Lpart2_abs /20);                     % inverted magnitude

% right IR phase & magnitude
Rpart_angle = angle(Rpart);
Rpart_abs = 10 .^(-Rpart_abs /20);                     % inverted magnitude

Rpart2_angle = angle(Rpart2);
Rpart2_abs = 10 .^(-Rpart2_abs /20);                     % inverted magnitude

% new Left complex values
% positive freqs
r = Rpart_abs .* cos(Lpart_angle);
imag = Rpart_abs .* sin(Lpart_angle);
FT_Lnew(startpt:endpt) = (r + imag .* sqrt(-1));         %/10;
% negative freqs
r = Rpart2_abs .* cos(Lpart2_angle);
imag = Rpart2_abs .* sin(Lpart2_angle);
FT_Lnew(Npts+1-endpt:Npts+1-startpt) = (r + imag .* sqrt(-1));            %/10;

% new Right complex values
% positive freqs
r = Lpart_abs .* cos(Rpart_angle);
imag = Lpart_abs .* sin(Rpart_angle);
FT_Rnew(startpt:endpt) = (r + imag .* sqrt(-1));             %/10;
% negative freqs
r = Lpart2_abs .* cos(Rpart2_angle);
imag = Lpart2_abs .* sin(Rpart2_angle);
FT_Rnew(Npts+1-endpt:Npts+1-startpt) = (r + imag .* sqrt(-1));            %/10;

% fade
if 0
    FT_Lnew(1:128) = fade(FT_Lnew(1:128),4,14,107,128);
    FT_Rnew(1:128) = fade(FT_Rnew(1:128),4,14,107,128);
    FT_Lnew(129:255) = fade(FT_Lnew(129:255),4,21,114,128);
    FT_Rnew(129:255) = fade(FT_Rnew(129:255),4,21,114,128);
    % fade
    FT_L(1:128) = fade(FT_L(1:128),4,14,107,128);
    FT_R(1:128) = fade(FT_R(1:128),4,14,107,128);
    FT_L(129:255) = fade(FT_L(129:255),4,21,114,128);
    FT_R(129:255) = fade(FT_R(129:255),4,21,114,128);
    
end

% new IRs
Lout = real(ifft(FT_Lnew,Npts));
Rout = real(ifft(FT_Rnew,Npts));

% rescale
Fs = 30000;							% sampling frequency
cF = round(1000*exp(([12:40]/12)*log(2)))'; n_cF = length(cF);
fcoefs = Make_ERBFiltA(Fs,cF);
tempLin = ERBFilterBankB(Lin, fcoefs);		% has dimensions n_cF x length(noi)
tempRin = ERBFilterBankB(Rin, fcoefs);
tempLout = ERBFilterBankB(Lout, fcoefs);		% has dimensions n_cF x length(noi)
tempRout = ERBFilterBankB(Rout, fcoefs);
[ILDin, ABLin] = calclevel_time(tempLin,tempRin, cF);
[ILDout, ABLout] = calclevel_time(tempLout,tempRout, cF);
factor = 10 .^((mean(ABLin) -  mean(ABLout)) /20);
Lout = Lout * factor;
Rout = Rout * factor;

% force ends to zero
Lout([1 255]) = [0 0];
Rout([1 255]) = [0 0];

if pflag
    Fs = 30000;							% sampling frequency
    cF = round(1000*exp(([12:40]/12)*log(2)))'; n_cF = length(cF);
    fcoefs = Make_ERBFiltA(Fs,cF);
    
    FT_Rout = fft(Rout,Npts);
    FT_Lout = fft(Lout,Npts);
    
    figure; hold on;
    plot(20*log10(abs(FT_Lout)) + 20*log10(abs(FT_R)),'g')
    plot(20*log10(abs(FT_Lout)),'r')
    plot(20*log10(abs(FT_R)),'b')
    title('MAGNITUDE SPECTRUM: out LEFT (R) + in RIGHT (B) === should be FLAT (G)')

    figure; hold on;
    plot(unwrap(angle(FT_Lout)),'r')
    plot(unwrap(angle(FT_L)),'b')
    title('PHASE SPECTRUM:  out LEFT (R) in RIGHT (B)')
    
    figure; hold on;
    plot(20*log10(abs(FT_Rout)) + 20*log10(abs(FT_L)),'g')
    plot(20*log10(abs(FT_Rout)),'r')
    plot(20*log10(abs(FT_L)),'b')
    title('MAGNITUDE SPECTRUM: out RIGHT (R) + in LEFT (B) === should be FLAT (G)')
    
    figure; hold on;
    plot(unwrap(angle(FT_Rout)),'r')
    plot(unwrap(angle(FT_R)),'b')
    title('PHASE SPECTRUM:  out RIGHT (R) in LEFT (B)')
    
    figure; hold on;
    tempLin = ERBFilterBankB(Lin, fcoefs);
    tempRin = ERBFilterBankB(Rin, fcoefs);
    [ITDin, time] = calcitd(tempLin,tempRin, cF, Fs, ones(size(cF)));
    tempLout = ERBFilterBankB(Lout, fcoefs);
    tempRout = ERBFilterBankB(Rout, fcoefs);
    [ITDout, time] = calcitd(tempLout,tempRout, cF, Fs, ones(size(cF)));
    plot(time, ITDin)
    plot(time, ITDout,'r')
    title('xcorr: OUT (R) IN (B)')
    
    figure; hold on
    ILDin = (20*log10(abs(FT_R)) - 20*log10(abs(FT_L)));
    ILDout = (20*log10(abs(FT_Rout)) - 20*log10(abs(FT_Lout)));
    plot(ILDin)
    plot(ILDout,'r')
    title('FFT ILD: OUT (R) IN (B)')
    
    figure; hold on
    [ILDin, ABLin] = calclevel_time(tempLin,tempRin, cF);
    [ILDout, ABLout] = calclevel_time(tempLout,tempRout, cF);
    plot(cF,ILDin,'b')
    plot(cF,ILDout,'r')
    title('FILTERBANK ILD: OUT (R) IN (B)')
    
    figure; hold on
    plot(cF,ABLin,'b')
    plot(cF,ABLout,'r')
    plot(cF,ones(size(cF))*mean(ABLin))
    plot(cF,ones(size(cF))*mean(ABLout),'r')
    title('FILTERBANK ABL: OUT (R) IN (B)')
    
end