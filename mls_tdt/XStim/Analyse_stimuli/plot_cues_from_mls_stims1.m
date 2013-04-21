function [ILD,ITD,ABL,cF,time] = plot_cues_from_mls_stims1(frame, pflag)

% function [ILD,ITD,ABL,cF,time] = plot_cues_from_mls_stims1(frame, pflag)
% frame [start stop] allows choosing where in stimulus to analyze
% pflag controls plotting
% asks for location, then calcs ILD, ITD and ABL and plots each

if nargin < 2
    pflag = 1;
end

[Lnoi, Rnoi] = read_mls_stim;
nPts = length(Lnoi);
if nargin < 1
    frame = [1 nPts];
end
Lnoi = Lnoi(frame(1):frame(2));
Rnoi = Rnoi(frame(1):frame(2));
nPts = length(Lnoi);

Fs = 30000;
% calc cF (center frequencies for filterbank) on 1/12th octave scale
cF = round(1000*exp(([12:40]/12)*log(2)))';
n_cF = length(cF);
fcoefs = Make_ERBFiltA(Fs,cF);

S = zeros(1024*16,1);
S(1024*8) = 1;
temp = ERBFilterBankB(S, fcoefs);
for icF = 1:n_cF
    M(icF) = mom(temp(icF,:),2);
end
Factor = max1(M)./M;
Factormat = repmat(Factor',1,nPts);
clear S tempL M

% gammatone filter each stimulus  
disp('filtering stimuli')
tempL = ERBFilterBankB(Lnoi, fcoefs).* Factormat;		% has dimensions n_cF x length(noi)
tempR = ERBFilterBankB(Rnoi, fcoefs).* Factormat;
disp('calculating binaural cues')
[ITD, time, IPD] = calcitd(tempL,tempR, cF, Fs, ones(size(cF)));
ABL = (20*log10(std(tempR')) + 20*log10(std(tempL')))/2;
ILD = 20*log10(std(tempR')) - 20*log10(std(tempL'));

% plot on new figures
if pflag
    figure; hold on
    plot(cF,ILD,'b')
    title('ILD')
    figure; hold on
    plot(cF,ABL,'b')
    title('ABL')
    figure; hold on
    plot(time,ITD,'b')
    title('ITD')
end
    