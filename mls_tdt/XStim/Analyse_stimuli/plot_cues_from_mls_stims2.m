function [ILD,ITD,ABL,cF,time,dIPD] = plot_cues_from_mls_stims2(frame, win, pflag, best_IPD, ind_best_IPD)

% function [ILD,ITD,ABL,cF,time,dIPD] = plot_cues_from_mls_stims2(frame, win, pflag)
% frame [start stop] allows choosing where in stimulus to analyze 
% (to avoid zero-padding use 3001:6000)
% pflag controls plotting
% asks for location, then calcs ILD, ITD and ABL and plots each
% VERSION to plot time-varying cues

if nargin < 4
    best_IPD = [];
    ind_best_IPD = [];
end
if nargin < 3 | isempty(pflag)
    pflag = 1;
end

if nargin < 2 | isempty(win)
    win = 30;
end

[Lnoi, Rnoi] = read_mls_stim;
nPts = length(Lnoi);
if nargin < 1 | isempty(frame)
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
Factormat = repmat(Factor',1,win);
clear S tempL M

nWin = floor((diff(frame)+1)/win);
for iwin = 1:nWin
    % gammatone filter each stimulus  
    disp(['filtering stimuli win# ' num2str(iwin) ' of ' num2str(nWin)])
    tempL = ERBFilterBankB(Lnoi((iwin-1)*win+1:iwin*win), fcoefs).* Factormat;		% has dimensions n_cF x length(noi)
    tempR = ERBFilterBankB(Rnoi((iwin-1)*win+1:iwin*win), fcoefs).* Factormat;
    
    [ITD(iwin,:), time, IPD] = calcitd(tempL,tempR, cF, Fs, ones(size(cF)));
    ABL(iwin,:) = (20*log10(std(tempR')) + 20*log10(std(tempL')))/2;
    ILD(iwin,:) = 20*log10(std(tempR')) - 20*log10(std(tempL'));
    if ~isempty(best_IPD) 
        for icF = 1:n_cF
            dIPD(iwin,icF) = (1+best_IPD(icF)) - (1+IPD(icF,ind_best_IPD(icF)));
        end
    end
end

% plot on new figures
if pflag
    figure
    pcolor((1:nWin),cF,ILD')
    shading flat
    title('ILD')
    figure; hold on
    plot((1:nWin),mean(ABL,2),'b')
    title('ABL')
    figure
    pcolor((1:nWin),time,ITD')
    shading flat
    title('ITD')
    if ~isempty(best_IPD) 
        figure
        pcolor((1:nWin),cF,dIPD')
        shading flat
        title('dIPD')
    end
end
