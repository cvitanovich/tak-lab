function [ILD,ITD,ABL,cF,time,dIPD,dILD,Env_dIPD,Env_dILD,Env_ABL,FT_Env_dIPD,FT_Env_dILD,FT_Env_ABL] = ...
    plot_cues_from_mls_stims3(frame, win, pflag, best_IPD, ind_best_IPD, best_ILD,el,az,stim_dir)

% function [ILD,ITD,ABL,cF,time,dIPD,dILD,Env_dIPD,Env_dILD,Env_ABL,FT_Env_dIPD,FT_Env_dILD,FT_Env_ABL] = ...
%   plot_cues_from_mls_stims3(frame, win, pflag, best_IPD, ind_best_IPD, best_ILD,el,az)
%
% frame [start stop] allows choosing where in stimulus to analyze 
% (to avoid zero-padding use 3001:6000)
% pflag controls plotting
% asks for location, then calcs ILD, ITD and ABL and plots each
% VERSION to plot time-varying cues
% to calc FFTs of envelope (etc)


if nargin < 9 | isempty(stim_dir)
    stim_dir = 'e:\kip\temp\ITA2\';
end
if nargin < 8 | isempty(az)
    az = [];
end
if nargin < 7 | isempty(el)
    el = [];
end
if nargin < 6 | isempty(best_ILD)
    best_ILD = [];
end
if nargin < 4 | isempty(best_IPD) | isempty(ind_best_IPD)
    best_IPD = [];
    ind_best_IPD = [];
end
if nargin < 3 | isempty(pflag)
    pflag = 1;
end
if nargin < 2 | isempty(win)
    win = 30;
end

fftpts = 512;

[Lnoi, Rnoi] = read_mls_stim(el,az, stim_dir);
if nargin < 1 | isempty(frame)
    nPts = length(Lnoi);
    if nPts >= 9300
        frame = [3301 9300];
    else
        frame = [3301 nPts-6000];
    end
end
Lnoi = Lnoi(frame(1):frame(2));
Rnoi = Rnoi(frame(1):frame(2));
nPts = length(Lnoi);

Fs = 30000;
% calc cF (center frequencies for filterbank) on 1/12th octave scale
cF = round(1000*exp(([12:40]/12)*log(2)))';
n_cF = length(cF);
% calc ERBfilters
[fcoefs,Factor] = calc_ERBfilters;
Factormat = repmat(Factor',1,win);

nWin = floor((diff(frame)+1)/win);

% allocate memory
ABL = zeros(nWin,n_cF);
ILD = zeros(nWin,n_cF);
ITD = zeros(nWin,511);
dIPD = zeros(nWin,n_cF);
dILD = zeros(nWin,n_cF);

warning off
for iwin = 1:nWin
    % gammatone filter each stimulus  
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
    if ~isempty(best_ILD) 
        dILD(iwin,:) = get_ILDdist2(ILD(iwin,:)', cF, best_ILD,[cF cF],1)';
    end
end
warning on

% calc envelopes and their FFTs for each frequency band
Env_dIPD = zeros(nWin,n_cF);
FT_Env_dIPD = zeros(fftpts/2,n_cF);
Env_dILD = zeros(nWin,n_cF);
FT_Env_dILD = zeros(fftpts/2,n_cF);
Env_ABL = zeros(nWin,n_cF);
FT_Env_ABL = zeros(fftpts/2,n_cF);

for icF = 1:n_cF
    ind = find(~isnan(dIPD(:,icF)) & ~isnan(dILD(:,icF)) & ~isnan(ABL(:,icF)));
    if ~isempty(best_IPD) 
        Env_dIPD(ind,icF) = abs(hilbert(dIPD(ind,icF)));
        temp = fft(Env_dIPD(:,icF),fftpts);
        FT_Env_dIPD(:,icF) = temp(1:fftpts/2);
    end
    if ~isempty(best_ILD) 
        Env_dILD(ind,icF) = abs(hilbert(dILD(ind,icF)));
        temp = fft(Env_dILD(:,icF),fftpts);
        FT_Env_dILD(:,icF) = temp(1:fftpts/2);
    end
    Env_ABL(ind,icF) = abs(hilbert(ABL(ind,icF)));
    temp = fft(Env_ABL(:,icF),fftpts);
    FT_Env_ABL(:,icF) = temp(1:fftpts/2);
end    


% plot cues on new figures
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
