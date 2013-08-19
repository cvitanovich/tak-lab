function [weights, IAresponse] = plot_MLSweights2(finalspikematrix, locind, locmax)

% function [weights, IAresponse] = plot_MLSweights2(finalspikematrix,locind, locmax)
% argins available from Engage_Space
% finalspikematrix is arranged as per locind
% locmax is a two-element vector (el az) for the best location (default is to find its own locmaxspikes)
% VERSION to use .mat-style HRTF files

global FN
global H

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
Factormat = repmat(Factor',1,255);
clear S tempL M

%%%%%%% read IA HRTFs & calc ILDs organized as locind
IAnLocs = size(locind,2);
dir = 0;
eval(['load ' FN.ILA_path FN.ildalone ' -mat']);
IAHRTF_dir = dir(:,locind);
ILD = zeros(n_cF,IAnLocs);
disp('calculating ILDs for each location');
for loc = 1:IAnLocs
    tempL = ERBFilterBankB(TF1(locind(loc),:), fcoefs).* Factormat;		% has dimensions n_cF x length(noi)
    tempR = ERBFilterBankB(TF2(locind(loc),:), fcoefs).* Factormat;
    ILD(:,loc) = (20*log10(std(tempR')) - 20*log10(std(tempL')))';	% returns level as function of frequency
end

% find index to maxspikes
if nargin <3
    locmaxspikes = min1(find(finalspikematrix == max(finalspikematrix)));
else
    locmaxspikes = find(IAHRTF_dir(1,:) == locmax(1) & IAHRTF_dir(2,:) == locmax(2));
end

%%%%%%% replot IA response as freq x ILD surface
disp('calculating IA response surface as freq x ILD');
IAresponse = zeros(n_cF,61);
for ifreq = 1:n_cF
    for iLoc = 1:IAnLocs
        bin = find(histc(ILD(ifreq,iLoc),-30:30) ==1);
        IAresponse(ifreq,bin) = IAresponse(ifreq,bin) + finalspikematrix(iLoc);
    end
end
H.IAresponse_fig = figure;
pcolor(cF,-30:30,IAresponse')
shading flat
title('IA response surface')
xlabel('frequency')
ylabel('ILD')
hold on
plot(cF,ILD(:,locmaxspikes),'w','Linewidth',4);
plot(cF,ILD(:,locmaxspikes),'k','Linewidth',2);

%%%%%%% get ILD distance (organized as locind)
dILD = get_ILDdist2(ILD, cF, ILD(:,locmaxspikes),[cF cF],1);

% get activity vs dILD curves for each band
disp('calculating activity curves for each frequency')
warning('off')
for ib = 1:n_cF
   [temp_dILD, ind_dILD(ib,:)] = sort(dILD(ib,:));
   temp_spikes = finalspikematrix(ind_dILD(ib,:));
   [activity(ib,ind_dILD(ib,:))] = locregress1(temp_dILD, temp_spikes', 0.05);
end 
warning('on')

% get weights using Michael's PC-regression
disp('calculating weights')
weights = get_dILDweights(dILD, finalspikematrix, cF, IAHRTF_dir, 1);

H.wts_fig = figure;
plot(cF,weights,'.','MarkerSize',5)
title('PCregression analysis of dILD')
xlabel('frequency')
ylabel('weight')