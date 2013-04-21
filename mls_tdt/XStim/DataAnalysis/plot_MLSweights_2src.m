function [weights, IAresponse] = plot_MLSweights_2src(finalspikematrix, locind, locind1, locind2, best_ILA)

% function [weights, IAresponse] = plot_MLSweights_2src(finalspikematrix,locind, locind1, locind2, best_ILA)
% argins available from Engage_Space
% finalspikematrix is arranged as per locind
% locmax is a two-element vector (el az) for the best location (default is to find its own locmaxspikes)

global FN
global H

Fs = 30000;

% calc cF (center frequencies for filterbank) on 1/12th octave scale
cF = round(1000*exp(([12:40]/12)*log(2)))';
n_cF = length(cF);
%Factor to adjust loudness of filters
maxFactor = .00003764*cF(n_cF)+.6236;
Factor = maxFactor ./ (.00003764*cF+.6236);

%%%%%%% read IA HRTFs & calc ILDs organized as locind
IAnLocs = size(locind,2);
IAHRTF_dir = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone]));
IAHRTF_dir1 = IAHRTF_dir(:,locind1);
IAHRTF_dir2 = IAHRTF_dir(:,locind2);
IAHRTF_dir = IAHRTF_dir(:,locind);
ILD = zeros(n_cF,IAnLocs);
disp('calculating ILDs for each location');
for loc = 1:IAnLocs
    L = mtlrch([FN.ILA_path FN.ildalone],locind1(loc)*2-1) + mtlrch([FN.ILA_path FN.ildalone],locind2(loc)*2-1);
    R = mtlrch([FN.ILA_path FN.ildalone],locind1(loc)*2) + mtlrch([FN.ILA_path FN.ildalone],locind2(loc)*2);
    tempL = ERBFilterBankA(L, cF, Fs) .* repmat(Factor,1,255);		% has dimensions n_cF x length(noi)
    tempR = ERBFilterBankA(R, cF, Fs) .* repmat(Factor,1,255);
    ILD(:,loc) = calclevel_time(tempL, tempR, cF)';	% returns level as function of frequency
end

% find index to maxspikes
    locmaxspikes = min1(find(finalspikematrix == max(finalspikematrix)));

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
plot(cF,best_ILA,'w','Linewidth',4);
plot(cF,best_ILA,'k','Linewidth',2);

%%%%%%% get ILD distance (organized as locind)
dILD = get_ILDdist2(ILD, cF, best_ILA,[cF cF],1);

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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ILD, ABI] = calclevel_time(tempL, tempR, cF);
% calculate levels in the time domain
n_cF = length(cF);

for i = 1:n_cF
   L = tempL(i,:);
   R = tempR(i,:);
   ABI(i) = (20*log10(std(R)) + 20*log10(std(L)))/2;
   ILD(i) = 20*log10(std(R)) - 20*log10(std(L));
end
