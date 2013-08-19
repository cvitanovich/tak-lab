function [locind, finalspikematrix, weights, IAresponse] = prep4_MLSweights_2src(dataFN1,dataFN2,dataPATH)

%function [locind, finalspikematrix, weights, IAresponse] = prep4_MLSweights_2src(dataFN1,dataFN2,dataPATH)
% reads data mat and provides argins for plot_MLSweights_2src
% also allows input of best location for plot_MLSweights_2src

global H
global FN

if nargin < 3
    dataPATH = input('Enter data PATH   ','s');
    if isempty(dataPATH)    dataPATH = 'e:\kip\datastor\';  end
end
if nargin < 2
    dataFN2 = input('Enter data FN2   ','s');
end

if nargin < 1
    dataFN1 = input('Enter data FN1   ','s');
end


Fs = 30000;

% calc cF (center frequencies for filterbank) on 1/12th octave scale
cF = round(1000*exp(([12:40]/12)*log(2)))';
n_cF = length(cF);
%Factor to adjust loudness of filters
maxFactor = .00003764*cF(n_cF)+.6236;
Factor = maxFactor ./ (.00003764*cF+.6236);

% display FN1 to get best ILA
read_data3_minus_spont(dataFN1,dataPATH);
H.tempdata_fig = gcf;
%%%%%% MAY NEED TO CHANGE FN.ildalone = 'out16ad'; and FN.ILA_path =
%%%%%% 'd:\hrtfdata\904\'; here, for example

%%% check for desired best location
locmax = input('Enter best el and az (as a 2-element vector), if desired   ');
dir = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone]));
ind = find(dir(1,:) == locmax(1) & dir(2,:) == locmax(2));
    
L = mtlrch([FN.ILA_path FN.ildalone],ind*2-1);
R = mtlrch([FN.ILA_path FN.ildalone],ind*2);
tempL = ERBFilterBankA(L, cF, Fs) .* repmat(Factor,1,255);		% has dimensions n_cF x length(L)
tempR = ERBFilterBankA(R, cF, Fs) .* repmat(Factor,1,255);
best_ILA = (20*log10(std(tempR,0,2)) - 20*log10(std(tempL,0,2)))';

%%% begin processing of two-source file
[DATA, params] = read_data3_minus_spont(dataFN2,dataPATH);
H.tempdata2_fig = gcf;

Nlocs = size(params.locations,2);
Nreps = params.numreps;
Ntrials = size(DATA,1)-1;         
Nparams = DATA(1,1);
maxspikes = max1(DATA(:,Nparams+1));
%%%%%% MAY NEED TO CHANGE FN.ildalone = 'out16ad'; and FN.ILA_path =
%%%%%% 'd:\hrtfdata\904\'; here, for example
dir = sph2dbl(mtlrdir([FN.ILA_path FN.ildalone]));

finalspikematrix = zeros(1,Nlocs);
for itrial = 2:Ntrials+1
    rep = DATA(itrial,Nparams);
    ind = find(params.locations(1,:) == DATA(itrial,1) & params.locations(2,:) == DATA(itrial,2));
    ind0 = find((DATA(itrial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(itrial,Nparams+2:Nparams + 1 + maxspikes) > 0));
    ind1 = find(DATA(itrial,Nparams+2:Nparams + 1 + maxspikes) > params.silence_lead);
    finalspikematrix(ind) = finalspikematrix(ind) + length(ind1) - length(ind0);
end

% params.locations contains the 'center-of-array' locations
% to find the actual speaker locations add or subtract the offsets
for iLoc = 1:Nlocs
    locind(iLoc) = find(dir(1,:) == params.locations(1,iLoc) & dir(2,:) == params.locations(2,iLoc));
    locind1(iLoc) = find(dir(1,:) == params.locations(1,iLoc) + params.offset_el/2 & dir(2,:) == params.locations(2,iLoc)+ params.offset_az/2);
    locind2(iLoc) = find(dir(1,:) == params.locations(1,iLoc) - params.offset_el/2 & dir(2,:) == params.locations(2,iLoc)- params.offset_az/2);
end

%%% move on to calcing and plotting weights (...)
[weights, IAresponse] = plot_MLSweights_2src(finalspikematrix, locind, locind1, locind2, best_ILA);

input('<CR> to close figures and continue...');
if exist1('H.tempdata_fig')     close(H.tempdata_fig);      end
if exist1('H.wts_fig')          close(H.wts_fig);           end
if exist1('H.IAresponse_fig')   close(H.IAresponse_fig);    end
H.tempdata_fig = [];
H.wts_fig = [];      
H.IAresponse_fig=[];        