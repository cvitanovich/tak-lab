function [locind, finalspikematrix, weights, IAresponse] = prep4_MLSweights2(dataFN,dataPATH)

%function [locind, finalspikematrix, weights, IAresponse] = prep4_MLSweights2(dataFN,dataPATH)
% reads data mat and provides argins for plot_MLSweights
% also allows input of best location for plot_MLSweights
% VERSION to use .mat style HRTF files

global H
global FN

if nargin < 2
    dataPATH = input('Enter data PATH   ','s');
    if isempty(dataPATH)    dataPATH = 'e:\kip\datastor\';  end
end

if nargin < 1
    dataFN = input('Enter data FN   ','s');
end

[DATA, params, FN] = read_data3_minus_spont(dataFN,dataPATH);
H.tempdata_fig = gcf;

Nlocs = size(params.locations,2);
Nreps = params.numreps;
Ntrials = size(DATA,1)-1;         
Nparams = DATA(1,1);
maxspikes = max1(DATA(:,Nparams+1));
%%%%%% MAY NEED TO CHANGE 
%FN.ildalone = 'out16ad'; and FN.ILA_path =
%%%%%% 'd:\hrtfdata\904\'; here, for example
dir = 0;
eval(['load ' FN.ILA_path FN.ildalone ' -mat']);

finalspikematrix = zeros(1,Nlocs);
for itrial = 2:Ntrials+1
    rep = DATA(itrial,Nparams);
    ind = find(params.locations(1,:) == DATA(itrial,1) & params.locations(2,:) == DATA(itrial,2));
    ind0 = find((DATA(itrial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(itrial,Nparams+2:Nparams + 1 + maxspikes) > 0));
    ind1 = find(DATA(itrial,Nparams+2:Nparams + 1 + maxspikes) > params.silence_lead);
    finalspikematrix(ind) = finalspikematrix(ind) + length(ind1) - length(ind0);
end

for iLoc = 1:Nlocs
    locind(iLoc) = min1(find(dir(1,:) == params.locations(1,iLoc) & dir(2,:) == params.locations(2,iLoc)));
end

%%% check for desired best location
locmax = input('Enter best el and az (as a 2-element vector), if desired   ');
if isempty(locmax)
    [weights, IAresponse] = plot_MLSweights2(finalspikematrix, locind);
else
    [weights, IAresponse] = plot_MLSweights2(finalspikematrix, locind, locmax);
end

input('<CR> to close figures and continue...');
if exist1('H.tempdata_fig')     close(H.tempdata_fig);      end
if exist1('H.wts_fig')          close(H.wts_fig);           end
if exist1('H.IAresponse_fig')   close(H.IAresponse_fig);    end
H.tempdata_fig = [];
H.wts_fig = [];      
H.IAresponse_fig=[];        