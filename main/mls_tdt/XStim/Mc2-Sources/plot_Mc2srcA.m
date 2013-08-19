function [RLFall, n1, SPLx, RLFm] = plot_Mc2srcA(spikes_trial, SPLseq1, RLFall, n1, repnum);

% function [RLFall, n1, SPLx, RLFm] = plot_Mc2srcA(spikes_trial, SPLseq1, RLFall, n1, repnum);
%
% plots RLFs rep-by-rep for Mc2src
%
% argin:    spikestrial: list of spiketimes for this trial (msec)
%           SPLseq1: list of SPLs for each epoch for SPL1
%           RLFall: rate for each instance of each level tested
%           n1: number of tests for each level
%           repnum

global XStimParams
global H
global FN

% parameters:
SPLbins = 0:40;
SPLx = 1.08*SPLbins + 33.75 - abs(XStimParams.curr_ABL);
nSPLs = length(SPLx);

onset = 1;                    %1+#epochs to skip at onset
offset = 0;                   % 0 or last epoch# to include

Factor_TDT_MII = 1.00165;       % time adjustment between TDT and MII clocks
Factor_HRIR = 0;                %3.1;          % msec delay for HRIRs
latency = 10;                   % msec latency of first spike (to assign spikes to proper epoch)

epoch_duration = XStimParams.epoch_duration(1);
epochBins = 0:epoch_duration:XStimParams.curr_stimdur;
nSecs = epoch_duration/1000;
nEpochs = XStimParams.curr_stimdur / epoch_duration;
if offset
    nEpochs = min1([offset nEpochs]);
end

% set up to loop through reps
if nargin <3 | isempty(RLFall)
    RLFall = cell(size(SPLbins'));
end
if nargin<4 | isempty(n1)
    n1 = zeros(size(SPLbins'));
end

% get RLF
n1 = n1 + histc(SPLseq1(onset:nEpochs),SPLbins)';      % # epochs for each SPL
tempSpikeTimes = (spikes_trial-Factor_HRIR)*Factor_TDT_MII;
ind0 = find(tempSpikeTimes & tempSpikeTimes <= XStimParams.silence_lead );
ind1 = find(tempSpikeTimes > XStimParams.silence_lead );
ind2 = find(tempSpikeTimes > (XStimParams.silence_lead + XStimParams.curr_stimdur));

N = zeros(nEpochs,1);
if ~isempty(ind1)
    spiketimes = tempSpikeTimes(setdiff(ind1,ind2)) - XStimParams.silence_lead - latency;
    N = histc(spiketimes ,epochBins);      % #spikes in each epoch
end
for iEpoch = onset:nEpochs
    RLFall{SPLseq1(iEpoch)+1} = [RLFall{SPLseq1(iEpoch) + 1}  N(iEpoch)];
end

RLFm = zeros(nSPLs,1);
RLFsem = zeros(nSPLs,1);
for iSPL = 1:nSPLs
    x = RLFall{iSPL}/nSecs;
    n = length(x);
    if n
        RLFm(iSPL) = mean(x);
        RLFsem(iSPL) = std(x) / sqrt(n);
    end
end

% plotting
figure(H.Mc2SourceRLFfig)
[temp ind0] = intersect(SPLbins,10+(XStimParams.focalSPL(1):XStimParams.focalSPL(2)));     % index into focal SPLs
[temp ind1] = setdiff(SPLbins(11:end),temp);            % index to non-focal SPLs
plot(SPLx(ind0) , RLFm(ind0),'o-','color', 'b', 'MarkerFaceColor', 'b');
hold on
for iind = 1:length(ind0)
    plot([SPLx(ind0(iind)) SPLx(ind0(iind))] , [RLFsem(ind0(iind)) -RLFsem(ind0(iind))]+RLFm(ind0(iind)),'-','color', 'b');
end
plot(SPLx(ind1+10) , RLFm(ind1+9),'o','color', 'b');
hold off
xlabel('SPL (dB)')
ylabel('spike rate (sp/s)')
title([FN.data ' RLF-response to probe after ' num2str(repnum) ' reps'])
