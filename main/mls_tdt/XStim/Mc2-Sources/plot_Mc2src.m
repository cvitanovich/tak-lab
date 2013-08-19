% script to plot Mc2src
% called by xstim

% load dataFN
eval(['load ' FN.data_path FN.data])
Nparams = DATA(1,1);
maxspikes = max1(DATA(:,Nparams+1));
nReps = length(params.numreps);
duration = params.curr_stimdur;

% parameters:
SPLbins = 0:40;
SPLx = 1.08*SPLbins + 33.75 - abs(params.curr_ABL);
nSPLs = length(SPLx);

onset = 1;                    %1+#epochs to skip at onset
offset = 0;                   % 0 or last epoch# to include

Factor_TDT_MII = 1.00165;       % time adjustment between TDT and MII clocks
Factor_HRIR = 0;                %3.1;          % msec delay for HRIRs
latency = 10;                   % msec latency of first spike (to assign spikes to proper epoch)

epoch_duration = params.epoch_duration(1);
epochspersec = 1000/epoch_duration;
epochBins = 0:epoch_duration:duration;
nSecs = epoch_duration/1000;
nEpochs = duration / epoch_duration;
if offset
    nEpochs = min1([offset nEpochs]);
end

if ~iscell(param3)
    SPLseq1 = param3(1,:);
end

% set up to loop through reps
spikenum = 0;
spont = 0;
nreps = 0;
RLF1 = zeros(size(SPLbins'));
RLF1all = cell(size(SPLbins'));
n1 = zeros(size(SPLbins'));

for iRep = params.numreps
    rep_ind = find(DATA(:,3)==iRep);
    N = zeros(nEpochs,1);
    for iRep_ind = 1:length(rep_ind)
        nreps=nreps+1;
        if iscell(param3)
            SPLseq1 = param3{rep_ind(iRep_ind)-1,1};
        end
        n1 = n1+histc(SPLseq1(onset:nEpochs),SPLbins)';      % # epochs for each SPL

        tempSpikeTimes = (DATA(rep_ind(iRep_ind),Nparams+2:end)-Factor_HRIR)*Factor_TDT_MII;
        ind0 = find(tempSpikeTimes & tempSpikeTimes <= params.silence_lead );
        ind1 = find(tempSpikeTimes > params.silence_lead );
        ind2 = find(tempSpikeTimes > (params.silence_lead + params.curr_stimdur));

        spont = spont + length(ind0);

        if ~isempty(ind1)
            spiketimes = tempSpikeTimes(setdiff(ind1,ind2)) - params.silence_lead - latency;
            N = histc(spiketimes ,epochBins);      % #spikes in each epoch
        end
        for iEpoch = onset:nEpochs
            RLF1(SPLseq1(iEpoch)+1) = RLF1(SPLseq1(iEpoch) + 1) + N(iEpoch);
            RLF1all{SPLseq1(iEpoch)+1} = [RLF1all{SPLseq1(iEpoch) + 1}  N(iEpoch)];
        end
    end
end

warning off MATLAB:divideByZero
RLF1 = RLF1./n1;            % now #spikes per epoch
RLF1 = RLF1/nSecs;          % per second
spont = spont/nreps * (1000/params.silence_lead);
warning on MATLAB:divideByZero

RLF1sem = zeros(nSPLs,1);
for iSPL = 1:nSPLs
    x = RLF1all{iSPL}/nSecs;
    n = length(x);
    if n
        RLF1sem(iSPL) = std(x) / sqrt(n);
    end
end

[temp ind0] = intersect(SPLbins,10+(params.focalSPL(1):params.focalSPL(2)));     % index into focal SPLs
[temp ind1] = setdiff(SPLbins,temp);            % index in non-focal SPLs
if ~isempty(ind1)   
    ind1 = ind1(2:end); 
end
plot(SPLx(ind0) , RLF1(ind0),'o-','color', 'b', 'MarkerFaceColor', 'b');
for iind = 1:length(ind0)
    plot([SPLx(ind0(iind)) SPLx(ind0(iind))] , [RLF1sem(ind0(iind)) -RLF1sem(ind0(iind))]+RLF1(ind0(iind)),'-','color', 'b');
end
plot([min(SPLx) max(SPLx)] , [1 1] *spont,'-','color', 'b');
plot(SPLx(ind1) , RLF1(ind1),'o','color', 'b');
xlabel('SPL (dB)')
ylabel('spike rate (sp/s)')
title([FN.data ' response to probe as RLF'])