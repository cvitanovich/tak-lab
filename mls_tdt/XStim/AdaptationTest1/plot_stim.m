% plot stim

% what's this stim???


DUR = XStimParams.curr_stimdur;
nPts = DUR*round(TDT.Fs/1000);
nEpochmsecs = XStimParams.epoch_duration(1);       % duration of each epoch at one SPL
nEpochPts = round(nEpochmsecs * TDT.Fs/1000);
nEpochs = round(nPts / nEpochPts);
nPts = nEpochs * nEpochPts;
DUR = nPts / round(TDT.Fs/1000);

SPLs_hi = 10+(XStimParams.focalSPL(1) : XStimParams.focalSPL(2));
SPLs_lo = setdiff([10:40], SPLs_hi);

%Randomize the stimuli
numtrials = size(XStimParams.locations,2);
randseq = randperm(numtrials);
trialnum = 1;

% make SPL sequences for each source
% source2 will jitter around in SPL, but fall to 0 SPL whenever source1 comes on
temp_ones = ones(nEpochPts,1);

% make SPLseq
randEpochs = rand(nEpochs,1);        
randEpochs2 = round(rand(nEpochs,1)*length(SPLs_hi)+.49999999999);  
randEpochs3 = round(rand(nEpochs,1)*length(SPLs_lo)+.49999999999); 
if isempty(SPLs_lo)
    XStimParams.focalProb = 1;
end
nEpochs1 = ceil(nEpochs/XStimParams.source1_ISI)*2;   % approx # epochs for probe at src1
%EpochJitter = round(2*(rand(nEpochs1,1)-.5))+1;
EpochJitter = ones(nEpochs1,1);

% determine which epochs will be on
if  XStimParams.repeatmode
    onEpochs = zeros(nEpochs,1); 
    onEpochs(1) = 1;
    for iepoch = 1:(XStimParams.onEpochs + XStimParams.offEpochs):nEpochs
        onEpochs(iepoch:iepoch+XStimParams.onEpochs-1) = ones(size(XStimParams.onEpochs));
    end
else
    onEpochs = ones(nEpochs,1); 
end

seq1 = zeros(nEpochs,1);                             % SPL for each epoch
seq2 = zeros(nEpochs,1);                             % SPL for each epoch
seqSPL1 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
seqSPL2 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
probenum = 0;
for iepoch = 1:nEpochs
    if onEpochs(iepoch)
        if iepoch == probenum*XStimParams.source1_ISI + EpochJitter(probenum+1)
            probenum = probenum+1;
            seq1(iepoch) = 10+XStimParams.source1_SPL;
            seq2(iepoch) = 0;
        else
            if randEpochs(iepoch) <= XStimParams.focalProb
                seq2(iepoch) = SPLs_hi(randEpochs2(iepoch));
            else
                seq2(iepoch) = SPLs_lo(randEpochs3(iepoch));
            end
        end
        seqSPL1((1+nEpochPts*(iepoch-1)):(nEpochPts*iepoch)) = seq1(iepoch) * temp_ones;
        seqSPL2((1+nEpochPts*(iepoch-1)):(nEpochPts*iepoch)) = seq2(iepoch) * temp_ones;
    else
        if iepoch == probenum*XStimParams.source1_ISI + EpochJitter(probenum)
            probenum = probenum+1;
        end
    end
end

% plot
H0 = figure; hold on
set(gcf,'position',[50 400 1200 300]);
set(gca,'position',[.05 .05 .9 .9]);

X = [0:length(seqSPL1)/(nEpochPts/2)-1; 0:length(seqSPL1)/(nEpochPts/2)-1; 1:length(seqSPL1)/(nEpochPts/2); 1:length(seqSPL1)/(nEpochPts/2)];
Y = [zeros(1,nEpochs*2); seqSPL1(1:nEpochPts*.5:end); seqSPL1(1:nEpochPts*.5:end); zeros(1,nEpochs*2)];
patch(X,Y,'r','facealpha',.3,'edgecolor','none')
X = [0:length(seqSPL1)/(nEpochPts/2)-1; 0:length(seqSPL1)/(nEpochPts/2)-1; 1:length(seqSPL2)/(nEpochPts/2); 1:length(seqSPL2)/(nEpochPts/2)];
Y = [zeros(1,nEpochs*2); seqSPL2(1:nEpochPts*.5:end); seqSPL2(1:nEpochPts*.5:end); zeros(1,nEpochs*2)];
patch(X,Y,'g','facealpha',.3,'edgecolor','none')
