% plot stim

% what's this stim???
if strcmp(XStimParams.test_type,'Mc2SourceReversed')
    flag_reverse = 1;
else
    flag_reverse = 0;
end


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

seq1 = zeros(nEpochs,1);                             % SPL for each epoch
seq2 = zeros(nEpochs,1);                             % SPL for each epoch
seqSPL1 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt
seqSPL2 = zeros(1,nEpochs * nEpochPts);              % SPL for each time pt


if XStimParams.repeatmode
    dur(1) = XStimParams.onEpochs + XStimParams.offEpochs;        % total epochs of each rep
    n(1) = XStimParams.onEpochs/(XStimParams.source1_ISI);                           % #probes during adapter ON
    n(2) = XStimParams.offEpochs/(XStimParams.source1_ISI);                           % #probes during adapter OFF
else
    dur(1) = nEpochs;        % total epochs of each rep
    n(1) = round(nEpochs/(XStimParams.source1_ISI));        % #probes during adapter ON    
end

if XStimParams.probeFlag(1)
    if flag_reverse
        if randEpochs(1) <= XStimParams.focalProb
            seq1(1) = SPLs_hi(randEpochs2(1));
        else
            seq1(1) = SPLs_lo(randEpochs3(1));
        end
    else
        seq1(1) = 10+XStimParams.source1_SPL;
    end
    ind0 = 1:nEpochPts;
    seqSPL1(ind0) = seq1(1) * ones(size(ind0));
end
for iepoch = 2:dur(1):nEpochs
    % probes
    ind1 = [];
    if XStimParams.probeFlag(1)
        ind1 = [(iepoch-1):XStimParams.source1_ISI:iepoch+XStimParams.source1_ISI*(n(1))];
    end
    if XStimParams.repeatmode
        if XStimParams.probeFlag(2)
            ind1 = [ind1 (iepoch+XStimParams.source1_ISI*(1+n(1))-1):XStimParams.source1_ISI:iepoch+XStimParams.source1_ISI*(n(1)+n(2))];
        end
    end
    
    for jepoch = 1:length(ind1)
        ind3 = ind1(jepoch);
        if flag_reverse     % in reverse mode, probes vary in level, adapter is constant
            if randEpochs(jepoch+1) <= XStimParams.focalProb
                seq1(ind3) = SPLs_hi(randEpochs2(jepoch+1));
            else
                seq1(ind3) = SPLs_lo(randEpochs3(jepoch+1));
            end
        else
            seq1(ind3) = 10+XStimParams.source1_SPL;
        end
        ind0 = (ind3-1)*nEpochPts+1:ind3*nEpochPts;
        seqSPL1(ind0) = seq1(ind1(jepoch)) * ones(size(ind0));
    end
    
    % adapter on during onEpochs when probe is not on
    if XStimParams.repeatmode
        ind2 = setdiff([iepoch:iepoch+XStimParams.onEpochs-1],ind1);
    else
        ind2 = setdiff([1:dur(1)],ind1);
    end
    for jepoch = 1:length(ind2) 
        ind3 = ind2(jepoch);
        if flag_reverse     % in reverse mode, probes vary in level, adapter is constant
            seq2(ind3) =  10+XStimParams.source1_SPL;
        else
            if randEpochs(ind3) <= XStimParams.focalProb
                seq2(ind3) = SPLs_hi(randEpochs2(ind3));
            else
                seq2(ind3) = SPLs_lo(randEpochs3(ind3));
            end
        end
        ind0 = (ind3-1)*nEpochPts+1:ind3*nEpochPts;
        seqSPL2(ind0) = seq2(ind2(jepoch)) * ones(size(ind0));
    end
end   

clear ind0 ind1 ind2 ind3 n dur

% plot
H0 = figure; hold on
set(gcf,'position',[50 400 1200 300]);
set(gca,'position',[.05 .05 .9 .9]);

X = [0:length(seqSPL1)/(nEpochPts/2)-1; 0:length(seqSPL1)/(nEpochPts/2)-1; 1:length(seqSPL1)/(nEpochPts/2); 1:length(seqSPL1)/(nEpochPts/2)];
Y = [zeros(1,nEpochs*2); seqSPL1(1:nEpochPts*.5:end); seqSPL1(1:nEpochPts*.5:end); zeros(1,nEpochs*2)];
patch(X,Y,'r','facealpha',.3,'edgecolor','none')
X1 = X(end);    Y1 = max1(Y);
X = [0:length(seqSPL1)/(nEpochPts/2)-1; 0:length(seqSPL1)/(nEpochPts/2)-1; 1:length(seqSPL2)/(nEpochPts/2); 1:length(seqSPL2)/(nEpochPts/2)];
Y = [zeros(1,nEpochs*2); seqSPL2(1:nEpochPts*.5:end); seqSPL2(1:nEpochPts*.5:end); zeros(1,nEpochs*2)];
patch(X,Y,'g','facealpha',.3,'edgecolor','none')
axis([0 max1([X1 X(end)]) 0 max1([Y(:); Y1])])
