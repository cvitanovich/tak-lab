% adapt1_plotstim

randStimFlag = get(H.Adapt1_randStimFlag,'value');
nEpochmsecs = XStimParams.epoch_duration;       % durations for epochs1 and 2
SPLs_1 = 10 + XStimParams.SPL1;                 % SPLs for epoch1 and epoch2
SPLs_2 = 10 + XStimParams.SPL2;

silence_lead_pts = (XStimParams.silence_lead * round(TDT.Fs/1000));
silence_mid_pts = (XStimParams.silence_mid * round(TDT.Fs/1000));
silence_trail_pts = (XStimParams.silence_trail * round(TDT.Fs/1000));

%%%%%%%%% make SPLseq, load buffer, while this one plays
nEpochPts = round(nEpochmsecs * TDT.Fs/1000);
ones1 = ones(nEpochPts(1),1);
ones2 = ones(nEpochPts(2),1);
zeros1= zeros(silence_mid_pts,1);
if randStimFlag
    tempACseq1 = SPLs_1(randperm(length(SPLs_1)));
    tempACseq2 = SPLs_2(randperm(length(SPLs_2)));
else
    tempACseq2 = SPLs_2;    
    tempACseq1 = SPLs_1;    
end
start = 0; iseg = 0;
for j = 1:length(SPLs_1)
    for i = 1:length(SPLs_2)
        iseg=iseg+1;
        seqSPL(start+1:start+nEpochPts(1)) = tempACseq1(j) *ones1;
        start = start + nEpochPts(1);
        seqSPL(start+1:start+nEpochPts(2)) = tempACseq2(i) *ones2;    
        start = start+nEpochPts(2);
        ACseq1(1,iseg) = tempACseq1(j);
        ACseq2(1,iseg) = tempACseq2(i);
    end
    if j<length(SPLs_1)
        seqSPL(start+1:start+silence_mid_pts) = zeros1;
        start = start+silence_mid_pts;
    end
end

H0 = figure; hold on
set(gcf,'position',[50 400 1200 300]);
set(gca,'position',[.05 .05 .9 .9]);
temp_fact = round((min1(XStimParams.epoch_duration)/2));
temp_pts = length(seqSPL)/temp_fact;
X = [0:temp_pts-1; 0:temp_pts-1; 1:temp_pts; 1:temp_pts];
Y = [zeros(1,temp_pts); seqSPL(1:temp_fact:end); seqSPL(1:temp_fact:end); zeros(1,temp_pts)];
patch(X,Y,'g','facealpha',.3,'edgecolor','none')
if length(SPLs_1)==1    plot([0 temp_pts],[SPLs_1 SPLs_1],'k--');    end
axis([0 max1(X) 0 max1(Y)])