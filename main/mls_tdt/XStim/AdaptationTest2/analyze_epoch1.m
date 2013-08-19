function [H0] = analyze_epoch1(datamatrix,param3,param4);

% plots all epochs from adapt1
% [H0] = analyze_epoch1(datamatrix,param3,param4);
%       datamatrix is the un-reorganized list of [param1 param2 rep# spiketimes]
%       param3 is seqSPL1
%       param4 is seqSPL2
%   returns the figure handle

global XStimParams


i = 1; Ntrials = 0; Nspikes = [];
while i < size(datamatrix,1)
    Nspikes(Ntrials+1) = datamatrix(i,1);
    i = i + Nspikes(Ntrials+1);
    Ntrials = Ntrials+1;
end
maxSpikes = max1(Nspikes);

Nparams = size(datamatrix,2)-2;
DATA = zeros(Ntrials+1,maxSpikes+Nparams+1);
DATA(1,1) = Nparams;

i = 1; Ntrial = 1;
while i < size(datamatrix,1)
    nexti = i + Nspikes(Ntrial);
    DATA(Ntrial+1,1:Nspikes(Ntrial) + Nparams+1) = [...		% first Ntrial reserved
            datamatrix(i,4:2+Nparams) ...		% params
            datamatrix(i,3) ...					% repnum
            Nspikes(Ntrial) ...						% Nspikes
            datamatrix(i:nexti-1,2)'];				% spiketimes
    i = nexti;
    Ntrial = Ntrial+1;   
end
%%%%%%%%%%%%%%%%%%%%%5    
% parameters:
atten = abs(XStimParams.curr_ABL);
nReps = max1(DATA(:,3));
binwidth = 1;

% make epochbins, E_dur, SPL
nSPL1 = size(XStimParams.SPL1,2);
nSPL2 = size(XStimParams.SPL2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
iseg = 1;
E_dur(iseg) = XStimParams.silence_lead;
epochbins = [0 E_dur];
SPL(iseg)=0;
% calc epoch duration and epochbins for all epochs
for j = 1:nSPL1
    for i = 1:nSPL2
        iseg=iseg+1;
        E_dur(iseg) = XStimParams.epoch_duration(1);
        epochbins = [epochbins epochbins(end)+E_dur(iseg)];
        SPL(iseg)=XStimParams.SPL1(j);
        
        iseg=iseg+1;
        E_dur(iseg) = XStimParams.epoch_duration(2);
        epochbins = [epochbins epochbins(end)+E_dur(iseg)];
        SPL(iseg)=XStimParams.SPL2(i);
    end
    if j<nSPL1
        iseg=iseg+1;
        E_dur(iseg) = XStimParams.silence_mid;
        epochbins = [epochbins epochbins(end)+E_dur(iseg)];
        SPL(iseg)=0;
    end
end
iseg = iseg+1;
E_dur(iseg) = XStimParams.silence_trail;
epochbins = [epochbins epochbins(end)+E_dur(iseg)];
nEpochs = iseg;
SPL(iseg)=0;

% maxspikes for any rep
maxspikes = max1(DATA(:,Nparams+1));
% estimate maxEspikes   (max spikes in any epoch)
ind1 = find(DATA(2,Nparams+2:Nparams + 1 + maxspikes) > XStimParams.silence_lead );
ind2 = find(DATA(2,Nparams+2:Nparams + 1 + maxspikes) > (XStimParams.silence_lead + XStimParams.curr_stimdur + XStimParams.silence_trail));
if ~isempty(ind1)
    spiketimes = (DATA(2,Nparams+1+setdiff(ind1,ind2)))*1.0017 - XStimParams.silence_lead;
    N = histc(spiketimes,epochbins);      % #spikes in each epoch
end
maxEspikes = roundn(max1(N),25)+25;                         % estimated max spikes in any epoch

% make Spikes array
iseg = 1;
Spikes{1} = struct( 'SPL_lead',   0, ...
    'SPL_current',   0,...
    'dur_lead',   0,...
    'dur_current',   [E_dur(iseg)],...
    'spiketimes',   zeros(nReps,maxEspikes));

for j = 1:nSPL1
    for i = 1:nSPL2
        iseg=iseg+1;
        Spikes{iseg} = struct( 'SPL_lead',   Spikes{iseg-1}.SPL_current, ...
            'SPL_current',   SPL(iseg),...
            'dur_lead',   [E_dur(iseg-1)],...
            'dur_current',   [E_dur(iseg)],...
            'spiketimes',   zeros(nReps,maxEspikes));
        
        iseg=iseg+1;
        Spikes{iseg} = struct( 'SPL_lead',   Spikes{iseg-1}.SPL_current, ...
            'SPL_current',   SPL(iseg),...
            'dur_lead',   [E_dur(iseg-1)],...
            'dur_current',   [E_dur(iseg)],...
            'spiketimes',   zeros(nReps,maxEspikes));
    end
    if j<nSPL1
        iseg=iseg+1;
        Spikes{iseg} = struct( 'SPL_lead',   Spikes{iseg-1}.SPL_current, ...
            'SPL_current',   0,...
            'dur_lead',   [E_dur(iseg-1)],...
            'dur_current',   [E_dur(iseg)],...
            'spiketimes',   zeros(nReps,maxEspikes));
    end
end
iseg = iseg+1;
Spikes{iseg} = struct( 'SPL_lead',   Spikes{iseg-1}.SPL_current, ...
    'SPL_current',   0,...
    'dur_lead',   [E_dur(iseg-1)],...
    'dur_current',   [E_dur(iseg)],...
    'spiketimes',   zeros(nReps,maxEspikes));

% have to do each rep separately as sequence differs between reps...
% allow choice of which reps to use
for iRep = 1:nReps
    % set up SPL for this rep
    tSPL = zeros(1,nEpochs);
    iseg = 1;
    for j = 1:nSPL1
        for i = 1:nSPL2
            iseg=iseg+1;
            tSPL(iseg) = param3(iRep,j)-10;
            iseg=iseg+1;
            tSPL(iseg) = param4(iRep,i)-10;
        end
        if j<nSPL1
            iseg=iseg+1;
        end
    end
    
    rep_ind = find(DATA(:,3)==iRep);
    for iRep_ind = 1:length(rep_ind)
        % index to all spikes
        ind1 = Nparams + 1 + find(DATA(rep_ind(iRep_ind),Nparams+2:Nparams + 1 + maxspikes) > 0 );
        if ~isempty(ind1)
            % index to spikes after trailing zeros
            ind2 = find(DATA(rep_ind(iRep_ind),ind1) > (XStimParams.silence_lead + XStimParams.curr_stimdur + XStimParams.silence_trail));
            ind1 = setdiff(ind1,ind1(ind2));
            % new time factor as of Apr 20, 2007
            spiketimes = DATA(rep_ind(iRep_ind),ind1) * 1.0017;
            % first epoch
            ind3 = find(spiketimes>epochbins(1)  & spiketimes<epochbins(2));
            temp_times = spiketimes(ind3);
            nEspikes = length(temp_times);
            Spikes{1}.spiketimes(iRep,1:nEspikes) = temp_times;
            
            % all other epochs
            for iEpoch = 2:nEpochs
                % determine which segment in Spikes array
                indEpoch = intersect(find(SPL==tSPL(iEpoch-1))+1,setdiff(find(SPL==tSPL(iEpoch)),1));
                if length(indEpoch)>1
                    indEpoch = intersect(find(SPL==tSPL(iEpoch-1))+2,indEpoch);
                end
                if isempty(indEpoch)
                    indEpoch = nEpochs+1;
                    Existflag = exist1('Spikes{indEpoch}');
                    while Existflag
                        if ~(Spikes{indEpoch}.SPL_current == tSPL(iEpoch) &...
                                Spikes{indEpoch}.SPL_lead == tSPL(iEpoch-1))
                            indEpoch = indEpoch+1;
                            Existflag = exist1('Spikes{indEpoch}');
                        else
                            Existflag =0;
                        end
                    end
                    if ~exist1('Spikes{indEpoch}')        
                        Spikes{indEpoch} = struct( 'SPL_lead',   tSPL(iEpoch-1), ...
                            'SPL_current',   tSPL(iEpoch),...
                            'dur_lead',   [E_dur(iEpoch-1)],...
                            'dur_current',   [E_dur(iEpoch)],...
                            'spiketimes',   zeros(nReps,maxEspikes));
                    end
                end 
                
                ind3 = find(spiketimes>epochbins(iEpoch)  & spiketimes<epochbins(iEpoch+1));
                if ~isempty(ind3)
                    temp_times = spiketimes(ind3)-epochbins(iEpoch);
                    nEspikes = length(temp_times);
                    if nEspikes > maxEspikes
                        for i = 1:nReps
                            Spikes{indEpoch}.spiketimes(i,:) = [Spikes{indEpoch}.spiketimes(i,:) zeros(1,nEspikes-maxEspikes)];
                        end
                        maxEspikes = nEspikes;
                    end
                    Spikes{indEpoch}.spiketimes(iRep,1:nEspikes) = temp_times;
                end % isempty(ind3)
            end     % iEpoch
        end         % isempty(ind1)
    end             % irep_ind
end                 % irep

% plot composite figure of all epochs
H0 = figure; hold on
set(H0,'position',[30 250 950 300])
last_time = 0; maxN = 0;
for iEpoch = 1:length(Spikes)
    timebin = 0:binwidth:Spikes{iEpoch}.dur_current;
    ind = find(Spikes{iEpoch}.spiketimes>0);
    if ~isempty(ind)
        N = histc(Spikes{iEpoch}.spiketimes(ind),timebin);      % #spikes in each epoch
        temp = smooth(N,3)/nReps * (1000/binwidth);
        plot(last_time+timebin(1:end-2),temp(1:end-2),'b')
        maxN1 = max1(temp);
        maxN = max1([maxN1 maxN]);
        if iEpoch<nEpochs
            text(last_time+Spikes{iEpoch}.dur_current*.2,max(temp)*1.1,[num2str(Spikes{iEpoch}.SPL_current)])
        else
            text(last_time+Spikes{iEpoch}.dur_current*.2,max(temp)*1.1,[num2str(Spikes{iEpoch}.SPL_lead) '->' num2str(Spikes{iEpoch}.SPL_current)])
        end
    end     % isempty(ind)
    last_time = last_time+timebin(end)+1;
end         % iEpoch
axis([0 last_time 0 maxN*1.25]);
set(gca,'position',[.05 .08 .94 .9])
title([ 'response to all epochs'])