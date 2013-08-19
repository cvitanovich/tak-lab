function [H0] = plotMcStim(McTestType,datamatrix,param,FigParams);

% based upon Engage_processMcStim
% [H0] = plotMcStim(McTestType,datamatrix,param3,param4,FigParams);
%       McTestType: 'SPLtest ' or 'AMdepth ' or 'AMperiod'
%       datamatrix is the un-reorganized list of [param1 param2 rep# spiketimes]
%       param is ACseq
%       FigParams.handle
%       FigParams.marker    (marker color, symbol as one string)
%   returns the figure handle
global XStimParams
global TDT

dt = 1000/TDT.Fs;

if isempty(FigParams.handle)
    H0 = figure; hold on
else
    H0 = figure(FigParams.handle);
end

i = 1; Ntrials = 0; Nspikes = [];
while i < size(datamatrix,1)
    Nspikes(Ntrials+1) = datamatrix(i,1);
    i = i + Nspikes(Ntrials+1);
    Ntrials = Ntrials+1;
end
maxspikes = max1(Nspikes);
if maxspikes
    
    Nparams = size(datamatrix,2)-2;
    DATA = zeros(Ntrials+1,maxspikes+Nparams+1);
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
    %%%%%%%%%%%%%%%%%%%%%
    onset = 2;                                  %1+#epochs to skip at onset
    latency = 0;                                % msec latency of first spike (to assign spikes to proper epoch)
    
    if strcmp(McTestType,'SPLtest  ')
        nSecs = XStimParams.epoch_duration(1)/1000;
        % parameters:
        SPLbins = 10:40;
        n = zeros(length(SPLbins),1);               % to hold #SPLs played from each bin
        Abscissa = 1.08*SPLbins + 33.75 - abs(XStimParams.curr_ABL);
        
        epochBins = 0:XStimParams.epoch_duration(1):XStimParams.curr_stimdur;
        Ordinate = zeros(size(SPLbins'));
        nEpochs = length(epochBins)-1;
        
        % have to do each rep separately as sequence differs between reps...
        
        for iRep = 1:XStimParams.numreps
            SPLseq = param(iRep,:);
            n = n+histc(SPLseq,SPLbins)';      % # epochs for each SPL
            lowSPL = min1(SPLseq);
            rep_ind = find(DATA(:,3)==iRep);
            for iRep_ind = 1:length(rep_ind)
                ind1 = find(DATA(rep_ind(iRep_ind),Nparams+2:Nparams + 1 + maxspikes) > XStimParams.silence_lead );
                ind2 = find(DATA(rep_ind(iRep_ind),Nparams+2:Nparams + 1 + maxspikes) > (XStimParams.silence_lead + XStimParams.curr_stimdur));
                if ~isempty(ind1)
                    spiketimes = (DATA(rep_ind(iRep_ind),Nparams+1+setdiff(ind1,ind2))-4.1)*1.0016 - XStimParams.silence_lead;
                    N = histc(spiketimes - latency,epochBins);      % #spikes in each epoch
                end
                for iEpoch = onset:nEpochs
                    Ordinate(SPLseq(iEpoch)+1-lowSPL) = Ordinate(SPLseq(iEpoch) + 1-lowSPL) + N(iEpoch);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    elseif strcmp(McTestType,'AMdepth  ') | strcmp(McTestType,'AMdepthdB')
        nSecs = 1/XStimParams.mod_freq(1);
        nCycles = XStimParams.mod_freq(1) * XStimParams.curr_stimdur/1000;    
        DEPTHbins = .1:.1:1;
        n = zeros(length(DEPTHbins),1);               % to hold #DEPTHs played from each bin
        Ordinate = zeros(size(DEPTHbins'));
        cycleBins = 0:nCycles;                        % to hold #spikes in each AMcycle
        modPer = TDT.Fs/XStimParams.mod_freq(1);
        Vstr = zeros(length(DEPTHbins),5);             % will hold nspikes, summed sin(phase), summed cos(phase), Vstr, Theta
        
        for iRep = 1:XStimParams.numreps
            DEPTHseq = roundn(param(iRep,:),.1);
            n = n+histc(DEPTHseq,DEPTHbins)';      % # epochs for each AMdepth
            rep_ind = find(DATA(:,3)==iRep);
            for iRep_ind = 1:length(rep_ind)
                ind1 = find(DATA(rep_ind(iRep_ind),Nparams+2:Nparams + 1 + maxspikes) > XStimParams.silence_lead );
                ind2 = find(DATA(rep_ind(iRep_ind),Nparams+2:Nparams + 1 + maxspikes) > (XStimParams.silence_lead + XStimParams.curr_stimdur));
                if ~isempty(ind1)
                    spiketimes = (DATA(rep_ind(iRep_ind),Nparams+1+setdiff(ind1,ind2))-4.1)*1.0016 - XStimParams.silence_lead;
                    N = histc(spiketimes - latency,cycleBins);      % #spikes in each epoch
                    st0 = 1;        %starting pt of current cycle
                    for icycle = onset:nCycles
                        indbin = round(DEPTHseq(icycle)*10);
                        Ordinate(indbin) = Ordinate(indbin) + N(icycle);
                        st1 = st0+modPer;            % starting pt of next cycle  
                        ind0 = find(spiketimes>=st0*dt & spiketimes<st1*dt);
                        for ispike = 1:length(ind0)
                            phas = (spiketimes(ind0(ispike))/dt - st0+1)/modPer * 2* pi;
                            Vstr(indbin,1) = Vstr(indbin,1)+1;
                            Vstr(indbin,2) = Vstr(indbin,2)+sin(phas);
                            Vstr(indbin,3) = Vstr(indbin,3)+cos(phas);
                        end
                        st0=st1;
                    end
                end
            end
        end
        Abscissa = DEPTHbins;    
        ind = find(Vstr(:,1));
        ycoord = (Vstr(ind,2)./Vstr(ind,1));
        xcoord = (Vstr(ind,3)./Vstr(ind,1));
        Vstr(ind,4) = sqrt( ycoord.^2 + xcoord.^2);     % Vstr
        Vstr(ind,5) = atan2(ycoord,xcoord);             % Theta
        ind = find(Vstr(:,5)<0);
        Vstr(ind,5) = Vstr(ind,5)+2*pi;
        Ordinate2 = Vstr(:,4);
       %%%%%%%%%%%%%%%%%%%%%%%%    
    elseif strcmp(McTestType,'AMperiod ') | strcmp(McTestType,'AMperiodB')
        PERIODbins = [150 158 167 176 187 200 214 231 250 273 300 333 375 429 500 600 750 1000 1500];        % in units of sampling pts
        nSecs = PERIODbins/TDT.Fs;
        % Mfreq = TDT.Fs./PERIODbins
        n = zeros(length(PERIODbins),1);              % to hold #PERIODs played from each bin
        Ordinate = zeros(size(PERIODbins'));
        Vstr = zeros(length(PERIODbins),5);             % will hold nspikes, summed sin(phase), summed cos(phase), Vstr, Theta
        
        for iRep = 1:XStimParams.numreps
            PERIODseq = roundn(param(iRep,:),.1);
            nCycles = length(PERIODseq);    
            n = n+histc(PERIODseq,PERIODbins)';      % # epochs for each AMdepth
            rep_ind = find(DATA(:,3)==iRep);
            for iRep_ind = 1:length(rep_ind)
                ind1 = find(DATA(rep_ind(iRep_ind),Nparams+2:Nparams + 1 + maxspikes) > XStimParams.silence_lead );
                ind2 = find(DATA(rep_ind(iRep_ind),Nparams+2:Nparams + 1 + maxspikes) > (XStimParams.silence_lead + XStimParams.curr_stimdur));
                if ~isempty(ind1)
                    spiketimes = (DATA(rep_ind(iRep_ind),Nparams+1+setdiff(ind1,ind2))-4.1)*1.0016 - XStimParams.silence_lead;
                    st0 = 1;        %starting pt of current cycle
                    for icycle = 1:nCycles
                        st1 = st0+PERIODseq(icycle);            % starting pt of next cycle  
                        ind0 = find(spiketimes>=st0*dt & spiketimes<st1*dt);
                        for ispike = 1:length(ind0)
                            indbin = find(PERIODbins==PERIODseq(icycle));
                            Ordinate(indbin) = Ordinate(indbin)+1;
                            phas = (spiketimes(ind0(ispike))/dt - st0+1)/PERIODseq(icycle) * 2* pi;
                            Vstr(indbin,1) = Vstr(indbin,1)+1;
                            Vstr(indbin,2) = Vstr(indbin,2)+sin(phas);
                            Vstr(indbin,3) = Vstr(indbin,3)+cos(phas);
                        end
                        st0=st1;
                    end     % icycle
                end     % if ~isempty(ind1)
            end         % iRep_ind
        end             % iRep
        Abscissa = PERIODbins;
        ind = find(Vstr(:,1));
        ycoord = (Vstr(ind,2)./Vstr(ind,1));
        xcoord = (Vstr(ind,3)./Vstr(ind,1));
        Vstr(ind,4) = sqrt( ycoord.^2 + xcoord.^2);         % Vstr
        Vstr(ind,5) = atan2(ycoord,xcoord);                 % Theta
        ind = find(Vstr(:,5)<0);
        Vstr(ind,5) = Vstr(ind,5)+2*pi;
        Ordinate2 = Vstr(:,4);
    end
    
    ind = find(n);
    Ordinate(ind) = Ordinate(ind)./n(ind);   % now #spikes per epoch
    Ordinate = Ordinate(:)./nSecs(:);              % per second
    
    if isempty(FigParams.handle)
        plot(Abscissa, Ordinate, 'bo-');
        if exist('Ordinate2') 
            plot(Abscissa, Ordinate2, 'b-');
        end
    else
        eval(['plot(Abscissa, Ordinate, ' '''' FigParams.marker '''' ')']);
        if exist('Ordinate2') 
            eval(['plot(Abscissa, Ordinate2, ' '''' FigParams.marker(1) '''' ')']);
        end
    end
    
end