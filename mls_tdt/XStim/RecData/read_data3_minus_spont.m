function [DATA, params, FN, nSpikes, Vstr, Zscores,Vstr_list,Zscores_list] = read_data3_minus_spont(dataFN,dataPATH, frame, pflag, flag1, repflag)

% function [DATA, params, FN, nSpikes, Vstr, Zscores,Vstr_list,Zscores_list] = read_data3_minus_spont(dataFN,dataPATH,frame,pflag,flag1,repflag);
%
% to read data created by record_data3 in Xstim testing module
%% data was saved as:
%	param1	param2	param3 ...	repnum	Nspikes	spiketimes...
% first line is a dummy with the number of params (including repnum) written to param1
% argin frame: [mintime maxtime] in msec to specify which spikes to count
% (defaults to all evoked spikes)
% argin pflag: 1) plot (default)
% argin flag1: allows Vstr plots in 2source and should be a vector of frequencies of interest
% argin rep: allows testing of one or more reps (e.g. [1 3]; default is all reps)
% plot and nSpikes returned as spikes/rep

%%%%%% check out %%%%% below to (use) first spike

if nargin < 5 | isempty(flag1)   
    flag1 = 0;
    %Vstr = [];  Zscores = [];   Vstr_list = []; Zscores_list = [];
end                         % Vstrength calc
if nargin < 4 | isempty(pflag)   pflag = 1;     end                         % plotting
if nargin < 2 | isempty(dataPATH)   dataPATH = 'e:\kip\datastor\';     end

marker_color= 'bgcmywbgcmyk';

if flag1
    modPer = 1000000 ./ [1:120];            % modulation period in usec
end
eval(['load ' dataPATH dataFN]);
if nargin < 3 | isempty(frame)
    frame = [0 params.curr_stimdur];    
end     % time frame for counting spikes

spont_factor = (frame(2) - frame(1))/params.curr_stimdur;   % part of spont estimate to subtract from evoked

Nparams = DATA(1,1);

if nargin < 6 | isempty(repflag)    repflag = 1:params.numreps;    end       % use all reps by default
Nreps = length(repflag);
maxspikes = max1(DATA(:,Nparams+1));

if ~exist1('params.silence_lead')
    params.silence_lead = 100;
end
if ~exist1('params.silence_trail')
    params.silence_trail = 50;
end

switch params.test_type
    case 'ABL'
        inc = (params.hiabl - params.loabl)/(params.numabls-1);
        ABL = params.loabl : inc : params.hiabl;
        ABLdata = zeros(params.numabls,1);
        for trial = 2:size(DATA,1)
            ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
            ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(1)));
            ind2 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(2)));
            ABLdata(DATA(trial,1)) = ABLdata(DATA(trial,1)) + length(ind1) - length(ind2) - length(ind0) * spont_factor;
        end
        nSpikes = ABLdata;
        if pflag    H.datafig = figure;     plot(ABL,ABLdata);  end
        
    case {'ITD'}  
        inc = (params.hiitd - params.loitd)/(params.numitds-1);
        ITD = params.loitd : inc : params.hiitd;
        ITDdata = zeros(params.numitds,1);
        for trial = 2:size(DATA,1)
            ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
            ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(1)));
            ind2 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(2)));
            ITDdata(DATA(trial,1)) = ITDdata(DATA(trial,1)) + length(ind1) - length(ind2) - length(ind0) * spont_factor;
        end
        if pflag    H.datafig = figure;     plot(ITD,ITDdata);   end
        nSpikes = ITDdata;
    case {'ITD_decorr'}  % params: itd#  iHiFreq  iLoFreq   repnum   Nspikes
        nTests = length(find(params.stims_to_play));
        inc = (params.hiitd - params.loitd)/(params.numitds-1);
        ITD = params.loitd : inc : params.hiitd;
        ITDdata = zeros(nTests,params.numitds);
        for trial = 2:size(DATA,1)
            testtype = find(find(params.stims_to_play)==sub2ind([4 4],DATA(trial,3),DATA(trial,2)));
            ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
            ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(1)));
            ind2 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(2)));
            ITDdata(testtype,DATA(trial,1)) = ITDdata(testtype,DATA(trial,1)) + length(ind1) - length(ind2) - length(ind0) * spont_factor;
        end
        if pflag
            % testtypes in DATA file assigned by XStimParams.stims_to_play down column 1, down column2,...
            % so that e.g. if three highest priority FNs played, 
            % 1)Lo (plotted blue)  2) Hi (green)  3) Hi & Lo (cyan)
            H.datafig = figure; hold on
            for itest = 1:nTests
                plot(ITD,ITDdata(itest,:),'*-','color',marker_color(itest));
            end
            legh = legend(params.legendstr); 
            set(legh,'position',[.6 (1 - length(params.legendstr)*.04) .4 length(params.legendstr)*.04]);
            set(legh,'color',[.9 .9 .9]);
        end
        nSpikes = ITDdata;
    case 'ILD'
        inc = (params.hiild - params.loild)/(params.numilds-1);
        ILD = params.loild : inc : params.hiild;
        ILDdata = zeros(params.numilds,1);
        for trial = 2:size(DATA,1)
            ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
            ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(1)));
            ind2 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(2)));
            ILDdata(DATA(trial,1)) = ILDdata(DATA(trial,1)) + length(ind1) - length(ind2) - length(ind0) * spont_factor;
        end
        if  pflag   H.datafig = figure;    plot(ILD,ILDdata);       end
        nSpikes = ILDdata;
        
    case 'FREQ'
        inc = (params.hifreq - params.lofreq)/(params.numfreqs-1);
        Freq = params.lofreq : inc : params.hifreq;
        Freqdata = zeros(params.numfreqs,1);
        for trial = 2:size(DATA,1)
            ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
            ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(1)));
            ind2 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(2)));
            Freqdata(DATA(trial,1)) = Freqdata(DATA(trial,1)) + length(ind1) - length(ind2) - length(ind0) * spont_factor;
        end
        if pflag   H.datafig = figure; plot(Freq,Freqdata);      end
        nSpikes = Freqdata;
        
    case {'Space', 'Space3', 'Space2', '2Source', 'AltIR', 'Space FC', 'Space ILA', 'Space ITA', ...
            'Space2 FC', 'Space2 ILA', 'Space2 ITA', 'AltIR FC', 'AltIR ILA', 'AltIR ITA'}
        ele = params.locations(1,:); 
        ele = min1(ele):5:max1(ele);   nEl = length(ele);
        azi = params.locations(2,:);
        azi = min1(azi):5:max1(azi);   nAz = length(azi);
        spacedata = zeros(nEl,nAz);
        index = spacedata;
        nSpikes = zeros(size(params.locations,2),1);
        
        if flag1    
            Vstr = zeros(nEl,nAz, length(modPer));  
            Vstr_list = zeros(size(params.locations,2),length(modPer));
            Zscores_list = zeros(size(params.locations,2),length(modPer));
        end
        %if pflag    figure; hold on;     end           % for raster
        
        for trial = 2:size(DATA,1)
            indEl = find(ele == DATA(trial,1));
            indAz = find(azi == DATA(trial,2));
            index(indEl,indAz) = 1;
            rep = DATA(trial,Nparams);
            if any(ismember(rep,repflag))
                ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
                ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(1)));
                ind2 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > (params.silence_lead + frame(2)));
                spacedata(indEl,indAz) = spacedata(indEl,indAz) + length(ind1) - length(ind2) - length(ind0) * spont_factor;
                indLoc = min1(find(params.locations(1,:) == DATA(trial,1) & params.locations(2,:) == DATA(trial,2)));
                nSpikes(indLoc) = nSpikes(indLoc) + length(ind1) - length(ind2) - length(ind0) * spont_factor;
                
                % plot raster
                if 0
                    plot(DATA(trial,Nparams+2:Nparams + 1 + maxspikes),trial* ones(maxspikes,1),'.')
                    plot(DATA(trial,Nparams+1+setdiff(ind1,ind2)),trial* ones(length(ind1)-length(ind2),1),'r.')
                end
                
                if flag1 & rep == repflag(1)          % Vstrength calc
                    spiketimes = DATA(trial,Nparams+1+ind1);
                    if Nreps > 1
                        for irep = 2:Nreps
                            % find trials for this AZ, EL and rep
                            ind = find(DATA(:,1) == DATA(trial,1) & DATA(:,2) == DATA(trial,2) & DATA(:,Nparams) == repflag(irep));
                            if ~isempty(ind)    % add these spiketimes onto spiketime vector
                                for i = 1:length(ind)
                                %%%%%%% trying to skip first spike
                                %%%%%%% Nparams+2: replaced with Nparams+3:
                                ind1 = find(DATA(ind(i),Nparams+3:Nparams + 1 + maxspikes) > params.silence_lead);
                                temp = DATA(ind(i),Nparams + 1 + ind1);
                                spiketimes = [spiketimes temp(:)'];
                                end
                            end
                        end   % end of ireps
                    end   % end of if Nreps
                    
                    % adjust spiketimes for sound's convolution with HRTFs
                    spiketimes = spiketimes - 4.1;
                    % adjust spiketimes for difference in MII and TDT clocks
                    spiketimes = spiketimes .* 1.0016;
                    
                    ind = find(spiketimes <= (params.curr_stimdur + params.silence_lead));
                    spiketimes = spiketimes(ind) - params.silence_lead;
                    % spiketimes > onset
                    if length(spiketimes) >= (2* Nreps)
                        Vstr(indEl,indAz,:) = Vstrength(spiketimes' * 1000, modPer);
                    else
                        Vstr(indEl,indAz,:) = zeros(1, length(modPer));
                    end   % end of spiketimes >=
                    Vstr(indEl,indAz,:) = Vstr(indEl,indAz,:) ./ Nreps;
                    Vstr_list(indLoc,:) = reshape(Vstr(indEl,indAz,:),size(modPer));       % list of Vstr
                    temp = Vstr_list(indLoc,:);
                    if std(temp) >0
                        Zscores_list(indLoc,:) = (temp - mean(temp))./ std(temp);
                    else
                        Zscores_list(indLoc,:) = zeros(size(modPer));
                    end
                end      % end of flag1 (Vstr calc)
                
            end          % end of repflag for inclusion of only chosen reps
        end         % end of trial
        % added on Feb 11, 2003 to handle locations that are tested but yield
        % no spikes and are therefore not listed
        for iloc = 1:size(params.locations,2)
            indEl = find(ele == params.locations(1,iloc));
            indAz = find(azi == params.locations(2,iloc));
            index(indEl,indAz) = 1;
        end
        nSpikes = nSpikes/Nreps;
        spacedata = spacedata/Nreps;
        warning off
        spacedata = spacedata ./ index;     % makes non-tested loci into NaNs
        warning on
        if pflag    plotdiam(azi, ele, spacedata); H.datafig = gcf;     end
        
        if flag1
            for iEL = 1:nEl      % calc Z-scores of Vstr
                for iAZ = 1:nAz
                    temp = reshape(Vstr(iEL,iAZ,:),1,size(Vstr,3));
                    if std(temp) >0
                        Zscores(iEL,iAZ,:) = (temp - mean(temp))./ std(temp);
                    else
                        Zscores(iEL,iAZ,:) = zeros(1,1,size(modPer,2));
                    end
                end
            end
            warning off;
            for ifig = 1:length(flag1)
                Zscores(:,:,flag1(ifig)) = Zscores(:,:,flag1(ifig)) ./ index;     % makes non-tested loci into NaNs
            end
            warning on;
            if pflag 
                cmin = min1(Zscores(:,:,flag1));
                cmax = max1(Zscores(:,:,flag1));
                for ifig = 1:length(flag1)
                    % plot Zvstr for important freqs
                    plotdiam(azi, ele, reshape(Zscores(:,:,flag1(ifig)),nEl,nAz));
                    if params.ildalone_flag
                        str = 'ILA';
                    elseif params.itdalone_flag
                        str = 'ITA';
                    elseif params.ABLalone_flag
                        str = 'ABL';
                    else
                        str = 'FC';
                    end
                    if strcmp(params.test_type,'2source')
                        str = [str '2'];
                    else
                        str = [str '1'];
                    end

                    title([dataFN '  ' str ' Vstr zscore for ' num2str(flag1(ifig)) ' Hz'])
                    caxis([cmin cmax]);
                    axis([-90 90 -90 90])
                    set(gca,'xtick',-80:10:80);
                    set(gca,'ytick',-80:10:80);
                    set(gca,'xticklabel',[]);
                    set(gca,'yticklabel',[]);
                end
            end
            Zscores = Zscores(:,:,flag1);
            Vstr = Vstr(:,:,flag1);
            Vstr_list = Vstr_list';
            Zscores_list = Zscores_list';
        end
    otherwise
        disp('This test Not implemented')
        return;
    end
    
    if pflag
        figure(H.datafig)
        tith = title([dataFN '  ' params.test_type]);
    end
    
    if ~exist('FN')
        FN = [];
    end