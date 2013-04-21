function [DATA, params, FN] = read_data3_minus_spontXXX(dataFN,dataPATH, flag1)

% function [DATA, params, FN] = read_data3_minus_spontXXX(dataFN,dataPATH,flag1);
%
% to read data created by record_data3 in Xstim testing module
%% data was saved as:
%	param1	param2	param3 ...	repnum	Nspikes	spiketimes...
% first line is a dummy with the number of params (including repnum) written to param1
% argin flag1: allows Vstr plots in 2source and should be a vector of frequencies of interest
%
%         %%%% KLUGE to find locs for 2 src data prior to 3/4/03

if nargin < 3
    flag1 = 0;
end
if nargin < 2
    dataPATH = 'e:\kip\datastor\';
end

modPer = 1000000 ./ [1:100];            % modulation period in usec

eval(['load ' dataPATH dataFN]);

Nparams = DATA(1,1);
Nreps = params.numreps;
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
      ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > params.silence_lead);
      ABLdata(DATA(trial,1)) = ABLdata(DATA(trial,1)) + length(ind1) - length(ind0);
    end
    H.datafig = figure; plot(ABL,ABLdata);
            
case 'ITD'
    inc = (params.hiitd - params.loitd)/(params.numitds-1);
    ITD = params.loitd : inc : params.hiitd;
    ITDdata = zeros(params.numitds,1);
    for trial = 2:size(DATA,1)
       ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
       ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > params.silence_lead);
       ITDdata(DATA(trial,1)) = ITDdata(DATA(trial,1)) + length(ind1) - length(ind0);
    end
    H.datafig = figure; plot(ITD,ITDdata);
            
case 'ILD'
    inc = (params.hiild - params.loild)/(params.numilds-1);
    ILD = params.loild : inc : params.hiild;
    ILDdata = zeros(params.numilds,1);
    for trial = 2:size(DATA,1)
       ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
       ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > params.silence_lead);
       ILDdata(DATA(trial,1)) = ILDdata(DATA(trial,1)) + length(ind1) - length(ind0);
    end
    H.datafig = figure; plot(ILD,ILDdata);
    
case 'FREQ'
    inc = (params.hifreq - params.lofreq)/(params.numfreqs-1);
    Freq = params.lofreq : inc : params.hifreq;
    Freqdata = zeros(params.numfreqs,1);
    for trial = 2:size(DATA,1)
      ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
      ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > params.silence_lead);
      Freqdata(DATA(trial,1)) = Freqdata(DATA(trial,1)) + length(ind1) - length(ind0);
    end
    H.datafig = figure; plot(Freq,Freqdata);

case {'Space', 'Space3', '2Source', 'AltIR'}
        %%%% kluge to find locs for 2 src data prior to 3/4/03
    ele = params.locations(1,:);    nEl = length(ele);  ele = ele(1:ceil(nEl/2));
    ele = min1(ele):5:max1(ele);   nEl = length(ele);
    azi = params.locations(2,:);    nAz = length(azi);  azi = azi(1:ceil(nAz/2));
    azi = min1(azi):5:max1(azi);   nAz = length(azi);
    spacedata = zeros(nEl,nAz);
    index = spacedata;
    if flag1    Vstr = zeros(nEl,nAz, length(modPer));  end
    for trial = 2:size(DATA,1)
        %%%% kluge to find locs for 2 src data prior to 3/4/03
       ind = ceil(find(params.locations(1,:) == DATA(trial,1) & params.locations(2,:) == DATA(trial,2))/2);
       indEl = find(ele == params.locations(1,ind));
       indAz = find(azi == params.locations(2,ind));
       index(indEl,indAz) = 1;
       rep = DATA(trial,Nparams);
       ind0 = find((DATA(trial,Nparams+2:Nparams + 1 + maxspikes) < params.silence_lead) & (DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > 0));
       ind1 = find(DATA(trial,Nparams+2:Nparams + 1 + maxspikes) > params.silence_lead);
       spacedata(indEl,indAz) = spacedata(indEl,indAz) + length(ind1) - length(ind0);
       if flag1 & rep == 1
          spiketimes = DATA(trial,Nparams+1+ind1);
          if Nreps > 1
          for irep = 2:Nreps
              ind = find(DATA(:,1) == DATA(trial,1) & DATA(:,2) == DATA(trial,2) & DATA(:,Nparams) == irep);
              if ~isempty(ind)
                  ind1 = find(DATA(ind,Nparams+2:Nparams + 1 + maxspikes) > params.silence_lead);
                  spiketimes = [spiketimes DATA(ind,Nparams + 1 + ind1)];
              end
          end   % end of ireps
          end   % end of if Nreps
          ind = find(spiketimes <= (params.curr_stimdur + params.silence_lead));
          spiketimes = spiketimes(ind);
          if length(spiketimes) >= (5* Nreps)
             Vstr(indEl,indAz,:) = Vstrength(spiketimes' * 1000, modPer);
          else
             Vstr(indEl,indAz,:) = zeros(1, length(modPer));
          end   % end of spiketimes >=4
       end      % end of flag1
    end         % end of trial
    % added on Feb 11, 2003 to handle locations that are tested but yield
    % no spikes and are therefore not listed
    for iloc = 1:length(params.locations)
       indEl = find(ele == params.locations(1,iloc));
       indAz = find(azi == params.locations(2,iloc));
       index(indEl,indAz) = 1;
   end

    warning off
    spacedata = spacedata ./ index;     % makes non-tested loci into NaNs
    warning on
    H.datafig = figure; plotdiam1(azi, ele, spacedata);

    if flag1
      for ifig = 1:length(flag1)
        H.Vstrfig(ifig) = figure; plotdiam(azi, ele, Vstr(:,:,flag1(ifig)));
        title([dataFN '  ' params.test_type ' Vstr for ' num2str(flag1(ifig)) ' Hz'])
      end
    end
otherwise
   disp('This test Not implemented')
end


figure(H.datafig)
title([dataFN '  ' params.test_type])
