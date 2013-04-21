function [DATA, params, FN] = read_data3(dataFN,dataPATH)

% function [DATA, params, FN] = read_data3(dataFN,dataPATH);
%
% to read data created by record_data3 in Xstim testing module
%% data was saved as:
%	param1	param2	param3 ...	repnum	Nspikes	spiketimes...
% first line is a dummy with the number of params (including repnum) written to param1

if nargin < 2
    dataPATH = 'e:\kip\datastor\';
end

eval(['load ' dataPATH dataFN]);

Nparams = DATA(1,1);
Nreps = params.numreps;

switch params.test_type
case 'ABL'
    inc = (params.hiabl - params.loabl)/(params.numabls-1);
    ABL = params.loabl : inc : params.hiabl;
    ABLdata = zeros(params.numabls,1);
    for trial = 2:size(DATA,1)
      ABLdata(DATA(trial,1)) = ABLdata(DATA(trial,1)) + DATA(trial,3);
    end
    figure; plot(ABL,ABLdata);
            
case 'ITD'
    inc = (params.hiitd - params.loitd)/(params.numitds-1);
    ITD = params.loitd : inc : params.hiitd;
    ITDdata = zeros(params.numitds,1);
    for trial = 2:size(DATA,1)
       ITDdata(DATA(trial,1)) = ITDdata(DATA(trial,1)) + DATA(trial,3);
    end
    figure; plot(ITD,ITDdata);
            
case 'ILD'
    inc = (params.hiild - params.loild)/(params.numilds-1);
    ILD = params.loild : inc : params.hiild;
    ILDdata = zeros(params.numilds,1);
    for trial = 2:size(DATA,1)
       ILDdata(DATA(trial,1)) = ILDdata(DATA(trial,1)) + DATA(trial,3);
    end
    figure; plot(ILD,ILDdata);
    
case 'FREQ'
    inc = (params.hifreq - params.lofreq)/(params.numfreqs-1);
    Freq = params.lofreq : inc : params.hifreq;
    Freqdata = zeros(params.numfreqs,1);
    for trial = 2:size(DATA,1)
       Freqdata(DATA(trial,1)) = Freqdata(DATA(trial,1)) + DATA(trial,3);
    end
    figure; plot(Freq,Freqdata);

case {'Space', 'Space3', '2Source', 'AltIR'}
    ele = params.locations(1,:); 
    ele = min1(ele):5:max1(ele);   nEl = length(ele);
    azi = params.locations(2,:);
    azi = min1(azi):5:max1(azi);   nAz = length(azi);
    spacedata = zeros(nEl,nAz);
    for trial = 2:size(DATA,1)
       indEl = find(ele == DATA(trial,1));
       indAz = find(azi == DATA(trial,2));
       spacedata(indEl,indAz) = spacedata(indEl,indAz) + DATA(trial,4);
    end
    figure;plotdiam(azi, ele, spacedata);
    
otherwise,
   disp('This test Not implemented')
end

title([dataFN '  ' params.test_type])