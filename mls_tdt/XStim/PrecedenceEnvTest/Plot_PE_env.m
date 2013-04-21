function [] = Plot_PE_env(datamatrix,param3)

%function [] = Plot_PE_env(datamatrix,param3)
%
% datamatrix from Engage_PE_env is converted into DATA herein
% param3 is the listing for each test/rep of the 12 variables (below)
%
% Follows closely on makeAMparams to generate sounds in the same order
%
% [params, n_trial] = makeAMparams;
% makes parameters for Precedence Effect AM test, [#trials x 12]
%
% 1) seed for noise carrier
% 2) seed for decorrelating noise carrier
% 3) stimdur (ms)
% 4) delay (ms) (if delay==0, then uses decorrelated carriers)
% 5) ramp_dur (ms)
% 6) seed for AM envelope
% 7) gammatone cF (Hz)
% 8) envelope decorrelation index (0-1)
% 9) seed for decorrelating envelope#1
% 10) seed for decorrelating envelope#2
% 11) index for sound1 (1 or 2)
% 12) index for sound2 (1 or 2)
%

global XStimParams
global H

% sort and combine datamatrix into DATA 
% (for compatability with later DATA analysis)

% find Ntrials and maxspikes
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

% recover nSpikes for each rep of each test
% params (XStimParams):
Factor_TDT_MII = 1.001654;
Factor_HRIR = 3.1;
delay = XStimParams.delay;
n_delay = length(delay);
cF = XStimParams.cF;
n_cF = length(cF);
nReps = XStimParams.numreps;
nTests = size(param3{1},1);

nSpikes = zeros(nTests,nReps);
for itest=1:nTests
    ind0 = find(DATA(2:end,1) == itest);
    if ~isempty(ind0)
        spiketimes = DATA(ind0+1,DATA(1,1)+2:end)*Factor_TDT_MII;
        reps = DATA(1+ind0, DATA(1,1));
        for jrep = 1:length(reps)
            irep = reps(jrep);
            ind0 = find(spiketimes(jrep,:)>XStimParams.silence_lead & spiketimes(jrep,:)<=XStimParams.silence_lead + XStimParams.curr_stimdur);
            if ~isempty(ind0)
                nSpikes(itest,irep) = length(ind0);
            end
        end
    end
end


% plot DATA
% open figure
if exist1('H.PE_env_fig')
    close(H.PE_env_fig);
    H.PE_env_fig = [];
end
H.PE_env_fig = figure('Units','pixels',...
    'Position',[400 500 72*10 72*5],...
    'Name','PE envelope Data plot',...
    'MenuBar','figure',...
    'ToolBar','figure',...
    'Tag','PE envelope Data plot',...
    'NumberTitle','off',...
    'Color',[1 1 1]*.9 );
xlabel('delay (ms)')
ylabel(['spike rate (/' num2str(XStimParams.curr_stimdur) ' ms)'])
set(gca,'color',[1 1 1]*.8);
hold on

% tests for standard target lead and target lag for each delay and each gammatone cF
linestyle{1} = '-'; linestyle{2} = '--';
itest = 0;
cmap = hsv(length(XStimParams.cF)*2);
tempLegend = cell(n_cF,1);
for iLoc=1:2
    for icF = 1:n_cF
        temp = zeros(n_delay,3);
        for idelay = 1:n_delay
            itest = itest+1;
            M = mom(nSpikes(itest,1:irep));
            temp(idelay,1) = M(1);
            %temp(idelay,2) = M(1) - M(2);
            %temp(idelay,3) = M(1) + M(2);
        end
        % Hplot = errorbar( delay, temp(:,1), temp(:,2), temp(:,3));
        Hplot = plot(delay, temp(:,1));
        set(Hplot,'color', cmap(icF+1,:),'Marker','.','markerSize',18,'linestyle',linestyle{iLoc},'linewidth',2);
        tempLegend{icF} = num2str(cF(icF));
        Hlegend(icF) = Hplot(1);
    end
    if iLoc==1
        Hlegend = legend(Hlegend, tempLegend);
    end
end
   
text(.85*max(get(gca,'xlim')),.7*max(get(gca,'ylim')),' - - lag')
text(.85*max(get(gca,'xlim')),.65*max(get(gca,'ylim')),' \_\_ lead')

% decorr env, decorr carrier (useful for getting spikes vs deltaI and
% spikes vs I/dt curves)
for icF = 1:n_cF
    itest = itest+1;
end

% identical env, decorr carrier (estimate of minimum rate(s))
for icF = 1:n_cF
    itest = itest+1;
    temp = mean(nSpikes(itest,1:irep));
    line([0 delay(end)], [1 1]*temp, 'color', cmap(icF+1,:));
end

% identical env, identical carrier one source in RF (estimate of maximum rate(s))
for icF = 1:n_cF
    itest = itest+1;
    temp = mean(nSpikes(itest,1:irep));
    line([0 delay(end)], [1 1]*temp, 'color', cmap(icF+1,:));
end

% identical env, identical carrier one source in anti-RF (no response?)
for icF = 1:n_cF
    itest = itest+1;    
end
