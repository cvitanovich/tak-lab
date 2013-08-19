function [phase, gain, mod_per] = pl_longSAMs(DATA, modFreq, overplot)

%function [phase, gain, mod_per] = pl_longSAMs(DATA, modFreq, overplot)
%
% argins:   DATA
%           modFreq
%           overplot 0) no, 1) yes
%
% argouts:  phase: 1st spike latency, Vphase (from Vstr calc), modalPhase
%           gain: firing rate (sp/sec), Vstr
%           mod_per (msec)

global XStimParams
global H

if nargin < 1 || isempty(DATA)
    error('pl_longSAMs requires DATA input')
end

if nargin < 2 || isempty(modFreq)
    modFreq = XStimParams.mod_freq(1);
end

if nargin < 3 || isempty(overplot)
    overplot = 0;
end

spikelatency = 10;

nparams = DATA(1,1);
nreps = max1(DATA(:,nparams));

mod_per = 1000/modFreq;
bins = 0:mod_per/360:mod_per - mod_per/360;
perhist = zeros(size(bins));
nper = XStimParams.curr_stimdur/mod_per;
%onsetNpers = 0;        % # of modpers to ignore at start
onsetNpers = max1([1 ceil(1000/mod_per)]);

factorMII = 1.001654;   % TDT and MII clocks differ
factorHRTF = 0;       % sound comes on ca 1.5 msec after spikeclock starts

phas = [];
Nspikes = 0; latency = 0;
for irep = 1:nreps
    spiketimes = DATA(1+irep,nparams+2:end) *factorMII + factorHRTF - XStimParams.silence_lead - spikelatency;
    spiketimes = spiketimes(spiketimes>0);
    latency = latency+spiketimes(1);
    spiketimes = spiketimes(spiketimes>mod_per*onsetNpers);
    Nspikes = Nspikes + length(spiketimes);
    pertime = mod(spiketimes,mod_per);
    phas = [phas 2 * pi * pertime/mod_per];       % in radians
    perhist = perhist + histc(pertime,bins);
end

%calc first spikelatency
phase(1) = latency/nreps;           % argout1

% calc firing rate
secsperrep = XStimParams.curr_stimdur/1000;
if onsetNpers
    secsperrep = secsperrep - (1/(modFreq*onsetNpers));
end
firingRate = Nspikes / (secsperrep * nreps);
gain(1) = firingRate;               % argout2

% calc vector strength
ycoord = mean(sin(phas(:)));
xcoord = mean(cos(phas(:)));
Vstr = roundn(sqrt(ycoord ^2 + xcoord^2),.01);
gain(2) = Vstr;                     % argout2

% calc best Vphase
temp = atan2(ycoord,xcoord);
if temp < 0
    temp = temp+ 2*pi;
end
ThetaRad = temp - (pi/2);
ThetaDeg = round(ThetaRad/(2*pi) * 360);
phase(2) = ThetaRad;                % argout1

temp = smooth([perhist perhist perhist]/max1(perhist),4);
perhist_sm = temp(361:2*360);

if(exist1('H.LongSAMsfig') & ~exist1('H.LongSAMs_PERIODfig'))
    H.LongSAMs_PERIODfig = figure('Position',[700 20 550 500],...
        'Name','LongSAMs Period Histogram',...
        'NumberTitle','off');
    H.LongSAMs_PERIODaxes = axes;
end
figure(H.LongSAMs_PERIODfig); 

if overplot
    hold on
    plot(perhist/max1(perhist),'g')
    plot(perhist_sm,'y','linewidth',2)
    plot([1 1]*mod(ThetaDeg+90,360),[-.2 1],'g')
    text(10,.65,['Vstr:  ' num2str(Vstr) '  mod period = ' num2str(mod_per) ' msec'],'color','g','fontsize',12)
    text(10,.85,['Vphase: ' num2str(ThetaDeg) ' deg    == ' num2str(roundn(ThetaDeg/360,.01)) ' cycles'],'color','g','fontsize',12)
else
    hold off
    plot(perhist/max1(perhist),'r')
    hold on
    plot(perhist_sm,'m','linewidth',2)
    plot([1 1]*mod(ThetaDeg+90,360),[-.2 1],'r')
    text(10,.7,['Vstr:  ' num2str(Vstr) '  mod period = ' num2str(mod_per) ' msec'],'color','r','fontsize',12)
    text(10,.9,['Vphase: ' num2str(ThetaDeg) ' deg    == ' num2str(roundn(ThetaDeg/360,.01)) ' cycles'],'color','r','fontsize',12)    
end

set(H.LongSAMs_PERIODaxes,'Color','black');
xlabel('period (deg)'); ylabel('normalized firing rate');

plot((1+sin(2*pi*bins/mod_per - pi/2))/2,'b','linewidth',2)
axis([0 360 0 1])
title(['period histogram after ' num2str(nreps) ' reps;  Nspikes = ' num2str(Nspikes)])