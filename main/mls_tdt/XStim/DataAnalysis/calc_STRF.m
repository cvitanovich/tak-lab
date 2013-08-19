function [outmat,time,cF] = calc_STRF(spikedata, stimFN)

%function [outmat,time,cF] = calc_STRF(spikedata, stimFN)
% spikedata is a vector of spiketimes in msec
% stim with Fs = 30000
% NO CHECKING FOR EMPTY SPIKEDATA

global FN
global XStimParams
global TDT

if nargin < 2
    disp(['Come on man, ya gotta have spikes and a stim ...']);
    return;
end

% params
Fs = 30000;
cF = 2000:100:10000;
ncF = length(cF);
outputpts= 1024;
ticspermsec = 1000 / Fs;
outputdur = outputpts * ticspermsec;
time = 1000*(-outputpts/Fs:1/Fs:(outputpts-1)/Fs);

lastspike = XStimParams.silence_lead * round(TDT.Fs/1000) *ticspermsec; 					%ms

% remove spike lag due to HRTF convolution
%spikedata = spikedata - 4.1;
% adjust for differing clock rates between MII and TDT (at Fs = 30000)
spikedata = spikedata .* 1.0016;
Nspikes = length(spikedata);

% load stim
stim = r_noi(FN.stim_path, stimFN,'int16');

% condition stimulus
stim = stim - mean(stim);
stim = stim / max(stim);
stimpts = length(stim);
stimdur = stimpts * ticspermsec;

% scaling factor relative to highest freq
maxF = cF(length(cF));
maxFactor = ones(size(cF)) * .00003764 * maxF + .6236;
Factor = maxFactor ./ (.00003764 * cF + .6236);

% filter stim and get envelopes
stim_filt = ERBfilterbankA(stim,cF,Fs);
stim_hilb = abs(hilbert(stim_filt'))';
clear stim_filt

outmat = zeros(ncF,outputpts*2);

% calc STRF
lastevent = stimdur - 5 - outputdur;
goodspikes = 0;
for i = 1:Nspikes
   spiketime = spikedata(i);
   if spiketime <= lastevent & (spiketime >= lastspike + .5)
      tic0 = round((spiketime/ticspermsec) - outputpts);
      data = stim_hilb(:,tic0:(tic0 + (outputpts*2)-1));
      outmat = outmat + data;
      lastspike = spiketime;
      goodspikes = goodspikes+1;
   end  	% <= lastevent
end		% 1:numspikes
outmat = outmat/goodspikes;

