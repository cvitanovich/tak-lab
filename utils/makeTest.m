function stim = makeTest(seedval,dur,minband,maxband,samprate,ramp,base_rms)
% makes ramped test sounds using specific parameters

rampPts = ceil((ramp/1000)*samprate);
[temp,stim] = seeded_whnoise(minband,maxband,samprate,dur/1000,seedval);
stim = stim';  %Stimuli are scaled between 1 and -1 in whnoise and not rescaled here; 
%MUST BE A ROW VECTOR!!!!
npts = length(stim);
% put envelope on it
env = [0:1/(rampPts-1):1 ones(1,npts - 2*rampPts) 1:-1/(rampPts-1):0];
stim = stim .* env;
% scale to desired rms:
rms = mean(stim.^2)/sqrt(2);
ratio = rms/base_rms;
stim = stim/sqrt(ratio);