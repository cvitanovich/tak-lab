function [Noise1, Noise2, Env1, Env2] = filteredstimuliAM_02(params, itrial)
%[Noise1, Noise2, Env1, Env2] = filteredstimuliAM_02(params, itrial);

% stimuli to be spatialized separately and then summed:
% Noise1 = lead stimulus
% Noise2 = lag stimulus
%
% their envelopes:
% Env1 = lead envelope
% Env2 = lag envelope

% params = stimulus parameters
% itrial = trial number

% trying somethings Jun 2012 - more BN-like sounds (not gammatoned)
% ATTN: added flatten of Noise1 and Noise2 (Jun26_12)

Fs = 30000; % may need to adjust derivative function if Fs is changed
% set separately in function SpikeProbabilitiesAM !!!

% Stimulus variables from params
% 1) seed for noise carrier
% 2) seed for decorrelating noise carrier
% 3) stimulus duration (ms)
% 4) delay (ms)
% 5) ramp duration (ms)
% 6) seed for envelope
% 7) center freq for gammatone (Hz)
% 8) 0) envelopes correlated >0) envelopes decorrelated by this
% 9) seed for 1st decorrelating envelope
% 10) seed for 2nd decorrelating envelope

% Increase duration by delay so that the stimuli can be gated/windowed
Delay = params(itrial,4);                       % delay (ms)
ramp_dur = params(itrial,5);                    % ramp duration (ms)
StimDur = params(itrial,3);
StimPnts = round((StimDur/1000)*Fs);
DelayPnts = round((Delay/1000)*Fs);
AMBandwidth = params(itrial,7);
AMDepth = 100;

% Make Envelope #1
Env1 = GetNoise_BB(params(itrial,6), Fs, 300, 1, AMBandwidth);
Env1 = Env1(1:StimPnts);
Env1 = Env1 - min(Env1); Env1 =  Env1 / max(Env1);  % normalize
Env1 = Env1 .* (AMDepth/100); % AM depth (as a percent)
Env1 = Env1 + 1 - (AMDepth/100);

if params(itrial,8) > 0
    DecorrNoise1 = GetNoise_BB(params(itrial,9), Fs, 300, 1, AMBandwidth);
    DecorrNoise1 = DecorrNoise1(1:StimPnts);
    DecorrNoise1 = DecorrNoise1 - min(DecorrNoise1);
    DecorrNoise1 =  DecorrNoise1 / max(DecorrNoise1);  % normalize
    
    DecorrNoise2 = GetNoise_BB(params(itrial,10), Fs, 300, 1, AMBandwidth);
    DecorrNoise2 = DecorrNoise2(1:StimPnts);
    DecorrNoise2 = DecorrNoise2 - min(DecorrNoise2);
    DecorrNoise2 =  DecorrNoise2 / max(DecorrNoise2);  % normalize
    
    Env1 = ((1-params(itrial,8)).*Env1)+((params(itrial,8)).*DecorrNoise1);
    Env2 = ((1-params(itrial,8)).*Env1)+((params(itrial,8)).*DecorrNoise2);
    
else
    Env2 = Env1;
end

% carriers
Noise1 = GetNoise_BB(params(itrial,1), Fs, StimDur, 2000, 11000);
Noise1 = flatten(Noise1);

% option for Decorrated carriers
if Delay == 0 % carriers are completely uncorrelated
    Noise2 = GetNoise_BB(params(itrial,2), Fs, StimDur, 2000, 11000);
    Noise2 = flatten(Noise2);
else
    Noise2 = Noise1;
end

% ramp envelopes on/off
Env1 = ramp_sound_onoff(Env1, Fs, ramp_dur, ramp_dur); % stim envelope
Env2 = ramp_sound_onoff(Env2, Fs, ramp_dur, ramp_dur); % stim envelope

%// normalize final envelopes
maxEnv = max([max(Env1) max(Env2)]);
Env1 = Env1 * 1/maxEnv;
Env2 = Env2 * 1/maxEnv;

% Delay to make "lead" and "lag"
Env1 = [Env1 zeros(1,DelayPnts)];
Env2 = [zeros(1,DelayPnts) Env2];

% multiply with envelopes
Noise1 = [Noise1 zeros(1,DelayPnts)];
Noise2 = [zeros(1,DelayPnts) Noise2];
Noise1 = Noise1 .* Env1;
Noise2 = Noise2 .* Env2;

Noise1 = Noise1 - mean(Noise1);
Noise2 = Noise2 - mean(Noise2);

function [bbnoise] = GetNoise_BB(state, Fs, dur, minfreq, maxfreq)
%GetNoise_BB:	Create a BroadBand Noise (2-11 kHz)
%state: rand state
%Fs: Sampling rate (Hz)
%dur: Stimulus duration (ms)
dur = dur/1000;
len = round(dur*Fs);
minfreq = round(((minfreq + 1)/Fs) * len);
maxfreq = round(((maxfreq + 1)/Fs) * len);
range = maxfreq-minfreq+1;
% mag spectrum = 1 between set frequencies:
mag = zeros(len,1);
mag(minfreq:maxfreq) = ones(range,1);
% random phase spectrum between set frequencies:
rand('state',state); % use state
phi = (rand(len,1) - 0.5) * (2*pi);
% combine phase and magnitude:
X = mag .* ( (cos(phi)) + (i .* sin(phi)) );
% convert to time domain:
bbnoise = real(ifft(X));
% scale to RMS = 0.23591 !!!!!!!!
rms = norm(bbnoise)/sqrt(length(bbnoise))';
bbnoise = bbnoise * (0.23591/rms);          % 0.23591 equals old code
bbnoise = bbnoise';


function [ramped_sound] = ramp_sound_onoff(vals,Fs,ramp_time_on, ramp_time_off)
% Ramp_Sound, [ramped_sound] = ramp_sound_onoff(vals,Fs,ramp_time)
% vals: sound to be ramped
% Fs: sampling rate (Hz)
% ramp_time_on: onset time in ms
% ramp_time_off: offset time in ms
len = length(vals);
ramp = ones(size(vals));
samp_per = 1000/Fs; %in ms
num_pts_on = round(ramp_time_on/samp_per);
num_pts_off = round(ramp_time_off/samp_per);
ramp(1:num_pts_on) = 0:1/num_pts_on:1-(1/num_pts_on);
ramp(len - num_pts_off+1:len) = 1-(1/num_pts_off):-1/num_pts_off:0;
ramped_sound = vals .* ramp;




