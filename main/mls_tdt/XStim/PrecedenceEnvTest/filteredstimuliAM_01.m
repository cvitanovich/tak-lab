function [Noise1, Noise2, Env1, Env2] = filteredstimuliAM_01(params, itrial)
    %[Noise1, Noise2, Env1, Env2] = filteredstimuliAM_01(params, itrial);
    
    % stimuli to be spatialized separately and then summed:
    % Noise1 = lead stimulus
    % Noise2 = lag stimulus
    %
    % their envelopes:
    % Env1 = lead envelope
    % Env2 = lag envelope
    
    % params = stimulus parameters
    % itrial = trial number 
    % ATTN: added flatten of Noise1 and Noise2 (Jun26_12)

    Fs = 30000; % may need to adjust derivative function if Fs is changed
                % 48828 30000
                % set separately in function SpikeProbabilitiesAM !!!
    
    % Stimulus variables from params  
    %NoiseState = params(itrial,1);                  % seed for noise carrier
    %DecorrNoiseState = params(itrial,2);            % seed for decorrelating noise carrier
    stim_dur = params(itrial,3);                    % stimulus duration (ms)
    stim_pts = stim_dur * Fs/1000;                  % total points in final stimulus
    Delay = params(itrial,4);                       % delay (ms)
    ramp_dur = params(itrial,5);                    % ramp duration (ms)
    %AMState = params(itrial,6);                     % seed for (AM) envelope
    %CF = params(itrial,7);                          % center freq for gammatone (Hz)
    %DecorrIndex = params(itrial,8);                 % 0) envelopes correlated >0) envelopes decorrelated by this
    %DecorrEnvState1 = params(itrial,9);             % seed for 1st decorrelating envelope
    %DecorrEnvState2 = params(itrial,10);            % seed for 2nd decorrelating envelope
    
    FilterRampBuffer = 30; % ms, to avoid filter effects at onset & offset
    FilterRampBufferPnts = round((FilterRampBuffer/1000)*Fs);
        
    % Increase duration by delay so that the stimuli can be gated/windowed
    StimDur = stim_dur + Delay + (FilterRampBuffer*2); % add delay
    StimPnts = round((StimDur/1000)*Fs);
    
    DelayPnts = round((Delay/1000)*Fs);
              
    %// make noise for envelope
    randn('state',params(itrial,6)); % use state, randn method
    EnvNoise1 = randn(1,StimPnts);
    %// normalize the carrier's amplitude
    %/ final envelope amplitude's can be normalized *after* decorrelation
    rms = norm(EnvNoise1)/sqrt(length(EnvNoise1));
    EnvNoise1 = EnvNoise1 * 1/rms;

    if params(itrial,8) > 0
        
        randn('state',params(itrial,9)); % use state, randn method
        DecorrNoise1 = randn(1,StimPnts);
        rms = norm(DecorrNoise1)/sqrt(length(DecorrNoise1));
        DecorrNoise1 = DecorrNoise1 * 1/rms;
       
        randn('state',params(itrial,10)); % use state, randn method
        DecorrNoise2 = randn(1,StimPnts);
        rms = norm(DecorrNoise2)/sqrt(length(DecorrNoise2));
        DecorrNoise2 = DecorrNoise2 * 1/rms;
       
        EnvNoise1 = ((1-params(itrial,8)).*EnvNoise1) + ((params(itrial,8)).*DecorrNoise1);
        EnvNoise2 = ((1-params(itrial,8)).*EnvNoise1) + ((params(itrial,8)).*DecorrNoise2);
                
    end
    
    EnvNoise1 = [zeros(1,500) EnvNoise1 zeros(1,2000)]; % pad with zeros
    EnvNoise1 = ERBFilter(Fs, params(itrial,7), EnvNoise1); % Filter (Hz, Hz, array)
    EnvNoise1 = EnvNoise1(500:end); % remove front zero pad
    Env1 = hilbert(EnvNoise1); % hilbert transform
    Env1 = sqrt(Env1.*conj(Env1)); % amplitude
    %// remove DC
    Env1 = Env1 - Env1(1); % adjust DC so that the first point is zero
    Env1 = max(Env1, 0); % zero-out any negative points near the end of the sound
    
    if params(itrial,8) > 0
        EnvNoise2 = [zeros(1,500) EnvNoise2 zeros(1,2000)]; % pad with zeros
        EnvNoise2 = ERBFilter(Fs, params(itrial,7), EnvNoise2); % Filter (Hz, Hz, array)
        EnvNoise2 = EnvNoise2(500:size(EnvNoise2,2)); % remove front zero pad
        Env2 = hilbert(EnvNoise2); % hilbert transform
        Env2 = sqrt(Env2.*conj(Env2)); % amplitude
        %// remove DC
        Env2 = Env2 - Env2(1); % adjust DC so that the first point is zero
        Env2 = max(Env2, 0); % zero-out any negative points near the end of the sound
        
    else
        Env2 = Env1;
        EnvNoise2 = EnvNoise1;
    end
    
    % High frequency carriers
    %// make noise carriers,  iFFT method
    minHz = 2000;
    maxHz = 11000;
    NoiseDur = stim_dur+Delay;
    NoisePnts = round((NoiseDur/1000)*Fs);
    Noise1 = GetNoise_BB(params(itrial,1), Fs, NoiseDur, minHz, maxHz);
    Noise1 = flatten(Noise1);
    
    % option for Decorrated carriers
    if Delay == 0 % carriers are completely uncorrelated
        Noise2 = GetNoise_BB(params(itrial,2), Fs, NoiseDur, minHz, maxHz);
    else
        Noise2 = Noise1;
        Noise2 = flatten(Noise2);
    end
    
    % Delay to make "lead" and "lag"
    Env1 = [Env1 zeros(1,DelayPnts)];
    Env2 = [zeros(1,DelayPnts) Env2];
    EnvNoise1 = [EnvNoise1 zeros(1,DelayPnts)];
    EnvNoise2 = [zeros(1,DelayPnts) EnvNoise2];
    
    Noise1 = [Noise1 zeros(1,DelayPnts)];
    Noise2 = [zeros(1,DelayPnts) Noise2];
    
    DelayPnts = max(1,DelayPnts);
        
    % window (i.e., gate) out onset/offset time disparity
    % AND filter "ramps"
    Env1 = Env1(DelayPnts+FilterRampBufferPnts:StimPnts-FilterRampBufferPnts-1);
    Env2 = Env2(DelayPnts+FilterRampBufferPnts:StimPnts-FilterRampBufferPnts-1);
    EnvNoise1 = EnvNoise1(DelayPnts+FilterRampBufferPnts:StimPnts-FilterRampBufferPnts-1);
    EnvNoise2 = EnvNoise2(DelayPnts+FilterRampBufferPnts:StimPnts-FilterRampBufferPnts-1);
        
    % window (i.e., gate) out onset/offset time disparity
    Noise1 = Noise1(DelayPnts:NoisePnts-1);
    Noise2 = Noise2(DelayPnts:NoisePnts-1);
    
    % make sure the two different signal generation methods didn't create
    % signal having different lengths -- due to roundoff errors
    if length(Noise1) > length(Env1)
        Noise1 = Noise1(1:length(Noise1)-1);
        Noise2 = Noise2(1:length(Noise2)-1);
    end
    if length(Env1) > length(Noise1)
        Env1 = Env1(1:length(Env1)-1);
        Env2 = Env2(1:length(Env2)-1);
        EnvNoise1 = EnvNoise1(1:length(EnvNoise1)-1);
        EnvNoise2 = EnvNoise2(1:length(EnvNoise2)-1);
    end
       
    % ramp envelopes on/off
    Env1 = ramp_sound_onoff(Env1, Fs, ramp_dur, ramp_dur); % stim envelope
    Env2 = ramp_sound_onoff(Env2, Fs, ramp_dur, ramp_dur); % stim envelope
    
    %// normalize final envelopes
    maxEnv = max([max(Env1) max(Env2)]);
    EnvNoise1 = EnvNoise1 * 1/maxEnv;
    EnvNoise2 = EnvNoise2 * 1/maxEnv;
    Env1 = Env1 * 1/maxEnv;
    Env2 = Env2 * 1/maxEnv;
        
    % multiply with envelopes
    Noise1 = Noise1 .* Env1;
    Noise2 = Noise2 .* Env2;
    
    EnvNoise1 = EnvNoise1 .* Env1;
    EnvNoise2 = EnvNoise2 .* Env2;
    
    Noise1 = Noise1 - mean(Noise1);
    Noise2 = Noise2 - mean(Noise2);
  
    % make sure stim_pts is correct
    temp = length(Env1) - stim_pts;
    if temp>0
        Env1=Env1(1:stim_pts);
        Env2=Env2(1:stim_pts);
        Noise1=Noise1(1:stim_pts);
        Noise2=Noise2(1:stim_pts);
    elseif temp<0
        Env1=[Env1 zeros(-temp,1)];
        Env2=[Env2 zeros(-temp,1)];
        Noise1=[Noise1 zeros(-temp,1)];
        Noise2=[Noise2 zeros(-temp,1)];
    end        
    
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
    
 

    
