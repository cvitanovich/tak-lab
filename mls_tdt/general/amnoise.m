function [td,time_matrix] = amnoise(minfreq_noise,maxfreq_noise,am_freq,mod_depth,Fs,dur,phas,seed)
%Script amnoise to produce amplitude-modulated noise, modulated at a specific frequency
%[td,time_matrix] = amnoise(minfreq_noise,maxfreq_noise,am_freq,mod_depth,Fs,dur,phas,seed)
%minfreq_noise: minimum frequency for noise
%maxfreq_noise: maximum frequency for noise
%amfreq: amplitude modulation frequency
%mod_depth: modulation depth (0 to 1)
%Fs: sampling frequency (Hz)
%dur: duration in seconds
%phase: 0-2*pi
%seed: seed for rand in whnoise

if nargin<7 | isempty(phas) phas=0; end
if nargin <8 | isempty(seed) seed = round(rand(1)*2^30);  end

timeVec = 0:dur/(Fs-1):dur;

%Noise signal
[fd,noisesig] = whnoise(minfreq_noise,maxfreq_noise,Fs,dur,seed);

%Modulation signal
modsig = sin(2*am_freq*timeVec*pi + phas);

%Modulate the noise
modenv = (1 - mod_depth/2) + (mod_depth/2)*modsig;
td = noisesig(:) .* modenv(:);