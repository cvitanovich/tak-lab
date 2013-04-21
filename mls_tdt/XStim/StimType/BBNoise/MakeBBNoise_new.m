function [bbnoise] = MakeBBNoise_new(Fs,dur)
%MakeBBNoise:	Create a BroadBand Noise (2-11 kHz)
%Fs:		Sampling rate (Hz)
%dur:		Stimulus duration (ms)

dur = dur/1000;

minfreq_noise = 2000;
maxfreq_noise = 11000;
am_freq = 55;

mod_depth = 0;
bbnoise = amnoise(minfreq_noise,maxfreq_noise,am_freq,mod_depth,Fs,dur);

return;