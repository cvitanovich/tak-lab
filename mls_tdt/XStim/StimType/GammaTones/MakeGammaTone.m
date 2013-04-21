function [gammatone] = MakeGammaTone(Fs,freq,dur)
%[gammatone] = MakeGammaTone(Fs,freq,dur)
%MakeGammaTone: 	Create a gammatonal stimulus
%Fs: 		sampling rate in Hz
%freq: 		frequency in Hz
%dur:		duration in ms

if(nargin < 3)
   error('Not enough input arguments');
end
dur = dur/1000;

numpts = Fs * dur;
gammatone = use1_ERBfilt(rand(1,numpts),Fs,freq);
return;