function [tone] = MakeTone(Fs,freq,dur,phi)
%MakeTone: 	Create a tonal stimulus with random (or specified) phase
%Fs: 		sampling rate in Hz
%freq: 		frequency in Hz
%dur:		duration in ms

if(nargin < 3)
   error('Not enough input arguments');
end
samp_period = 1/(Fs-1);
dur = dur/1000;
t = 0:samp_period:dur;
rand('state',sum(100*clock));
if nargin < 4 | isempty(phi)
    phi = (rand(1) - 0.5) * (2*pi);
end

tone = sin(2*pi*freq*t + phi);
tone = tone/mom(tone,2);