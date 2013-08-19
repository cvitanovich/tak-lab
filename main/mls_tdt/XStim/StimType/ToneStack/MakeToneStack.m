function [tonestackL,tonestackR] = MakeToneStack(Fs,tonestackfreqs,stim_dur,targ_freq,ILD)
%[tonestackl,tonestackr] = MakeToneStack(Fs,tonestackfreqs,stim_dur,targ_freq,ILD)
%Fs, sampling rate in Hz
%tonestackfreqs, a vector containing the frequencies to be included
%stim_dur, the stimulus duration in ms

stim_dur = stim_dur/1000;
stim_len = round(stim_dur*Fs);
mag = zeros(1,stim_len);

% mag spectrum = 1 at specified frequencies:
for freqnum = 1:length(tonestackfreqs)
   randit = rand(1);
   if(randit <= 0.5) randit = -randit; end
   freqpos = round(((tonestackfreqs(freqnum) + (10*randit) + 1)/Fs) * stim_len);
   %freqpos = round(((tonestackfreqs(freqnum) + 1)/Fs) * stim_len);
   mag(freqpos) = 10000;
end

if(nargin > 3)
   magL = mag; magR = mag;
   freqpos = round(((targ_freq + 1)/Fs) * stim_len);
   %Make tone at specified frequency and ild - NOT GETTING ILD by ATTENUATORS
   %Lright - Lleft = 20*log10(Ampright/Ampleft), Hartmann p. 29
   mult_fact = 10^(ILD/20);
   magL(freqpos) = 10000/sqrt(mult_fact);
   magR(freqpos) = 10000*sqrt(mult_fact);
   
   % random phase spectrum between set frequencies:
   rand('state',sum(100*clock));
   phi = (rand(1,stim_len) - 0.5) * (2*pi);
   % combine phase and magnitude:
   XL = magL .* ( (cos(phi)) + (i .* sin(phi)) );
   XR = magR .* ( (cos(phi)) + (i .* sin(phi)) );
   
   % convert to time domain:
   tonestackL = real(ifft(XL));
   tonestackR = real(ifft(XR));
else
   % random phase spectrum between set frequencies:
   rand('state',sum(100*clock));
   phi = (rand(1,stim_len) - 0.5) * (2*pi);
   % combine phase and magnitude:
   X = mag .* ( (cos(phi)) + (i .* sin(phi)) );
   
   % convert to time domain:
   tonestack = real(ifft(X));
   tonestackL = tonestack/max(abs(tonestack));
   tonestackR = tonestack/max(abs(tonestack));
end

