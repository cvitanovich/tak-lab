function [tonestackL,tonestackR] = MakeToneStackmissing(Fs,tonestackfreqs,stim_dur,targ_freq)
%[tonestackl,tonestackr] = MakeToneStack(Fs,tonestackfreqs,stim_dur,targ_freq,ILD)
%Fs, sampling rate in Hz
%tonestackfreqs, a vector containing the frequencies to be included
%stim_dur, the stimulus duration in ms

stim_dur = stim_dur/1000;
stim_len = round(stim_dur*Fs);
mag = zeros(1,stim_len);

tonestackfreqs = tonestackfreqs(find(tonestackfreqs >= 2000));

% mag spectrum = 1 at specified frequencies:
for freqnum = 1:length(tonestackfreqs)
   if(tonestackfreqs(freqnum) ~= targ_freq)
      randit = rand(1);
      if(randit <= 0.5) randit = -randit; end
      freqpos = round(((tonestackfreqs(freqnum) + (10*randit) + 1)/Fs) * stim_len);
      %freqpos = round(((tonestackfreqs(freqnum) + 1)/Fs) * stim_len);
      mag(freqpos) = 10000;
   end
end

% random phase spectrum between set frequencies:
rand('state',sum(100*clock));
phi = (rand(1,stim_len) - 0.5) * (2*pi);
% combine phase and magnitude:
X = mag .* ( (cos(phi)) + (i .* sin(phi)) );

% convert to time domain:
tonestack = real(ifft(X));
tonestackL = tonestack/max(abs(tonestack));
tonestackR = tonestack/max(abs(tonestack));

