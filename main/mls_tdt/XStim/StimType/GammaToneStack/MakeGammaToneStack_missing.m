function [gtonestackL,gtonestackR] = MakeToneStack_missing(Fs,tonestackfreqs,stim_dur,targ_freq)
%[tonestackl,tonestackr] = MakeGammaToneStack_missing(Fs,tonestackfreqs,stim_dur,targ_freq)
%Fs, sampling rate in Hz
%tonestackfreqs, a vector containing the frequencies to be included
%stim_dur, the stimulus duration in ms

dur = stim_dur/1000;
numpts = dur*Fs;

for freqnum = 1:length(tonestackfreqs)
   if(tonestackfreqs(freqnum) ~= targ_freq)
      randit = rand(1);
      if(randit <= 0.5) randit = -randit; end
      newfreqs(freqnum) = round(((tonestackfreqs(freqnum) + (10*randit) + 1)));
   end
end


for freqnum = 1:length(newfreqs)
   gammatone(freqnum,:) = use1_ERBfilt(rand(1,numpts),Fs,newfreqs(freqnum));
   gammatone(freqnum,:) = gammatone(freqnum,:)/max(abs(gammatone(freqnum,:)));
end

gtonestackL = sum(gammatone,1);
gtonestackR = sum(gammatone,1);


return;
   
