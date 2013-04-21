function [gtonestackL,gtonestackR] = MakeToneStack(Fs,tonestackfreqs,stim_dur,targ_freq,ILD)
%[tonestackl,tonestackr] = MakeGammaToneStack(Fs,tonestackfreqs,stim_dur,targ_freq,ILD)
%Fs, sampling rate in Hz
%tonestackfreqs, a vector containing the frequencies to be included
%stim_dur, the stimulus duration in ms

dur = stim_dur/1000;
numpts = dur*Fs;

for freqnum = 1:length(tonestackfreqs)
   randit = rand(1);
   if(randit <= 0.5) randit = -randit; end
   newfreqs(freqnum) = round(((tonestackfreqs(freqnum) + (10*randit) + 1)));
end


for freqnum = 1:length(newfreqs)
   gammatone(freqnum,:) = use1_ERBfilt(rand(1,numpts),Fs,newfreqs(freqnum));
   gammatone(freqnum,:) = gammatone(freqnum,:)/max(abs(gammatone(freqnum,:)));
end

if(nargin > 3)
   mult_fact = 10^(ILD/20);
   gammatone_targ = use1_ERBfilt(rand(1,numpts),Fs,targ_freq);
   gammatone_targ = gammatone_targ/max(abs(gammatone_targ));
   gammatone_left = gammatone_targ/sqrt(mult_fact);
   gammatone_right = sqrt(mult_fact) * gammatone_targ;
   gtonestackL = sum([gammatone;gammatone_left],1);
   gtonestackR = sum([gammatone;gammatone_right],1);
else
   gtonestackL = sum(gammatone,1);
   gtonestackR = sum(gammatone,1);
end

return;
   
