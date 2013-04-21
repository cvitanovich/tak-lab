function [ramped_sound] = ramp_sound(var_name,ramptime,Fs)
%Ramp_Sound:	Function to ramp a sound using APOS
%var_name:		Must be a 1-D array
%ramptime:		Ramp time in ms

if(~isreal(var_name))
   error('Sound variable must contain a real array');
else
   [r,c] = size(var_name);
   if(r > 1 & c > 1)
      error('Sound variable must be a 1-D arrary');
   end
end

samp_per = 1e6/Fs;
while (S232('APactive')) pause(0); end
S232('pushf',var_name,length(var_name));
while (S232('APactive')) pause(0); end
S232('qwind',ramptime,samp_per);
while (S232('APactive')) pause(0); end
ramped_sound = S232('popf');

return