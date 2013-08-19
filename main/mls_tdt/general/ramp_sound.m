function [ramped_sound] = ramp_sound(var_name,Fs,ramp_time)
%Ramp_Sound, [ramped_sound] = ramp_sound(var_name,Fs,ramp_time)
%var_name:		sound to be ramped
%Fs:			sampling rate (Hz)
%ramp_time:		time in ms

if(~isreal(var_name))
   error('Sound variable must contain a real array');
else
   [r,c] = size(var_name);
   if(r > 1 & c > 1)
      error('Sound variable must be a 1-D arrary');
   end
end

lenvarname = length(var_name);
ramp = ones(size(var_name));
samp_per = 1000/Fs; %in ms
num_pts = round(ramp_time/samp_per);
ramp(1:num_pts) = 0:1/num_pts:1-(1/num_pts);
ramp(lenvarname - num_pts+1:lenvarname) = 1-(1/num_pts):-1/num_pts:0;

ramped_sound = var_name .* ramp;

return