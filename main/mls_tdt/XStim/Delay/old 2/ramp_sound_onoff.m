function [ramped_sound] = ramp_sound_onoff(var_name,Fs,ramp_time_on, ramp_time_off)
%Ramp_Sound, [ramped_sound] = ramp_sound_onoff(var_name,Fs,ramp_time)
%var_name:		sound to be ramped
%Fs:			sampling rate (Hz)
%ramp_time_on:		onset time in ms
%ramp_time_off:		offset time in ms

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
num_pts_on = round(ramp_time_on/samp_per);
num_pts_off = round(ramp_time_off/samp_per);
ramp(1:num_pts_on) = 0:1/num_pts_on:1-(1/num_pts_on);
ramp(lenvarname - num_pts_off+1:lenvarname) = 1-(1/num_pts_off):-1/num_pts_off:0;

ramped_sound = var_name .* ramp;

return