function ramped_snd = rampMySound(snd,ramplen,Fs)
ramplen = ramplen*Fs;
len = length(snd);
ramp = 0:(1/ramplen):1;
flip_ramp = fliplr(ramp);
env = [ramp ones(1,len - 2*length(ramp)) flip_ramp];
ramped_snd = env.*snd;
