close all
tonestackfreqs = 2000:1000:11000;
Fs = 30000;
stim_dur = 2000;
stim_len = round((stim_dur/1000)*Fs);
minfreq = 2000; maxfreq = 11000;
minfreq = round(((minfreq + 1)/Fs) * stim_len);
maxfreq = round(((maxfreq + 1)/Fs) * stim_len);
range = maxfreq-minfreq+1;

% mag spectrum = 1 between set frequencies:
fftlen = 4096;
mag = zeros(1,fftlen);
minfreq = 2000; maxfreq = 11000;
minfreq = round(((minfreq + 1)/(Fs)) * fftlen);
maxfreq = round(((maxfreq + 1)/(Fs)) * fftlen);
range = maxfreq-minfreq+1;
mag(minfreq:maxfreq) = ones(1,range);

tonestack = MakeToneStack(Fs,tonestackfreqs,stim_dur);

psd(tonestack,256,Fs);
figure
plot(tonestack)

tonestackflat = flatten(tonestack,0);
load e:\spezio\matlab\scripts\tdt\XStim\HRTFfilts\bp2to11
tonestacknew = conv(tonestackflat,filt1.tf.num);

figure
psd(tonestackflat,256,Fs);
figure
plot(tonestackflat)


figure
psd(tonestacknew,2048,Fs);
figure
plot(tonestacknew)
