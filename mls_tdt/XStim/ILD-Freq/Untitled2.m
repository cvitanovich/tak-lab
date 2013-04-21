clear
Fs = 30000;
t = 0:1/(Fs-1):0.2;
tone1 = sin(2*pi*2000*t + (2*pi*rand(1)));
amp1 = sqrt(max(psd(tone1,256,Fs)));
tone2 = sin(2*pi*2000*t + (2*pi*rand(1)));
amp2 = sqrt(max(psd(tone2,256,Fs)));

ild = -10;
mult_fact = 10^(ild/20);
tone2new = sqrt(mult_fact) * tone2;
amp2new = sqrt(max(psd(tone2new,256,Fs)));
tone1new = tone1/sqrt(mult_fact);
amp1new = sqrt(max(psd(tone1new,256,Fs)));

newild = 20*log10(amp2new/amp1new);