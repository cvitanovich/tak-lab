%
% TGRPDEMO shows the influence of noise to the computation
%	   of the group delay.

fftlen=2048;
samplerate=44100;
h=[0.1 0.2 0.3 0.2 0.1 0.0 -0.1 -0.2 0.1 0];
h=[h,zeros(1,fftlen-length(h))];
n=1E-4*randn(1,fftlen);
hn=h+n;

En=energy(n);
Eh=energy(h);
disp(sprintf('SNR=%gdB',20*log10(Eh/En)));
subplot(2,1,1);
plot(h(1:50),'y');
subplot(2,1,2);
plot(hn(1:50),'r');

H=fft(h,fftlen); H=H(1:fftlen/2);
Hn=fft(hn,fftlen); Hn=Hn(1:fftlen/2);
tgrpplot([Hn(:),H(:)],samplerate);



