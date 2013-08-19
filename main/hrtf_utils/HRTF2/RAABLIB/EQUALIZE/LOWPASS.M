%
% LOWPASS.M
% Design of a special linear phase lowpass filter.
% passband: f < 15000 Hz
% stopband: f > 18500 Hz
% Frequency values for 44100 Hz sampling rate

L=2048;
fa=samplfrq;
ymin=1E-7;
ymax=1.00;
N=46;
pass=[0,15000];
stop=[18500,22100];

s='-- LINEAR PHASE EQUIRIPPLE LOWPASS FILTER --';
disp(s);
disp(' ');
disp(sprintf('passband: f<%g Hz, stopband: f>%g Hz, N=%i',pass(2),stop(1),N));

mstop=ymin*ones(size(stop));
mpass=ymax*ones(size(pass));

f=[pass stop]*2/fa;
m=[mpass mstop];

b=remez(N,f,m);

H=freqz(b,1,L);

f1=0+1;
f2=round(2*L*pass(2)/fa)+1;
f3=round(2*L*stop(1)/fa)+1;
f4=L;;
max_pb=max(abs(H(f1:f2)));
min_pb=min(abs(H(f1:f2)));
max_sb=max(abs(H(f3:f4)));
disp(sprintf('max. magnitude passband: %g dB',20*log10(max_pb)));
disp(sprintf('min. magnitude passband: %g dB',20*log10(min_pb)));
disp(sprintf('max. magnitude stopband: %g dB',20*log10(max_sb)));

magplot(H,samplfrq,s,'r');
grid;

info='LOWPASS FIR FILTER COEFFICIENTS';
save lp15k b info
