function remdemo(f);

% 300 <= f < 15000

% f = 300 => n = 344;
% f = 15000 => n = 15;

Fs = 44100;

selected_frq= f;



%df=df_erb(f)
df = df_zwick(f)
% 3 dB corner frequency x=1.6746 passband

mp = 1;
rp = 0.5;
fp = 1.6746/4*df



% 50 dB corner frequency x=1.6746 stopband

ms = 0;
rs = 50;
fs = 14.236/4*df



f=[-fs -fp fp fs]
f=f+selected_frq

if any(f>Fs/2) f(f>Fs/2)=Fs/2; end;
if any(f<0) f(f<0)=10; end;

rp=(10^(rp/20)-1)/(10^(rp/20)+1);
rs=10^(-rs/20);

dev = [ rs rp  rs];
m= [ms mp ms];

[n,f0,m0,w]=remezord(f,m,dev,Fs)

b=remez(n,f0,m0,w);

[h,f]=freqz(b,1,1024,Fs);

%plot(f,20*log10(abs(h)))

magplot(h,Fs);
