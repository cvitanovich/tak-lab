function mtlextrem(H)

% function mtlextrem(H)

subplot(2,1,1);

i=1:1:length(H);
li=logspace(0,log10(length(H)),1000);

Hl=20*log10(abs(H));
Hl=interp1(i,Hl,li,'spline');
plot(Hl);
n=30
Hl=conv(Hl,boxcar(n)/n);
subplot(2,1,2);

plot(Hl(n:1000+n-1));
d_Hl=del2(Hl(n:1000+n-1));

figure;
plot(d_Hl(30:800));
e=find(d_Hl==0)
