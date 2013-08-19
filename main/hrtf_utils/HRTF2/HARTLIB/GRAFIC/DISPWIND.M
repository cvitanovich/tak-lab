function []=dispwind(filename,t);
% DISPWIND displays the measured impulse responses and transfer-
%          functions.
%          DISPDATA(filename,t);
%          estimated length of impulse response
%          PRELIMINARY VERSION !!

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
%     RUB
% (p) A. Raab, K. Hartung, 09.12.1994, 8.00 pm
%  V1.2 K. Hartung 13.12.94 2.30 pm

if (nargin < 2) t = 2.5 ;end;
samplerate=24000;
epsilon=1E-8;
n=t*samplerate/1000.0;

[Hl,Hr,azi,ele]=readsc2(filename);
Hl=Hl(:);
Hr=Hr(:);

load 802B_DTF
[hl,hr] = luxusfft(Hl,Hr,n,Mmean);


gtitle=['File: ' filename 'Azimuth: '];
gtitle=sprintf('file:%s, azimuth=%i, elevation=%i',filename,azi,ele);
irplot([hl(:) hr(:)],24000,gtitle,'gr');
legend('g','left','r','right',-1);
grcntrl
prcntrl

Hl = rfft(hl,2048);
Hr = rfft(hr,2048);


magplot([Hl(:) Hr(:)],samplerate,gtitle,'g-r-.',3000,10000);
grcntrl
prcntrl

ax=axis;
ax(3)=-20; ax(4)=10; 
axis(ax);

gtitle=sprintf('file:%s, azimuth=%i, elevation=%i (interaural)',filename,azi,ele);
Hinteraural=Hl(:)./Hr(:);
Hinteraural(1:171)=abs(real(Hinteraural(1:171)));
magplot(Hinteraural,samplerate,gtitle,'c',3000,10000);
grcntrl
prcntrl

ax=axis;
ax(3)=-20; ax(4)=20; 
axis(ax);


tphplot(Hinteraural,samplerate,gtitle,'c',3000,10000);
grcntrl
prcntrl
ax=axis;
ax(3)=-0.3E-3; ax(4)=0.3E-3; 
axis(ax);

