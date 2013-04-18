function []=dispdata(filename);
% DISPDATA displays the measured impulse responses and transfer-
%          functions.
%          DISPDATA(filename);
%          PRELIMINARY VERSION !!

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
%     RUB
% (p) A. Raab, K. Hartung, 09.12.1994, 8.00 pm

samplerate=24000;
epsilon=1E-8;

[Hl,Hr,azi,ele]=readsc2(filename);
Hl=Hl(:);
Hr=Hr(:);

Hl=fade(Hl,171,271,839,939);
Hr=fade(Hr,171,271,839,939);
index=find(Hr==0); Hr(index)=ones(size(index))*epsilon;
index=find(Hl==0); Hl(index)=ones(size(index))*epsilon;

hl=rifft(Hl);
hr=rifft(Hr);


%hl=hl(1060:length(hl));
%hr=hr(1060:length(hr));
%hl=hl(1:150);
%hr=hr(1:150);

gtitle=['File: ' filename 'Azimuth: '];
gtitle=sprintf('file:%s, azimuth=%i, elevation=%i',filename,azi,ele);
irplot([hl(:) hr(:)],24000,gtitle,'gr');
legend('g','left','r','right',-1);
grcntrl

Hl = rfft(hl,2048);

Hr = rfft(hr,2048);

%load 802B_DTF
%Hr=Hr./Mmean(:,2);
%Hl=Hl./Mmean(:,1);

magplot([Hl(:) Hr(:)],samplerate,gtitle,'gr',3000,10000);
grcntrl
phasplot([Hl(:) Hr(:)],samplerate,gtitle,'gr',3000,10000);
grcntrl
