function [hl, hr]=luxusfft(Hl,Hr,n,Mmean)

%function [hl, hr]=luxusfft(H,n,Mmean)

epsilon=1E-8;

Hl=fade(Hl,171,271,839,939);
Hr=fade(Hr,171,271,839,939);
index=find(Hr==0); Hr(index)=ones(size(index))*epsilon;
index=find(Hl==0); Hl(index)=ones(size(index))*epsilon;

hl=rifft(Hl);
hr=rifft(Hr);

[dl,tl]=h_max_en(hl,n);
[dr,tr]=h_max_en(hr,n);

t_start=min(tl,tr);

hl=hl(t_start:t_start+n);
hr=hr(t_start:t_start+n);

hil=iresponse(Mmean(:,1));
hir=iresponse(Mmean(:,2));

hl=conv(hl,hil);
hr=conv(hr,hir);

[dl,tl]=h_max_en(hl,n);
[dr,tr]=h_max_en(hr,n);
t_start=min(tl,tr);
hl=hl(t_start:t_start+n);
hr=hr(t_start:t_start+n);

return;

