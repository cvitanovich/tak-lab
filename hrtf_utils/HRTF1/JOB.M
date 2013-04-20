%function [hil,hir]=iresponse(Ml,Mr);

epsilon=1E-2;
load 802B_DTF
X=ones(size(Mmean))./Mmean;

for j=1:2
  %X(:,j)=fade(X(:,j),171,271,839,939);
X(:,j)=fade(X(:,j),220,271,839,880);

  index=find(X(:,j)==0); X(index,j)=ones(size(index))*epsilon;
  h=rifft(X(:,j));
  h=fftshift(h);
  plot(h(900:1100));
  H=fft(h,2048); H=H(1:1024);
  magplot(H(:),24000);
  pause;
end;

