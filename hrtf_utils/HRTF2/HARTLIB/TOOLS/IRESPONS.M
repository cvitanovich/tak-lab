function [hi]=iresponse(M);

epsilon=1E-2;
X=ones(size(M))./M;

X=fade(X,171,271,839,880);
index=find(X==0); X(index)=ones(size(index))*epsilon;
hi=rifft(X);
hi=fftshift(hi);
