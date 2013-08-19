
function X = rfft(x,len)

%function X = rfft(x,len)

if nargin<2 len=length(x); end;
X=fft(x,len);
X=X(1:length(X)/2);

return;
