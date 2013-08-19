function x = rifft(X)

%function x = rifft(X)
%

len = length(X);
X=[X ; imag(X(1)) ; conj(flipud(X(2:len)))];
x=real(ifft(X));

return;
