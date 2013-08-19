function y=shift(x,n)

% SHIFT cyclic shift of vector x for n taps
%
% function y=shift(x,n)
% 
% positive values shift vector to the right
% negative values shift vcetor to the lift 
% 
% example:
%
% x=[ 1 2 3 4 5];
% y=shift(x,1);
% y
% y =
%
%     5     1     2     3     4
%
% y=shift(x,-1);
% y
% y =
%
%     2     3     4     5     1

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl für allgemeine Elektrotechnik und Akustik
% Ruhr-Universität Bochum
% Germany
%
% 18.09.1995


if (nargin == 1) n=1; end;


% x=x(:);
l=length(x);
n=-n;

if (n<0) n=rem(n,l); n=n+l; end;

if (n>=l) n=rem(n,l); end;


k(1:(l-n))=(1+n):l;
if (n~=0) k(l-n+1:l)=1:n; end;

y=x(k);

return;
