function y=smooth(x,M);
%SMOOTH y=SMOOTH(x,M) 
%       Smooth data vector x using M values on the left and 
%	right of each point.

% (p) 10.07.1994 A. Raab

LX=length(x);
y=zeros(size(x));
y(1:M+1)=x(1:M+1); y(LX-M+1:LX)=x(LX-M+1:LX);
for k=M+1:LX-M
  y(k)=sum(x(k-M:k+M))/(2*M+1);
end;
  
  