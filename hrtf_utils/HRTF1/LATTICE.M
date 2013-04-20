function [K,ec] = lattice(A)
% LATTICE [K,ec]=LATTICE(A) calculates lattice reflection coefficients 
%         K from polynomial coefficents given in A. 
%         The error condition flag ec is set to 1 if any of the 
%         reflection coefficients is equal to one.

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 03.07.1994 A. Raab

ec=0;
if min(size(A))>1 error('A must be a vector'); end;
if size(A,1)>1 A=A'; end;
K=zeros(size(A));
A=A./A(1);              % normalize polynomial coefficients
for m=length(A):-1:2
  K(m) = A(m);
  if K(m)==1 ec=1; break; end;  
  A=(A-K(m)*fliplr(A))./(1-K(m).^2);
  A=A(1:(m-1));
end;
K=fliplr(K(2:length(K)));


