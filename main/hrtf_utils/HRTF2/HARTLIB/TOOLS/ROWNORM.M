function y=rownorm(x)

% ROWNORM Normalizes each row of a matrix
%
% y=rownorm(x)

% (C) Klaus Hartung
% Ruhr-Universit�t Bochum
% Lehrstuhl f�r allgemeine Eletrotechnik und Akustik
% 17.8.1995

[r,c]=size(x);

y=zeros(r,c);

[m,p]=max(abs(x'));
for n=1:r
  y(n,:)=x(n,:)/m(n);
end

return
end
