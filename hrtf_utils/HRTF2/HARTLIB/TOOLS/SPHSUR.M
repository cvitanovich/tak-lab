function y=sphsur(n,m,phi,theta)

% SPHSUR calculates the spherical surface function? Kugelflaechenfunktion
%
% function y=sphsur(n,m,phi,theta)
% 
% 

p=legendre(n,cos(theta))
m=abs(m);
p(m+1)
y = sqrt((2*n+1)/(4*pi)*fac(n-m)/fac(n+m))*p(m+1)*exp(i*m*phi);

return
end

