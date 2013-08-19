function b = stable(a);
%STABLE Polynomial stabilization.
%       STABLE(A), where A is a vector of polynomial coefficients,
%       stabilizes the polynomial with respect to the unit circle;
%       roots whose magnitudes are greater than one are reflected
%       inside the unit circle.

%       Copyright (c) 1985-1988 by the MathWorks, Inc.
%       Revised 7-25-89 JNL (Handles roots at zero)

%       Changed 22.05.1994 by A. Raab
%       in order to keep "magnitude" of the polynomial
%       (Multiply polynomial coefficients by product of the
%       magnitude all roots outside the unit circle)

if length(a) == 1, b = a; return, end

v = roots(a);
% section added by A.Raab
i = find(abs(v)>1);
w = prod(abs(v(i)));
%
i = find(v~=0);
vs = 0.5*(sign(abs(v(i))-1)+1);
v(i) = (1-vs).*v(i) + vs./conj(v(i));
% section changed by A.Raab
b = w*a(1)*poly(v);

% Return only real coefficients if input was real:
if ~any(imag(a))
  b = real(b);
end
