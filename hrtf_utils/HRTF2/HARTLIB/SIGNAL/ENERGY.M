function ex = energy(x);
%ENERGY Energy of vector x.
%       ENERGY(x) returns the energy of vector x.
%       (the summed square of all elements).

% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 27.06.1994 A. Raab

x = abs(x(:));
ex = x'*x;
