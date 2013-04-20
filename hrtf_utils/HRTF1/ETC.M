function etch = etc(h)
%ETC Energy time curve
%    ETC(h) returns the energy time curve of the (real) vector h.
 
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 27.06.1994 A. Raab

etch=abs(hilbert(h));
