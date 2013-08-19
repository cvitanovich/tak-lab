function f=checkfft(H);
% CHECKFFT Check complex transfer function
%          f=CHECKFFT(H);
%          The variable H is interpreted as a part of the FFT
%          of a real valued sequence. CHECKFFT tries to determine
%          whether H contains samples in the range 0...pi
%          (returns 1) or in the range 0...2*pi (returns 0)

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 30.08.1994 A. Raab

L=length(H)/2;
H=H(:);
f=~(all(H(2:L)==conj(flipud(H(L+2:2*L)))) ...
    & imag(H(1))==0 ...
    & imag(H(L+1))==0);
