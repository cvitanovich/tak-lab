function bark=hz2bark(f);
%HZ2BARK Convert Hz to Bark.
%        The frequency values in vector f are converted to
%        the corresponding Bark-values.

% Reference: Zwicker, E. and Terhard, E.
%            Analytical expression for critical-band rate
%            and critical bandwidth as a function of frequency.
%            J. Acoust. Soc. Am., 1980, Vol.68 (5),
%            p.1523-1525

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 05.07.1994 A. Raab

bark=13*atan(0.76*f/1000) + 3.5*(atan(f/7500)).^2;
