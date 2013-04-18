function signal = abs2db(samples)

% function signal = abs2db(samples)
%
% ABS2DB changes the complex elements or elements with absolute magnitude to
%        magnitude with logarithm scaling in [dB]. To provide logarithm of
%        zero, magnitudes of zero are changed to 1e-6 before logarithm.
%        - samples: matrix or vector of complex or real datas with absolute
%                   magnitude
%        - signal:  matrix or vector of datas with magnitude in [dB]
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 09.09.1994 Joerg Schmatter

signal = abs(samples);   % magnitude of the samples
[i,j] = find(signal==0);
for k=1:length(i)
  signal(i(k),j(k)) = 1e-6;
end;
signal = 20*log10(signal);