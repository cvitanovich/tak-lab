function fa=samplfrq()
% SAMPLEFRQ returns the default sampling frequency.
%           This function is used from nearly all other 
%           functions, so the samplerate could be changed 
%           easily. If a global variable SAMPLERATE is 
%           defined, SAMPLEFRQ returns the value of this 
%           variable.

% (c) Lehrstuhl fuer allgemeine Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 02.10.1994 A. Raab

if exist('SAMPLERATE')==1
  fa=SAMPLERATE;
else
  fa=44100;
end;ÿ
