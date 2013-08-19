function Out = lp_filt(In,fs,tau);
% Low pass 
%   calculation with 'filter'
%
% Parameter:
%   In:  input signal
%   fs:  sampling frequency (in Hz)
%   tau: time constant (in ms)
%
% Result:
%   Out: output signal

% Author: Karsten Lehn
% Creation: 1.12.95
% Modification: 13.2.96

% File: lp_filt.m

dt = 1/fs;
c = exp(-dt/(tau / 1000));  
b = [1-c 0];
a = [1 -c];

Out = filter(b,a,In);
