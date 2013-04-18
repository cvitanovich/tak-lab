function pd = power(x);

% POWER     pseudopower of vector 
% 
% power(x) returns pseudo-power of x 

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektortechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
% Date 25.09.1995
% Modified: 7.7.96
%        improved help
%  

% Calls		toolbox
%
% no calls to other toolboxes



if isreal(x) 
   pd = x .* x; 
  else
   pd = x .* conj(x);
end

% end of file power.m