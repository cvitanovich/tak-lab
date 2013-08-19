function result = gauss(value)

% function result = gauss(value)
%
% GAUSS returns the next integer value that is higher or greater than the
%       given number 'value'.
%
% (c) Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% (p) 24.09.1994 Joerg Schmatter


difference = fix(value) - value;
if (difference~=0)
  result = fix(value) + 1;
else
  result = fix(value);
end;