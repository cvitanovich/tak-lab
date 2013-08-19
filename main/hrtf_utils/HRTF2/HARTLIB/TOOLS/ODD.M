function a=odd(x)

% ODD(X) is odd
%
% odd(x) returns 1 for odd integers and 0 for even integers.
% Zero is an even number.

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%
% Date 25.09.1995
% Modified: 8.7.96 Klaus Hartung
%     improved help

% Calls 	toolbox
% 
% no calls

x = round(x);

if ( floor(abs(x/2))*2 ~= x ) a=1;
  else a=0;
end;


% end of file odd.m
