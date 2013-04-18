function a=even(x)

% EVEN(X) 	is even 
%
% 	 even(x) returns 1 for even integers and 0 for odd integers.
%        Zero is an even number.

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektrotechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%
% 
% Date 25.09.1995
% Modified: 8.7.96 Klaus Hartung
%     improved help

% Calls 	Toolbox
%
% No calls

if (nargin<1) error('No Input defined'); end;

x = round(x);

if ( floor(abs(x/2))*2 ~= x ) a=0;
  else a=1;
end;


% end of file even.m
