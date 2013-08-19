function a=fexist(filename)

% FEXIST  check if a file exists
%
%         fexist(filename) checks if a file exists 
%         and returns 1 if true and 0 if not.
%

% Klaus Hartung (hartung@aea.ruhr-uni-bochum.de)
% Lehrstuhl fuer allg. Elektortechnik und Akustik
% Ruhr-Universitaet Bochum
% 44780 Bochum 
%
% Date: 6.7.96  

% Calls		toolbox
%
% no calls to other toolboxes

if (nargin < 1) error('Insufficient number of arguments'); end;
fid = fopen(filename,'r');


if (fid>0)
  a=1;
  fclose(fid);
else
  a=0;
end;



end; 


% end of fexist

